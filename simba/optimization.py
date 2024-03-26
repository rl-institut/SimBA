""" Collection of procedures optimizing arbitrary parameters of a bus schedule or infrastructure.
"""
from copy import deepcopy
import datetime
import logging

from simba.rotation import Rotation
from simba.trip import Trip


# defaults
DEFAULT_DEPOT_DISTANCE = 5  # km
DEFAULT_MEAN_SPEED = 30  # km/h


def service_optimization(schedule, scenario, args):
    """ Optimize rotations based on feasability.
        Try to find sets of rotations that produce no negative SoC

    :param schedule: Schedule to be optimized
    :type schedule: simba.Schedule
    :param scenario: Simulated schedule
    :type scenario: spice_ev.Scenario
    :param args: Command line arguments
    :type args: argparse.Namespace
    :return: original and most optimized scenario (highest electrification rate)
    :rtype: dict
    """
    common_stations = schedule.get_common_stations(only_opps=True)

    original = (deepcopy(schedule), deepcopy(scenario))
    optimal = (None, None)

    # single out negative rotations. Try to run these with common non-negative rotations
    negative_rotations = schedule.get_negative_rotations(scenario)
    logging.info(f"Initially, rotations {sorted(negative_rotations)} have neg. SoC.")

    if not negative_rotations:
        return {
            "original": original,
            "optimized": original,
        }

    negative_sets = {}
    for rot_key in negative_rotations:
        logging.debug(f"Looking at negative rotation {rot_key}")
        # remove negative rotation from initial schedule
        rotation = schedule.rotations.pop(rot_key)
        if rotation.charging_type != "oppb":
            # only oppb rotations are optimized -> skip others
            logging.warning(f"Rotation {rot_key} should be optimized, "
                            f"but is of type {rotation.charging_type}.")
            continue
        # oppb: build non-interfering sets of negative rotations
        # (these include the dependent non-negative rotations)
        s = {rot_key}
        vid = rotation.vehicle_id
        last_neg_soc_time = scenario.negative_soc_tracker[vid][-1]
        last_neg_soc_time = datetime.datetime.fromisoformat(last_neg_soc_time)
        dependent_rotations = {r: t for r, t in common_stations[rot_key].items()
                               if t <= last_neg_soc_time}
        logging.debug(f"Dependent stations: {dependent_rotations}")
        while dependent_rotations:
            # get next rotation ID and last common time
            r, t = dependent_rotations.popitem()
            logging.debug(f"\tRotation {r} @ {t}:")
            if r not in negative_rotations:
                # r not negative: add to set of dependent rotations
                s.add(r)
                # add dependencies of r, avoid circular dependency
                dependencies = {r2: t2 for r2, t2 in common_stations[r].items()
                                if t2 <= t and r2 not in s}
                dependent_rotations.update(dependencies)
                logging.debug(f"\t{dependencies}")
            elif original.rotations[r].charging_type != "obbp":
                logging.warning(f"Rotation {rot_key} depends on negative non-oppb rotation")

        negative_sets[rot_key] = s
        logging.debug(f"Negative sets: {negative_sets}")

    # run scenario with non-negative rotations only
    if schedule.rotations:
        scenario = schedule.run(args)
        optimal = (deepcopy(schedule), deepcopy(scenario))

    # run singled-out negative rotations
    ignored = []
    for i, (rot, s) in enumerate(negative_sets.items()):
        schedule.rotations = {r: original[0].rotations[r] for r in s}
        logging.debug(f"{i+1} / {len(negative_sets)} negative schedules: {rot}")
        scenario = schedule.run(args)
        if scenario.negative_soc_tracker:
            # still fail: try just the negative rotation
            schedule.rotations = {rot: original[0].rotations[rot]}
            scenario = schedule.run(args)
            if scenario.negative_soc_tracker:
                # no hope, this just won't work
                logging.info(f"Rotation {rot} will stay negative")
            else:
                # works alone with other non-negative rotations
                logging.info(f"Rotation {rot} works alone")
            ignored.append(rot)

    negative_sets = {k: v for k, v in negative_sets.items() if k not in ignored}

    # now, in negative_sets are just rotations that work with others still
    # try to combine them
    possible = [(a, b) for a in negative_sets for b in negative_sets if a != b]
    while possible:
        logging.debug(f"{len(possible)} combinations remain")
        r1, r2 = possible.pop()
        combined = negative_sets[r1].union(negative_sets[r2])
        schedule.rotations = {r: original[0].rotations[r] for r in combined}
        scenario = schedule.run(args)
        logging.debug(f"{r1} + {r2}: {not scenario.negative_soc_tracker}")
        if not scenario.negative_soc_tracker:
            # compatible (don't interfere): keep union, remove r2
            negative_sets[r1] = combined
            # simplified. What about triangle? (r1+r2, r1+r3, r2!+r3)?
            possible = [t for t in possible if r2 not in t]
            logging.debug(possible)

            # save scenario with highest electrification rate
            if optimal is None or len(scenario[0].rotations) > len(optimal[0].rotations):
                optimal = (deepcopy(schedule), deepcopy(scenario))

    logging.info("Service optimization finished")

    return {
        "original": original,
        "optimized": optimal,
    }


def split_trips(schedule):
    """
    Split all rotations into their individual trips.

    Preparation for recombination. Extracts all trips from or to a depot into their own dictionary.
    :param schedule: Schedule to be optimized
    :type schedule: simba.Schedule
    :return: individual trips between stations (no depot trips), dict of depot trips
    :rtype: tuple
    """
    trips = []
    depot_trips = {}  # take note of trips from/to depots. station -> depot -> trip
    # iterate over all trips of schedule
    for rot in schedule.rotations.values():
        for trip_idx, trip in enumerate(rot.trips):
            if trip.delta_soc > 1:
                logging.info(f"Trip {trip_idx} of rotation {rot.id} too long, discarded")
                continue
            # check if trip starts or ends in depot
            departure_station = schedule.stations.get(trip.departure_name)
            arrival_station = schedule.stations.get(trip.arrival_name)
            depot_station = None
            trip_station = None
            if departure_station is not None and departure_station["type"] == "deps":
                depot_station = trip.departure_name
                trip_station = trip.arrival_name
            elif arrival_station is not None and arrival_station["type"] == "deps":
                depot_station = trip.arrival_name
                trip_station = trip.departure_name
            if depot_station is None:
                # not a depot trip: add to trips for recombination
                trips.append(trip)
            else:
                # trip starts or ends in depot: save trip for later, but do not add to trips
                if trip_station not in depot_trips:
                    depot_trips[trip_station] = {}
                try:
                    old_trip = depot_trips[trip_station][depot_station]
                    # trip to depot already known:
                    # must have same distance and duration
                    # consumption can be different, as it depends on environment
                    if trip.distance != old_trip.distance:
                        logging.error(f"Depot trips {depot_station} - {trip_station} "
                                      f"differ in distance: {old_trip.distance} / {trip.distance}")
                    t1 = (trip.arrival_time - trip.departure_time).total_seconds() // 60
                    t2 = (old_trip.arrival_time - old_trip.departure_time).total_seconds() // 60
                    if t1 != t2:
                        logging.error(f"Depot trips {depot_station} - {trip_station} "
                                      f"differ in duration: {t1} / {t2}")
                except KeyError:
                    # trip to depot not yet known
                    if depot_trips[trip_station]:
                        # but there is a trip to another depot
                        logging.info(f"Multiple depot trips for {trip_station}")
                    depot_trips[trip_station][depot_station] = trip
    return (trips, depot_trips)


def generate_depot_trip(station, depot_trips):
    """
    Look up closest depot to a station.

    Assumes default if no trip between this station and a depot exists in the original dataset.
    :param station: starting location
    :type station: str
    :param depot_trips: look-up-table for depot trips
    :type depot_trips: dict
    :return: name and distance to closest depot
    :rtype: dict
    """
    try:
        # get depot with shortest distance
        depots = sorted(depot_trips[station].items(), key=lambda x: x[1].distance)
        depot_trip = depots[0]
        return {
            "name": depot_trip[0],
            "distance": depot_trip[1].distance,
        }
    except KeyError:
        # no trip from depot to departure station known: assume defaults
        assert len(depot_trips) > 0, "No depot trips known"
        return {
            "name": list(depot_trips)[0],  # get first depot
            "distance": DEFAULT_DEPOT_DISTANCE * 1000,
        }


def recombination(schedule, args, trips, depot_trips):
    """
    Recombine individual trips to new rotations.

    Gathers rotations greedy from given trips. Checks if rotation can be done at all.
    :param schedule: Schedule to be recombined
    :type schedule: simba.Schedule
    :param args: Command line arguments
    :type args: argparse.Namespace
    :param trips: individual trips
    :type trips: list of simba.Trip
    :param depot_trips: look-up-table for depot trips
    :type depot_trips: dict
    :return: recombined schedule
    :rtype: simba.Schedule
    """
    rotations = {}
    # sort trips by departure time
    trips.sort(key=lambda t: t.departure_time)

    while trips:
        # still trips to be done: new rotation, get earliest trip
        trip = trips.pop(0)
        rotation = Rotation(
            id=f"recombination_{len(rotations) + 1}",
            vehicle_type=trip.rotation.vehicle_type,
            schedule=schedule,
        )
        vehicle_type = rotation.vehicle_type
        charging_type = trip.rotation.charging_type
        rotation.set_charging_type(charging_type)

        # begin rotation in depot: add initial trip
        depot_trip = generate_depot_trip(trip.departure_name, depot_trips)
        depot_trip_time = depot_trip["distance"] / (DEFAULT_MEAN_SPEED * 1000)
        rotation.add_trip({
            "departure_time": trip.departure_time - datetime.timedelta(hours=depot_trip_time),
            "departure_name": depot_trip["name"],
            "arrival_time": trip.departure_time,
            "arrival_name": trip.departure_name,
            "distance": depot_trip["distance"],
            "line": trip.line,
            "charging_type": charging_type,
            "temperature": trip.temperature,
            # "height_difference": None,  # compute from station data
            "level_of_loading": 0,  # no passengers from depot
            "mean_speed": DEFAULT_MEAN_SPEED,
            "station_data": schedule.station_data,
        })

        # add first trip to rotation
        # use info from original trip
        trip_dict = vars(trip)
        del trip_dict["rotation"]
        trip_dict["charging_type"] = charging_type
        rotation.add_trip(trip_dict)

        # initial check: can bus reach depot again?
        soc = args.desired_soc_deps
        # trip.calculate_consumption does not return delta_soc, just consumption:
        # use Trip.consumption.calculate_consumption directly (should exist, because of initial sim)
        for trip in rotation.trips:
            _, delta_soc = Trip.consumption.calculate_consumption(
                trip.arrival_time, trip.distance, vehicle_type, charging_type,
                trip.temperature, trip.height_diff, trip.level_of_loading, trip.mean_speed)
            soc += delta_soc
        # soc now contains battery level after first proper trip
        # it will be updated with every trip added to the new rotation
        # probe return trip (same depot)
        try:
            height_diff = (
                schedule.station_data[trip.arrival_name]["elevation"]
                - schedule.station_data[depot_trip["name"]]["elevation"])
        except (KeyError, ValueError):
            height_diff = 0
        _, delta_soc = Trip.consumption.calculate_consumption(
            trip.departure_time, depot_trip["distance"], vehicle_type, charging_type,
            trip.temperature, height_diff, level_of_loading=0, mean_speed=DEFAULT_MEAN_SPEED)
        if -delta_soc > soc:
            logging.info(f"Trip of line {trip.line} between {trip.departure_name} "
                         f"and {trip.arrival_name} too long with added depot trips, discarded")
            # ignore new rotation
            continue

        # initial trip can come back to depot: make rotation official
        rotations[rotation.id] = rotation
        last_depot_trip = depot_trip

        # add further trips
        trip_idx = 0
        arrival_time = trip.arrival_time
        current_name = trip.arrival_name
        while trip_idx < len(trips):
            trip = trips[trip_idx]
            # find next trip from same station with same vehicle type
            if (
                    trip.departure_name != current_name or
                    trip.rotation.vehicle_type != vehicle_type or
                    arrival_time > trip.departure_time
            ):
                trip_idx += 1
                continue

            # can bus come back from this station?
            trip = trips[trip_idx]
            _, delta_soc_trip = Trip.consumption.calculate_consumption(
                trip.departure_time, trip.distance, vehicle_type, charging_type,
                trip.temperature, trip.height_diff, trip.level_of_loading, trip.mean_speed)

            depot_trip = generate_depot_trip(trip.arrival_name, depot_trips)
            depot_trip_time = depot_trip["distance"] / (DEFAULT_MEAN_SPEED * 1000)
            try:
                height_diff = (
                    schedule.station_data[trip.arrival_name]["elevation"]
                    - schedule.station_data[depot_trip["name"]]["elevation"])
            except (KeyError, ValueError):
                height_diff = 0
            _, delta_soc_depot = Trip.consumption.calculate_consumption(
                trip.departure_time, depot_trip["distance"], vehicle_type, charging_type,
                trip.temperature, height_diff, level_of_loading=0, mean_speed=DEFAULT_MEAN_SPEED)
            if soc < -(delta_soc_trip + delta_soc_depot):
                # can't return from this station: skip
                trip_idx += 1
                continue

            # add trip to rotation
            trip = trips.pop(trip_idx)
            trip_dict = vars(trip)
            del trip_dict["rotation"]
            trip_dict["charging_type"] = charging_type
            rotation.add_trip(trip_dict)
            soc += delta_soc_trip
            last_depot_trip = depot_trip

            arrival_time = trip.arrival_time
            current_name = trip.arrival_name
        else:
            # add final depot trip after no more possible trips remain
            last_trip = rotation.trips[-1]
            depot_trip_time = last_depot_trip["distance"] / (DEFAULT_MEAN_SPEED * 1000)
            rotation.add_trip({
                "departure_time": arrival_time,
                "departure_name": current_name,
                "arrival_time": arrival_time + datetime.timedelta(hours=depot_trip_time),
                "arrival_name": last_depot_trip["name"],
                "distance": last_depot_trip["distance"],
                "line": last_trip.line,
                "charging_type": charging_type,
                "temperature": last_trip.temperature,
                "level_of_loading": 0,  # no passengers to depot
                "mean_speed": DEFAULT_MEAN_SPEED,
                "station_data": schedule.station_data,
            })
    schedule.rotations = rotations
    schedule.calculate_consumption()
    return schedule
