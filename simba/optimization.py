""" Collection of procedures optimizing arbitrary parameters of a bus schedule or infrastructure """
from copy import deepcopy
import datetime
import logging

from simba.rotation import Rotation
from simba.trip import Trip


def service_optimization(schedule, scenario, args):
    """ Optimize rotations based on feasability.

    Try to find sets of rotations that produce no negative SoC.

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


def prepare_trips(schedule, negative_rotations=None):
    """ Prepares schedule for recombination.

    Find trips to/from depots, prepare lists of trips to recombine.

    :param schedule: Schedule to be optimized
    :type schedule: simba.Schedule
    :param negative_rotations: ID of negative rotations. May be None (recombine all depb rotations)
    :type negative_rotations: iterable
    :return: rotation ID -> trip list (without depot trips)
    :rtype: dict
    """
    trips = {}  # structure: rotation ID -> trip list
    depot_trips = {}  # take note of trips from/to depots. station -> depot -> trip
    for rot_id, rot in schedule.rotations.items():
        rotation_trips = sorted(rot.trips, key=lambda t: t.departure_time)
        # make note of trips from negative depot rotations (except initial/final trip)
        qualified_id = negative_rotations is None or rot_id in negative_rotations
        if rot.charging_type == "depb" and qualified_id:
            trips[rot_id] = rotation_trips[1:-1]

        # assumption: rotation ends at same depot where it starts at
        if rotation_trips[0].departure_name != rotation_trips[-1].arrival_name:
            logging.warning(f"Rotation {rot_id} has different start and end depots")

        # check all trips of rotation: gather all trips from/to depots for later reference
        for trip in rotation_trips:
            # does trip start or end in depot?
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
            if depot_station is not None:
                # trip starts or ends in depot: save trip in depot_trips
                add_depot_trip(trip_station, depot_station, trip, depot_trips)
                # dummy trip from depot to same depot for lookup
                depot_trip = Trip(
                    rotation=rot,
                    departure_time=trip.arrival_time,
                    departure_name=depot_station,
                    arrival_time=trip.arrival_time,
                    arrival_name=depot_station,
                    distance=0.1,  # must not be zero?
                    line=None,
                    temperature=trip.temperature,
                    height_difference=0,
                    level_of_loading=0,
                    mean_speed=18,
                )
                add_depot_trip(depot_station, depot_station, depot_trip, depot_trips)
        # all trips done
    return (trips, depot_trips)


def add_depot_trip(station_name, depot_name, trip, depot_trips):
    if station_name in depot_trips:
        if depot_name in depot_trips[station_name]:
            # trip from station to depot already known: use trip with lesser distance
            if trip.distance < depot_trips[station_name][depot_name].distance:
                depot_trips[station_name][depot_name] = trip
        else:
            # trip from station to depot unknown: save trip as reference
            depot_trips[station_name][depot_name] = trip
    else:
        # no trip to any depot known for this station: save trip as reference
        depot_trips[station_name] = {depot_name: trip}
    return depot_trips


def generate_depot_trip_data_dict(
        station, depot, depot_trips,
        default_depot_distance, default_mean_speed):
    """
    Look up trip from station to depot and give info about that trip.

    Returns name of depot, distance in km, mean speed in km/h and travel time as timedelta.
    Assumes defaults if no trip between this station and the depot exists in the original dataset.
    :param station: starting location
    :type station: str
    :param depot: depot name to drive from/to
    :type depot: str
    :param depot_trips: look-up-table for depot trips
    :type depot_trips: dict
    :param default_depot_distance: default assumed average distance from any station to a depot [km]
    :type default_depot_distance: numeric
    :param default_mean_speed: default assumed average speed for busses [km/h]
    :type default_mean_speed: numeric
    :return: distance, speed and travel time to given depot
    :rtype: dict
    """
    try:
        # get depot trip if it exists
        depot_trip = depot_trips[station][depot]
        return {
            "distance": depot_trip.distance,  # in meters
            "mean_speed": depot_trip.mean_speed,  # in km/h
            "travel_time": depot_trip.arrival_time - depot_trip.departure_time,  # timedelta
        }
    except KeyError:
        # no trip from station to depot known: assume defaults
        return {
            "distance": default_depot_distance * 1000,
            "mean_speed": default_mean_speed,
            "travel_time": datetime.timedelta(hours=default_depot_distance/default_mean_speed),
        }


def recombination(schedule, args, trips, depot_trips):
    """
    Recombine individual trips to new rotations.

    Gathers rotations greedy from given trips.
    Rotations are not mixed, keep trips together if possible.
    Checks if rotation can be done at all.

    :param schedule: Schedule to be recombined
    :type schedule: simba.Schedule
    :param args: Command line arguments
    :type args: argparse.Namespace
    :param trips: look-up-table depb rotation ID -> trip list
    :type trips: dict
    :param depot_trips: look-up-table for depot trips
    :type depot_trips: dict
    :return: recombined schedule
    :rtype: simba.Schedule
    """
    for rot_id, trip_list in trips.items():
        # remove original rotation to recombine
        original_rotation = schedule.rotations.pop(rot_id)
        # end rotation at same depot again
        depot_name = original_rotation.trips[0].departure_name

        # examine single rotation: which trips can be made safely with a single battery charge?
        # generate initial trip from depot to first station,
        # then add trips from original rotation until soc is too low to reach a depot again

        # sort trips by departure time
        trip_list.sort(key=lambda t: t.departure_time)

        # initial values
        rotation = None  # new rotation is generated during first trip
        rot_name = f"{rot_id}_r"  # new rotation name, if it needs to be split a counter is added
        rot_counter = 0  # into how many rotations has the original already been split?
        soc = None  # remaining soc of new rotation
        last_depot_trip = None  # keep track of last possible depot trip

        while trip_list:
            # probe next trip
            trip = trip_list.pop(0)
            if rotation is None:
                # first trip of new rotation: generate new rotation
                rotation = Rotation(
                    id=rot_name,
                    vehicle_type=trip.rotation.vehicle_type,
                    schedule=schedule,
                )
                charging_type = trip.rotation.charging_type
                rotation.set_charging_type(charging_type)

                # begin rotation in depot: add initial trip
                depot_trip = generate_depot_trip_data_dict(
                    trip.departure_name, depot_name, depot_trips,
                    args.default_depot_distance, args.default_mean_speed)
                rotation.add_trip({
                    "departure_time": trip.departure_time - depot_trip["travel_time"],
                    "departure_name": depot_name,
                    "arrival_time": trip.departure_time,
                    "arrival_name": trip.departure_name,
                    "distance": depot_trip["distance"],
                    "line": trip.line,
                    "charging_type": charging_type,
                    "temperature": trip.temperature,
                    # "height_difference": None,  # compute from station data
                    "level_of_loading": 0,  # no passengers from depot
                    "mean_speed": depot_trip["mean_speed"],
                    "station_data": schedule.station_data,
                })

                # calculate consumption for initial trip
                soc = args.desired_soc_deps  # vehicle leaves depot with this soc
                rotation.calculate_consumption()
                soc += rotation.trips[0].delta_soc  # new soc after initial trip

                if rot_counter > 0:
                    logging.info(
                        f'Rotation {rot_id}: New Einsetzfahrt '
                        f'{trip.departure_time - depot_trip["travel_time"]} '
                        f'{depot_name} - {trip.departure_name} '
                        f'({depot_trip["travel_time"].total_seconds() / 60} min, '
                        f'{depot_trip["distance"] / 1000} km)')

            # can trip be completed while having enough energy to reach depot again?
            # probe return trip
            depot_trip = generate_depot_trip_data_dict(
                trip.arrival_name, depot_name, depot_trips,
                args.default_depot_distance, args.default_mean_speed)
            depot_trip = {
                "departure_time": trip.arrival_time,
                "departure_name": trip.arrival_name,
                "arrival_time": trip.arrival_time + depot_trip["travel_time"],
                "arrival_name": depot_name,
                "distance": depot_trip["distance"],
                "line": trip.line,
                "charging_type": charging_type,
                "temperature": trip.temperature,
                "level_of_loading": 0,
                "mean_speed": depot_trip["mean_speed"],
                "station_data": schedule.station_data,
            }
            # rotation.add_trip needs dict, but consumption calculation is better done on Trip obj:
            # create temporary depot trip object for consumption calculation
            tmp_trip = Trip(rotation, **depot_trip)
            tmp_trip.consumption, tmp_trip.delta_soc = tmp_trip.calculate_consumption()
            if soc >= -(trip.delta_soc + tmp_trip.delta_soc):
                # next trip is possible: add trip, use info from original trip
                trip_dict = vars(trip)
                del trip_dict["rotation"]
                trip_dict["charging_type"] = charging_type
                rotation.add_trip(trip_dict)
                last_depot_trip = depot_trip
                soc += trip.delta_soc
            else:
                # next trip is not safely possible
                if last_depot_trip is None:
                    # this was initial trip: discard rotation
                    logging.info(f"Trip of line {trip.line} between {trip.departure_name} and "
                                 f"{trip.arrival_name} too long with added depot trips, discarded")
                else:
                    # not initial trip: add last possible depot trip
                    rotation.add_trip(last_depot_trip)
                    assert last_depot_trip["arrival_name"] == depot_name
                    travel_time = last_depot_trip["arrival_time"]-last_depot_trip["departure_time"]
                    logging.info(
                        f'Rotation {rot_id}: New Aussetzfahrt '
                        f'{last_depot_trip["departure_time"]} '
                        f'{last_depot_trip["departure_name"]} - {last_depot_trip["arrival_name"]} '
                        f'({travel_time.total_seconds() / 60} min, '
                        f'{last_depot_trip["distance"] / 1000} km)')
                    schedule.rotations[rotation.id] = rotation
                    last_depot_trip = None
                    rot_counter += 1
                    rot_name = f"{rot_id}_r_{rot_counter}"
                    # re-add current trip to trip list for next Einsetzfahrt
                    trip_list = [trip] + trip_list
                # prepare for new rotation with next trip
                rotation = None
            # end of trip
        # rotation finished (trip list empty): add final depot trip, if needed
        if last_depot_trip is not None:
            rotation.add_trip(last_depot_trip)
            schedule.rotations[rotation.id] = rotation

    schedule.calculate_consumption()
    return schedule
