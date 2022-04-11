import csv
import datetime
from datetime import timedelta
import json
from os import path
import random

from ebus_toolbox.rotation import Rotation


class Schedule:

    def __init__(self, vehicle_types) -> None:
        """Constructs Schedule object from CSV file containing all trips of schedule

        :param vehicle_types: Collection of vehicle types and their properties.
        :type vehicle_types: dict
        """
        # Check if all bus types have both an opp and depot version
        # Also make sure that both versions have the same mileage
        vehicle_type_names = list(vehicle_types.keys())
        for name in vehicle_type_names:
            try:
                base, ct = name.rsplit('_', 1)
            except ValueError:
                continue
            if f"{base}_oppb" in vehicle_types and ct == 'dep':
                assert vehicle_types[name]["mileage"] == vehicle_types[f"{base}_oppb"]["mileage"]
            elif f"{base}_depb" in vehicle_types and ct == 'opp':
                assert vehicle_types[name]["mileage"] == vehicle_types[f"{base}_depb"]["mileage"]
        self.vehicle_types = vehicle_types

        self.rotations = {}
        self.consumption = 0

    @classmethod
    def from_csv(cls, path_to_csv, vehicle_types):
        """Constructs Schedule object from CSV file containing all trips of schedule.

        :param path_to_csv: Path to csv file containing trip data
        :type path_to_csv: str
        :param vehicle_types: Collection of vehicle types and their properties.
        :type vehicle_types: dict
        :return: Returns a new instance of Schedule with all trips from csv loaded.
        :rtype: Schedule
        """
        schedule = cls(vehicle_types)

        with open(path_to_csv, 'r') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                rotation_id = trip['rotation_id']
                if rotation_id not in schedule.rotations.keys():
                    schedule.rotations.update({
                        rotation_id: Rotation(id=rotation_id, vehicle_type=trip['vehicle_type'])})
                schedule.rotations[rotation_id].add_trip(trip)

        return schedule

    def remove_rotation_by_id(self, id):
        """ Removes a rotation from schedule

        :param id: Rotation ID to be removed
        :type id: int
        """
        del self.rotations[id]

    def filter_rotations(self):
        """Based on a given filter definition (tbd), rotations will be dropped from schedule."""
        pass

    def set_charging_type(self, preferred_ct, args, rotation_ids=None):
        """ Change charging type of either all or specified rotations. Adjust minimum standing time
            at depot after completion of rotation.

        :param preferred_ct: Choose this charging type wheneever possible. Either 'depb' or 'oppb'.
        :type preferred_ct: str
        :param args: Command line arguments and/or arguments from config file.
        :type args: argparse.Namespace
        :param rotation_ids: IDs of rotations for which to set charging type. If None set charging
                             charging type for all rotations.
        :type rotation_ids: list
        """
        assert preferred_ct in ["oppb", "depb"], f"Invalid charging type: {preferred_ct}"
        if rotation_ids is None:
            rotation_ids = self.rotations.keys()

        for id in rotation_ids:
            rot = self.rotations[id]
            vehicle_type = self.vehicle_types[f"{rot.vehicle_type}_{rot.charging_type}"]
            capacity = vehicle_type["capacity"]
            if preferred_ct == "oppb" or capacity < rot.consumption:
                self.rotations[id].charging_type = "oppb"
                min_standing_time = \
                    (capacity / args.cs_power_deps_oppb) * args.min_recharge_deps_oppb
            else:
                self.rotations[id].charging_type = "depb"
                min_standing_time = (rot.consumption / args.cs_power_deps_depb)
                desired_max_standing_time = \
                    (capacity / args.cs_power_deps_depb) * args.min_recharge_deps_oppb
                if min_standing_time > desired_max_standing_time:
                    min_standing_time = desired_max_standing_time

            rot.earliest_departure_next_rot = rot.arrival_time + timedelta(hours=min_standing_time)

    def assign_vehicles(self):
        """ Assign vehicle IDs to rotations. A FIFO approach is used.
        For every rotation it is checked whether vehicles with matching type are idle, in which
        case the one with longest standing time since last rotation is used.
        If no vehicle is available a new vehicle ID is generated.
        """
        rotations_in_progress = []
        idle_vehicles = []
        vehicle_type_counts = {vehicle_type: 0 for vehicle_type in self.vehicle_types.keys()}

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)

        for rot in rotations:
            # find vehicles that have completed rotation and stood for a minimum staning time
            # mark those vehicle as idle
            for r in rotations_in_progress:
                # calculate min_standing_time deps
                if rot.departure_time > r.earliest_departure_next_rot:
                    idle_vehicles.append(r.vehicle_id)
                    rotations_in_progress.pop(0)
                else:
                    break

            # find idle vehicle for rotation if exists
            # else generate new vehicle id
            vt_ct = f"{rot.vehicle_type}_{rot.charging_type}"
            vehicle_id = next((id for id in idle_vehicles if vt_ct in id), None)
            if vehicle_id is None:
                vehicle_type_counts[vt_ct] += 1
                vehicle_id = f"{vt_ct}_{vehicle_type_counts[vt_ct]}"
            else:
                idle_vehicles.remove(vehicle_id)

            rot.vehicle_id = vehicle_id

            # keep list of rotations in progress sorted
            i = 0
            for i, r in enumerate(rotations_in_progress):
                # go through rotations in order, stop at same or higher departure
                if r.earliest_departure_next_rot >= rot.earliest_departure_next_rot:
                    break
            # insert at calculated index
            rotations_in_progress.insert(i, rot)

        self.vehicle_type_counts = vehicle_type_counts

    def calculate_consumption(self):
        """ Computes consumption for all trips of all rotations.

        :return: Total consumption for entire schedule [kWh]
        :rtype: float
        """
        self.consumption = 0
        for rot in self.rotations.values():
            self.consumption += rot.calculate_consumption()

        return self.consumption

    def delta_soc_all_trips(self):
        for rot in self.rotations.values():
            rot.delta_soc_all_trips()

    def get_departure_of_first_trip(self):
        """ Finds earliest departure time among all rotations.

        :return: Date and time of earliest departure of schedule.
        :rtype: datetime.datetime
        """
        sorted_rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)
        return sorted_rotations[0].departure_time

    def get_arrival_of_last_trip(self):
        """Finds latest arrival time among all rotations.

        :return: Date and time of latest arrival of schedule.
        :rtype: datetime.datetime
        """
        sorted_rotations = sorted(self.rotations.values(), key=lambda rot: rot.arrival_time)
        return sorted_rotations[-1].arrival_time

    def readjust_charging_type(self, args):
        """
        Loads rotations with negative soc from spice_ev results and adjusts charging type from
        depb to oppb and vice versa.
        # todo: it would maybe make sense to change each rotations charging type at a time and not
        # todo: all together, as they may influence one another
        :param args: Command line arguments and/or arguments from config file.
        :type args: argparse.Namespace
        """
        negative_rotations = self.get_negative_rotations(args)

        print(f"Rotations {self.get_negative_rotations(args)} have negative SoC.")
        print("Adjust charging types for rotations with negative soc.")

        for rot in negative_rotations:
            if self.rotations[rot].charging_type == "depb":
                self.set_charging_type("oppb", args, [rot])
                # todo: actually this case should not happen, but it still does happen.. why?
            else:
                self.set_charging_type("depb", args, [rot])

    def get_negative_rotations(self, args):
        """
        Get rotations with negative soc from spice_ev outputs

        :param args: Command line arguments and/or arguments from config file.
        :type args: argparse.Namespace
        :return: list of negative rotation_id's
        :rtype: list

        :raises TypeError: If args.save_results is not set.
        """

        # load any json output file of sice_ev
        gcID = list(self.scenario["constants"]["grid_connectors"].keys())[0]
        try:
            ext = path.splitext(args.save_results)
        except TypeError:
            raise TypeError("In order to get negative totations from spice_ev results, please "
                            "specify 'save_results' in your input arguments.")
        filename = f"{ext[0]}_{gcID}{ext[-1]}"
        f = open(filename)
        results = json.load(f)

        # get dict of vehicles with negative soc's
        negative_vehicles = results["vehicles with negative soc"]
        # get matching rotations
        negative_rotations = []
        for v_id, time in negative_vehicles.items():
            time = datetime.datetime.fromisoformat(time)
            rides = {k: v for k, v in self.rotations.items() if v.vehicle_id == v_id}
            negative_rotations.append([k for k, v in rides.items() if time >= v.departure_time and
                                       time <= v.arrival_time][0])
        return negative_rotations

    def generate_scenario_json(self, args):
        """ Generate scenario.json for spiceEV

        :param args: Command line arguments and/or arguments from config file.
        :type args: argparse.Namespace
        """
        # load stations file
        if args.electrified_stations is None:
            args.electrified_stations = "examples/electrified_stations.json"
        ext = args.electrified_stations.split('.')[-1]
        if ext != "json":
            print("File extension mismatch: electrified_stations file should be .json")
        with open(args.electrified_stations) as json_file:
            stations_dict = json.load(json_file)

        interval = datetime.timedelta(minutes=args.interval)

        vehicles = {}
        batteries = {}
        charging_stations = {}
        grid_connectors = {}
        events = {
            "grid_operator_signals": [],
            "external_load": {},
            "energy_feed_in": {},
            "vehicle_events": []
        }

        # add vehicle events
        for vehicle_id in {rot.vehicle_id for rot in self.rotations.values()}:
            v_name = vehicle_id
            vt = vehicle_id.split("_")[0]
            ct = vehicle_id.split("_")[1]
            # define start conditions
            vehicles[v_name] = {
                "connected_charging_station": None,
                "estimated_time_of_departure": None,
                "desired_soc": None,
                "soc": args.desired_soc,
                "vehicle_type": vt + "_" + ct
            }
            # filter all rides for that bus
            v_id = {k: v for k, v in self.rotations.items() if v.vehicle_id == v_name}
            # sort events for their departure time, so that the matching departure time of an
            # arrival event can be read out of the next element in vid_list
            v_id = {k: v for k, v in sorted(v_id.items(), key=lambda x: x[1].departure_time)}
            key_list = list(v_id.keys())
            for i, v in enumerate(key_list):
                departure_event_in_input = True
                # create events for all trips of one rotation
                for j, trip in enumerate(v_id[v].trips):
                    cs_name = "{}_{}".format(v_name, trip.arrival_name)
                    gc_name = trip.arrival_name
                    arrival = trip.arrival_time
                    try:
                        departure = v_id[v].trips[j + 1].departure_time
                        next_arrival = v_id[v].trips[j + 1].arrival_time
                    except IndexError:
                        # get departure of the first trip of the next rotation
                        try:
                            departure = v_id[key_list[i + 1]].departure_time
                            next_arrival = v_id[key_list[i + 1]].trips[0].arrival_time
                        except IndexError:
                            departure_event_in_input = False
                            departure = arrival + datetime.timedelta(hours=8)
                            # no more rotations

                    # connect cs and add gc if station is electrified
                    connected_charging_station = None
                    if gc_name in stations_dict["depot_stations"]:
                        station_type = "deps"
                        desired_soc = args.desired_soc
                    elif gc_name in stations_dict["opp_stations"] and ct == "oppb":
                        station_type = "opps"
                        desired_soc = 1
                    else:
                        # either station has no charger or current bus cannot charge at this station
                        station_type = None
                        desired_soc = 0

                    # 1. if departure - arrival shorter than min_charging_time,
                    # do not connect charging station
                    # 2. if current station has no charger or a depot bus arrives at opp charger,
                    # do not connect charging station either
                    if (((departure - arrival).seconds / 60 >= args.min_charging_time_opps) and
                            station_type is not None):

                        cs_name_and_type = cs_name + "_" + station_type
                        connected_charging_station = cs_name_and_type

                        if cs_name not in charging_stations or gc_name not in grid_connectors:
                            if station_type == "deps":
                                number_cs = stations_dict["depot_stations"][gc_name]
                                cs_power = args.cs_power_deps_oppb if ct == 'oppb' \
                                    else args.cs_power_deps_depb
                                gc_power = args.gc_power_deps
                            elif station_type == "opps":
                                number_cs = stations_dict["opp_stations"][gc_name]
                                cs_power = args.cs_power_opps
                                gc_power = args.gc_power_opps

                            # gc power is not set in config
                            if gc_power is None:
                                if number_cs != "None":
                                    gc_power = number_cs * cs_power
                                else:
                                    # add a really large number
                                    gc_power = 100 * cs_power

                            # add one charging station for each bus at bus station
                            charging_stations[cs_name_and_type] = {
                                "max_power": cs_power,
                                "min_power": 0.1 * cs_power,
                                "parent": gc_name
                            }
                            # add one grid connector for each bus station
                            number_cs = None if number_cs == 'None' else number_cs
                            grid_connectors[gc_name] = {
                                "max_power": gc_power,
                                "cost": {"type": "fixed", "value": 0.3},
                                "number_cs": number_cs
                            }

                    # create arrival events
                    events["vehicle_events"].append({
                        "signal_time": (arrival + datetime.timedelta(minutes=-args.signal_time_dif)
                                        ).isoformat(),
                        "start_time": arrival.isoformat(),
                        "vehicle_id": v_name,
                        "event_type": "arrival",
                        "update": {
                            "connected_charging_station": connected_charging_station,
                            "estimated_time_of_departure": departure.isoformat(),
                            "soc_delta": trip.delta_soc,
                            "desired_soc": desired_soc
                        }
                    })
                    # create departure events
                    if departure_event_in_input:
                        events["vehicle_events"].append({
                            "signal_time": (departure + datetime.timedelta(
                                            minutes=-args.signal_time_dif)).isoformat(),
                            "start_time": departure.isoformat(),
                            "vehicle_id": v_name,
                            "event_type": "departure",
                            "update": {
                                "estimated_time_of_arrival": next_arrival.isoformat()
                            }
                        })

        # ######## END OF VEHICLE EVENTS ########## #

        # define start and stop times
        start = self.get_departure_of_first_trip() - timedelta(minutes=args.signal_time_dif)
        stop = self.get_arrival_of_last_trip() + interval
        if args.days is not None:
            stop = min(stop, start + datetime.timedelta(days=args.days))
        daily = datetime.timedelta(days=1)
        # price events
        for key in grid_connectors.keys():
            if not args.include_price_csv:
                now = start - daily
                while now < stop + 2 * daily:
                    now += daily
                    for v_id, v in vehicles.items():
                        if now >= stop:
                            # after end of scenario: keep generating trips, but don't include in
                            # scenario
                            continue

                    # generate prices for the day
                    if now < stop:
                        morning = now + datetime.timedelta(hours=6)
                        evening_by_month = now + datetime.timedelta(
                            hours=22 - abs(6 - now.month))
                        events['grid_operator_signals'] += [{
                            # day (6-evening): 15ct
                            "signal_time": max(start, now - daily).isoformat(),
                            "grid_connector_id": key,
                            "start_time": morning.isoformat(),
                            "cost": {
                                "type": "fixed",
                                "value": 0.15 + random.gauss(0, 0.05)
                            }
                        }, {
                            # night (depending on month - 6): 5ct
                            "signal_time": max(start, now - daily).isoformat(),
                            "grid_connector_id": key,
                            "start_time": evening_by_month.isoformat(),
                            "cost": {
                                "type": "fixed",
                                "value": 0.05 + random.gauss(0, 0.03)
                            }
                        }]

        # add timeseries from csv
        # save path and options for CSV timeseries

        if args.include_ext_load_csv:
            for filename, gc_name in args.include_ext_load_csv:
                basename = path.splitext(path.basename(filename))[0]
                options = {
                    "csv_file": filename,
                    "start_time": start.isoformat(),
                    "step_duration_s": 900,  # 15 minutes
                    "grid_connector_id": gc_name,
                    "column": "energy"
                }
                if args.include_ext_csv_option:
                    for key, value in args.include_ext_csv_option:
                        if key == "step_duration_s":
                            value = int(value)
                        options[key] = value
                events['external_load'][basename] = options
                # check if CSV file exists
                ext_csv_path = path.join(args.output_directory, filename)
                if not path.exists(ext_csv_path):
                    print("Warning: external csv file '{}' does not exist yet".format(ext_csv_path))

        if args.include_feed_in_csv:
            for filename, gc_name in args.include_feed_in_csv:
                basename = path.splitext(path.basename(filename))[0]
                options = {
                    "csv_file": filename,
                    "start_time": start.isoformat(),
                    "step_duration_s": 3600,  # 60 minutes
                    "grid_connector_id": gc_name,
                    "column": "energy"
                }
                if args.include_feed_in_csv_option:
                    for key, value in args.include_feed_in_csv_option:
                        if key == "step_duration_s":
                            value = int(value)
                        options[key] = value
                events['energy_feed_in'][basename] = options
                feed_in_path = path.join(args.output_directory, filename)
                if not path.exists(feed_in_path):
                    print("Warning: feed-in csv file '{}' does not exist yet".format(feed_in_path))

        if args.include_price_csv:
            for filename, gc_name in args.include_price_csv:
                options = {
                    "csv_file": filename,
                    "start_time": start.isoformat(),
                    "step_duration_s": 3600,  # 60 minutes
                    "grid_connector_id": gc_name,
                    "column": "price [ct/kWh]"
                }
                for key, value in args.include_price_csv_option:
                    if key == "step_duration_s":
                        value = int(value)
                    options[key] = value
                events['energy_price_from_csv'] = options
                price_csv_path = path.join(args.output_directory, filename)
                if not path.exists(price_csv_path):
                    print("Warning: price csv file '{}' does not exist yet".format(price_csv_path))

        if args.battery:
            for idx, (capacity, c_rate, gc) in enumerate(args.battery):
                if capacity > 0:
                    max_power = c_rate * capacity
                else:
                    # unlimited battery: set power directly
                    max_power = c_rate
                batteries["BAT{}".format(idx + 1)] = {
                    "parent": gc,
                    "capacity": capacity,
                    "charging_curve": [[0, max_power], [1, max_power]]
                }
        # create final dict
        j = {
            "scenario": {
                "start_time": start.isoformat(),
                "interval": interval.days * 24 * 60 + interval.seconds // 60,
                "n_intervals": (stop - start) // interval
            },
            "constants": {
                "vehicle_types": self.vehicle_types,
                "vehicles": vehicles,
                "grid_connectors": grid_connectors,
                "charging_stations": charging_stations,
                "batteries": batteries
            },
            "events": events
        }
        # Write JSON
        self.scenario = j
        with open(args.input, 'w+') as f:
            json.dump(j, f, indent=2)
