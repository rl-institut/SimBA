import csv
import json
import random
import datetime
from os import path
import warnings

from ebus_toolbox.rotation import Rotation
from src.scenario import Scenario


class Schedule:

    def __init__(self, vehicle_types, stations_file, **kwargs):
        """Constructs Schedule object from CSV file containing all trips of schedule

        :param vehicle_types: Collection of vehicle types and their properties.
        :type vehicle_types: dict
        :param stations_file: json of electrified stations
        :type stations_file: string
        :param kwargs: Command line arguments
        :type kwargs: dict
        """
        self.vehicle_types = vehicle_types

        # load stations file
        if stations_file is None:
            stations_file = "examples/electrified_stations.json"
        ext = stations_file.split('.')[-1]
        if ext != "json":
            print("File extension mismatch: electrified_stations file should be .json")
        with open(stations_file) as f:
            self.stations = json.load(f)

        self.rotations = {}
        self.original_rotations = None
        self.consumption = 0

        # mandatory config parameters
        mandatory_options = [
            "min_recharge_deps_oppb",
            "min_recharge_deps_depb",
            "gc_power_opps",
            "gc_power_deps",
            "cs_power_opps",
            "cs_power_deps_depb",
            "cs_power_deps_oppb",
        ]
        for opt in mandatory_options:
            assert opt in kwargs, f"Missing config paramter: {opt}"
            setattr(self, opt, kwargs.get(opt))

    @classmethod
    def from_csv(cls, path_to_csv, vehicle_types, stations, **kwargs):
        """Constructs Schedule object from CSV file containing all trips of schedule.

        :param path_to_csv: Path to csv file containing trip data
        :type path_to_csv: str
        :param vehicle_types: Collection of vehicle types and their properties.
        :type vehicle_types: dict
        :param stations: json of electrified stations
        :type stations: string
        :param kwargs: Command line arguments
        :type kwargs: dict
        :return: Returns a new instance of Schedule with all trips from csv loaded.
        :rtype: Schedule
        """
        schedule = cls(vehicle_types, stations, **kwargs)

        with open(path_to_csv, 'r') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                rotation_id = trip['rotation_id']
                if rotation_id not in schedule.rotations.keys():
                    schedule.rotations.update({
                        rotation_id: Rotation(id=rotation_id,
                                              vehicle_type=trip['vehicle_type'],
                                              schedule=schedule)})
                schedule.rotations[rotation_id].add_trip(trip)

        return schedule

    def run(self, args):
        # (re)calculate the change in SoC for every trip
        # charging types may have changed which may impact battery capacity
        # while mileage is assumed to stay constant
        self.delta_soc_all_trips()

        # each rotation is assigned a vehicle ID
        self.assign_vehicles()

        scenario = self.generate_scenario(args)

        print("Running Spice EV...")
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            scenario.run('distributed', vars(args).copy())

        return scenario

    def set_charging_type(self, ct, rotation_ids=None):
        """ Change charging type of either all or specified rotations. Adjust minimum standing time
        at depot after completion of rotation.
        :param ct: Choose this charging type wheneever possible. Either 'depb' or 'oppb'.
        :type ct: str
        :param rotation_ids: IDs of rotations for which to set charging type. If None set charging
                            charging type for all rotations.
        :type rotation_ids: list
        """
        assert ct in ["oppb", "depb"], f"Invalid charging type: {ct}"

        for id, rot in self.rotations.items():
            if rotation_ids is None or id in rotation_ids:
                rot.set_charging_type(ct)

    def assign_vehicles(self):
        """ Assign vehicle IDs to rotations. A FIFO approach is used.
            For every rotation it is checked whether vehicles with matching type are idle, in which
            case the one with longest standing time since last rotation is used.
            If no vehicle is available a new vehicle ID is generated.
        """
        rotations_in_progress = []
        idle_vehicles = []
        # TODO: create vehicle type counts dict with ct and vt values of all types
        vehicle_type_counts = {vehicle_type: 0 for vehicle_type in self.vehicle_types.keys()}

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)

        for rot in rotations:
            # find vehicles that have completed rotation and stood for a minimum staning time
            # mark those vehicle as idle
            while rotations_in_progress:
                # calculate min_standing_time deps
                r = rotations_in_progress.pop(0)
                if rot.departure_time > r.earliest_departure_next_rot:
                    idle_vehicles.append((r.vehicle_id, r.arrival_name))
                else:
                    rotations_in_progress.insert(0, r)
                    break

            # find idle vehicle for rotation if exists
            # else generate new vehicle id
            vt_ct = f"{rot.vehicle_type}_{rot.charging_type}"
            for item in idle_vehicles:
                id, deps = item
                if vt_ct in id and deps == rot.departure_name:
                    rot.vehicle_id = id
                    idle_vehicles.remove(item)
                    break
            else:
                # no vehicle available for dispatch, generate new one
                vehicle_type_counts[vt_ct] += 1
                rot.vehicle_id = f"{vt_ct}_{vehicle_type_counts[vt_ct]}"

            # keep list of rotations in progress sorted
            i = 0
            for i, r in enumerate(rotations_in_progress):
                # go through rotations in order, stop at same or higher departure
                if r.earliest_departure_next_rot >= rot.earliest_departure_next_rot:
                    break
            else:
                # highest departure, insert at
                i = i + 1
            # insert at calculated index
            rotations_in_progress.insert(i, rot)

        self.vehicle_type_counts = vehicle_type_counts

    def calculate_consumption(self):
        """ Computes consumption for all trips of all rotations.
            Depends on vehicle type only, not on charging type.

        :return: Total consumption for entire schedule [kWh]
        :rtype: float
        """
        self.consumption = 0
        for rot in self.rotations.values():
            self.consumption += rot.calculate_consumption()

        return self.consumption

    def delta_soc_all_trips(self):
        """ Computes delta SOC for all trips of all rotations.
            Depends on vehicle type and on charging type, since
            busses of the same vehicle type may have different
            battery sizes for different charging types."""

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

    def get_common_stations(self, only_opps=True):
        """
        for each rotation key, return set of rotations
            that share a station during any trip (with time info)
        :param only_opps: only search for opps stations
        :type only_opps: boolean
        :return: dictionary of rotations
        """

        # rot -> stations with timings
        rotations = {}
        for rot_key, rotation in self.rotations.items():
            rotations[rot_key] = {}
            for t in rotation.trips:
                arr_station = self.stations.get(t.arrival_name)
                if not only_opps or arr_station is not None and arr_station["type"] == "opps":
                    # dependent station: add to rotation with arrival time
                    if t.arrival_name in rotations[rot_key]:
                        rotations[rot_key][t.arrival_name].append([t.arrival_time, None])
                    else:
                        rotations[rot_key][t.arrival_name] = [[t.arrival_time, None]]
                if t.departure_name in rotations[rot_key]:
                    rotations[rot_key][t.departure_name][-1][1] = t.departure_time

        # check stations for overlaps
        rot_set = {}
        for rot_key, stations in rotations.items():
            rot_set[rot_key] = {}
            for r, alt_stations in rotations.items():
                if rot_key == r:
                    continue
                for station_name, times in stations.items():
                    if station_name in alt_stations:
                        # both at same station. Check for same time
                        for t in times:
                            for alt_t in alt_stations[station_name]:
                                if t[1] is None or alt_t[1] is None:
                                    # no departure: ignore
                                    continue
                                if max(t[0], alt_t[0]) < min(t[1], alt_t[1]):
                                    # overlap
                                    rot_set[rot_key][r] = min(t[0], alt_t[0])
                                    break
                            if r in rot_set[rot_key]:
                                break
        return rot_set

    def get_negative_rotations(self, scenario):
        """
        Get rotations with negative soc from spice_ev outputs

        :param scenario: Simulation scenario containing simulation results
                         including the SoC of all vehicles over time
        :type scenario: spice_ev.Scenario
        :return: list of negative rotation_id's
        :rtype: list
        """

        # get dict of vehicles with negative soc's
        try:
            negative_vehicles = scenario.negative_soc_tracker
        except AttributeError:
            return []
        # get matching rotations
        negative_rotations = set()
        for v_id, times in negative_vehicles.items():
            for t_str in times:
                time = datetime.datetime.fromisoformat(t_str)
                rides = {k: v for k, v in self.rotations.items() if v.vehicle_id == v_id}
                negative_rotations.add([k for k, v in rides.items()
                                        if time >= v.departure_time and time <= v.arrival_time][0])
        return list(negative_rotations)

    def generate_scenario(self, args):
        """ Generate scenario.json for spiceEV

        :param args: Command line arguments and/or arguments from config file.
        :type args: argparse.Namespace
        :return: A spiceEV Scenario instance that can be run and also collects all
                 simulation outputs.
        :rtype:  spice_ev.src.Scenario
        """

        random.seed(args.seed)

        # load stations file
        if args.electrified_stations is None:
            args.electrified_stations = "examples/electrified_stations.json"
        ext = args.electrified_stations.split('.')[-1]
        if ext != "json":
            print("File extension mismatch: electrified_stations file should be .json")
        # with open(args.electrified_stations) as json_file:
            # stations_dict = json.load(json_file)

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

            # filter all rides for that bus
            vehicle_rotations = {k: v for k, v in self.rotations.items() if v.vehicle_id == v_name}
            # sort events for their departure time, so that the matching departure time of an
            # arrival event can be read out of the next element in vid_list
            vehicle_rotations = {k: v for k, v in sorted(
                                 vehicle_rotations.items(), key=lambda x: x[1].departure_time)}
            rotation_ids = list(vehicle_rotations.keys())

            # define start conditions
            first_rotation = list(vehicle_rotations.values())[0]
            departure = first_rotation.departure_time
            # trips list is sorted by time
            arrival = first_rotation.trips[0].arrival_time
            vehicles[v_name] = {
                "connected_charging_station": f'{v_name}_{first_rotation.departure_name}_deps',
                "estimated_time_of_departure": departure.isoformat(),
                "desired_soc": None,
                "soc": args.desired_soc,
                "vehicle_type": vt + "_" + ct
            }
            # initial departure event from starting depot
            events["vehicle_events"].append({
                "signal_time": (departure + datetime.timedelta(
                                minutes=-args.signal_time_dif)).isoformat(),
                "start_time": departure.isoformat(),
                "vehicle_id": v_name,
                "event_type": "departure",
                "update": {
                    "estimated_time_of_arrival": arrival.isoformat()
                }
            })

            for i, v in enumerate(rotation_ids):
                departure_event_in_input = True
                # create events for all trips of one rotation
                for j, trip in enumerate(vehicle_rotations[v].trips):
                    cs_name = "{}_{}".format(v_name, trip.arrival_name)
                    gc_name = trip.arrival_name
                    arrival = trip.arrival_time
                    try:
                        departure = vehicle_rotations[v].trips[j + 1].departure_time
                        next_arrival = vehicle_rotations[v].trips[j + 1].arrival_time
                    except IndexError:
                        # get departure of the first trip of the next rotation
                        try:
                            departure = vehicle_rotations[rotation_ids[i + 1]].departure_time
                            next_arrival = \
                                vehicle_rotations[rotation_ids[i + 1]].trips[0].arrival_time
                        except IndexError:
                            departure_event_in_input = False
                            departure = arrival + datetime.timedelta(hours=8)
                            # no more rotations

                    # calculate total minutes spend at station
                    standing_time = (departure - arrival).seconds / 60
                    # get buffer time from user configuration
                    # buffer time resembles amount of time deducted off of the planned standing
                    # time. It may resemble things like delays and/or docking procedures
                    # use buffer time from electrified stations JSON or in case none is
                    # provided use global default from config file
                    buffer_time = self.stations.get(trip.arrival_name, {})\
                        .get('buffer_time', args.default_buffer_time_opps)

                    # distinct buffer times depending on time of day can be provided
                    # in that case buffer time is of type dict instead of int
                    if isinstance(buffer_time, dict):
                        # sort dict to make sure 'else' key is last key
                        buffer_time = {key: buffer_time[key] for key in sorted(buffer_time)}
                        current_hour = arrival.hour
                        for time_range, buffer in buffer_time.items():
                            if time_range == 'else':
                                buffer_time = buffer
                                break
                            else:
                                start_hour, end_hour = [int(t) for t in time_range.split('-')]
                                if end_hour < start_hour:
                                    if current_hour >= start_hour or current_hour < end_hour:
                                        buffer_time = buffer
                                        break
                                else:
                                    if start_hour <= current_hour < end_hour:
                                        buffer_time = buffer
                                        break
                        else:
                            # buffer time not specified for hour of current stop
                            buffer_time = args.default_buffer_time_opps

                    # connect cs and add gc if station is electrified
                    connected_charging_station = None
                    desired_soc = 0
                    try:
                        # assume electrified station
                        station = self.stations[gc_name]
                        station_type = station["type"]
                        if station_type == 'opps' and vehicle_rotations[v].charging_type == 'depb':
                            # a depot bus cannot charge at an opp station
                            station_type = None
                        else:
                            # at opp station, always try to charge as much as you can
                            desired_soc = 1 if station_type == "opps" else args.desired_soc
                    except KeyError:
                        # non-electrified station
                        station_type = None

                    # 1. if standing time - buffer time shorter than min_charging_time,
                    # do not connect charging station
                    # 2. if current station has no charger or a depot bus arrives at opp charger,
                    # do not connect charging station either
                    if (station_type is not None and
                            (standing_time - buffer_time >= args.min_charging_time_opps)):

                        cs_name_and_type = cs_name + "_" + station_type
                        connected_charging_station = cs_name_and_type

                        if cs_name not in charging_stations or gc_name not in grid_connectors:
                            number_cs = station["n_charging_stations"]
                            if station_type == "deps":
                                cs_power = args.cs_power_deps_oppb if ct == 'oppb' \
                                    else args.cs_power_deps_depb
                                gc_power = args.gc_power_deps
                            elif station_type == "opps":
                                # delay the arrival time to account for docking procedure
                                arrival = arrival + datetime.timedelta(minutes=buffer_time)
                                cs_power = args.cs_power_opps
                                gc_power = args.gc_power_opps

                            # gc power is not set in config
                            if gc_power is None:
                                if number_cs is not None:
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
        start = \
            self.get_departure_of_first_trip() - datetime.timedelta(minutes=args.signal_time_dif)
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
            for idx, bat in enumerate(args.battery, 1):
                capacity, c_rate, gc = bat[:3]
                if gc not in grid_connectors:
                    # no vehicle charges here
                    continue
                if capacity > 0:
                    max_power = c_rate * capacity
                else:
                    # unlimited battery: set power directly
                    max_power = c_rate
                batteries[f"BAT{idx}"] = {
                    "parent": gc,
                    "capacity": capacity,
                    "charging_curve": [[0, max_power], [1, max_power]],
                }
                if len(bat) == 4:
                    # discharge curve can be passed as optional 4th param
                    batteries[f"BAT{idx}"].update({
                        "discharge_curve": bat[3]
                    })

        # TODO: restructure vehicle types for SpiceEV

        # create final dict
        self.scenario = {
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

        return Scenario(self.scenario, args.output_directory)

    def generate_rotations_overview(self, scenario, args):
        rotation_infos = []

        negative_rotations = self.get_negative_rotations(scenario)

        interval = datetime.timedelta(minutes=args.interval)
        sim_start_time = \
            self.get_departure_of_first_trip() - datetime.timedelta(minutes=args.signal_time_dif)

        for id, rotation in self.rotations.items():
            # get SOC timeseries for this rotation
            vehicle_id = rotation.vehicle_id

            # get soc timeseries for current rotation
            vehicle_soc = scenario.vehicle_socs[vehicle_id]
            start_idx = (rotation.departure_time - sim_start_time) // interval
            end_idx = start_idx + ((rotation.arrival_time-rotation.departure_time) // interval)
            if end_idx > len(vehicle_soc):
                # SpiceEV stopped before rotation was fully simulated
                continue
            rotation_soc_ts = vehicle_soc[start_idx:end_idx]

            rotation_info = {
                                "rotation_id": id,
                                "start_time": rotation.departure_time.isoformat(),
                                "end_time": rotation.arrival_time.isoformat(),
                                "vehicle_type": rotation.vehicle_type,
                                "vehicle_id": rotation.vehicle_id,
                                "depot_name": rotation.departure_name,
                                "lines": ':'.join(rotation.lines),
                                "total_consumption_[kWh]": rotation.consumption,
                                "distance": rotation.distance,
                                "charging_type": rotation.charging_type,
                                "SOC_at_arrival": rotation_soc_ts[-1],
                                "Minumum_SOC": min(rotation_soc_ts),
                                "Negative_SOC": 1 if id in negative_rotations else 0
                             }
            rotation_infos.append(rotation_info)

        with open(path.join(args.output_directory, "rotations.csv"), "w+") as f:
            csv_writer = csv.DictWriter(f, list(rotation_infos[0].keys()))
            csv_writer.writeheader()
            csv_writer.writerows(rotation_infos)
