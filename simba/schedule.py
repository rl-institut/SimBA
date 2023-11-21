import csv
import datetime
import json
import logging
from pathlib import Path
import random
import warnings

from simba import util
from simba.rotation import Rotation
from spice_ev.scenario import Scenario


class Schedule:

    def __init__(self, vehicle_types, stations, **kwargs):
        """Constructs Schedule object from CSV file containing all trips of schedule

        :param vehicle_types: Collection of vehicle types and their properties.
        :type vehicle_types: dict
        :param stations: electrified stations
        :type stations: dict

        :raises Exception: In case not all mandatory options are provided

        :param kwargs: Command line arguments
        :type kwargs: dict
        """

        self.stations = stations
        self.rotations = {}
        self.consumption = 0
        self.vehicle_types = vehicle_types
        self.original_rotations = None
        self.station_data = None

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
        missing = [opt for opt in mandatory_options if kwargs.get(opt) is None]
        if missing:
            raise Exception("The following arguments are required: {}".format(", ".join(missing)))
        else:
            for opt in mandatory_options:
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

        station_data = dict()
        station_path = kwargs.get("station_data_path")

        # find the temperature and elevation of the stations by reading the .csv file.
        # this data is stored in the schedule and passed to the trips, which use the information
        # for consumption calculation. Missing station data is handled with default values.
        if station_path is not None:
            try:
                with open(station_path, "r", encoding='utf-8') as f:
                    delim = util.get_csv_delim(station_path)
                    reader = csv.DictReader(f, delimiter=delim)
                    for row in reader:
                        station_data.update({str(row['Endhaltestelle']):
                                            {"elevation": float(row['elevation'])}})
            except FileNotFoundError or KeyError:
                warnings.warn("Warning: external csv file '{}' not found or not named properly "
                              "(Needed column names are 'Endhaltestelle' and 'elevation')".
                              format(station_path),
                              stacklevel=100)
            except ValueError:
                warnings.warn("Warning: external csv file '{}' does not contain numeric "
                              "values in the column 'elevation'. Station data is discarded.".
                              format(station_path),
                              stacklevel=100)
            schedule.station_data = station_data

        with open(path_to_csv, 'r', encoding='utf-8') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                rotation_id = trip['rotation_id']
                # trip gets reference to station data and calculates height diff during trip
                # initialization. Could also get the height difference from here on
                trip["station_data"] = station_data
                if rotation_id not in schedule.rotations.keys():
                    schedule.rotations.update({
                        rotation_id: Rotation(id=rotation_id,
                                              vehicle_type=trip['vehicle_type'],
                                              schedule=schedule)})
                schedule.rotations[rotation_id].add_trip(trip)

        # set charging type for all rotations without explicitly specified charging type
        # charging type may have been set above if a trip of a rotation has a specified
        # charging type
        for rot in schedule.rotations.values():
            if rot.charging_type is None:
                rot.set_charging_type(ct=kwargs.get('preferred_charging_type', 'oppb'))

        if kwargs.get("check_rotation_consistency"):
            # check rotation expectations
            inconsistent_rotations = cls.check_consistency(schedule)
            if inconsistent_rotations:
                # write errors to file
                filepath = kwargs["output_directory"] / "inconsistent_rotations.csv"
                with open(filepath, "w", encoding='utf-8') as f:
                    for rot_id, e in inconsistent_rotations.items():
                        f.write(f"Rotation {rot_id}: {e}\n")
                        logging.error(f"Rotation {rot_id}: {e}")
                        if kwargs.get("skip_inconsistent_rotations"):
                            # remove this rotation from schedule
                            del schedule.rotations[rot_id]
        elif kwargs.get("skip_inconsistent_rotations"):
            warnings.warn("Option skip_inconsistent_rotations ignored, "
                          "as check_rotation_consistency is not set to 'true'")

        return schedule

    @classmethod
    def check_consistency(cls, schedule):
        """
        Check rotation expectations, such as
        - the rotation starts and ends at the same station
        - every trip within a rotation starts where the previous trip ended
        - trips are chronologically sorted
        - trips have positive breaks in between
        - trips have positive times between departure and arrival

        :param schedule: the schedule to check
        :type schedule: dict
        :return: inconsistent rotations. Dict of rotation ID -> error message
        :rtype: dict
        """
        inconsistent_rotations = {}
        for rot_id, rotation in schedule.rotations.items():
            # iterate over trips, looking for initial and final stations
            prev_trip = None
            try:
                for trip in rotation.trips:
                    # positive duration for a trip
                    assert trip.arrival_time >= trip.departure_time, "Trip time is negative"
                    if not prev_trip:
                        prev_trip = trip
                        continue
                    # positive break time in between trips
                    assert trip.departure_time >= prev_trip.arrival_time, "Break time is negative"
                    assert trip.departure_name == prev_trip.arrival_name, "Trips are not sequential"

                    prev_trip = trip

                rot_departure_trip = list(rotation.trips)[0]
                rot_arrival_trip = list(rotation.trips)[-1]
                dep_name = rot_departure_trip.departure_name
                arr_name = rot_arrival_trip.arrival_name
                assert dep_name == arr_name, "Start and end of rotation differ"

                error = "Rotation data differs from trips data"
                assert dep_name == rotation.departure_name, error
                assert arr_name == rotation.arrival_name, error
                assert rot_arrival_trip.arrival_time == rotation.arrival_time, error
                assert rot_departure_trip.departure_time == rotation.departure_time, error
            except AssertionError as e:
                # some assumption is violated
                # save error text
                inconsistent_rotations[rot_id] = str(e)
        return inconsistent_rotations

    def run(self, args):
        # each rotation is assigned a vehicle ID
        self.assign_vehicles()

        scenario = self.generate_scenario(args)

        logging.info("Running SpiceEV...")
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            scenario.run('distributed', vars(args).copy())
        assert scenario.step_i == scenario.n_intervals, \
            'SpiceEV simulation aborted, see above for details'
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
        # count number of vehicles per type
        # used for unique vehicle id e.g. vehicletype_chargingtype_id
        vehicle_type_counts = {f'{vehicle_type}_{charging_type}': 0
                               for vehicle_type, charging_types in self.vehicle_types.items()
                               for charging_type in charging_types.keys()}

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)

        for rot in rotations:
            # find vehicles that have completed rotation and stood for a minimum standing time
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

        # update vehicle ID for nice natural sorting
        for rot in rotations:
            # get old vehicle id (vehicle type + charging type, sequential number)
            vt_ct, old_num = rot.vehicle_id.rsplit('_', 1)
            # how many vehicles of this type?
            vt_cnt = vehicle_type_counts[vt_ct]
            # easy log10 to get max needed number of digits
            digits = len(str(vt_cnt))
            # how many digits have to be added for this vehicle ID?
            missing = digits - len(str(old_num))
            # assign new zero-padded ID (all of same vehicle type have same length)
            rot.vehicle_id = f"{vt_ct}_{'0'*missing}{old_num}"

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

    def get_departure_of_first_trip(self):
        """ Finds earliest departure time among all rotations.

        :return: Date and time of earliest departure of schedule. None if rotations are empty.
        :rtype: datetime.datetime
        """
        if not self.rotations:
            return None
        sorted_rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)
        return sorted_rotations[0].departure_time

    def get_arrival_of_last_trip(self):
        """Finds latest arrival time among all rotations.

        :return: Date and time of latest arrival of schedule. None if rotations are empty.
        :rtype: datetime.datetime
        """
        if not self.rotations:
            return None
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
        Get rotations with negative soc from SpiceEV outputs

        :param scenario: Simulation scenario containing simulation results
                         including the SoC of all vehicles over time
        :type scenario: spice_ev.Scenario
        :return: list of negative rotation_id's
        :rtype: list
        """

        # get dict of vehicles with negative SOCs
        try:
            negative_vehicles = scenario.negative_soc_tracker
        except AttributeError:
            return []
        # get matching rotations
        negative_rotations = set()
        for v_id, times in negative_vehicles.items():
            for t_str in times:
                time = datetime.datetime.fromisoformat(t_str)
                # filter rotations by vehicle ID
                rides = {k: v for k, v in self.rotations.items() if v.vehicle_id == v_id}
                # sort rotations by departure time
                rides = sorted(rides.items(), key=lambda r: r[1].departure_time)
                # find rotation before negative arrival
                last_rot_id = None
                for k, v in rides:
                    if v.departure_time >= time:
                        # departure after negative arrival found: stop search
                        break
                    last_rot_id = k
                # end of iteration (next departure found / end of list): add last rotation ID
                if last_rot_id is not None:
                    negative_rotations.add(last_rot_id)

        return list(negative_rotations)

    def rotation_filter(self, args, rf_list=[]):
        """Edits rotations according to args.rotation_filter_variable.

        :param args: used arguments are rotation_filter, path to rotation ids,
                     and rotation_filter_variable that sets mode (options: include, exclude)
        :type args: argparse.Namespace
        :param rf_list: rotation filter list with strings of rotation ids (default is None)
        :type rf_list: list
        """
        if args.rotation_filter_variable is None:
            # filtering disabled
            return
        # cast rotations in filter to string
        rf_list = [str(i) for i in rf_list]

        if args.rotation_filter is None and not rf_list:
            warnings.warn("Rotation filter variable is enabled but file and list are not used.")
            return

        if args.rotation_filter:
            # read out rotations from file (one rotation ID per line)
            try:
                with open(args.rotation_filter, encoding='utf-8') as f:
                    for line in f:
                        rf_list.append(line.strip())
            except FileNotFoundError:
                warnings.warn(f"Path to rotation filter {args.rotation_filter} is invalid.")
                # no file, no change
                return
        # filter out rotations in self.rotations
        if args.rotation_filter_variable == "exclude":
            self.rotations = {k: v for k, v in self.rotations.items() if k not in rf_list}
        elif args.rotation_filter_variable == "include":
            self.rotations = {k: v for k, v in self.rotations.items() if k in rf_list}

    def generate_scenario(self, args):
        """ Generate scenario.json for SpiceEV

        :param args: Command line arguments and/or arguments from config file.
        :type args: argparse.Namespace
        :return: A SpiceEV Scenario instance that can be run and also collects all
                 simulation outputs.
        :rtype:  SpiceEV.Scenario
        """

        interval = datetime.timedelta(minutes=args.interval)

        vehicles = {}
        batteries = {}
        photovoltaics = {}
        charging_stations = {}
        grid_connectors = {}
        events = {
            "grid_operator_signals": [],
            "external_load": {},
            "energy_feed_in": {},
            "vehicle_events": []
        }

        # define start and stop times
        if self.rotations:
            start_simulation = self.get_departure_of_first_trip()
            start_simulation -= datetime.timedelta(minutes=args.signal_time_dif)
            stop_simulation = self.get_arrival_of_last_trip() + interval
            if args.days is not None:
                stop_simulation = min(
                    stop_simulation, start_simulation + datetime.timedelta(days=args.days))
        else:
            # no rotations to base times off: use any datetime
            start_simulation = datetime.datetime.fromtimestamp(0)
            stop_simulation = datetime.datetime.fromtimestamp(0)

        # add vehicle events
        for vehicle_id in sorted({rot.vehicle_id for rot in self.rotations.values()}):
            # filter all rides for that bus
            vehicle_rotations = {k: v for k, v in self.rotations.items()
                                 if v.vehicle_id == vehicle_id}
            # get sorted list of all trips for current vehicle
            # sort rotations by time leveraging the fact that the
            # list of trips per rotation is always sorted
            vehicle_rotations = {k: v for k, v in sorted(
                vehicle_rotations.items(), key=lambda x: x[1].departure_time)}
            vehicle_trips = [t for rot in vehicle_rotations.values() for t in rot.trips]

            for i, trip in enumerate(vehicle_trips):
                # don't generate events that start after simulation has stopped
                if trip.departure_time >= stop_simulation:
                    break

                cs_name = f"{vehicle_id}_{trip.arrival_name}"
                gc_name = trip.arrival_name

                # departure time of next trip for standing time calculation
                try:
                    next_departure_time = vehicle_trips[i+1].departure_time
                except IndexError:
                    # last trip
                    next_departure_time = trip.arrival_time + datetime.timedelta(hours=8)

                # get current station if electrified
                connected_charging_station = None
                desired_soc = 0
                try:
                    # assume electrified station
                    station = self.stations[gc_name]
                    station_type = station["type"]
                    if station_type == 'opps' and trip.rotation.charging_type == 'depb':
                        # a depot bus cannot charge at an opp station
                        station_type = None
                    else:
                        # get desired soc by station type
                        desired_soc = vars(args).get("desired_soc_" + station_type)
                except KeyError:
                    # non-electrified station
                    station_type = None

                # get buffer time from user configuration
                # buffer_time is an abstraction of delays like docking procedures and
                # is added to the planned arrival time
                # ignore buffer time for end of last trip to make sure vehicles arrive
                # before simulation ends
                if i < len(vehicle_trips) - 1:
                    buffer_time = util.get_buffer_time(trip=trip,
                                                       default=args.default_buffer_time_opps)
                else:
                    buffer_time = 0
                # arrival event must occur no later than next departure and
                # one step before simulation terminates for arrival event to be taken into account
                arrival_time = min(trip.arrival_time + datetime.timedelta(minutes=buffer_time),
                                   next_departure_time)

                # total minutes spend at station
                standing_time = (next_departure_time - arrival_time).total_seconds() / 60

                # connect to charging station
                # generate gc and cs if they dont exist
                if station_type is not None and standing_time >= args.min_charging_time:
                    # vehicle connects to charging station
                    connected_charging_station = f"{cs_name}_{station_type}"
                    # create charging station and grid connector if necessary
                    if cs_name not in charging_stations:
                        # get CS and GC power from stations file,
                        # default back to input arguments from config
                        # trip type only relevant to depot stations
                        trip_type = trip.rotation.charging_type
                        trip_type = '_' + trip_type if station_type == "deps" else ""
                        cs_power_type = f"cs_power_{station_type}{trip_type}"
                        cs_power = station.get(cs_power_type, vars(args)[cs_power_type])
                        gc_power = station.get("gc_power", vars(args)[f"gc_power_{station_type}"])
                        # add one charging station for each bus at bus station
                        charging_stations[connected_charging_station] = {
                            "type": station_type,
                            "max_power": cs_power,
                            "min_power": 0.1 * cs_power,
                            "parent": gc_name
                        }
                    if gc_name not in grid_connectors:
                        # add one grid connector for each bus station
                        grid_connectors[gc_name] = {
                            "max_power": gc_power,
                            "cost": {"type": "fixed", "value": 0.3},
                            "number_cs": station["n_charging_stations"],
                            "grid_operator": station.get("grid_operator", "default_grid_operator"),
                            "voltage_level":
                                station.get("voltage_level", args.default_voltage_level)
                        }
                        # check for stationary battery
                        battery = station.get("battery")
                        if battery is not None:
                            # add stationary battery at this station/GC
                            battery["parent"] = gc_name
                            batteries[f"{gc_name} storage"] = battery

                        # add feed-in name and power at grid connector if exists
                        feed_in = station.get("energy_feed_in")
                        if feed_in:
                            feed_in_path = Path(feed_in["csv_file"])
                            if not feed_in_path.exists():
                                warnings.warn("feed-in csv file '{}' does not exist".format(
                                    feed_in_path), category=UserWarning)
                            feed_in["grid_connector_id"] = gc_name
                            feed_in["csv_file"] = str(feed_in_path)
                            events["energy_feed_in"][gc_name + " feed-in"] = feed_in
                            # add PV component
                            photovoltaics[gc_name] = {
                                "parent": gc_name,
                                "nominal_power": feed_in.get("nominal_power", 0)
                            }

                        # add external load if exists
                        ext_load = station.get("external_load")
                        if ext_load:
                            ext_load_path = Path(ext_load["csv_file"])
                            if not ext_load_path.exists():
                                warnings.warn("external load csv file '{}' does not exist".format(
                                    ext_load_path), category=UserWarning)
                            ext_load["grid_connector_id"] = gc_name
                            ext_load["csv_file"] = str(ext_load_path)
                            events["external_load"][gc_name + " ext. load"] = ext_load

                # initial condition of vehicle
                if i == 0:
                    vehicles[vehicle_id] = {
                        "connected_charging_station": None,
                        "estimated_time_of_departure": trip.departure_time.isoformat(),
                        "desired_soc": None,
                        "soc": args.desired_soc_deps,
                        "vehicle_type":
                            f"{trip.rotation.vehicle_type}_{trip.rotation.charging_type}"
                    }

                # create departure event
                events["vehicle_events"].append({
                    "signal_time": (trip.departure_time
                                    - datetime.timedelta(minutes=args.signal_time_dif)
                                    ).isoformat(),
                    "start_time": trip.departure_time.isoformat(),
                    "vehicle_id": vehicle_id,
                    "event_type": "departure",
                    "update": {
                        "estimated_time_of_arrival": arrival_time.isoformat()
                    }
                })
                assert trip.delta_soc is not None, (
                    "Trip delta_soc is None. Did you forget to calculate consumption?")

                # create arrival event
                events["vehicle_events"].append({
                    "signal_time": (arrival_time
                                    - datetime.timedelta(minutes=args.signal_time_dif)
                                    ).isoformat(),
                    "start_time": arrival_time.isoformat(),
                    "vehicle_id": vehicle_id,
                    "event_type": "arrival",
                    "update": {
                        "connected_charging_station": connected_charging_station,
                        "estimated_time_of_departure": next_departure_time.isoformat(),
                        "soc_delta": trip.delta_soc,
                        "desired_soc": desired_soc
                    }
                })

        # ######## END OF VEHICLE EVENTS ########## #

        daily = datetime.timedelta(days=1)
        # price events
        if not args.include_price_csv:
            random.seed(args.seed)
            for key in grid_connectors.keys():
                now = start_simulation - daily
                while now < stop_simulation + 2 * daily:
                    now += daily
                    for v_id, v in vehicles.items():
                        if now >= stop_simulation:
                            # after end of scenario: keep generating trips, but don't include in
                            # scenario
                            continue

                    # generate prices for the day
                    if now < stop_simulation:
                        morning = now + datetime.timedelta(hours=6)
                        evening_by_month = now + datetime.timedelta(
                            hours=22 - abs(6 - now.month))
                        events['grid_operator_signals'] += [{
                            # day (6-evening): 15ct
                            "signal_time": max(start_simulation, now - daily).isoformat(),
                            "grid_connector_id": key,
                            "start_time": morning.isoformat(),
                            "cost": {
                                "type": "fixed",
                                "value": 0.15 + random.gauss(0, 0.05)
                            }
                        }, {
                            # night (depending on month - 6): 5ct
                            "signal_time": max(start_simulation, now - daily).isoformat(),
                            "grid_connector_id": key,
                            "start_time": evening_by_month.isoformat(),
                            "cost": {
                                "type": "fixed",
                                "value": 0.05 + random.gauss(0, 0.03)
                            }
                        }]

        # add timeseries from csv
        # save path and options for CSV timeseries
        if args.include_price_csv:
            for filename, gc_name in args.include_price_csv:
                options = {
                    "csv_file": filename,
                    "start_time": start_simulation.isoformat(),
                    "step_duration_s": 3600,  # 60 minutes
                    "grid_connector_id": gc_name,
                    "column": "price [ct/kWh]"
                }
                for key, value in args.include_price_csv_option:
                    if key == "step_duration_s":
                        value = int(value)
                    options[key] = value
                events['energy_price_from_csv'] = options
                price_csv_path = args.output_directory / filename
                if not price_csv_path.exists():
                    logging.warning(f"Price csv file '{price_csv_path}' does not exist yet")

        # reformat vehicle types for SpiceEV
        vehicle_types_spiceev = {
            f'{vehicle_type}_{charging_type}': body
            for vehicle_type, subtypes in self.vehicle_types.items()
            for charging_type, body in subtypes.items()
        }

        # create final dict
        self.scenario = {
            "scenario": {
                "start_time": start_simulation.isoformat(),
                "interval": interval.days * 24 * 60 + interval.seconds // 60,
                "n_intervals": (stop_simulation - start_simulation) // interval
            },
            "components": {
                "vehicle_types": vehicle_types_spiceev,
                "vehicles": vehicles,
                "grid_connectors": grid_connectors,
                "charging_stations": charging_stations,
                "batteries": batteries,
                "photovoltaics": photovoltaics,
            },
            "events": events
        }
        if vars(args).get("create_scenario_file"):
            p = args.output_directory_input / args.create_scenario_file
            with p.open('w', encoding='utf-8') as f:
                json.dump(self.scenario, f, indent=2)
        return Scenario(self.scenario, Path())
