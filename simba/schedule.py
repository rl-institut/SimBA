import csv # noqa
import datetime
import json
import logging
from pathlib import Path
import random
import warnings
from simba import util, optimizer_util
from simba.rotation import Rotation
from spice_ev.battery import Battery
from spice_ev.scenario import Scenario
import spice_ev.util as spice_ev_util


class Schedule:

    def __init__(self, vehicle_types, stations, **kwargs):
        """ Constructs Schedule object from CSV file containing all trips of schedule

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
        """ Constructs Schedule object from CSV file containing all trips of schedule.

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
        """ Check rotation expectations.

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
        # self.assign_vehicles()
        self.new_assign_vehicles(args)

        scenario = self.generate_scenario(args)

        logging.info("Running SpiceEV...")
        with warnings.catch_warnings():
            if logging.root.level > logging.DEBUG:
                warnings.simplefilter('ignore', UserWarning)
            scenario.run('distributed', vars(args).copy())
        assert scenario.step_i == scenario.n_intervals, \
            'SpiceEV simulation aborted, see above for details'
        return scenario

    def set_charging_type(self, ct, rotation_ids=None):
        """ Change charging type of either all or specified rotations.

        Adjust minimum standing time at depot after completion of rotation.

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
            case the one with the longest standing time since last rotation is used.
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
            rot.vehicle_id = f"{vt_ct}_{'0' * missing}{old_num}"

        self.vehicle_type_counts = vehicle_type_counts

    def new_assign_vehicles(self, args):
        """ Assign vehicle IDs to rotations.

        A FIFO approach is used.
        For every rotation it is checked whether vehicles with matching type are idle,
        in which case the one with longest standing time since last rotation is used.
        If no vehicle is available, a new vehicle ID is generated.
        :param args: arguments
        :type args: Namespace
        """
        rotations_in_progress = []
        standing_vehicles = []
        vehicle_data = dict()
        print("Assigning")
        initial_soc = args.desired_soc_deps

        # count number of vehicles per type
        # used for unique vehicle id e.g. vehicletype_chargingtype_id
        vehicle_type_counts = {f'{vehicle_type}_{charging_type}': 0
                               for vehicle_type, charging_types in self.vehicle_types.items()
                               for charging_type in charging_types.keys()}

        charge_levels = {station.get(param) for name, station in self.stations.items()
                         for param in ["cs_power_deps_oppb", "cs_power_deps_depb"]
                         if station["type"] == "deps"}
        # Remove None as charge level
        charge_levels = charge_levels.difference([None])

        charge_levels.add(args.cs_power_deps_depb)
        charge_levels.add(args.cs_power_deps_oppb)

        # Calculates numeric charge curve
        charge_curves = self.get_charge_curves(charge_levels, time_step=1)

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)
        for k, rot in enumerate(rotations):

            # find vehicles that have completed rotation and stood for a minimum standing time
            # mark those vehicle as idle
            while rotations_in_progress:
                # calculate min_standing_time deps
                r = rotations_in_progress.pop(0)
                if rot.departure_time > r.arrival_time:
                    standing_vehicles.append((r.vehicle_id, r.arrival_name))
                else:
                    rotations_in_progress.insert(0, r)
                    break

            # find idle vehicle for rotation if exists
            # else generate new vehicle id
            vt = rot.vehicle_type
            ct = rot.charging_type
            vt_ct = f"{vt}_{ct}"
            for item in standing_vehicles:
                vehicle_id, deps = item
                if vt_ct in vehicle_id and deps == rot.departure_name:
                    start_soc = vehicle_data[vehicle_id]["soc"]

                    station_power = self.stations[deps].get(f"cs_power_deps_{ct}",
                                                            vars(args).get(f"cs_power_deps_{ct}"))
                    buffer_time = datetime.timedelta(minutes=util.get_buffer_time(rot.trips[0]))
                    duration_in_m = (rot.departure_time - vehicle_data[vehicle_id][
                        "arrival_time"] - buffer_time) / datetime.timedelta(minutes=1)
                    # When arriving the vehicle has this soc
                    # This soc is reached when leaving for the current rotation
                    charge_delta_soc = get_charge_delta_soc(charge_curves, vt, ct, station_power,
                                                            duration_in_m,
                                                            start_soc, max_soc=initial_soc)
                    # Assign vehicles
                    charged_soc = max(start_soc, 0) + charge_delta_soc
                    delta_soc = rot.calculate_consumption() / self.vehicle_types[vt][ct]["capacity"]
                    end_soc = charged_soc - delta_soc
                    if end_soc > 0 or charged_soc >= initial_soc:
                        print(delta_soc, charged_soc)
                        rot.vehicle_id = vehicle_id
                        standing_vehicles.remove(item)
                        vehicle_data[vehicle_id] = {"soc": end_soc,
                                                    "arrival_time": rot.arrival_time}
                        break
            else:
                # no vehicle available for dispatch, generate new one
                vehicle_type_counts[vt_ct] += 1
                vehicle_id = f"{vt_ct}_{vehicle_type_counts[vt_ct]}"
                rot.vehicle_id = vehicle_id
                start_soc = initial_soc
                delta_soc = rot.calculate_consumption() / self.vehicle_types[vt][ct]["capacity"]
                end_soc = start_soc - delta_soc
                vehicle_data[vehicle_id] = {"soc": end_soc,
                                            "arrival_time": rot.arrival_time}

            print(f"assigned {vehicle_id} to {rot.id} with {start_soc=} and {end_soc=}")

            # keep list of rotations in progress sorted
            def sort_by_soc_at_next_rot_departure(rot):
                vt = rot.vehicle_type
                ct = rot.charging_type
                start_soc = vehicle_data[rot.vehicle_id]["soc"]
                station_power = self.stations[rot.arrival_name].get(f"cs_power_deps_{ct}",
                                                                    vars(args).get(
                                                                        f"cs_power_deps_{ct}"))
                buffer_time = datetime.timedelta(minutes=util.get_buffer_time(rot.trips[0]))
                duration_in_m = (next_rot.departure_time - vehicle_data[vehicle_id][
                    "arrival_time"] - buffer_time) / datetime.timedelta(minutes=1)
                charge_delta_soc = get_charge_delta_soc(charge_curves, vt, ct, station_power,
                                                        duration_in_m,
                                                        start_soc, max_soc=initial_soc)
                return max(0, start_soc) + charge_delta_soc

            rotations_in_progress.append(rot)
            try:
                next_rot = rotations[k + 1]
            except IndexError:
                pass
            else:
                rotations_in_progress = sorted(rotations_in_progress,
                                               key=sort_by_soc_at_next_rot_departure)

        # update vehicle ID for nice natural sorting
        for rot in rotations:
            # get old vehicle id (vehicle type + charging type, sequential number)
            vt_ct, old_num = rot.vehicle_id.rsplit('_', 1)
            # how many vehicles of this type?
            vt_cnt = vehicle_type_counts[vt_ct]
            # easy log10 to get max needed number of digits
            digits = len(str(vt_cnt))
            # assign new zero-padded ID (all of same vehicle type have same length)
            rot.vehicle_id = f"{vt_ct}_{int(old_num):0{digits}}"
        self.vehicle_type_counts = vehicle_type_counts

    def get_charge_curves(self, charge_levels, time_step):
        charge_curves = dict()
        for vehicle_name, vehicle_type in self.vehicle_types.items():
            charge_curves[vehicle_name] = dict()
            for ct_name, v_info in vehicle_type.items():
                charge_curves[vehicle_name][ct_name] = dict()
                default_eff = Battery(1, [[0, 1], [1, 1]], 0).efficiency
                eff = v_info.get("battery_efficiency", default_eff)
                for charge_level in charge_levels:
                    if ct_name == "depb":
                        final_value = self.min_recharge_deps_depb
                    else:
                        assert ct_name == "oppb"
                        final_value = self.min_recharge_deps_oppb
                    curve = optimizer_util.charging_curve_to_soc_over_time(
                        v_info["charging_curve"],
                        v_info["capacity"],
                        final_value,
                        charge_level,
                        time_step=time_step,
                        efficiency=eff,
                        logger=logging.getLogger())
                    charge_curves[vehicle_name][ct_name][charge_level] = curve

        return charge_curves

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
        """ Finds latest arrival time among all rotations.

        :return: Date and time of latest arrival of schedule. None if rotations are empty.
        :rtype: datetime.datetime
        """
        if not self.rotations:
            return None
        sorted_rotations = sorted(self.rotations.values(), key=lambda rot: rot.arrival_time)
        return sorted_rotations[-1].arrival_time

    def get_common_stations(self, only_opps=True):
        """ For each rotation key, return set of rotations that share a station during any trip.

        :param only_opps: only search for opps stations
        :type only_opps: boolean
        :return: rotations with time info
        :rtype: dict
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
        """ Get rotations with negative SoC from SpiceEV outputs.

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
        """ Edits rotations according to args.rotation_filter_variable.

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
        """ Generate SpiceEV Scenario.

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
            "fixed_load": {},
            "local_generation": {},
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

        # read in time windows file
        time_windows = None
        if vars(args).get('time_windows'):
            time_windows_path = Path(args.time_windows)
            if time_windows_path.exists():
                with time_windows_path.open('r', encoding='utf-8') as f:
                    time_windows = util.uncomment_json_file(f)
                # convert time window strings to date/times
                for grid_operator, grid_operator_seasons in time_windows.items():
                    for season, info in grid_operator_seasons.items():
                        info["start"] = datetime.date.fromisoformat(info["start"])
                        info["end"] = datetime.date.fromisoformat(info["end"])
                        for level, windows in info.get("windows", {}).items():
                            info["windows"][level] = [
                                (datetime.time.fromisoformat(t[0]),
                                 datetime.time.fromisoformat(t[1]))
                                for t in windows]
                        time_windows[grid_operator][season] = info
                logging.debug("time windows read: " + str(time_windows))
            else:
                warnings.warn(f"Time windows file {args.time_windows} does not exist")

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
                        voltage_level = station.get('voltage_level', args.default_voltage_level)
                        grid_operator = station.get("grid_operator", "default_grid_operator")
                        grid_connectors[gc_name] = {
                            "max_power": gc_power,
                            "cost": {"type": "fixed", "value": 0.3},
                            "number_cs": station["n_charging_stations"],
                            "grid_operator": grid_operator,
                            "voltage_level": voltage_level,
                        }
                        # check for stationary battery
                        battery = station.get("battery")
                        if battery is not None:
                            # add stationary battery at this station/GC
                            battery["parent"] = gc_name
                            batteries[f"{gc_name} storage"] = battery

                        # add feed-in name and power at grid connector if exists
                        local_generation = station.get("energy_feed_in")
                        if local_generation:
                            local_generation = update_csv_file_info(local_generation, gc_name)
                            events["local_generation"][gc_name + " feed-in"] = local_generation
                            # add PV component
                            photovoltaics[gc_name] = {
                                "parent": gc_name,
                                "nominal_power": local_generation.get("nominal_power", 0)
                            }

                        # add external load if exists
                        fixed_load = station.get("external_load")
                        if fixed_load:
                            fixed_load = update_csv_file_info(fixed_load, gc_name)
                            events["fixed_load"][gc_name + " ext. load"] = fixed_load

                        # temporary lowering of grid connector max power during peak load windows
                        if time_windows is not None:
                            # get relevant peak load windows
                            windows = time_windows.get(grid_operator)
                            reduced_power = station.get(
                                'peak_load_window_power',
                                vars(args).get('peak_load_window_power_' + station_type))
                            events["grid_operator_signals"] += generate_time_window_event_list(
                                windows, gc_name, voltage_level, (gc_power, reduced_power),
                                (start_simulation, stop_simulation, interval))

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

        # price events
        if args.include_price_csv:
            # SpiceEV's include_price_csv supports only one grid connector
            # declare price input files for each station in electrified stations instead
            logging.warning("Unsupported SpiceEV option: include_price_csv. Use at your own risk!")
        random.seed(args.seed)
        for gc_name in grid_connectors.keys():
            price_csv = self.stations[gc_name].get("price_csv")
            if price_csv is None:
                # generate pseudo-random price events
                events["grid_operator_signals"] += generate_random_price_list(
                    gc_name, start_simulation, stop_simulation)
            else:
                # read prices from CSV, convert to events
                prices = get_price_list_from_csv(price_csv)
                events["grid_operator_signals"] += generate_event_list_from_prices(
                    prices, gc_name, start_simulation, stop_simulation,
                    price_csv.get('start_time'), price_csv.get('step_duration_s'))

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


def update_csv_file_info(file_info, gc_name):
    """
    add infos to csv information dictionary from electrified station

    - set grid_connector_id
    - update csv_file path
    - set start_time and step_duration_s from CSV information if not given
    :param file_info: csv information from electrified station
    :type file_info: dict
    :param gc_name: station name
    :type gc_name: string
    :return: updated file_info
    :rtype: dict
    """
    file_path = Path(file_info["csv_file"])
    start_time = file_info.get("start_time")
    if not file_path.exists():
        warnings.warn("csv file '{}' does not exist".format(file_info))
    elif start_time is None:
        # start_time is optional (read from first line), but file must exist
        with file_path.open('r', encoding='utf-8') as f:
            # skip header
            f.readline()
            # start timestep in first column, second row
            row = f.readline()
            ts = row.split(',')[0]
            start_time = spice_ev_util.datetime_from_isoformat(ts)
            # force timezone-unaware (like all SimBA datetimes)
            start_time = start_time.replace(tzinfo=None)
            # read next timestep to get interval length
            row = f.readline()
            ts = row.split(',')[0]
            time = spice_ev_util.datetime_from_isoformat(ts)
            time = time.replace(tzinfo=None)
            csv_interval = (time - start_time).total_seconds()
        file_info["start_time"] = start_time.isoformat()
        file_info["step_duration_s"] = round(csv_interval)

    file_info["grid_connector_id"] = gc_name
    file_info["csv_file"] = str(file_path)
    return file_info


def generate_time_window_event_list(time_windows, gc_name, voltage_level, power_levels, time_info):
    """
    generate grid operator signals to lower station power during peak load windows

    :param time_windows: grid operator time window information
    :type time_windows: dict
    :param gc_name: station name
    :type gc_name: string
    :param voltage_level: station voltage level
    :type voltage_level: string
    :param power_levels: (normal power, reduced power)
    :type power_levels: tuple
    :param time_info: (start of simulation, end of simulation, interval)
    :type time_info: triple
    :return: grid operator signals
    :rtype: list
    """
    if time_windows is None:
        warnings.warn(f'{gc_name}: grid operator not in time windows')
        return []
    logging.debug(f"{gc_name} time windows: {time_windows}")

    # check reduced power
    # if higher than normal max power, don't generate events
    if power_levels[1] > power_levels[0]:
        warnings.warn(f"{gc_name}: maximum power is not reduced during peak load windows")
        return []

    events = []
    in_window = False
    day = datetime.timedelta(days=1)
    start_simulation, stop_simulation, interval = time_info
    # can't directly iterate over datetimes:
    # get number of simulation timesteps (round up)
    n_intervals = -((start_simulation - stop_simulation) // interval)
    for t_idx in range(n_intervals):
        time = start_simulation + t_idx * interval
        cur_window = spice_ev_util.datetime_within_time_window(time, time_windows, voltage_level)
        if cur_window ^ in_window:
            # window change: generate event
            events.append({
                # known one day in advance
                "signal_time": (time - day).isoformat(),
                "start_time": time.isoformat(),
                "grid_connector_id": gc_name,
                "max_power": power_levels[cur_window],
            })
            in_window = cur_window
    return events


def generate_random_price_list(gc_name, start_simulation, stop_simulation):
    """ Generate random price events.

    :param gc_name: grid connector ID
    :type gc_name: string
    :param start_simulation: start of simulation
    :type start_simulation: datetime
    :param stop_simulation: end of simulation
    :type stop_simulation: datetime
    :return: newly generated grid operator signals
    :rtype: list
    """
    day = datetime.timedelta(days=1)
    events = []
    now = start_simulation - day
    while now < stop_simulation + 2 * day:
        now += day
        # generate prices for the day
        if now < stop_simulation:
            morning = now + datetime.timedelta(hours=6)
            evening_by_month = now + datetime.timedelta(
                hours=22 - abs(6 - now.month))
            events += [{
                # day (6 to evening): 15ct
                "signal_time": max(start_simulation, now - day).isoformat(),
                "grid_connector_id": gc_name,
                "start_time": morning.isoformat(),
                "cost": {
                    "type": "fixed",
                    "value": 0.15 + random.gauss(0, 0.05)
                }
            }, {
                # night (evening to 6 ): 5ct
                "signal_time": max(start_simulation, now - day).isoformat(),
                "grid_connector_id": gc_name,
                "start_time": evening_by_month.isoformat(),
                "cost": {
                    "type": "fixed",
                    "value": 0.05 + random.gauss(0, 0.03)
                }
            }]
    return events


def get_price_list_from_csv(price_csv_dict):
    """ Read out price CSV.

    :param price_csv_dict: price CSV info
    :type price_csv_dict: dict
    :return: price timestamp (if in first column of CSV) and values
    :rtype: list of tuples
    """
    csv_path = Path(price_csv_dict["csv_file"])
    if not csv_path.exists():
        logging.error(f"price csv file {csv_path} does not exist, skipping price generation")
        return []

    prices = []
    column = price_csv_dict.get("column")
    factor = price_csv_dict.get("factor", 1)
    with csv_path.open('r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        for idx, row in enumerate(reader):
            row_values = list(row.values())
            # read value from given column or last column
            value_str = row[column] if column else row_values[-1]
            value = float(value_str) * factor
            # add entry: first cell (may contain datetime), value
            prices.append((row_values[0], value))
    return prices


def generate_event_list_from_prices(
        prices, gc_name, start_simulation, stop_simulation,
        start_events=None, price_interval_s=None):
    """ Generate grid operator signals from price list.

    :param prices: price timestamp and values
    :type prices: list of tuples
    :param gc_name: grid connector ID
    :type gc_name: string
    :param start_simulation: start of simulation
    :type start_simulation: datetime.datetime
    :param stop_simulation: end of simulation
    :type stop_simulation: datetime.datetime
    :param start_events: timestamp of first event in list (optional). Read from list if not given
    :type start_events: string (ISO-format)
    :param price_interval_s: interval between list entries in seconds (only with start_events)
    :type price_interval_s: number
    :return: grid operator signals
    :rtype: list
    """
    if start_events is None:
        assert price_interval_s is None, (
            f"{gc_name} price CSV: must have start_time if price_interval is given")
        logging.debug(f"{gc_name} price CSV: using timestamps from file")
    else:
        start_events = spice_ev_util.datetime_from_isoformat(start_events)
        price_interval = datetime.timedelta(seconds=price_interval_s)

    day = datetime.timedelta(days=1)
    event = None
    events = []
    stop = None

    for idx, price in enumerate(prices):
        if price_interval_s is None:
            # default: use CSV timestamps in first column
            event_time = spice_ev_util.datetime_from_isoformat(price[0])
            if start_events is None:
                # keep track of first event time
                start_events = event_time
        else:
            # compute event timestamp from given arguments
            event_time = idx * price_interval + start_events
        if stop is None or stop < event_time:
            # keep track of last event time
            stop = event_time
        if event_time > stop_simulation:
            # don't generate prices after stop of scenario
            break
        if event_time >= start_simulation and len(events) == 0 and event is not None:
            # first event within scenario time, but there were events before: add prior event
            events.append(event)
        # price events known one day in advance
        signal_time = event_time - day
        # read value from given column or last column
        event = {
            "start_time": event_time.isoformat(),
            "signal_time": signal_time.isoformat(),
            "grid_connector_id": gc_name,
            "cost": {"type": "fixed", "value": price[1]},
        }
        if event_time >= start_simulation:
            # event within scenario time: add to list
            events.append(event)
    # post check: events before scenario start?
    if len(events) == 0 and event is not None:
        # prior events: add last
        events.append(event)

    if events:
        # check if CSV timestamps cover simulation time
        if start_events > start_simulation or stop < stop_simulation:
            logging.info(f"{gc_name} price csv does not cover simulation time")
    return events


def get_charge_delta_soc(charge_curves: dict, vt: str, ct: str, max_power: float,
                         time_delta: datetime.timedelta, start_soc: float, max_soc: float) -> float:
    charge_curve = charge_curves[vt][ct][max_power]
    d_soc = optimizer_util.get_delta_soc(charge_curve, start_soc, time_delta=time_delta,
                                         max_soc=max_soc)
    return d_soc
