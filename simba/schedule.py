import csv
import datetime
import json
import logging
from pathlib import Path
import random
import warnings
from typing import Dict, Type, Iterable

import simba.rotation
from simba.data_container import DataContainer
from simba import util, optimizer_util
from simba.rotation import Rotation

from spice_ev.components import VehicleType
from spice_ev.scenario import Scenario
import spice_ev.util as spice_ev_util


class SocDispatcher:
    """Dispatches the right initial SoC for every vehicle id at scenario generation.

    Used for specific vehicle initialization for example when coupling tools."""

    def __init__(self,
                 default_soc_deps: float,
                 default_soc_opps: float,
                 # vehicle_socs stores the departure soc of a rotation as a dict of the previous
                 # trip, since this is how the SpiceEV scenario is generated.
                 # The first trip of a vehicle has no previous trip and therefore is None
                 vehicle_socs: Dict[str, Type[Dict["simba.trip.Trip", float]]] = None):
        self.default_soc_deps = default_soc_deps
        self.default_soc_opps = default_soc_opps
        self.vehicle_socs = {}
        if vehicle_socs is not None:
            self.vehicle_socs = vehicle_socs

    def get_soc(self, vehicle_id: str, trip: "simba.trip.Trip",  station_type: str = None):
        try:
            v_socs = self.vehicle_socs[vehicle_id]
            return v_socs[trip]
        except KeyError:
            return vars(self).get("default_soc_" + station_type)


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
        self.soc_dispatcher: SocDispatcher = None
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
    def from_datacontainer(cls, data_container: DataContainer, args):

        schedule = cls(data_container.vehicle_types_data, data_container.stations_data, **args)
        schedule.station_data = data_container.station_geo_data

        for trip in data_container.trip_data:
            rotation_id = trip['rotation_id']
            # trip gets reference to station data and calculates height diff during trip
            # initialization. Could also get the height difference from here on
            # get average hour of trip if level of loading or temperature has to be read from
            # auxiliary tabular data

            # get average hour of trip and parse to string, since tabular data has strings
            # as keys
            hour = (trip["departure_time"] +
                    (trip["arrival_time"] - trip["departure_time"]) / 2).hour
            # Get height difference from station_data

            trip["height_difference"] = schedule.get_height_difference(
                trip["departure_name", trip["arrival_name"]])

            # Get level of loading from trips.csv or from file
            try:
                # Clip level of loading to [0,1]
                lol = max(0, min(float(trip["level_of_loading"]), 1))
            # In case of empty temperature column or no column at all
            except (KeyError, ValueError):
                lol = data_container.level_of_loading_data[hour]

            trip["level_of_loading"] = lol

            # Get temperature from trips.csv or from file
            try:
                # Cast temperature to float
                temperature = float(trip["temperature"])
            # In case of empty temperature column or no column at all
            except (KeyError, ValueError):
                temperature = data_container.temperature_data[hour]
            trip["temperature"] = temperature
            if rotation_id not in schedule.rotations.keys():
                schedule.rotations.update({
                    rotation_id: Rotation(id=rotation_id,
                                          vehicle_type=trip['vehicle_type'],
                                          schedule=schedule)})
            schedule.rotations[rotation_id].add_trip(trip)

        # set charging type for all rotations without explicitly specified charging type.
        # charging type may have been set above if a trip of a rotation has a specified
        # charging type
        for rot in schedule.rotations.values():
            if rot.charging_type is None:
                rot.set_charging_type(ct=args.get('preferred_charging_type', 'oppb'))

        if args.get("check_rotation_consistency"):
            # check rotation expectations
            inconsistent_rotations = cls.check_consistency(schedule)
            if inconsistent_rotations:
                # write errors to file
                filepath = args["output_directory"] / "inconsistent_rotations.csv"
                with open(filepath, "w", encoding='utf-8') as f:
                    for rot_id, e in inconsistent_rotations.items():
                        f.write(f"Rotation {rot_id}: {e}\n")
                        logging.error(f"Rotation {rot_id}: {e}")
                        if args.get("skip_inconsistent_rotations"):
                            # remove this rotation from schedule
                            del schedule.rotations[rot_id]
        elif args.get("skip_inconsistent_rotations"):
            logging.warning("Option skip_inconsistent_rotations ignored, "
                          "as check_rotation_consistency is not set to 'true'")

        return schedule



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
        :raises NotImplementedError: if stations is neither a str,Path nor dictionary.
        :return: Returns a new instance of Schedule with all trips from csv loaded.
        :rtype: Schedule
        """

        # Check station type
        if isinstance(stations, (str, Path)):
            with open(Path(stations), "r") as f:
                stations_dict = util.uncomment_json_file(f)
        elif isinstance(stations, dict):
            stations_dict = stations
        else:
            raise NotImplementedError

        schedule = cls(vehicle_types, stations_dict, **kwargs)

        station_data = dict()
        station_path = kwargs.get("station_data_path")

        level_of_loading_path = kwargs.get("level_of_loading_over_day_path", None)

        if level_of_loading_path is not None:
            index = "hour"
            column = "level_of_loading"
            level_of_loading_dict = cls.get_dict_from_csv(column, level_of_loading_path, index)

        temperature_path = kwargs.get("outside_temperature_over_day_path", None)
        if temperature_path is not None:
            index = "hour"
            column = "temperature"
            temperature_data_dict = cls.get_dict_from_csv(column, temperature_path, index)

        # find the temperature and elevation of the stations by reading the .csv file.
        # this data is stored in the schedule and passed to the trips, which use the information
        # for consumption calculation. Missing station data is handled with default values.
        if station_path is not None:
            try:
                with open(station_path, "r", encoding='utf-8') as f:
                    delim = util.get_csv_delim(station_path)
                    reader = csv.DictReader(f, delimiter=delim)
                    for row in reader:
                        station_data.update({str(row['Endhaltestelle']): {
                            "elevation": float(row['elevation']), "lat": float(row.get('lat', 0)),
                            "long": float(row.get('long', 0))}
                                             })
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
                # get average hour of trip if level of loading or temperature has to be read from
                # auxiliary tabular data
                arr_time = datetime.datetime.fromisoformat(trip["arrival_time"])
                dep_time = datetime.datetime.fromisoformat(trip["departure_time"])

                # get average hour of trip and parse to string, since tabular data has strings
                # as keys
                hour = (dep_time + (arr_time - dep_time) / 2).hour
                # Get height difference from station_data
                try:
                    height_diff = station_data[trip["arrival_name"]]["elevation"] \
                                  - station_data[trip["departure_name"]]["elevation"]
                except (KeyError, TypeError):
                    height_diff = 0
                trip["height_difference"] = height_diff

                # Get level of loading from trips.csv or from file
                try:
                    # Clip level of loading to [0,1]
                    lol = max(0, min(float(trip["level_of_loading"]), 1))
                # In case of empty temperature column or no column at all
                except (KeyError, ValueError):
                    lol = level_of_loading_dict[hour]
                trip["level_of_loading"] = lol

                # Get temperature from trips.csv or from file
                try:
                    # Cast temperature to float
                    temperature = float(trip["temperature"])
                # In case of empty temperature column or no column at all
                except (KeyError, ValueError):
                    temperature = temperature_data_dict[hour]
                trip["temperature"] = temperature
                if rotation_id not in schedule.rotations.keys():
                    schedule.rotations.update({
                        rotation_id: Rotation(id=rotation_id,
                                              vehicle_type=trip['vehicle_type'],
                                              schedule=schedule)})
                schedule.rotations[rotation_id].add_trip(trip)

        # set charging type for all rotations without explicitly specified charging type.
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
        :type schedule: Schedule
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
        """Runs a schedule without assigning vehicles.

        For external usage the core run functionality is accessible through this function. It
        allows for defining a custom-made assign_vehicles method for the schedule.
        :param args: used arguments are rotation_filter, path to rotation ids,
            and rotation_filter_variable that sets mode (options: include, exclude)
        :type args: argparse.Namespace
        :return: scenario
        :rtype spice_ev.Scenario
        """
        # Make sure all rotations have an assigned vehicle
        assert all([rot.vehicle_id is not None for rot in self.rotations.values()])
        scenario = self.generate_scenario(args)

        logging.info("Running SpiceEV...")
        if logging.root.level > logging.DEBUG:
            # don't log SpiceEV warnings for log levels above debug
            # logging.root.level is lowest of console and file (if present)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                scenario.run('distributed', vars(args).copy())
        else:
            # debug: log SpiceEV warnings as well
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

    def assign_vehicles(self, args):
        """ Assign vehicles using the strategy given in the arguments
        :param args: Arguments with attribute assign_strategy
        :type args: Namespace
        :raises NotImplementedError: if args.assign_strategy has a no allowed value
        """
        assign_strategy = vars(args).get("assign_strategy") or "adaptive"

        if assign_strategy == "adaptive":
            self.assign_vehicles_w_adaptive_soc(args)
        elif assign_strategy == "fixed_recharge":
            self.assign_vehicles_w_min_recharge_soc()
        else:
            logging.error('Allowed values for assign_strategy are "adaptive" and "fixed_recharge"')
            raise NotImplementedError

    def assign_vehicles_w_min_recharge_soc(self):
        """ Assign vehicle IDs to rotations.

        A FIFO approach is used.
        A vehicle is added to the idle stack if the min_recharge_soc is roughly reached.
        For every rotation it is checked whether vehicles with matching type are idle,
        in which case the one with longest standing time since last rotation is used.
        If no vehicle is available, a new vehicle ID is generated.
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

    def assign_vehicles_for_django(self, eflips_output: Iterable[dict]):
        """Assign vehicles based on eflips outputs

        eflips couples vehicles and returns for every rotation the departure soc and vehicle id.
        This is included into simba by assigning new vehicles with the respective values. I.e. in
        simba every rotation gets a new vehicle.
        :param eflips_output: output from eflips meant for simba. Iterable contains
            rotation_id, vehicle_id and start_soc for each rotation
        :type eflips_output: iterable of dataclass "simba_input"
        :raises KeyError: If not every rotation has a vehicle assigned to it
        """
        eflips_rot_dict = {d["rot"]: {"v_id": d["v_id"], "soc": d["soc"]} for d in eflips_output}
        unique_vids = {d["v_id"] for d in eflips_output}
        vehicle_socs = {v_id: dict() for v_id in unique_vids}
        eflips_vid_dict = {v_id: sorted([d["rot"] for d in eflips_output
                                         if d["v_id"] == v_id],
                                        key=lambda r_id: self.rotations[r_id].departure_time)
                           for v_id in unique_vids}

        # Calculate vehicle counts
        # count number of vehicles per type
        # used for unique vehicle id e.g. vehicletype_chargingtype_id
        vehicle_type_counts = {f'{vehicle_type}_{charging_type}': 0
                               for vehicle_type, charging_types in self.vehicle_types.items()
                               for charging_type in charging_types.keys()}
        for vid in unique_vids:
            v_ls = vid.split("_")
            vehicle_type_counts[f'{v_ls[0]}_{v_ls[1]}'] += 1

        self.vehicle_type_counts = vehicle_type_counts

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)
        for rot in rotations:
            try:
                v_id = eflips_rot_dict[rot.id]["v_id"]
            except KeyError as exc:
                raise KeyError(f"SoC-data does not include the rotation with the id: {rot.id}. "
                               "Externally generated vehicles assignments need to include all "
                               "rotations") from exc
            rot.vehicle_id = v_id
            index = eflips_vid_dict[v_id].index(rot.id)
            # if this rotation is not the first rotation of the vehicle, find the previous trip
            if index != 0:
                prev_rot_id = eflips_vid_dict[v_id][index - 1]
                trip = self.rotations[prev_rot_id].trips[-1]
            else:
                # if the rotation has no previous trip, trip is set as None
                trip = None
            vehicle_socs[v_id][trip] = eflips_rot_dict[rot.id]["soc"]
        self.soc_dispatcher.vehicle_socs = vehicle_socs

    def init_soc_dispatcher(self, args):
        self.soc_dispatcher = SocDispatcher(default_soc_deps=args.desired_soc_deps,
                                            default_soc_opps=args.desired_soc_opps)

    def assign_only_new_vehicles(self):
        """ Assign new vehicle IDs to rotations
        """
        # count number of vehicles per type
        # used for unique vehicle id e.g. vehicletype_chargingtype_id
        vehicle_type_counts = {f'{vehicle_type}_{charging_type}': 0
                               for vehicle_type, charging_types in self.vehicle_types.items()
                               for charging_type in charging_types.keys()}
        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)
        for rot in rotations:
            vt_ct = f"{rot.vehicle_type}_{rot.charging_type}"
            # no vehicle available for dispatch, generate new one
            vehicle_type_counts[vt_ct] += 1
            rot.vehicle_id = f"{vt_ct}_{vehicle_type_counts[vt_ct]}"
        self.vehicle_type_counts = vehicle_type_counts

    def assign_vehicles_w_adaptive_soc(self, args):
        """ Assign vehicle IDs to rotations.

        A greedy approach is used to assign a vehicle to every rotation.
        If an existing vehicle of the same type and location has enough soc to service the rotation,
        not accounting for possible opportunity charging of the next rotation it is used.
        If multiple such vehicles exist, the one with the lowest soc is used.
        If no vehicle is available, a new vehicle ID is generated.
        :param args: arguments
        :type args: Namespace
        """
        rotations_in_progress = []
        all_standing_vehicles = []
        vehicle_data = dict()
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
        # Add default values from arguments
        charge_levels.add(self.cs_power_deps_depb)
        charge_levels.add(self.cs_power_deps_oppb)

        # Calculates numeric charge curves for each charge_level.
        # The charge levels might clip the charge curve of the vehicle.
        charge_curves = self.get_charge_curves(charge_levels, initial_soc, time_step_min=1)

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)
        for k, rot in enumerate(rotations):
            # find vehicles that have completed their rotation
            while rotations_in_progress:
                r = rotations_in_progress.pop(0)
                if rot.departure_time > r.arrival_time:
                    all_standing_vehicles.append((r.vehicle_id, r.arrival_name))
                else:
                    rotations_in_progress.insert(0, r)
                    break

            # find standing vehicle for rotation if exists
            # else generate new vehicle id
            vt = rot.vehicle_type
            ct = rot.charging_type
            vt_ct = f"{vt}_{ct}"

            station_is_electrified = rot.departure_name in self.stations
            if not station_is_electrified:
                logging.warning(f"Rotation {rot.id} ends at a non electrified station.")

            # filter vehicles of the same vehicle type and same location
            standing_vehicles = list(filter(lambda x: vt_ct in x[0] * (x[1] == rot.departure_name),
                                            all_standing_vehicles))

            # join standing vehicles with their expected soc and sort by soc
            socs = list(map(lambda v_id_deps: soc_at_departure_time(v_id_deps[0], v_id_deps[1],
                                                                    departure_time, vehicle_data,
                                                                    self.stations, charge_curves,
                                                                    args), standing_vehicles))

            # pair with socs
            standing_vehicles_w_soc = [(*v, rot_idx) for v, rot_idx in zip(standing_vehicles, socs)]

            standing_vehicles_w_soc = sorted(standing_vehicles_w_soc, key=lambda x: x[-1])
            consumption_soc = rot.calculate_consumption() / self.vehicle_types[vt][ct]["capacity"]
            for vehicle_id, depot, soc in standing_vehicles_w_soc:
                # end soc of vehicle if it services this rotation
                end_soc = soc - consumption_soc

                if end_soc > 0 or soc >= initial_soc or not station_is_electrified:
                    # Assign existing vehicle if rotation can be done or desired soc is reached.
                    # (Generating a new vehicle would not make a difference.)
                    # Special case non-electrified station:
                    # vehicles should not strand here even if rotation is not possible
                    rot.vehicle_id = vehicle_id
                    all_standing_vehicles.remove((vehicle_id, depot))
                    vehicle_data[vehicle_id] = {"soc": end_soc, "arrival_time": rot.arrival_time}
                    break
            else:
                # no vehicle available for dispatch, generate new one
                vehicle_type_counts[vt_ct] += 1
                vehicle_id = f"{vt_ct}_{vehicle_type_counts[vt_ct]}"
                rot.vehicle_id = vehicle_id
                end_soc = initial_soc - consumption_soc
                vehicle_data[vehicle_id] = {"soc": end_soc, "arrival_time": rot.arrival_time}

            rotations_in_progress.append(rot)
            try:
                # check if a next rotation exists to get the optional value of departure_time.
                # this is needed to calculate expected socs after greedy charging up to the
                # departure time
                next_rot = rotations[k + 1]
                departure_time = next_rot.departure_time
            except IndexError:
                # last rotation got a vehicle. Exit loop
                break
            # sort rotations in progress by their arrival time
            rotations_in_progress = sorted(rotations_in_progress, key=lambda x: x.arrival_time)

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

    def get_charge_curves(self, charge_levels, final_value: float, time_step_min: float) -> dict:
        """ Get the numeric charge curves.

        :param charge_levels: different power levels which clip the charge curves
        :type charge_levels: set[float]
        :param final_value: soc value at which the numeric calculation stops
        :type final_value: float
        :param time_step_min: time_step in minutes for the numeric calculation
        :type time_step_min: float
        :return: soc over time curves
        :rtype: dict
        """
        charge_curves = dict()
        for vehicle_name, vehicle_type in self.vehicle_types.items():
            charge_curves[vehicle_name] = dict()
            for ct_name, v_info in vehicle_type.items():
                charge_curves[vehicle_name][ct_name] = dict()
                obj = {"name": "default", "capacity": 1, "charging_curve": [[0, 1], [1, 1]]}
                default_eff = VehicleType(obj).battery_efficiency
                eff = v_info.get("battery_efficiency", default_eff)
                for charge_level in charge_levels:
                    curve = optimizer_util.charging_curve_to_soc_over_time(
                        v_info["charging_curve"],
                        v_info["capacity"],
                        final_value,
                        charge_level,
                        time_step=time_step_min,
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

    def get_height_difference(self, departure_name, arrival_name):
        """ Get the height difference of two stations.

        Defaults to 0 if height data is not found
        :param departure_name: Departure station
        :type departure_name: str
        :param arrival_name: Arrival station
        :type arrival_name: str
        :return: Height difference
        :rtype: float
        """
        if isinstance(self.station_data, dict):
            station = departure_name
            try:
                start_height = self.station_data[station]["elevation"]
                station = arrival_name
                end_height = self.station_data[arrival_name]["elevation"]
                return end_height - start_height
            except KeyError:
                logging.error(f"No elevation data found for {station}. Height Difference set to 0")
        else:
            logging.error("No Station Data found for schedule. Height Difference set to 0")
        return 0

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

    def get_total_distance(self):
        """ Calculate the total distance of all trips in the schedule.

        :return: total distance of schedule
        :rtype: float
        """
        total_distance = 0
        for rotation in self.rotations.values():
            total_distance += rotation.distance
        return total_distance

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
            arrival_of_last_trip = self.get_arrival_of_last_trip()
            stop_simulation = arrival_of_last_trip + interval
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
                    next_departure_time = vehicle_trips[i + 1].departure_time
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
                        # get desired soc by station type and trip
                        desired_soc = self.soc_dispatcher.get_soc(
                            vehicle_id=vehicle_id, trip=trip, station_type=station_type)
                except KeyError:
                    # non-electrified station
                    station_type = None

                # get buffer time from user configuration
                # buffer_time is an abstraction of delays like docking procedures and
                # is added to the planned arrival time
                # arrival + buffer time is clipped to the arrival of last trip and next departure
                if not station_type:
                    buffer_time = 0
                elif station_type == "deps":
                    buffer_time = util.get_buffer_time(trip=trip,
                                                       default=args.default_buffer_time_deps)
                else:
                    assert station_type == "opps"
                    buffer_time = util.get_buffer_time(trip=trip,
                                                       default=args.default_buffer_time_opps)
                # arrival event must occur no later than next departure and
                # one step before simulation terminates for arrival event to be taken into account
                arrival_time = min(trip.arrival_time + datetime.timedelta(minutes=buffer_time),
                                   next_departure_time, arrival_of_last_trip)

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
                    gc_name = trip.departure_name
                    try:
                        departure_station_type = self.stations[gc_name]["type"]
                    except KeyError:
                        departure_station_type = "deps"

                    vehicles[vehicle_id] = {
                        "connected_charging_station": None,
                        "estimated_time_of_departure": trip.departure_time.isoformat(),
                        "desired_soc": None,
                        "soc": self.soc_dispatcher.get_soc(vehicle_id=vehicle_id,
                                                           trip=None,
                                                           station_type=departure_station_type),
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

    @classmethod
    def get_dict_from_csv(cls, column, file_path, index):
        output = dict()
        with open(file_path, "r") as f:
            delim = util.get_csv_delim(file_path)
            reader = csv.DictReader(f, delimiter=delim)
            for row in reader:
                output[float(row[index])] = float(row[column])
        return output

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
    :type start_simulation: datetime.datetime
    :param stop_simulation: end of simulation
    :type stop_simulation: datetime.datetime
    :return: newly generated grid operator signals
    :rtype: list
    """
    day = datetime.timedelta(days=1)
    events = []
    MORNING_HOUR = 6
    # First price events on the day before the simulation.
    # This is needed if the simulation start time is between 00:00 and 06:00
    first_day = start_simulation - datetime.timedelta(hours=MORNING_HOUR)
    # Reset beginning to midnight,
    # so price events always start at 6:00 and get signaled at midnight each day
    first_day = first_day.replace(hour=0, minute=0, second=0)
    # Iterate over the simulation duration
    for current_time in util.daterange(first_day, stop_simulation, day):
        # create price events covering 24h from 6am onwards
        morning = current_time + datetime.timedelta(hours=MORNING_HOUR)
        evening_by_month = current_time + datetime.timedelta(
            hours=22 - abs(6 - current_time.month))
        events += [{
            # day (6 to evening): 15ct
            "signal_time": max(start_simulation, current_time - day).isoformat(),
            "grid_connector_id": gc_name,
            "start_time": morning.isoformat(),
            "cost": {
                "type": "fixed",
                "value": 0.15 + random.gauss(0, 0.05)
            }
        }, {
            # night (evening to 6 ): 5ct
            "signal_time": max(start_simulation, current_time - day).isoformat(),
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
                         duration_min: float, start_soc: float) -> float:
    """ Get the delta soc of a charge event for a given vehicle and charge type

    :param charge_curves: charge curves for a given vehicle_type, charge_type and max_power
    :type charge_curves: dict
    :param vt: vehicle type
    :type vt: str
    :param ct: charge type
    :type ct: str
    :param max_power: max power
    :type max_power: float
    :param duration_min: duration in minutes
    :type duration_min: float
    :param start_soc: start soc
    :type start_soc: float
    :return: delta soc of charge event
    :rtype: float
    """
    charge_curve = charge_curves[vt][ct][max_power]
    return optimizer_util.get_delta_soc(charge_curve, start_soc, duration_min=duration_min)


def soc_at_departure_time(v_id, deps, departure_time, vehicle_data, stations, charge_curves, args):
    """ Get the possible SoC of the vehicle at a specified departure_time

    :param v_id: vehicle_id
    :type v_id: str
    :param deps: depot name
    :type deps: str
    :param departure_time: time for which the soc is evaluated
    :type departure_time: datetime.datetime
    :param vehicle_data: data of used vehicles including soc and arrival_times
    :type vehicle_data: dict
    :param stations: station data from schedule
    :type stations: dict
    :param charge_curves: numeric charge_curves with soc over time in minutes
    :type charge_curves:  dict
    :param args: Arguments
    :type args: Namespace
    :return: soc at departure time for the given vehicle
    :rtype: float
    """

    vt = "_".join(v_id.split("_")[:-2])
    ct = v_id.split("_")[-2]

    start_soc = vehicle_data[v_id]["soc"]

    station_power = stations[deps].get(f"cs_power_deps_{ct}", vars(args).get(f"cs_power_deps_{ct}"))

    buffer_time = datetime.timedelta(minutes=vars(args).get("default_buffer_time_deps", 0))

    standing_duration = departure_time - vehicle_data[v_id]["arrival_time"] - buffer_time
    duration_in_m = standing_duration.total_seconds() / 60
    if ct == "oppb":
        # for opportunity chargers assume a minimal soc>=0
        start_soc = max(start_soc, 0)
    charge_delta_soc = \
        get_charge_delta_soc(charge_curves, vt, ct, station_power, duration_in_m,
                             start_soc)
    return min(start_soc + charge_delta_soc, args.desired_soc_deps)
