import csv
import random
import datetime
import warnings
from pathlib import Path

from ebus_toolbox import util
from ebus_toolbox.rotation import Rotation
from src.scenario import Scenario


class Schedule:

    def __init__(self, vehicle_types, stations, **kwargs):
        """Constructs Schedule object from CSV file containing all trips of schedule

        :param vehicle_types: Collection of vehicle types and their properties.
        :type vehicle_types: dict
        :param stations: electrified stations
        :type stations: dict

        :raises SystemExit: In case not all mandatory options are provided

        :param kwargs: Command line arguments
        :type kwargs: dict
        """

        self.stations = stations
        self.rotations = {}
        self.consumption = 0
        self.vehicle_types = vehicle_types
        self.original_rotations = None

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
            raise SystemExit("The following arguments are required: {}".format(", ".join(missing)))
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

        # define start and stop times
        start_simulation = \
            self.get_departure_of_first_trip() - datetime.timedelta(minutes=args.signal_time_dif)
        stop_simulation = self.get_arrival_of_last_trip() + interval
        if args.days is not None:
            stop_simulation = min(
                stop_simulation, start_simulation + datetime.timedelta(days=args.days))

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
                        # at opp station, always try to charge as much as you can
                        desired_soc = 1 if station_type == "opps" else args.desired_soc
                except KeyError:
                    # non-electrified station
                    station_type = None

                # get buffer time from user configuration
                # buffer time resembles amount of time deducted off of the planned standing
                # time. It may resemble things like delays and/or docking procedures
                # use buffer time from electrified stations JSON or in case none is
                # provided use global default from config file
                buffer_time = util.get_buffer_time(schedule=self,
                                                   trip=trip,
                                                   default=args.default_buffer_time_opps)
                arrival_time = min(trip.arrival_time + datetime.timedelta(minutes=buffer_time),
                                   next_departure_time)

                # total minutes spend at station
                standing_time = (next_departure_time - arrival_time).seconds / 60

                # connect to charging station
                # generate gc and cs if they dont exist
                if station_type is not None and standing_time >= args.min_charging_time:
                    # vehicle connects to charging station
                    connected_charging_station = f"{cs_name}_{station_type}"
                    # create charging station and grid connector if necessary
                    if cs_name not in charging_stations or gc_name not in grid_connectors:
                        number_cs = station["n_charging_stations"]
                        if station_type == "deps":
                            cs_power = args.cs_power_deps_oppb \
                                if trip.rotation.charging_type == 'oppb' \
                                else args.cs_power_deps_depb
                            gc_power = args.gc_power_deps
                        elif station_type == "opps":
                            cs_power = args.cs_power_opps
                            gc_power = args.gc_power_opps

                        # add one charging station for each bus at bus station
                        charging_stations[connected_charging_station] = {
                            "type": station_type,
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

                # initial condition of vehicle
                if i == 0:
                    vehicles[vehicle_id] = {
                        "connected_charging_station": None,
                        "estimated_time_of_departure": trip.departure_time.isoformat(),
                        "desired_soc": None,
                        "soc": args.desired_soc,
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
        for key in grid_connectors.keys():
            if not args.include_price_csv:
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

        if args.include_ext_load_csv:
            for filename, gc_name in args.include_ext_load_csv:
                options = {
                    "csv_file": filename,
                    "start_time": start_simulation.isoformat(),
                    "step_duration_s": 900,  # 15 minutes
                    "grid_connector_id": gc_name,
                    "column": "energy"
                }
                if args.include_ext_csv_option:
                    for key, value in args.include_ext_csv_option:
                        if key == "step_duration_s":
                            value = int(value)
                        options[key] = value
                events['external_load'][Path(filename).stem] = options
                # check if CSV file exists
                ext_csv_path = args.output_directory / filename
                if not ext_csv_path.exists():
                    print("Warning: external csv file '{}' does not exist yet".format(ext_csv_path))

        if args.include_feed_in_csv:
            for filename, gc_name in args.include_feed_in_csv:
                options = {
                    "csv_file": filename,
                    "start_time": start_simulation.isoformat(),
                    "step_duration_s": 3600,  # 60 minutes
                    "grid_connector_id": gc_name,
                    "column": "energy"
                }
                if args.include_feed_in_csv_option:
                    for key, value in args.include_feed_in_csv_option:
                        if key == "step_duration_s":
                            value = int(value)
                        options[key] = value
                events['energy_feed_in'][Path(filename).stem] = options
                feed_in_path = args.output_directory / filename
                if not feed_in_path.exists():
                    print("Warning: feed-in csv file '{}' does not exist yet".format(feed_in_path))

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
                "start_time": start_simulation.isoformat(),
                "interval": interval.days * 24 * 60 + interval.seconds // 60,
                "n_intervals": (stop_simulation - start_simulation) // interval
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
