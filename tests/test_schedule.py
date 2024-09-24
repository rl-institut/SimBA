import math
from argparse import Namespace
from copy import deepcopy
from datetime import timedelta, datetime # noqa
import pytest
import sys
import spice_ev.scenario as scenario

from simba.simulate import pre_simulation
from spice_ev.events import VehicleEvent
from tests.conftest import example_root, file_root
from tests.helpers import generate_basic_schedule
from simba import rotation, schedule, util
from simba.data_container import DataContainer

mandatory_args = {
    "min_recharge_deps_oppb": 1,
    "min_recharge_deps_depb": 1,
    "gc_power_opps": 1000,
    "gc_power_deps": 1000,
    "cs_power_opps": 100,
    "cs_power_deps_depb": 50,
    "cs_power_deps_oppb": 150,
}


class BasicSchedule:
    temperature_path = example_root / 'default_temp_winter.csv'
    lol_path = example_root / 'default_level_of_loading_over_day.csv'
    vehicle_types_path = example_root / "vehicle_types.json"
    electrified_stations_path = example_root / "electrified_stations.json"
    path_to_all_station_data = example_root / "all_stations.csv"

    @pytest.fixture
    def default_schedule_arguments(self):
        arguments = {"vehicle_types_path": self.vehicle_types_path,
                     "electrified_stations_path": self.electrified_stations_path,
                     "station_data_path": self.path_to_all_station_data,
                     "outside_temperature_over_day_path": self.temperature_path,
                     "level_of_loading_over_day_path": self.lol_path
                     }
        return arguments

    def basic_run(self):
        """Returns a schedule, scenario and args after running SimBA from config.
        :return: schedule, scenario, args
        """
        # set the system variables to imitate the console call with the config argument.
        # first element has to be set to something or error is thrown
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()

        data_container = DataContainer().fill_with_args(args)
        first_trip = min([t["arrival_time"] for t in data_container.trip_data])
        data_container.trip_data = [t for t in data_container.trip_data
                                    if t["arrival_time"] - first_trip < timedelta(hours=10)]
        sched, args = pre_simulation(args, data_container)
        scen = sched.run(args)
        return sched, scen, args


class TestSchedule(BasicSchedule):

    def test_optional_timeseries(self):
        # Test if simulation runs if level of loading and temperature timeseries is not given
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        data_container = DataContainer().fill_with_args(args)
        vehicle_mileage_path = None
        found = False
        for vt in data_container.vehicle_types_data.values():
            for ct in vt.values():
                if isinstance(ct["mileage"], str):
                    vehicle_mileage_path = ct["mileage"]
                    found = True
                    break
            if found:
                break

        sched, args = pre_simulation(args, data_container)

        # Make sure at least a single used vehicle type has mileage lookup.
        some_used_vehicle_type = next(iter(sched.rotations.values())).vehicle_type
        data_container.vehicle_types_data[some_used_vehicle_type][
            "oppb"]["mileage"] = vehicle_mileage_path
        data_container.vehicle_types_data[some_used_vehicle_type][
            "depb"]["mileage"] = vehicle_mileage_path

        # Delete the lol and temp data sources -> pre-simulation should fail,
        # since consumption cannot be calculated
        data_container.level_of_loading_data = {}
        data_container.temperature_data = {}
        for trip in data_container.trip_data:
            print(trip)
        with pytest.raises(Exception):
            sched, args = pre_simulation(args, data_container)

        # if all vehicle types have constant consumption, pre-simulation should work
        for vt in data_container.vehicle_types_data.values():
            for ct in vt.values():
                ct["mileage"] = 1
        _, _ = pre_simulation(args, data_container)

    def test_timestep(self):
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        data_container = DataContainer().fill_with_args(args)

        assert args.interval == 1

        # Reduce rotations to increase simulation / test speed
        some_rotation = data_container.trip_data[0]["rotation_id"]
        data_container.trip_data = [trip_dict for trip_dict in data_container.trip_data
                                    if trip_dict["rotation_id"] == some_rotation]

        first_trip = None
        last_trip = None
        for trip in data_container.trip_data:
            if first_trip is None or trip["departure_time"] < first_trip["departure_time"]:
                first_trip = trip
            if last_trip is None or trip["arrival_time"] > last_trip["arrival_time"]:
                last_trip = trip

        # Make sure the duration is not an integer multiple of interval
        first_trip["departure_time"] -= timedelta(minutes=1)
        first_trip["departure_time"] = first_trip["departure_time"].replace(second=0)
        last_trip["arrival_time"] += timedelta(minutes=1)
        last_trip["arrival_time"] = last_trip["arrival_time"].replace(second=30)
        duration_in_s = (last_trip["arrival_time"] - first_trip["departure_time"]).total_seconds()
        assert duration_in_s % (args.interval * 60) != 0

        n_steps = math.ceil((duration_in_s + args.signal_time_dif * 60) / (args.interval * 60)) + 1
        sched, args = pre_simulation(args, data_container)
        scenario = sched.generate_scenario(args)
        assert n_steps == scenario.n_intervals
        scenario.run("distributed", vars(args).copy())
        assert scenario.strat.current_time >= last_trip["arrival_time"]
        assert all(
            not isinstance(e, VehicleEvent) for e in scenario.strat.world_state.future_events)

        # Make sure that if the last step would not have been simulated,
        # a vehicle_event would still be up for simulation
        scenario = sched.generate_scenario(args)
        scenario.n_intervals -= 1
        scenario.run("distributed", vars(args).copy())
        assert scenario.strat.current_time < last_trip["arrival_time"]
        assert any(isinstance(e, VehicleEvent) for e in scenario.strat.world_state.future_events)

    def test_mandatory_options_exit(self):
        """
        Check if the schedule creation properly throws an error in case of missing mandatory options
        """
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        data_container = DataContainer().fill_with_args(args)

        for key in mandatory_args.keys():
            args = util.get_args()
            args.__delattr__(key)
            with pytest.raises(Exception):
                # schedule creation without mandatory arg
                schedule.Schedule.from_datacontainer(data_container, args)

    def test_station_data_reading(self, caplog):
        """ Test if the reading of the geo station data works and outputs warnings in
        case the data was problematic, e.g. not numeric or not existent

        :param caplog: pytest fixture to capture logging
        """
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        data_container = DataContainer().fill_with_args(args)

        generated_schedule, args = pre_simulation(args, data_container)
        assert generated_schedule.station_data is not None

        # check if reading a non valid station.csv throws warnings
        args.station_data_path = file_root / "not_existent_file"
        with pytest.raises(FileNotFoundError):
            data_container = DataContainer().fill_with_args(args)

        args.station_data_path = file_root / "not_numeric_stations.csv"
        with pytest.raises(ValueError):
            data_container = DataContainer().fill_with_args(args)

    def test_basic_run(self):
        """ Check if running a basic example works and if a scenario object is returned
        """
        sched, scen, args = self.basic_run()
        assert type(scen) is scenario.Scenario

    def test_assign_vehicles_fixed_recharge(self):
        """ Test if assigning vehicles works as intended using the fixed_recharge strategy
        """

        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.min_recharge_deps_oppb = 1
        args.min_recharge_deps_depb = 1
        args.schedule_path = file_root / "trips_assign_vehicles_extended.csv"
        data_container = DataContainer().fill_with_args(args)

        generated_schedule, args = pre_simulation(args, data_container)

        all_rotations = [r for r in generated_schedule.rotations]
        args.assign_strategy = "fixed_recharge"
        generated_schedule.assign_vehicles(args)
        gen_rotations = generated_schedule.rotations
        vehicle_ids = {rot.vehicle_id for key, rot in gen_rotations.items()}
        assert len(vehicle_ids) == 6

        # only check the first day
        for r in all_rotations:
            if "_2" in r or "_3" in r:
                del generated_schedule.rotations[r]

        args.assign_strategy = "fixed_recharge"
        generated_schedule.assign_vehicles(args)
        gen_rotations = generated_schedule.rotations
        vehicle_ids = {rot.vehicle_id for key, rot in gen_rotations.items()}

        assert len(vehicle_ids) == 4

        # Assertions based on the following output  / graphic
        # day 1
        # two opp rotations which would match but min recharge does not allow same vehicle
        # assignment
        assert gen_rotations["1_1"].vehicle_id != gen_rotations["2_1"].vehicle_id
        # longer rotation which cant be serviced by the previous
        assert gen_rotations["3_1a"].vehicle_id != gen_rotations["1_1"].vehicle_id
        assert gen_rotations["3_1a"].vehicle_id != gen_rotations["2_1"].vehicle_id

        assert gen_rotations["3_1b"].vehicle_id == gen_rotations["2_1"].vehicle_id
        assert gen_rotations["3_1c"].vehicle_id == gen_rotations["1_1"].vehicle_id
        assert gen_rotations["3_1d"].vehicle_id == gen_rotations["3_1a"].vehicle_id
        # depot Rotation
        assert gen_rotations["4_1"].vehicle_id == gen_rotations["5_1"].vehicle_id

    def test_assign_vehicles_adaptive(self):
        """ Test if assigning vehicles works as intended using the adaptive strategy
        """
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.schedule_path = file_root / "trips_assign_vehicles_extended.csv"
        data_container = DataContainer().fill_with_args(args)

        generated_schedule, args = pre_simulation(args, data_container)

        args = Namespace(**{"desired_soc_deps": 1})
        args.assign_strategy = None
        generated_schedule.assign_vehicles(args)
        gen_rotations = generated_schedule.rotations
        vehicle_ids = {rot.vehicle_id for key, rot in gen_rotations.items()}
        assert len(vehicle_ids) == 4

        # Assertions based on the following output  / graphic
        # https://github.com/rl-institut/SimBA/pull/177#issuecomment-2066447393
        # day 1
        # two opp rotations which match
        assert gen_rotations["1_1"].vehicle_id == gen_rotations["2_1"].vehicle_id
        # longer rotation which cant be serviced by the previous
        assert gen_rotations["3_1a"].vehicle_id != gen_rotations["1_1"].vehicle_id
        assert gen_rotations["3_1b"].vehicle_id == gen_rotations["1_1"].vehicle_id
        assert gen_rotations["3_1c"].vehicle_id == gen_rotations["1_1"].vehicle_id
        # first vehicle charged enough
        assert gen_rotations["3_1d"].vehicle_id == gen_rotations["3_1a"].vehicle_id
        # depot Rotation
        assert gen_rotations["4_1"].vehicle_id == gen_rotations["5_1"].vehicle_id
        # day 2
        # Identical rotation to 1_1 and 1_2 --> same assignments
        assert gen_rotations["1_2"].vehicle_id == gen_rotations["1_1"].vehicle_id
        assert gen_rotations["2_2"].vehicle_id == gen_rotations["1_1"].vehicle_id
        # 3_2a starts a little after 2_2 ended. But soc does not allow for the same vehicle
        assert gen_rotations["3_2a"].vehicle_id != gen_rotations["1_1"].vehicle_id
        # depot rotations which barely allow for the same vehicle with the final soc at almost 0%
        assert gen_rotations["4_2"].vehicle_id == gen_rotations["5_2"].vehicle_id
        # day 3
        # vehicle of 6_3 ends with enough soc and can also charge soc enough for 7_3a
        assert gen_rotations["6_3a"].vehicle_id == gen_rotations["7_3a"].vehicle_id
        # opportunity 6_3b ends in negative soc. The extra charge is bigger than the consumption of
        # 7_3b therefore it is used again. In this case this leads to a negative rotation 7_3b.
        # This is allowed since fixing 6_3b will also fix 7_3b.
        assert gen_rotations["6_3b"].vehicle_id == gen_rotations["7_3b"].vehicle_id

        # depot rotations are not assigned when their charged energy is above the next rotations
        # energy consumption, but only when they reach a soc which is enough or full
        assert gen_rotations["6_3c"].vehicle_id != gen_rotations["7_3c"].vehicle_id
        assert gen_rotations["6_3c"].vehicle_id == gen_rotations["8_3c"].vehicle_id

    def test_calculate_consumption(self, default_schedule_arguments):
        """ Test if calling the consumption calculation works
        :param default_schedule_arguments: basic arguments the schedule needs for creation
        """
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.schedule_path = file_root / "trips_assign_vehicles.csv"
        data_container = DataContainer().fill_with_args(args)

        generated_schedule, args = pre_simulation(args, data_container)

        # set mileage to a constant
        mileage = 10
        for v_typ in generated_schedule.vehicle_types.values():
            for charge_typ in v_typ.values():
                charge_typ['mileage'] = mileage
        calc_consumption = generated_schedule.calculate_consumption()
        # set mileage to half of the previous constant
        for v_typ in generated_schedule.vehicle_types.values():
            for charge_typ in v_typ.values():
                charge_typ['mileage'] = mileage / 2
        assert calc_consumption == generated_schedule.calculate_consumption() * 2

    def test_get_common_stations(self, default_schedule_arguments):
        """Test if getting common_stations works. Rotation 1 is on the first day, rotation 2 and 3
        on the second day. rotation 1 should not share any stations with other rotations and
        2 and 3 are almost simultaneous.

        :param default_schedule_arguments: basic arguments the schedule needs for creation
        """
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.schedule_path = file_root / "trips_assign_vehicles.csv"
        data_container = DataContainer().fill_with_args(args)
        generated_schedule, args = pre_simulation(args, data_container)

        common_stations = generated_schedule.get_common_stations(only_opps=False)
        assert len(common_stations["1"]) == 0
        assert len(common_stations["2"]) == 1
        assert len(common_stations["3"]) == 1
        assert '3' in (common_stations["2"])
        assert '2' in (common_stations["3"])

    def test_get_negative_rotations(self):
        """Check if rotation '11' with negative SOCs is found """
        # make use of the test_run() which has to return schedule and scenario object
        sched, scen, args = self.basic_run()
        for rot in sched.rotations.values():
            for t in rot.trips:
                t.distance = 0.01
        sched.rotations["1"].trips[-1].distance = 99_999
        sched.calculate_consumption()
        scen = sched.run(args)
        neg_rots = sched.get_negative_rotations(scen)
        assert ['1'] == neg_rots

    def test_rotation_filter(self, tmp_path):
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.schedule_path = file_root / "trips_assign_vehicles.csv"
        data_container = DataContainer().fill_with_args(args)

        s, args = pre_simulation(args, data_container)

        args = Namespace(**{
            "rotation_filter_variable": None,
            "rotation_filter_path": None,
        })
        # add dummy rotations
        s.rotations = {
            str(i): rotation.Rotation(id=str(i), vehicle_type="", schedule=None)
            for i in range(6)
        }
        s.original_rotations = deepcopy(s.rotations)
        # filtering disabled
        args.rotation_filter_variable = None
        s.rotation_filter(args)
        assert s.rotations.keys() == s.original_rotations.keys()

        # filtering not disabled, but neither file nor list given -> warning
        args.rotation_filter_variable = "include"
        args.rotation_filter_path = None
        with pytest.warns(UserWarning):
            s.rotation_filter(args)
        assert s.rotations.keys() == s.original_rotations.keys()

        # filter file not found -> warning
        args.rotation_filter_path = tmp_path / "filter.txt"
        with pytest.warns(UserWarning):
            s.rotation_filter(args)
        assert s.rotations.keys() == s.original_rotations.keys()

        # filter (include) from JSON file
        args.rotation_filter_path.write_text("3 \n 4\n16")
        s.rotation_filter(args)
        assert sorted(s.rotations.keys()) == ['3', '4']

        # filter (exclude) from given list
        args.rotation_filter_variable = "exclude"
        args.rotation_filter_path = None
        s.rotations = deepcopy(s.original_rotations)
        s.rotation_filter(args, rf_list=['3', '4'])
        assert sorted(s.rotations.keys()) == ['0', '1', '2', '5']

        # filter (include) from integer list
        s.rotations = deepcopy(s.original_rotations)
        args.rotation_filter_variable = "include"
        s.rotation_filter(args, rf_list=[3, 4])
        assert sorted(s.rotations.keys()) == ['3', '4']

        # filter nothing
        s.rotations = deepcopy(s.original_rotations)
        args.rotation_filter_path = tmp_path / "filter.txt"
        args.rotation_filter_path.write_text('')
        args.rotation_filter_variable = "exclude"
        s.rotation_filter(args, rf_list=[])
        assert s.rotations.keys() == s.original_rotations.keys()

        # filter all (is this intended?)
        args.rotation_filter_variable = "include"
        s.rotation_filter(args, rf_list=[])
        assert not s.rotations

    def test_scenario_with_feed_in(self, default_schedule_arguments):
        """ Check if running a example with an extended electrified stations file
         with feed in, external load and battery works and if a scenario object is returned

        :param default_schedule_arguments: basic arguments the schedule needs for creation
        """
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        data_container = DataContainer().fill_with_args(args)
        generated_schedule, args = pre_simulation(args, data_container)

        scen = generated_schedule.generate_scenario(args)
        assert "Station-0" in scen.components.photovoltaics
        assert "Station-3" in scen.components.photovoltaics
        assert "Station-0 storage" in scen.components.batteries
        assert scen.components.batteries["Station-0 storage"].capacity == 300
        assert scen.components.batteries["Station-0 storage"].efficiency == 0.95
        assert scen.components.batteries["Station-0 storage"].min_charging_power == 0

        scen = generated_schedule.run(args)
        assert type(scen) is scenario.Scenario

        # Change a station
        station = data_container.stations_data["Station-0"]
        station["energy_feed_in"]["csv_file"] = file_root / "notafile"
        station["external_load"]["csv_file"] = file_root / "notafile"

        generated_schedule, args = pre_simulation(args, data_container)

        # check that 2 user warnings are put out for missing files and an error is thrown
        with pytest.warns(Warning) as record:
            try:
                scen = generated_schedule.generate_scenario(args)
            except FileNotFoundError:
                user_warning_count = sum([1 for warning in record.list
                                          if warning.category == UserWarning])
                assert user_warning_count == 2
            else:
                assert 0, "No error despite wrong file paths"

    def test_schedule_from_datacontainer(self):
        generated_schedule, _ = generate_basic_schedule()
        assert len(generated_schedule.rotations) == 8
        assert type(generated_schedule) is schedule.Schedule

    def test_consistency(self):
        sched, _ = generate_basic_schedule()
        # check if no error is thrown in the basic case
        assert len(schedule.Schedule.check_consistency(sched)) == 0

        error = "Trip time is negative"
        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_trip = faulty_rot.trips[0]
        # create error through moving trip arrival 1 day before departure
        faulty_trip.arrival_time = faulty_trip.departure_time - timedelta(days=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Break time is negative"
        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_trip = faulty_rot.trips[1]
        # create error through moving trip departure before last arrival
        faulty_trip.departure_time = faulty_rot.trips[0].arrival_time-timedelta(minutes=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Trips are not sequential"
        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_rot.trips[1].arrival_name = "foo"
        faulty_rot.trips[0].departure_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Start and end of rotation differ"
        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        departure_trip = list(faulty_rot.trips)[0]
        departure_trip.departure_name = "foo"
        faulty_rot.trips[-1].arrival_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Rotation data differs from trips data"

        # check arrival data in rotation
        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_rot.trips[-1].arrival_name = "foo"
        faulty_rot.trips[0].departure_name = "foo"
        faulty_rot.arrival_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        arrival_trip = faulty_rot.trips[-1]
        faulty_rot.arrival_time = arrival_trip.arrival_time - timedelta(minutes=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        # check departure data in rotation
        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_rot.trips[-1].arrival_name = "foo"
        faulty_rot.trips[0].departure_name = "foo"
        faulty_rot.departure_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        sched, _ = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        departure_trip = faulty_rot.trips[0]
        faulty_rot.departure_time = departure_trip.departure_time - timedelta(minutes=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

    def test_peak_load_window(self):
        # generate events to lower GC max power during peak load windows
        # setup basic schedule (reuse during test)
        generated_schedule, args = generate_basic_schedule()
        generated_schedule.init_soc_dispatcher(args)
        # needed for timeseries, but default is false
        args.cost_calculation = True
        for station in generated_schedule.stations.values():
            station["gc_power"] = 1000
            station.pop("peak_load_window_power", None)

        def count_max_power_events(scenario):
            # count how many grid operator events (re)set GC max_power
            return sum([e.max_power is not None for e in scenario.events.grid_operator_signals])

        # no time windows
        args.time_windows = None
        scenario = generated_schedule.generate_scenario(args)
        assert count_max_power_events(scenario) == 0

        # wrong time windows file
        args.time_windows = "does-not-exist"
        with pytest.warns(UserWarning):
            generated_schedule.generate_scenario(args)  # warning, no events

        # example time windows, but no corresponding grid operators
        args.time_windows = example_root / "time_windows.json"
        for station in generated_schedule.stations.values():
            station["grid_operator"] = "wrong-operator"
        with pytest.warns(UserWarning):
            scenario = generated_schedule.generate_scenario(args)
        # reset grid operator
        for station in generated_schedule.stations.values():
            station["grid_operator"] = "default_grid_operator"
        assert count_max_power_events(scenario) == 0

        # reduced power too high
        args.peak_load_window_power_deps = 10000
        args.peak_load_window_power_opps = 10000
        with pytest.warns(UserWarning):
            scenario = generated_schedule.generate_scenario(args)
        assert count_max_power_events(scenario) == 0

        # test reduced power events
        args.peak_load_window_power_deps = 10
        args.peak_load_window_power_opps = 10
        scenario = generated_schedule.generate_scenario(args)
        assert count_max_power_events(scenario) == 8

        # test that max_power is actually reduced during simulation
        args.time_windows = None
        basic_run = generated_schedule.run(args)
        args.time_windows = example_root / "time_windows.json"
        args.peak_load_window_power_deps = 75
        args.peak_load_window_power_opps = 75
        reduced_run = generated_schedule.run(args)

        window_start = datetime(year=2022, month=3, day=8, hour=3, minute=0)
        window_end = datetime(year=2022, month=3, day=8, hour=4, minute=55)
        start_index = (window_start - scenario.start_time) // scenario.interval
        end_index = (window_end - scenario.start_time) // scenario.interval
        idx_slice = slice(start_index, end_index, 1)
        timeseries_no_reduction = getattr(basic_run, "Station-0_timeseries")
        sum_grid_power_no_red = -sum(timeseries_no_reduction["grid supply [kW]"][idx_slice])
        timeseries_with_reduction = getattr(reduced_run, "Station-0_timeseries")
        sum_grid_power_with_red = -sum(timeseries_with_reduction["grid supply [kW]"][idx_slice])
        assert sum_grid_power_no_red > sum_grid_power_with_red

    def test_generate_price_lists(self):
        # setup basic schedule
        generated_schedule, args = generate_basic_schedule()

        # only test individual price CSV and random price generation
        args.include_price_csv = None
        # Station-0: all options
        generated_schedule.stations["Station-0"]["price_csv"] = {
            "csv_file": example_root / "price_timeseries.csv",
            "start_time": "2022-03-07 00:00:00",
            "step_duration_s": 86400,
            "column": "price",
            "factor": 2
        }
        # Station-3: minimal options (only CSV path)
        generated_schedule.stations["Station-3"]["price_csv"] = {
            "csv_file": example_root / "price_timeseries.csv"
        }
        # Station-10: wrong option (wrong CSV file)
        generated_schedule.stations["Station-10"]["price_csv"] = {
            "csv_file": example_root / "does-not-exist.csv"
        }
        # Station-21: start after end of schedule
        generated_schedule.stations["Station-21"]["price_csv"] = {
            "csv_file": example_root / "price_timeseries.csv",
            "start_time": "3333-03-03 00:00:00",
            "step_duration_s": 3600
        }
        scenario = generated_schedule.generate_scenario(args)
        events = [e for e in scenario.events.grid_operator_signals if e.cost is not None]
        events_by_gc = {
            gc: [e for e in events if e.grid_connector_id == gc]
            for gc in generated_schedule.stations.keys()}
        assert len(events_by_gc["Station-0"]) > 0
        # first entry is 0.07995, factor is 2
        assert events_by_gc["Station-0"][0].cost["value"] == 0.1599
        assert len(events_by_gc["Station-3"]) > 0
        # default factor of 1, price at 20:00 (example scenario start)
        assert events_by_gc["Station-3"][0].cost["value"] == 0.2107
        # wrong file: no events
        assert len(events_by_gc["Station-10"]) == 0
        # after scenario: no events
        assert len(events_by_gc["Station-21"]) == 0

        # run schedule and check prices
        # example scenario covers 2022-03-07 20:15 - 2022-03-09 04:59
        scenario = generated_schedule.run(args)
        # Station-0: price change every 24h, starting midnight => two price changes at midnight
        assert set(scenario.prices["Station-0"]) == {0.1599, 0.18404, 0.23492}
        assert scenario.prices["Station-0"][224:226] == [0.1599, 0.18404]
        # Station-3: price change every hour, starting 04:00 (from csv timestamp)
        # => 32 price changes
        assert len(set(scenario.prices["Station-3"])) == 33
        # same price for last 59 minutes
        assert set(scenario.prices["Station-3"][-59:]) == {0.15501}

    def test_get_price_list_from_csv(self, tmp_path):
        # wrong file given
        assert len(schedule.get_price_list_from_csv({'csv_file': 'does-not-exist'})) == 0

        # generate tmp csv file
        with open(tmp_path / 'price.csv', 'w') as f:
            f.write('\n'.join([
                'date,price',
                '2022-03-07 00:00:00, 1',
                '2022-03-07 01:00:00, 2',
                '2022-03-07 02:00:00, 3']))

        # just file given
        prices = schedule.get_price_list_from_csv({'csv_file': tmp_path / 'price.csv'})
        assert len(prices) == 3
        assert prices[0] == ('2022-03-07 00:00:00', 1)

        # change column
        with pytest.raises(KeyError):
            schedule.get_price_list_from_csv({
                'csv_file': tmp_path / 'price.csv',
                'column': 'does-not-exist'
            })
        assert len(schedule.get_price_list_from_csv({
            'csv_file': tmp_path / 'price.csv',
            'column': 'price'
        })) == 3

        # change factor
        prices = schedule.get_price_list_from_csv({
            'csv_file': tmp_path / 'price.csv',
            'factor': 1.5
        })
        assert prices[0][1] == 1.5

    def test_generate_event_list_from_prices(self):
        prices = [
            ('2022-03-07 00:00:00', 1),
            ('2022-03-07 01:00:00', 2),
            ('2022-03-07 02:00:00', 3)]
        start = datetime.fromisoformat('2022-03-07')
        stop = datetime.fromisoformat('2022-03-08')

        # no prices: no events
        assert len(schedule.generate_event_list_from_prices([], 'GC', start, stop)) == 0

        # basic functionality: all events
        events = schedule.generate_event_list_from_prices(prices, 'GC', start, stop)
        assert len(events) == 3

        # change list start time
        # in-between: all events, different time
        events = schedule.generate_event_list_from_prices(
            prices, 'GC', start, stop,
            start_events='2022-03-07 01:30:00',
            price_interval_s=60
        )
        assert len(events) == 3
        assert events[-1]["start_time"] == '2022-03-07T01:32:00'

        # before: only last event
        events = schedule.generate_event_list_from_prices(
            prices, 'GC', start, stop,
            start_events='1970-01-01',
            price_interval_s=3600
        )
        assert len(events) == 1 and events[0]["cost"]["value"] == 3
        # change start date after price list: only last event
        start = datetime.fromisoformat('2022-03-07 12:00:00')
        events = schedule.generate_event_list_from_prices(prices, 'GC', start, stop)
        assert len(events) == 1 and events[0]["cost"]["value"] == 3

        # after: no events
        events = schedule.generate_event_list_from_prices(
            prices, 'GC', start, stop,
            start_events='3000-01-01',
            price_interval_s=3600
        )
        assert len(events) == 0

    def test_rotation_consumption_calc(self):
        s, args = generate_basic_schedule()
        rot_iter = iter(s.rotations.values())
        r = next(rot_iter)
        r.trips = []
        assert s.calculate_rotation_consumption(r) == 0

        some_rot = next(rot_iter)
        first_trip = some_rot.trips[0]
        del first_trip.rotation
        r.add_trip(vars(first_trip))
        r.trips[0].consumption = None
        s.calculate_rotation_consumption(r)
        assert r.trips[0].consumption is not None
        second_trip = some_rot.trips[1]
        del second_trip.rotation
        r.add_trip(vars(second_trip))
        for trip in r.trips:
            trip.consumption = None
        s.calculate_rotation_consumption(r)
        for trip in r.trips:
            assert trip.consumption is not None
