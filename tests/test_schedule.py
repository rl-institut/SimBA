from argparse import Namespace
from copy import deepcopy
from datetime import timedelta
import pytest
import sys
import spice_ev.scenario as scenario
from spice_ev.util import set_options_from_config

from tests.conftest import example_root, file_root
from tests.helpers import generate_basic_schedule
from simba import consumption, rotation, schedule, trip, util


mandatory_args = {
    "min_recharge_deps_oppb": 0,
    "min_recharge_deps_depb": 0,
    "gc_power_opps": 1000,
    "gc_power_deps": 1000,
    "cs_power_opps": 100,
    "cs_power_deps_depb": 50,
    "cs_power_deps_oppb": 150
}


class TestSchedule:
    temperature_path = example_root / 'default_temp_winter.csv'
    lol_path = example_root / 'default_level_of_loading_over_day.csv'

    with open(example_root / "electrified_stations.json", "r", encoding='utf-8') as file:
        electrified_stations = util.uncomment_json_file(file)

    with open(example_root / "vehicle_types.json", "r", encoding='utf-8') as file:
        vehicle_types = util.uncomment_json_file(file)

    def basic_run(self):
        """Returns a schedule and scenario after running SimBA.
        :return: schedule, scenario
        """
        path_to_trips = example_root / "trips_example.csv"
        # set the system variables to imitate the console call with the config argument.
        # first element has to be set to something or error is thrown
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.config = example_root / "simba.cfg"
        args.days = None
        args.seed = 5

        trip.Trip.consumption = consumption.Consumption(
            self.vehicle_types, outside_temperatures=self.temperature_path,
            level_of_loading_over_day=self.lol_path)

        path_to_all_station_data = example_root / "all_stations.csv"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args,
            station_data_path=path_to_all_station_data)

        set_options_from_config(args, verbose=False)
        args.ALLOW_NEGATIVE_SOC = True
        args.attach_vehicle_soc = True

        scen = generated_schedule.run(args)
        return generated_schedule, scen

    def test_mandatory_options_exit(self):
        """
        Check if the schedule creation properly throws an error in case of missing mandatory options
        """
        args = mandatory_args.copy()
        for key in mandatory_args.keys():
            value = args.pop(key)
            with pytest.raises(Exception):
                # schedule creation without mandatory arg
                schedule.Schedule(self.vehicle_types, self.electrified_stations, **args)
            args[key] = value

    def test_station_data_reading(self):
        """ Test if the reading of the geo station data works and outputs warnings in
        case the data was problematic, e.g. not numeric or not existent"""
        path_to_trips = example_root / "trips_example.csv"

        trip.Trip.consumption = consumption.Consumption(
            self.vehicle_types, outside_temperatures=self.temperature_path,
            level_of_loading_over_day=self.lol_path)

        path_to_all_station_data = example_root / "all_stations.csv"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args,
            station_data_path=path_to_all_station_data)
        assert generated_schedule.station_data is not None

        # check if reading a non valid station.csv throws warnings
        with pytest.warns(Warning) as record:
            path_to_all_station_data = file_root / "not_existent_file"
            schedule.Schedule.from_csv(
                path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args,
                station_data_path=path_to_all_station_data)
            assert len(record) == 1

            path_to_all_station_data = file_root / "not_numeric_stations.csv"
            schedule.Schedule.from_csv(
                path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args,
                station_data_path=path_to_all_station_data)
            assert len(record) == 2

    def test_basic_run(self):
        """ Check if running a basic example works and if a scenario object is returned
        """
        schedule, scen = self.basic_run()
        assert type(scen) is scenario.Scenario

    def test_assign_vehicles(self):
        """ Test if assigning vehicles works as intended.

        Use a trips csv with two rotations ("1","2") a day apart.
        SimBA should assign the same vehicle to both of them.
        Rotation "3" starts shortly after "2" and should be a new vehicle.
        """

        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_trips = file_root / "trips_assign_vehicles.csv"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args)
        generated_schedule.assign_vehicles()
        gen_rotations = generated_schedule.rotations
        assert gen_rotations["1"].vehicle_id == gen_rotations["2"].vehicle_id
        assert gen_rotations["1"].vehicle_id != gen_rotations["3"].vehicle_id

    def test_calculate_consumption(self):
        """ Test if calling the consumption calculation works
        """
        trip.Trip.consumption = consumption.Consumption(
            self.vehicle_types, outside_temperatures=self.temperature_path,
            level_of_loading_over_day=self.lol_path)

        path_to_trips = file_root / "trips_assign_vehicles.csv"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args)

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

    def test_get_common_stations(self):
        """Test if getting common_stations works. Rotation 1 is on the first day, rotation 2 and 3
        on the second day. rotation 1 should not share any stations with other rotations and
        2 and 3 are almost simultaneous.
        """
        trip.Trip.consumption = consumption.Consumption(
            self.vehicle_types, outside_temperatures=self.temperature_path,
            level_of_loading_over_day=self.lol_path)

        path_to_trips = file_root / "trips_assign_vehicles.csv"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, self.electrified_stations, **mandatory_args)

        common_stations = generated_schedule.get_common_stations(only_opps=False)
        assert len(common_stations["1"]) == 0
        assert len(common_stations["2"]) == 1
        assert len(common_stations["3"]) == 1
        assert '3' in (common_stations["2"])
        assert '2' in (common_stations["3"])

    def test_get_negative_rotations(self):
        """Check if the single rotation '1' with a negative soc is found """

        # make use of the test_run() which has to return schedule and scenario object
        sched, scen = self.basic_run()

        neg_rots = sched.get_negative_rotations(scen)
        assert '1' in neg_rots

    def test_rotation_filter(self, tmp_path):
        s = schedule.Schedule(self.vehicle_types, self.electrified_stations, **mandatory_args)
        args = Namespace(**{
            "rotation_filter_variable": None,
            "rotation_filter": None,
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
        args.rotation_filter = None
        with pytest.warns(UserWarning):
            s.rotation_filter(args)
        assert s.rotations.keys() == s.original_rotations.keys()

        # filter file not found -> warning
        args.rotation_filter = tmp_path / "filter.txt"
        with pytest.warns(UserWarning):
            s.rotation_filter(args)
        assert s.rotations.keys() == s.original_rotations.keys()

        # filter (include) from JSON file
        args.rotation_filter.write_text("3 \n 4\n16")
        s.rotation_filter(args)
        assert sorted(s.rotations.keys()) == ['3', '4']

        # filter (exclude) from given list
        args.rotation_filter_variable = "exclude"
        args.rotation_filter = None
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
        args.rotation_filter = tmp_path / "filter.txt"
        args.rotation_filter.write_text('')
        args.rotation_filter_variable = "exclude"
        s.rotation_filter(args, rf_list=[])
        assert s.rotations.keys() == s.original_rotations.keys()

        # filter all (is this intended?)
        args.rotation_filter_variable = "include"
        s.rotation_filter(args, rf_list=[])
        assert not s.rotations

    def test_scenario_with_feed_in(self):
        """ Check if running a example with an extended electrified stations file
         with feed in, external load and battery works and if a scenario object is returned"""

        path_to_trips = example_root / "trips_example.csv"
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.config = example_root / "simba.cfg"
        electrified_stations_path = example_root / "electrified_stations_with_feeds.json"
        args.electrified_stations = electrified_stations_path
        with open(electrified_stations_path, "r", encoding='utf-8') as file:
            electrified_stations = util.uncomment_json_file(file)

        args.days = None
        args.seed = 5

        trip.Trip.consumption = consumption.Consumption(
            self.vehicle_types, outside_temperatures=self.temperature_path,
            level_of_loading_over_day=self.lol_path)

        path_to_all_station_data = example_root / "all_stations.csv"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, electrified_stations, **mandatory_args,
            station_data_path=path_to_all_station_data)

        set_options_from_config(args, verbose=False)
        args.ALLOW_NEGATIVE_SOC = True
        args.attach_vehicle_soc = True
        scen = generated_schedule.generate_scenario(args)
        assert "Station-0" in scen.components.photovoltaics
        assert "Station-3" in scen.components.photovoltaics
        assert "Station-0" in scen.components.batteries
        assert scen.components.batteries["Station-0"].capacity == 300
        assert scen.components.batteries["Station-0"].efficiency == 0.95
        assert scen.components.batteries["Station-0"].min_charging_power == 0
        scen = generated_schedule.run(args)
        assert type(scen) is scenario.Scenario

        with open(electrified_stations_path, "r", encoding='utf-8') as file:
            electrified_stations = util.uncomment_json_file(file)

        electrified_stations["Station-0"]["energy_feed_in"]["csv_file"] = file_root / "notafile"
        electrified_stations["Station-0"]["external_load"]["csv_file"] = file_root / "notafile"
        generated_schedule = schedule.Schedule.from_csv(
            path_to_trips, self.vehicle_types, electrified_stations, **mandatory_args,
            station_data_path=path_to_all_station_data)

        set_options_from_config(args, verbose=False)

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

    def test_schedule_from_csv(self):
        generated_schedule = generate_basic_schedule()
        assert len(generated_schedule.rotations) == 4
        assert type(generated_schedule) is schedule.Schedule

    def test_consistency(self):
        sched = generate_basic_schedule()
        # check if no error is thrown in the basic case
        assert len(schedule.Schedule.check_consistency(sched)) == 0

        error = "Trip time is negative"
        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_trip = faulty_rot.trips[0]
        # create error through moving trip arrival 1 day before departure
        faulty_trip.arrival_time = faulty_trip.departure_time - timedelta(days=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Break time is negative"
        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_trip = faulty_rot.trips[1]
        # create error through moving trip departure before last arrival
        faulty_trip.departure_time = faulty_rot.trips[0].arrival_time-timedelta(minutes=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Trips are not sequential"
        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_rot.trips[1].arrival_name = "foo"
        faulty_rot.trips[0].departure_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Start and end of rotation differ"
        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        departure_trip = list(faulty_rot.trips)[0]
        departure_trip.departure_name = "foo"
        faulty_rot.trips[-1].arrival_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        error = "Rotation data differs from trips data"

        # check arrival data in rotation
        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_rot.trips[-1].arrival_name = "foo"
        faulty_rot.trips[0].departure_name = "foo"
        faulty_rot.arrival_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        arrival_trip = faulty_rot.trips[-1]
        faulty_rot.arrival_time = arrival_trip.arrival_time - timedelta(minutes=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        # check departure data in rotation
        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        faulty_rot.trips[-1].arrival_name = "foo"
        faulty_rot.trips[0].departure_name = "foo"
        faulty_rot.departure_name = "bar"
        assert schedule.Schedule.check_consistency(sched)["1"] == error

        sched = generate_basic_schedule()
        faulty_rot = list(sched.rotations.values())[0]
        departure_trip = faulty_rot.trips[0]
        faulty_rot.departure_time = departure_trip.departure_time - timedelta(minutes=1)
        assert schedule.Schedule.check_consistency(sched)["1"] == error
