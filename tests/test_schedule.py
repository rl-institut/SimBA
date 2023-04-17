import sys

import pytest
import spice_ev.scenario as scenario
from spice_ev.util import set_options_from_config

from tests.helpers import generate_basic_schedule
import pathlib
from ebus_toolbox import schedule, trip, consumption, util

test_root = pathlib.Path(__file__).parent
file_root = test_root / "test_input_files"
example_root = pathlib.Path(__file__).parent.parent / "data/examples"

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
        """Returns a schedule and scenario after running the Ebus-Toolbox.
        :return: schedule, scenario
        """
        path_to_trips = example_root / "trips_example.csv"
        # set the system variables to imitate the console call with the config argument.
        # first element has to be set to something or error is thrown
        sys.argv = ["foo", "--config", str(example_root / "ebus_toolbox.cfg")]
        args = util.get_args()
        args.config = example_root / "ebus_toolbox.cfg"
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
        """ Check if the schedule creation properly throws a SystemExit error in case of missing
        mandatory options"""

        # check if System exit is properly thrown when mandatory options are not given
        with pytest.raises(SystemExit):
            # schedule creation without mandatory args
            schedule.Schedule(self.vehicle_types, self.electrified_stations)

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
        assert type(scen) == scenario.Scenario

    def test_assign_vehicles(self):
        """ Test if assigning vehicles works as intended.

        Use a trips csv with two rotations ("1","2") a day apart. Ebus toolbox should assign
        the same vehicle to both of them. rotation "3", starts shortly after "2" and should be
        a new vehicle."""

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

    def test_scenario_with_feed_in(self):
        """ Check if running a example with an extended electrified stations file
         with feed in, external load and battery works and if a scenario object is returned"""

        path_to_trips = example_root / "trips_example.csv"
        sys.argv = ["foo", "--config", str(example_root / "ebus_toolbox.cfg")]
        args = util.get_args()
        args.config = example_root / "ebus_toolbox.cfg"
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
        assert type(scen) == scenario.Scenario

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
        assert type(generated_schedule) == schedule.Schedule
