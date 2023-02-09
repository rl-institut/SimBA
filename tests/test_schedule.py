"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""
import pytest

from tests.helpers import generate_basic_schedule
import json
import pathlib
from ebus_toolbox import schedule, trip, consumption

test_root = pathlib.Path(__file__).parent
file_root = test_root / "test_input_files"
path_to_trips = file_root / "trips.csv"

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
    temperature_path = file_root / 'default_temp_winter.csv'
    lol_path = file_root / 'default_level_of_loading_over_day.csv'

    with open(file_root / "electrified_stations.json", "r", encoding='utf-8') as file:
        electrified_stations = json.load(file)

    with open(file_root / "vehicle_types.json", "r", encoding='utf-8') as file:
        vehicle_types = json.load(file)

    def test_mandatory_options_exit(self):
        """ Check if the schedule creation properly throws a SystemExit error in case of missing
        mandatory options"""
        try:
            schedule.Schedule(self.vehicle_types, self.electrified_stations)
        except SystemExit:
            # exception was properly thrown
            pass
        else:
            raise Exception(
                "Schedule initialization did not get mandatory options but did not throw an"
                "exception")

    def test_station_data_reading(self):
        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_all_station_data = file_root / "all_stations_example.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                               self.electrified_stations, **mandatory_args,
                               station_data_pata=path_to_all_station_data)
        assert generated_schedule.stations

        path_to_all_station_data = file_root / "not_existent_file"
        with pytest.warns(Warning) as record:
            schedule.Schedule.from_csv(path_to_trips, self.vehicle_types, self.electrified_stations,
                                       **mandatory_args, station_data_path=path_to_all_station_data)
            assert len(record) == 1

            path_to_all_station_data = file_root / "not_numeric_stations.csv"
            schedule.Schedule.from_csv(path_to_trips, self.vehicle_types, self.electrified_stations,
                                       **mandatory_args, station_data_path=path_to_all_station_data)
            assert len(record) == 2
            assert 0


    def test_rotation_consistency(self):
        pass

    def test_run(self):
        pass

    def test_assign_vehicles(self):
        pass

    def test_calculate_consumption(self):
        pass

    def test_get_departure_of_first_trip(self):
        pass

    def test_get_common_stations(self):
        pass

    def test_generate_scenario(self):
        pass

    def test_get_negative_rotations(self):
        pass

    def test_set_charging_type(self):
        pass

    def test_schedule_from_csv(self):
        generated_schedule = generate_basic_schedule()
        assert len(generated_schedule.rotations) == 1
        assert type(generated_schedule) == schedule.Schedule

t =TestSchedule()
t.test_station_data_reading()