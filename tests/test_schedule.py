"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""
import pytest
import spice_ev.scenario as scenario
from tests.helpers import generate_basic_schedule
import json
import pathlib
from ebus_toolbox import schedule, trip, consumption, util

test_root = pathlib.Path(__file__).parent
file_root = test_root / "test_input_files"

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
            # schedule creation without mandatory args
            schedule.Schedule(self.vehicle_types, self.electrified_stations)
        except SystemExit:
            # exception was properly thrown
            pass
        else:
            # create assertion error since schedule creation did no
            assert 0, "Schedule initialization did not get mandatory options but did" \
                      "not throw an exception"

    def test_station_data_reading(self):
        """ Test if the reading of the geo station data works and outputs warnings in
        case the data was problematic, e.g. not numeric or not existent"""
        path_to_trips = file_root / "trips.csv"

        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_all_station_data = file_root / "all_stations_example.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args,
                                                        station_data_path=path_to_all_station_data)
        assert generated_schedule.station_data is not None

        # check if reading a non valid station.csv throws warnings
        with pytest.warns(Warning) as record:
            path_to_all_station_data = file_root / "not_existent_file"
            schedule.Schedule.from_csv(path_to_trips, self.vehicle_types, self.electrified_stations,
                                       **mandatory_args, station_data_path=path_to_all_station_data)
            assert len(record) == 1

            path_to_all_station_data = file_root / "not_numeric_stations.csv"
            schedule.Schedule.from_csv(path_to_trips, self.vehicle_types, self.electrified_stations,
                                       **mandatory_args, station_data_path=path_to_all_station_data)
            assert len(record) == 2

    def test_rotation_consistency(self):
        """ Test if the consistency check for rotations finds missing einsetz or ausetzfahrten,
        if it finds gaps between arrival and departure stations of successive trips and last
        it checks if the arguments for check_rotation_consistency and ignore_inconsistent_rotations
        were set correctly, i.e. if ignore_inconsistent_rotations is True, check_rotation_consistency
        has to be true as well. If not a warning has to be thrown"""

        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)
        output_dir = test_root / "test_output_files"
        rotation_consistency_dict = dict(check_rotation_consistency=True,
                                         ignore_inconsistent_rotations=True,
                                         output_directory=output_dir)

        # check if consistency check does not remove rotations of a "good" trips.csv
        path_to_trips = file_root / "trips.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args)

        generated_schedule2 = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                         self.electrified_stations,
                                                         **mandatory_args,
                                                         **rotation_consistency_dict)

        assert len(generated_schedule.rotations) == len(generated_schedule2.rotations)

        # check if consistency check does removes bad rotation of a "bad" trips.csv
        path_to_trips = file_root / "trips_no_einsetz.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args,
                                                        **rotation_consistency_dict)
        assert len(generated_schedule.rotations) == 0

        # check if consistency check does removes bad rotation of a "bad" trips.csv
        path_to_trips = file_root / "trips_no_aussetz.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args,
                                                        **rotation_consistency_dict)
        assert len(generated_schedule.rotations) == 0

        # check if consistency check does removes bad rotation of a "bad" trips.csv.
        # in this case a skipping of stations between trips
        path_to_trips = file_root / "inconsistent_stations_trips.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args,
                                                        **rotation_consistency_dict)
        assert len(generated_schedule.rotations) == 0

        # check if consistency throws warning for bad arguments
        with pytest.warns(Warning) as record:
            path_to_trips = file_root / "inconsistent_stations_trips.csv"
            generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                            self.electrified_stations,
                                                            **mandatory_args,
                                                            check_rotation_consistency=False,
                                                            ignore_inconsistent_rotations=True,
                                                            output_directory=output_dir)
            assert len(record) == 1

    def test_run(self):
        """ Check if running a basic example works and if a scenario object is returned"""
        path_to_trips = file_root / "trips.csv"
        parser = util.create_ArgumentParser_with_arguments()
        args = parser.parse_args(args="")
        args.config = file_root / "ebus_toolbox.cfg"
        args.days = None
        args.seed = 5

        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_all_station_data = file_root / "all_stations_example.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args,
                                                        station_data_path=path_to_all_station_data)

        util.set_options_from_config(args, check=False, verbose=False)
        args.ALLOW_NEGATIVE_SOC = True
        args.attach_vehicle_soc = True

        scen = generated_schedule.run(args)
        assert type(scen) == scenario.Scenario

        return generated_schedule, scen

    def test_assign_vehicles(self):
        """ Test if assigning vehicles works as intended.

        Use a trips csv with two rotations ("1","2") a day apart. Ebus toolbox should assign
        the same vehicle to both of them. rotation "3", starts shortly after "2" and should be
        a new vehicle."""

        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_trips = file_root / "trips_assign_vehicles.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args)
        generated_schedule.assign_vehicles()

        assert generated_schedule.rotations["1"].vehicle_id == \
               generated_schedule.rotations["2"].vehicle_id
        assert generated_schedule.rotations["1"].vehicle_id != \
               generated_schedule.rotations["3"].vehicle_id

    def test_calculate_consumption(self):
        """ Test if calling the consumption calculation works
        """
        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_trips = file_root / "trips_assign_vehicles.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args)
        # set mileage constant
        mileage = 10
        for v_typ in generated_schedule.vehicle_types.values():
            for charge_typ in v_typ.values():
                charge_typ['mileage'] = mileage
        calc_consumption = generated_schedule.calculate_consumption()
        # set mileage to half of the prev constant
        for v_typ in generated_schedule.vehicle_types.values():
            for charge_typ in v_typ.values():
                charge_typ['mileage'] = mileage / 2
        assert calc_consumption == generated_schedule.calculate_consumption() * 2


    def test_get_common_stations(self):
        """ Test if getting common_stations works. Rotation 1 is on the first day, rotation 2 and 3
        on the second day. rotation 1 should not share any stations with other rotations and
        2 and 3 are basically simulati
        """
        trip.Trip.consumption = consumption.Consumption(self.vehicle_types,
                                                        outside_temperatures=self.temperature_path,
                                                        level_of_loading_over_day=self.lol_path)

        path_to_trips = file_root / "trips_assign_vehicles.csv"
        generated_schedule = schedule.Schedule.from_csv(path_to_trips, self.vehicle_types,
                                                        self.electrified_stations, **mandatory_args)

        common_stations=generated_schedule.get_common_stations(only_opps=False)
        assert len(common_stations["1"]) == 0
        assert len(common_stations["2"]) == 1
        assert len(common_stations["3"]) == 1
        assert '3' in (common_stations["2"])
        assert '2' in (common_stations["3"])


    def test_get_negative_rotations(self):
        """Check if the single rotation '1' with a negative soc is found """

        # make use of the test_run() which has to return schedule and scenario object
        sched, scen = self.test_run()

        neg_rots= sched.get_negative_rotations(scen)
        assert '1' in neg_rots


    def test_set_charging_type(self):
        pass

    def test_schedule_from_csv(self):
        generated_schedule = generate_basic_schedule()
        assert len(generated_schedule.rotations) == 1
        assert type(generated_schedule) == schedule.Schedule


t = TestSchedule()
t.test_get_negative_rotations()
