import pytest
import pathlib
from tests.test_schedule import TestSchedule
from datetime import datetime

test_root = pathlib.Path(__file__).parent
file_root = test_root / "test_input_files"


class TestConsumption:
    """Class to test Consumption functionality"""
    consumption_path = file_root / "testing_energy_consumption.csv"

    def test_calculate_consumption(self):
        """ Various tests to trigger errors and check if behaviour is as expected"""
        schedule, scenario = TestSchedule().test_basic_run()
        trip = next(iter(schedule.rotations.values())).trips.pop(0)
        consumption = trip.__class__.consumption
        consumption.temperatures_by_hour = {hour: hour * 2 - 15 for hour in range(0, 24)}
        time = datetime(year=2023, month=1, day=1, hour=1)
        dist = 10
        vehicle = next(iter(consumption.vehicle_types.items()))
        vehicle_type = vehicle[0]
        charging_type = next(iter(vehicle[1].keys()))

        # check distance scaling
        def calc_c(distance):
            return consumption.calculate_consumption(time, distance, vehicle_type, charging_type,
                                                     temp=10, height_diff=0, level_of_loading=0,
                                                     mean_speed=18)[0]

        assert calc_c(dist) * 2 == calc_c(dist * 2)
        assert calc_c(dist) / 2 == calc_c(dist / 2)

        vehicle[1][charging_type]["mileage"] = str(self.consumption_path)

        # consumption is based on this formula
        def true_cons(lol, incline, speed, t_amb):
            return lol + incline + speed / 10 + abs(t_amb - 20) / 10

        consumption.vehicle_types[vehicle_type][charging_type] = vehicle[1][charging_type]

        lol = 0.5
        incline = 0
        speed = 15
        t_amb = 20
        distance = 1000  # 1000m =1km, since true_cons give the consumption per 1000 m

        # Check various inputs, which need interpolation
        assert true_cons(lol, incline, speed, t_amb) == \
               consumption.calculate_consumption(time, distance, vehicle_type, charging_type,
                                                 temp=t_amb, height_diff=incline * distance,
                                                 level_of_loading=lol,
                                                 mean_speed=speed)[0]

        incline = 0.05
        assert true_cons(lol, incline, speed, t_amb) == \
               consumption.calculate_consumption(time, distance, vehicle_type, charging_type,
                                                 temp=t_amb, height_diff=incline * distance,
                                                 level_of_loading=lol,
                                                 mean_speed=speed)[0]
        t_amb = 15
        assert true_cons(lol, incline, speed, t_amb) == \
               consumption.calculate_consumption(time, distance, vehicle_type, charging_type,
                                                 temp=t_amb, height_diff=incline * distance,
                                                 level_of_loading=lol,
                                                 mean_speed=speed)[0]
        lol = 0.1
        assert true_cons(lol, incline, speed, t_amb) == \
               consumption.calculate_consumption(time, distance, vehicle_type, charging_type,
                                                 temp=t_amb, height_diff=incline * distance,
                                                 level_of_loading=lol,
                                                 mean_speed=speed)[0]

        # check for out of bounds consumption. Max consumption in the table is 6.6.
        t_amb = -99999
        incline = 99999
        lol = 99999
        speed = 99999
        assert consumption.calculate_consumption(time, distance, vehicle_type, charging_type,
                                                 temp=t_amb, height_diff=incline * distance,
                                                 level_of_loading=lol,
                                                 mean_speed=speed)[0] < 6.7

        # check temperature default from temperature timeseries when temp is None
        _ = consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                              temp=None, height_diff=0, level_of_loading=0,
                                              mean_speed=18)[0]

        # check temperature default from temperature time series error throwing
        last_hour = 12
        consumption.temperatures_by_hour = {hour: hour * 2 - 15 for hour in range(0, last_hour)}
        time = datetime(year=2023, month=1, day=1, hour=last_hour + 2)
        with pytest.raises(KeyError):
            _ = consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                                  temp=None, height_diff=0, level_of_loading=0,
                                                  mean_speed=18)[0]
        del consumption.temperatures_by_hour
        with pytest.raises(AttributeError):
            _ = consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                                  temp=None, height_diff=0, level_of_loading=0,
                                                  mean_speed=18)[0]

        # reset temperature_by_hour
        consumption.temperatures_by_hour = {hour: hour * 2 - 15 for hour in range(0, 24)}

        # check level_of_loading default from level_of_loading time series error throwing
        last_hour = 12
        consumption.lol_by_hour = {hour: hour * 2 - 15 for hour in range(0, last_hour)}
        time = datetime(year=2023, month=1, day=1, hour=last_hour + 2)
        with pytest.raises(KeyError):
            _ = consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                                  temp=20, height_diff=0, level_of_loading=None,
                                                  mean_speed=18)[0]
        del consumption.lol_by_hour
        with pytest.raises(AttributeError):
            _ = consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                                  temp=20, height_diff=0, level_of_loading=None,
                                                  mean_speed=18)[0]