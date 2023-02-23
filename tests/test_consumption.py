import pytest
import pathlib
from tests.test_schedule import TestSchedule
from datetime import datetime

test_root = pathlib.Path(__file__).parent
file_root = test_root / "test_input_files"


class TestConsumption:
    consumption_path = file_root / "testing_energy_consumption.csv"

    def test_calculate_consumption(self):
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
        # def true_cons(lol,incline,speed,t_amb):
        #     return lol+incline+speed/10+abs(t_amb-20)/10

        consumption.vehicle_types[vehicle_type][charging_type] = vehicle[1][charging_type]
        # check temperature default from temperature timeseries
        consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                          temp=None, height_diff=0, level_of_loading=0,
                                          mean_speed=18)[0]

        # check temperature default from temperature time series error throwing
        last_hour = 12
        consumption.temperatures_by_hour = {hour: hour * 2 - 15 for hour in range(0, last_hour)}
        time = datetime(year=2023, month=1, day=1, hour=last_hour+2)
        with pytest.raises(KeyError):
            consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                              temp=None, height_diff=0, level_of_loading=0,
                                              mean_speed=18)[0]
        del consumption.temperatures_by_hour
        with pytest.raises(AttributeError):
            consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                              temp=None, height_diff=0, level_of_loading=0,
                                              mean_speed=18)[0]

        # reset temperature_by_hour
        consumption.temperatures_by_hour = {hour: hour * 2 - 15 for hour in range(0, 24)}

        # check level_of_loading default from level_of_loading time series error throwing
        last_hour = 12
        consumption.lol_by_hour = {hour: hour * 2 - 15 for hour in range(0, last_hour)}
        time = datetime(year=2023, month=1, day=1, hour=last_hour+2)
        with pytest.raises(KeyError):
            consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                              temp=20, height_diff=0, level_of_loading=None,
                                              mean_speed=18)[0]
        del consumption.lol_by_hour
        with pytest.raises(AttributeError):
            consumption.calculate_consumption(time, dist, vehicle_type, charging_type,
                                              temp=20, height_diff=0, level_of_loading=None,
                                              mean_speed=18)[0]
