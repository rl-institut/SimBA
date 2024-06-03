from tests.test_schedule import BasicSchedule
from tests.conftest import example_root
import pandas as pd


class TestConsumption:
    """Class to test Consumption functionality"""
    consumption_path = example_root / "energy_consumption_example.csv"

    def test_calculate_consumption(self, tmp_path):
        """Various tests to trigger errors and check if behaviour is as expected

        :param tmp_path: pytest fixture to create a temporary path
        """
        schedule, scenario, _ = BasicSchedule().basic_run()
        consumption_calculator = schedule.consumption_calculator
        dist = 10
        vehicle = next(iter(schedule.vehicle_types.items()))
        vehicle_type = vehicle[0]
        charging_type = next(iter(vehicle[1].keys()))
        vehicle_info = schedule.vehicle_types[vehicle_type][charging_type]

        # check distance scaling
        def calc_c(distance):
            return consumption_calculator(
                distance, vehicle_type, vehicle_info, temp=10, height_difference=0,
                level_of_loading=0, mean_speed=18)[0]

        assert calc_c(dist) * 2 == calc_c(dist * 2)
        assert calc_c(dist) / 2 == calc_c(dist / 2)

        # create a custom and well defined consumption file based on this formula

        def true_cons(lol, incline, speed, t_amb):
            return lol + incline + speed / 10 + (t_amb - 20) / 10

        # apply the formula on the consumption dataframe
        consumption_df = pd.read_csv(self.consumption_path)
        consumption_col = consumption_df["consumption_kwh_per_km"].view()
        lol_col = consumption_df["level_of_loading"]
        incline_col = consumption_df["incline"]
        speed_col = consumption_df["mean_speed_kmh"]
        temp_col = consumption_df["t_amb"]

        consumption_col[:] = true_cons(lol_col, incline_col, speed_col, temp_col)

        # save the file in a temp folder and use from now on
        vehicle[1][charging_type]["mileage"] = "new_consumption"
        consumption_calculator.set_consumption_interpolation(vehicle[1][charging_type]["mileage"],
                                                             consumption_df)

        lol = 0.5
        incline = 0
        speed = 15
        t_amb = 20
        distance = 1000  # 1000m =1km, since true_cons give the consumption per 1000 m

        # Check various inputs, which need interpolation. Inputs have to be inside the data, i.e.
        # not out of bounds
        assert true_cons(lol, incline, speed, t_amb) == consumption_calculator(
            distance, vehicle_type, vehicle_info, temp=t_amb,
            height_difference=incline * distance, level_of_loading=lol, mean_speed=speed)[0]

        incline = 0.02
        assert true_cons(lol, incline, speed, t_amb) == consumption_calculator(
            distance, vehicle_type, vehicle_info, temp=t_amb,
            height_difference=incline * distance, level_of_loading=lol, mean_speed=speed)[0]

        t_amb = 15
        assert true_cons(lol, incline, speed, t_amb) == consumption_calculator(
            distance, vehicle_type, vehicle_info, temp=t_amb,
            height_difference=incline * distance, level_of_loading=lol, mean_speed=speed)[0]

        lol = 0.1
        assert true_cons(lol, incline, speed, t_amb) == consumption_calculator(
            distance, vehicle_type, vehicle_info, temp=t_amb,
            height_difference=incline * distance, level_of_loading=lol, mean_speed=speed)[0]

        # check for out of bounds consumption. Max consumption in the table is 6.6.
        t_amb = 99999
        incline = 99999
        lol = 99999
        speed = 99999
        assert consumption_calculator(
            distance, vehicle_type, vehicle_info, temp=t_amb,
            height_difference=incline * distance, level_of_loading=lol, mean_speed=speed)[0] < 6.7
