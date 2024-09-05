from tests.test_schedule import BasicSchedule
from tests.conftest import example_root
from datetime import timedelta
from simba.schedule import get_idle_consumption
import pandas as pd


class TestConsumption:
    """Class to test Consumption functionality"""
    consumption_path = example_root / "energy_consumption_example.csv"

    def test_calculate_idle_consumption(self, tmp_path):
        """Various tests to trigger errors and check if behaviour is as expected

        :param tmp_path: pytest fixture to create a temporary path
        """
        schedule, scenario, args = BasicSchedule().basic_run()
        rotation = schedule.rotations["1"]
        first_trip = rotation.trips[0]
        second_trip = rotation.trips[1]

        # shift all trips except the first
        for t in schedule.rotations["1"].trips[1:]:
            t.arrival_time = t.arrival_time + timedelta(minutes=70)
            t.departure_time = t.departure_time + timedelta(minutes=70)

        vt, ct = schedule.rotations["1"].vehicle_type, rotation.charging_type
        v_info = schedule.vehicle_types[vt][ct]
        v_info["idle_consumption"] = 0

        # Make sure that there is a break duration.
        # By shifing all trips earlier we made sure that this "second_trip" stays the second trip
        second_trip.departure_time = first_trip.arrival_time + timedelta(minutes=60)
        idle_consumption, idle_delta_soc = get_idle_consumption(first_trip, second_trip, v_info)
        assert idle_consumption == 0
        v_info["idle_consumption"] = 1
        second_trip.departure_time = first_trip.arrival_time + timedelta(minutes=60)
        idle_consumption, idle_delta_soc = get_idle_consumption(first_trip, second_trip, v_info)
        assert idle_consumption == 1

        second_trip.departure_time = first_trip.arrival_time + timedelta(minutes=30)
        idle_consumption, idle_delta_soc = get_idle_consumption(first_trip, second_trip, v_info)
        assert idle_consumption == 0.5

        second_trip.departure_time = first_trip.arrival_time
        idle_consumption, idle_delta_soc = get_idle_consumption(first_trip, second_trip, v_info)
        assert idle_consumption == 0

        # make all rotations depb
        for r in schedule.rotations.values():
            r.set_charging_type("depb")
        # make all trips of all rotations consecutive
        last_trip = schedule.rotations["1"].trips[0]

        for r in schedule.rotations.values():
            r.departure_time = last_trip.arrival_time + timedelta(minutes=10)
            for t in r.trips:
                t.departure_time = last_trip.arrival_time + timedelta(minutes=10)
                t.arrival_time = t.departure_time + timedelta(minutes=1)
                last_trip = t
            r.arrival_time = t.arrival_time

        # Check that assignment of vehicles changes due to increased consumption. Only works
        # with adaptive_soc assignment
        for vt in schedule.vehicle_types.values():
            for ct in vt:
                vt[ct]["idle_consumption"] = 0
                vt[ct]["mileage"] = 0
        schedule.calculate_consumption()
        assert schedule.consumption == 0

        schedule.assign_vehicles_w_adaptive_soc(args)
        no_idle_consumption = schedule.vehicle_type_counts.copy()

        for vt in schedule.vehicle_types.values():
            for ct in vt:
                vt[ct]["idle_consumption"] = 9999
        schedule.calculate_consumption()
        schedule.assign_vehicles_w_adaptive_soc(args)
        # Without consumption, single vehicle can service all rotations.
        # With high idling, every rotation needs its own vehicle
        assert no_idle_consumption["AB_depb"] * 4 == schedule.vehicle_type_counts["AB_depb"]

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
