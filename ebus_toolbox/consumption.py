import csv
import pandas as pd
from ebus_toolbox import util


class Consumption:
    def __init__(self, vehicle_types, **kwargs) -> None:
        # load temperature of the day, now dummy winter day
        self.temperatures_by_hour = {}

        temperature_file_path = kwargs.get("outside_temperatures", None)
        # parsing the Temperature to a dict
        if temperature_file_path is not None:
            with open(temperature_file_path, encoding='utf-8') as f:
                delim = util.get_csv_delim(temperature_file_path)
                reader = csv.DictReader(f, delimiter=delim)
                for row in reader:
                    self.temperatures_by_hour.update({int(row['hour']): float(row['temperature'])})

        lol_file_path = kwargs.get("level_of_loading_over_day", None)
        # parsing the level of loading to a dict
        if lol_file_path is not None:
            with open(lol_file_path, encoding='utf-8') as f:
                delim = util.get_csv_delim(lol_file_path)
                reader = csv.DictReader(f, delimiter=delim)
                self.lol_by_hour = {}
                for row in reader:
                    self.lol_by_hour.update({int(row['hour']): float(row['level_of_loading'])})

        self.consumption_files = {}
        self.vehicle_types = vehicle_types

    def calculate_consumption(self, time, distance, vehicle_type, charging_type, temp=None,
                              height_diff=0, level_of_loading=None, mean_speed=18):
        """ Calculates consumed amount of energy for a given distance.
        :param time: The date and time at which the trip ends
        :type time: datetime.datetime
        :param distance: Distance travelled [m]
        :type distance: float
        :param vehicle_type: The vehicle type for which to calculate consumption
        :type vehicle_type: str
        :param charging_type: Charging type for the trip. Consumption differs between
                              distinct types.
        :type charging_type: str
        :param temp: Temperature outside of the bus in °Celsius
        :type temp: float
        :param height_diff: difference in height between stations in meters-
        :type height_diff: float
        :param level_of_loading: Level of loading of the bus between empty (=0) and
                                 completely full (=1.0). If None
                                 is provided, Level of loading will be interpolated from time series
        :type level_of_loading: float
        :param mean_speed: Mean speed between two stops in km/h
        :type mean_speed: float
        :return: Consumed energy [kWh] and delta SOC as tuple
        :rtype: (float, float)

        :raises IndexError: if there is missing data for temperature or lol data
        :raises AttributeError: if there is no path to temperature or lol data provided
        """

        assert self.vehicle_types.get(vehicle_type, {}).get(charging_type),\
            f"Combination of vehicle type {vehicle_type} and {charging_type} not defined."

        vehicle_info = self.vehicle_types[vehicle_type][charging_type]

        # in case a constant mileage is provided
        if isinstance(vehicle_info['mileage'], (int, float)):
            consumed_energy = vehicle_info['mileage'] * distance / 1000
            delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])
            return consumed_energy, delta_soc

        # if no specific Temperature is given, lookup temperature
        if temp is None:
            try:
                temp = self.temperatures_by_hour[time.hour]
            except AttributeError:
                print("Neither of these conditions is met:\n"
                      "1. Temperature data is available for every trip through the trips file "
                      "or a temperature over day file.\n"
                      f"2. A constant mileage for the vehicle "
                      f"{vehicle_info['mileage']} is provided.")
                raise AttributeError
            except IndexError:
                print(f"No temperature data for the hour {time.hour} is provided")
                raise IndexError

        # if no specific LoL is given, lookup temperature
        if level_of_loading is None:
            try:
                level_of_loading = self.lol_by_hour[time.hour]
            except AttributeError:
                print("Neither of these conditions is met:\n"
                      "1. Level of loading data is available for every trip through the trips file "
                      "or a level of loading over day file.\n"
                      f"2. A constant mileage for the vehicle "
                      f"{vehicle_info['mileage']} is provided.")
                raise AttributeError
            except IndexError:
                print(f"No level of loading data for the hour {time.hour} is provided")
                raise IndexError

        # load consumption csv
        consumption_path = vehicle_info["mileage"]

        # consumption_files holds interpol functions of csv files which are called directly

        # try to use the interpol function. If it does not exist yet its created in except case.
        consumption_function = vehicle_type+"_from_"+consumption_path
        try:
            mileage = self.consumption_files[consumption_function](
                                                               this_incline=height_diff / distance,
                                                               this_temp=temp,
                                                               this_lol=level_of_loading,
                                                               this_speed=mean_speed)
        except KeyError:
            # creating the interpol function from csv file.
            delim = util.get_csv_delim(consumption_path)
            df = pd.read_csv(consumption_path, sep=delim)
            # create lookup table and make sure its in the same order as the input point
            # which will be the input for the nd lookup
            df = df[df["vehicle_type"] == vehicle_type]
            assert len(df) > 0, f"Vehicle type {vehicle_type} not found in {consumption_path}"
            inc_col = df["incline"]
            tmp_col = df["t_amb"]
            lol_col = df["level_of_loading"]
            speed_col = df["mean_speed_kmh"]
            cons_col = df["consumption_kwh_per_km"]
            data_table = list(zip(inc_col, tmp_col, lol_col, speed_col, cons_col))

            def interpol_function(this_incline, this_temp, this_lol, this_speed):
                input_point = (this_incline, this_temp, this_lol, this_speed)
                return util.nd_interp(input_point, data_table)

            self.consumption_files.update({consumption_function: interpol_function})

            mileage = self.consumption_files[consumption_function](
                                                               this_incline=height_diff / distance,
                                                               this_temp=temp,
                                                               this_lol=level_of_loading,
                                                               this_speed=mean_speed)

        consumed_energy = mileage * distance / 1000  # kWh
        delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])

        return consumed_energy, delta_soc
