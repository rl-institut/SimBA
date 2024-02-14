import csv
import warnings

import pandas as pd
from simba import util
from simba.ids import INCLINE, LEVEL_OF_LOADING, SPEED, T_AMB, CONSUMPTION, VEHICLE_TYPE


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
        """

        assert self.vehicle_types.get(vehicle_type, {}).get(charging_type), (
            f"Combination of vehicle type {vehicle_type} and {charging_type} not defined.")

        vehicle_info = self.vehicle_types[vehicle_type][charging_type]

        # in case a constant mileage is provided
        if isinstance(vehicle_info['mileage'], (int, float)):
            consumed_energy = vehicle_info['mileage'] * distance / 1000
            delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])
            return consumed_energy, delta_soc

        # if no specific Temperature is given, lookup temperature
        if temp is None:
            temp = self.get_temperature(time, vehicle_info)

        # if no specific LoL is given, lookup LoL
        if level_of_loading is None:
            level_of_loading = self.get_level_of_loading(time, vehicle_info)

        # load consumption csv
        consumption_path = str(vehicle_info["mileage"])

        # consumption_files holds interpol functions of csv files which are called directly
        # try to use the interpol function. If it does not exist yet its created in except case.
        consumption_lookup_name = self.get_consumption_lookup_name(consumption_path, vehicle_type)

        # This lookup includes the vehicle type. If the consumption data did not include vehicle
        # types this key will not be found. In this case use the file path instead
        try:
            interpol_function = self.consumption_files[consumption_lookup_name]
        except KeyError:
            interpol_function = self.consumption_files[consumption_path]

        mileage = interpol_function(
            this_incline=height_diff / distance, this_temp=temp,
            this_lol=level_of_loading, this_speed=mean_speed)

        consumed_energy = mileage * distance / 1000  # kWh
        delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])

        return consumed_energy, delta_soc

    def set_consumption_interpolation(self, consumption_lookup_name: str, df: pd.DataFrame):
        """
         Set interpolation function for consumption lookup.

         If dataframes contain vehicle_types the name of the vehicle_type will be added and the df
         filtered.

         :param consumption_lookup_name: Name for the consumption lookup.
         :type consumption_lookup_name: str
         :param df: DataFrame containing consumption data.
         :type df: pd.DataFrame
         """

        if VEHICLE_TYPE in df.columns:
            unique_vts = df[VEHICLE_TYPE].unique()
            for vt in unique_vts:
                mask = df[VEHICLE_TYPE] == vt
                df_vt = df.loc[mask, [INCLINE, SPEED, LEVEL_OF_LOADING, T_AMB, CONSUMPTION]]
                interpol_function = self.get_nd_interpolation(df_vt)
                vt_specific_name = self.get_consumption_lookup_name(consumption_lookup_name, vt)
                if vt_specific_name in self.consumption_files:
                    warnings.warn("Overwriting exising consumption function")
                self.consumption_files.update({vt_specific_name: interpol_function})
            return

        interpol_function = self.get_nd_interpolation(df)
        if consumption_lookup_name in self.consumption_files:
            warnings.warn("Overwriting exising consumption function")
        self.consumption_files.update({consumption_lookup_name: interpol_function})

    def get_nd_interpolation(self, df):
        """
        Get n-dimensional interpolation function.

        :param df: DataFrame containing consumption data.
        :type df: pd.DataFrame
        :return: Interpolation function.
        :rtype: function
        """

        inc_col = df[INCLINE]
        tmp_col = df[T_AMB]
        lol_col = df[LEVEL_OF_LOADING]
        speed_col = df[SPEED]
        cons_col = df[CONSUMPTION]
        data_table = list(zip(inc_col, tmp_col, lol_col, speed_col, cons_col))

        def interpol_function(this_incline, this_temp, this_lol, this_speed):
            input_point = (this_incline, this_temp, this_lol, this_speed)
            return util.nd_interp(input_point, data_table)

        return interpol_function

    def get_temperature(self, time, vehicle_info):
        """
        Get temperature for the given time.

        :param time: time of the lookup.
        :type time: datetime.datetime
        :param vehicle_info: Information about the vehicle.
        :type vehicle_info: dict
        :return: Temperature.
        :rtype: float
        :raises AttributeError: if temperature data is not available.
        :raises KeyError: if temperature data is not available for the given hour.
        """

        try:
            return self.temperatures_by_hour[time.hour]
        except AttributeError as e:
            raise AttributeError(
                "Neither of these conditions is met:\n"
                "1. Temperature data is available for every trip through the trips file "
                "or a temperature over day file.\n"
                f"2. A constant mileage for the vehicle "
                f"{vehicle_info['mileage']} is provided."
            ) from e
        except KeyError as e:
            raise KeyError(f"No temperature data for the hour {time.hour} is provided") from e

    def get_level_of_loading(self, time, vehicle_info):
        """
         Get level of loading for the given time.

         :param time: time of the lookup.
         :type time: datetime.datetime
         :param vehicle_info: Information about the vehicle.
         :type vehicle_info: dict
         :return: Level of loading.
         :rtype: float
         :raises AttributeError: if level of loading data is not available.
         :raises KeyError: if level of loading data is not available for the given hour.
         """
        try:
            return self.lol_by_hour[time.hour]
        except AttributeError as e:
            raise AttributeError(
                "Neither of these conditions is met:\n"
                "1. Level of loading data is available for every trip through the trips file "
                "or a level of loading over day file.\n"
                f"2. A constant mileage for the vehicle "
                f"{vehicle_info['mileage']} is provided."
            ) from e
        except KeyError as e:
            raise KeyError(f"No level of loading for the hour {time.hour} is provided") from e

    def get_consumption_lookup_name(self, consumption_path, vehicle_type):
        """
        Get name for the consumption lookup.

        :param consumption_path: Path to consumption data.
        :type consumption_path: str
        :param vehicle_type: Type of vehicle.
        :type vehicle_type: str
        :return: Name for the consumption lookup.
        :rtype: str
        """

        consumption_function = vehicle_type + "_from_" + consumption_path
        return consumption_function
