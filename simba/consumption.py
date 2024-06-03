import warnings

import pandas as pd
from simba import util
from simba.data_container import DataContainer
from simba.ids import INCLINE, LEVEL_OF_LOADING, SPEED, T_AMB, CONSUMPTION, VEHICLE_TYPE


class Consumption:
    def __init__(self) -> None:
        self.consumption_interpolation = {}

    def __call__(self, distance, vehicle_type, vehicle_info, temp,
                 height_difference, level_of_loading, mean_speed):
        """ Calculates consumed amount of energy for a given distance.

        :param distance: Distance travelled [m]
        :type distance: float
        :param vehicle_type: The vehicle type for which to calculate consumption
        :type vehicle_type: str
        :param vehicle_info: vehicle type information including mileage / consumption path
        :type vehicle_info: dict
        :param temp: Temperature outside of the bus in °Celsius
        :type temp: float
        :param height_difference: difference in height between stations in meters-
        :type height_difference: float
        :param level_of_loading: Level of loading of the bus between empty (=0) and
                                 completely full (=1.0). If None
                                 is provided, Level of loading will be interpolated from time series
        :type level_of_loading: float
        :param mean_speed: Mean speed between two stops in km/h
        :type mean_speed: float
        :return: Consumed energy [kWh] and delta SOC as tuple
        :rtype: (float, float)
        """

        # in case a constant mileage is provided
        if isinstance(vehicle_info['mileage'], (int, float)):
            consumed_energy = vehicle_info['mileage'] * distance / 1000
            delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])
            return consumed_energy, delta_soc

        # load consumption csv
        consumption_path = str(vehicle_info["mileage"])

        # consumption_interpolation holds interpol functions of csv files which are called directly
        # try to use the interpol function. If it does not exist yet its created in except case.
        consumption_lookup_name = self.get_consumption_lookup_name(consumption_path, vehicle_type)

        # This lookup includes the vehicle type. If the consumption data did not include vehicle
        # types this key will not be found. In this case use the file path instead
        try:
            interpol_function = self.consumption_interpolation[consumption_lookup_name]
        except KeyError:
            interpol_function = self.consumption_interpolation[consumption_path]

        mileage = interpol_function(
            this_incline=height_difference / distance, this_temp=temp,
            this_lol=level_of_loading, this_speed=mean_speed)

        consumed_energy = mileage * distance / 1000  # kWh
        delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])

        return consumed_energy, delta_soc

    @staticmethod
    def create_from_data_container(data_container: DataContainer):
        """Build a consumption instance from the stored data in the data container.
        :returns: Consumption instance
        :param data_container: container with consumption and vehicle type data
        :type data_container: DataContainer
        """
        # setup consumption calculator that can be accessed by all trips
        consumption = Consumption()
        for name, df in data_container.consumption_data.items():
            consumption.set_consumption_interpolation(name, df)
        return consumption

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
                if vt_specific_name in self.consumption_interpolation:
                    warnings.warn("Overwriting exising consumption function")
                self.consumption_interpolation.update({vt_specific_name: interpol_function})
            return

        interpol_function = self.get_nd_interpolation(df)
        if consumption_lookup_name in self.consumption_interpolation:
            warnings.warn("Overwriting exising consumption function")
        self.consumption_interpolation.update({consumption_lookup_name: interpol_function})

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

    def get_consumption_lookup_name(self, consumption_path, vehicle_type):
        """ Get name for the consumption lookup.

        :param consumption_path: Path to consumption data.
        :type consumption_path: str
        :param vehicle_type: Type of vehicle.
        :type vehicle_type: str
        :return: Name for the consumption lookup.
        :rtype: str
        """

        consumption_function = vehicle_type + "_from_" + consumption_path
        return consumption_function
