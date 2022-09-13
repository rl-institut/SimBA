import numpy as np
import csv
import pandas as pd
from scipy.interpolate import interpn


class Consumption:

    def __init__(self, vehicle_types, **kwargs) -> None:
        # load temperature of the day, now dummy winter day
        self.temperatures_by_hour = {}
        with open(kwargs.get("outside_temperatures", "data/examples/default_temp_winter.csv")) as f:
            reader = csv.DictReader(f)
            for row in reader:
                self.temperatures_by_hour.update({int(row['hour']): float(row['temperature'])})

        with open(kwargs.get("level_of_loading_over_day", "data/examples/default_level_of_loading_over_day.csv")) as f:
            reader = csv.DictReader(f)
            self.lol_by_hour = {}
            for row in reader:
                self.lol_by_hour.update({int(row['hour']): float(row['level_of_loading'])})

        self.consumption_files = {}
        self.vehicle_types = vehicle_types

    def calculate_consumption(self, time, distance, vehicle_type, charging_type,
                              height_difference=0, level_of_loading=None, mean_speed=18):
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
        :param height_difference: difference in height between stations in meters-
        :type height_difference: float
        :param level_of_loading: Level of loading of the bus between empty (=0) and completely full (=1.0). If None
        is provided, Level of loading will be interpolated from timeseries
        :type level_of_loading: float
        :param mean_speed: Mean speed between two stops in km/h
        :type mean_speed: float

        :return: Consumed energy [kWh]
        :rtype: float
        """

        # the charging type may not be set
        # picking one of the available charging types as only vehicle's mileage is
        # needed which does not depend on charging type
        if charging_type is None:
            vt_ct = next((t for t in self.vehicle_types.keys() if vehicle_type in t))
        else:
            vt_ct = f"{vehicle_type}_{charging_type}"

        # in case a constant mileage is provided
        if isinstance(self.vehicle_types[vt_ct]['mileage'], float):
            consumed_energy = self.vehicle_types[vt_ct]['mileage'] * distance / 1000
            delta_soc = -1 * (consumed_energy / self.vehicle_types[vt_ct]["capacity"])
            return consumed_energy, delta_soc

        temp = np.interp(time.hour,
                         list(self.temperatures_by_hour.keys()),
                         list(self.temperatures_by_hour.values()))

        # If no specific LoL is given, interpolate from demand time series.
        if level_of_loading is None:
            level_of_loading = np.interp(time.hour,
                                         list(self.lol_by_hour.keys()),
                                         list(self.lol_by_hour.values()))

        # load consumption csv
        consumption_file = self.vehicle_types[vt_ct]["mileage"]

        # Consumption_files holds interpol functions of csv files which are called directly and held in memory
        vehicle_type_nr = dict(SB=0, VDL=0, AB=1, CKB=1)[vehicle_type]
        try:
            mileage = self.consumption_files[consumption_file](vehicle_type=vehicle_type_nr,
                                                               incline=height_difference / distance,
                                                               temp=temp,
                                                               lol=level_of_loading,
                                                               speed=mean_speed)
            return mileage * distance
        except KeyError:
            # Creating the interpol function from csv file.
            df = pd.read_csv(consumption_file, sep=",")
            vehicle_type_u = list(df["vehicle_type"].unique())
            incline_u = list(df["incline"].unique())
            temp_u = list(df["temp"].unique())
            lol_u = list(df["level_of_loading"].unique())
            speed_u = list(df["mean_speed_kmh"].unique())
            points = (vehicle_type_u, incline_u, temp_u, lol_u, speed_u)
            values = np.zeros((len(vehicle_type_u), len(incline_u), len(temp_u), len(lol_u), len(speed_u)))

            # Nested loops for creating values in grid
            for i, vt in enumerate(vehicle_type_u):
                for ii, inc in enumerate(incline_u):
                    for iii, temp in enumerate(temp_u):
                        for iv, lol in enumerate(lol_u):
                            for v, speed in enumerate(speed_u):
                                values[i, ii, iii, iv, v] = df[(df["vehicle_type"] == vt) * (df["incline"] == inc) *
                                                               (df["temp"] == temp) * (df["level_of_loading"] == lol) *
                                                               (df["mean_speed_kmh"] == speed)][
                                    "consumption_kwh_per_km"]

            def interpol_function(this_vehicle_type, this_incline, this_temp, this_lol, this_speed):
                point = (this_vehicle_type, this_incline, this_temp, this_lol, this_speed)
                return interpn(points, values, point)[0]

            self.consumption_files.update({consumption_file: interpol_function})

        mileage = self.consumption_files[consumption_file](vehicle_type=vehicle_type_nr,
                                                           incline=height_difference / distance,
                                                           temp=temp,
                                                           lol=level_of_loading,
                                                           speed=mean_speed)
        consumed_energy = mileage * distance / 1000  # kWh

        delta_soc = -1 * (consumed_energy / self.vehicle_types[vt_ct]["capacity"])

        return consumed_energy, delta_soc