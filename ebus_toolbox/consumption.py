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

        with open(kwargs.get("level_of_loading_over_day",
                             "data/examples/default_level_of_loading_over_day.csv")) as f:
            reader = csv.DictReader(f)
            self.lol_by_hour = {}
            for row in reader:
                self.lol_by_hour.update({int(row['hour']): float(row['level_of_loading'])})

        self.consumption_files = {}
        self.vehicle_types = vehicle_types

    def calculate_consumption(self, time, distance, vehicle_type, charging_type,
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

        # the charging type may not be set
        # picking one of the available charging types as only vehicle's mileage is
        # needed which does not depend on charging type
        if charging_type is None:
            vt_ct = next((t for t in self.vehicle_types.keys() if vehicle_type in t))
        else:
            vt_ct = f"{vehicle_type}_{charging_type}"

        # in case a constant mileage is provided
        if isinstance(self.vehicle_types[vt_ct]['mileage'], (int, float)):

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
        consumption_func = self.vehicle_types[vt_ct]["mileage"]

        # Consumption_files holds interpol functions of csv files which are called directly
        vehicle_type_nr = dict(SB=0, VDL=0, AB=1, CKB=1)[vehicle_type]

        # Try to use the interpol function. If it does not exist yet its created in except case.
        try:
            mileage = self.consumption_files[consumption_func](this_vehicle_type=vehicle_type_nr,
                                                               this_incline=height_diff / distance,
                                                               this_temp=temp,
                                                               this_lol=level_of_loading,
                                                               this_speed=mean_speed)
        except KeyError:
            # Creating the interpol function from csv file.
            # Find the unique values of the grid for each dimension, e.g. unique temperatures
            # for which consumption values are present.
            df = pd.read_csv(consumption_func, sep=",")
            vehicle_type_u = df["vehicle_type"].unique()
            incline_u = df["incline"].unique()
            temp_u = df["temp"].unique()
            lol_u = df["level_of_loading"].unique()
            speed_u = df["mean_speed_kmh"].unique()
            points = (vehicle_type_u, incline_u, temp_u, lol_u, speed_u)
            values = np.zeros((len(vehicle_type_u), len(incline_u),
                               len(temp_u), len(lol_u), len(speed_u)))

            # Nested loops for creating values in grid
            # Each given dependence of consumption, e.g. temperature, speed, level of loading,
            # creates a dimension of the value matrix. Eg. 2*Temperatures,4*speeds and
            # 3* levels of loading  create a 3d matrix, with the shape (2,4,3). While one matrix
            # holds all the input values (see "points") another one has to hold all the output
            # consumption values. This matrix is filled with this nested for loop (loops in loops)
            # Later on the scypy interpn function can get created with the "points" from above and
            # the here created value grid.
            for i, vt in enumerate(vehicle_type_u):
                mask_vt = df["vehicle_type"] == vt
                for ii, inc in enumerate(incline_u):
                    mask_inc = df["incline"] == inc
                    for iii, temp in enumerate(temp_u):
                        mask_temp = df["temp"] == temp
                        for iv, lol in enumerate(lol_u):
                            mask_lol = df["level_of_loading"] == lol
                            for v, speed in enumerate(speed_u):
                                values[i, ii, iii, iv, v] = \
                                    df[mask_vt * mask_inc * mask_temp * mask_lol *
                                       (df["mean_speed_kmh"] == speed)]["consumption_kwh_per_km"]

            def interpol_function(this_vehicle_type, this_incline, this_temp, this_lol, this_speed):
                point = (this_vehicle_type, this_incline, this_temp, this_lol, this_speed)
                return interpn(points, values, point)[0]

            self.consumption_files.update({consumption_func: interpol_function})

            mileage = self.consumption_files[consumption_func](this_vehicle_type=vehicle_type_nr,
                                                               this_incline=height_diff / distance,
                                                               this_temp=temp,
                                                               this_lol=level_of_loading,
                                                               this_speed=mean_speed)

        consumed_energy = mileage * distance / 1000  # kWh
        delta_soc = -1 * (consumed_energy / self.vehicle_types[vt_ct]["capacity"])

        return consumed_energy, delta_soc
