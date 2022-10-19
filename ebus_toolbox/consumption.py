import numpy as np
import csv
import pandas as pd


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

        assert self.vehicle_types.get(vehicle_type, {}).get(charging_type), \
            f"Combination of vehicle type {vehicle_type} and {charging_type} not defined."

        vehicle_info = self.vehicle_types[vehicle_type][charging_type]

        # in case a constant mileage is provided
        if isinstance(vehicle_info['mileage'], (int, float)):
            consumed_energy = vehicle_info['mileage'] * distance / 1000
            delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])
            return consumed_energy, delta_soc

        # If no specific LoL is given, interpolate from demand time series.
        if temp is None:
            temp = np.interp(time.hour,
                             list(self.temperatures_by_hour.keys()),
                             list(self.temperatures_by_hour.values()))

        # If no specific LoL is given, interpolate from demand time series.
        if level_of_loading is None:
            level_of_loading = np.interp(time.hour,
                                         list(self.lol_by_hour.keys()),
                                         list(self.lol_by_hour.values()))

        # load consumption csv
        consumption_path = vehicle_info["mileage"]

        # Consumption_files holds interpol functions of csv files which are called directly
        vehicle_type_nr = dict(SB=0, VDL=0, AB=1, CKB=1)[vehicle_type]

        # Try to use the interpol function. If it does not exist yet its created in except case.
        try:
            mileage = self.consumption_files[consumption_path](this_vehicle_type=vehicle_type_nr,
                                                               this_incline=height_diff / distance,
                                                               this_temp=temp,
                                                               this_lol=level_of_loading,
                                                               this_speed=mean_speed)
        except KeyError:
            # Creating the interpol function from csv file.
            df = pd.read_csv(consumption_path, sep=",")
            # Create lookup table and make sure its in the same order as the input point
            # which will be the input for the nd lookup
            vt_col = df["vehicle_type"]
            inc_col = df["incline"]
            tmp_col = df["temp"]
            lol_col = df["level_of_loading"]
            speed_col = df["mean_speed_kmh"]
            cons_col = df["consumption_kwh_per_km"]
            data_table = list(zip(vt_col, inc_col, tmp_col, lol_col, speed_col, cons_col))

            def interpol_function(this_vehicle_type, this_incline, this_temp, this_lol, this_speed):
                input_point = (this_vehicle_type, this_incline, this_temp, this_lol, this_speed)
                return nd_interp(input_point, data_table)

            self.consumption_files.update({consumption_path: interpol_function})

            mileage = self.consumption_files[consumption_path](this_vehicle_type=vehicle_type_nr,
                                                               this_incline=height_diff / distance,
                                                               this_temp=temp,
                                                               this_lol=level_of_loading,
                                                               this_speed=mean_speed)

        consumed_energy = mileage * distance / 1000  # kWh
        delta_soc = -1 * (consumed_energy / vehicle_info["capacity"])

        return consumed_energy, delta_soc


def nd_interp(input_values, lookup_table):
    # find all unique values in table per column
    dim_sets = [set() for _ in input_values]
    for row in lookup_table:
        for i, v in enumerate(row[:-1]):
            dim_sets[i].add(v)
    dim_values = [sorted(s) for s in dim_sets]
    # find nearest value(s) per column
    # go through sorted column values until last less / first greater
    lower = [None] * len(input_values)
    upper = [None] * len(input_values)
    for i, v in enumerate(input_values):
        # initialize for out of bound values -> Constant value
        lower[i] = dim_values[i][0]
        upper[i] = dim_values[i][-1]
        for c in dim_values[i]:
            if v >= c:
                lower[i] = c
            if v <= c:
                upper[i] = c
                break
    # find rows in table made up of only lower or upper values
    points = []
    for row in lookup_table:
        for i, v in enumerate(row[:-1]):
            if lower[i] != v and upper[i] != v:
                break
        else:
            points.append(row)

    # interpolate between points that differ only in current dimension
    for i, x in enumerate(input_values):
        new_points = []
        # find points that differ in just that dimension
        for j, p1 in enumerate(points):
            for p2 in points[j + 1:]:
                for k in range(len(input_values)):
                    if p1[k] != p2[k] and i != k:
                        break
                else:
                    # differing row found
                    x1 = p1[i]
                    y1 = p1[-1]
                    x2 = p2[i]
                    y2 = p2[-1]
                    dx = x2 - x1
                    dy = y2 - y1
                    m = dy / dx
                    n = y1 - m * x1
                    y = m * x + n
                    # generate new point at interpolation
                    p = [v for v in p1]
                    p[i] = x
                    p[-1] = y
                    new_points.append(p)
                    # only couple
                    break
            else:
                # no matching row (singleton dimension?)
                new_points.append(p1)
        points = new_points

    return points[0][-1]
