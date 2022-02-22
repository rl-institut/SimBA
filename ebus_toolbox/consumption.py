import numpy as np
import csv
import json


class Consumption:

    def __init__(self) -> None:
        # load temperature of the day, now dummy winter day
        self.temperatures_by_hour = {}
        with open("./data/examples/default_temp_winter.csv") as f:
            reader = csv.DictReader(f)
            for row in reader:
                self.temperatures_by_hour.update({int(row['hour']): float(row['temperature'])})

        self.consumption_files = {}
        with open("./data/examples/vehicle_types.json") as f:
            self.vehicle_types = json.load(f)

    def calculate_consumption(self, time, distance, vehicle_type):
        """ Calculates consumed amount of energy for a given distance.

        :param time: The date and time at which the trip ends
        :type time: datetime.datetime
        :param distance: Distance travelled [km]
        :type distance: float
        :param vehicle_type: The vehicle type for which to calculate consumption
        :type vehicle_type: str
        :return: Consumed energy [kWh] and delta SOC as tuple
        :rtype: (float, float)
        """
        temp = np.interp(time.hour,
                         list(self.temperatures_by_hour.keys()),
                         list(self.temperatures_by_hour.values()))

        # load consumption csv
        consumption_file = self.vehicle_types[vehicle_type]["mileage"]
        try:
            consumption = self.consumption_files[consumption_file]
        except KeyError:
            consumption = {'temperature': [], 'consumption': []}
            with open(consumption_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    consumption['temperature'].append(float(row['Temp.']))
                    consumption['consumption'].append(float(row['Kat. B']))
            self.consumption_files.update({consumption_file: consumption})

        xp = self.consumption_files[consumption_file]['temperature']
        fp = self.consumption_files[consumption_file]['consumption']

        mileage = np.interp(temp, xp, fp)  # kWh / km ???
        consumed_energy = mileage * distance  # kWh
        delta_SOC = consumed_energy / self.vehicle_types[vehicle_type]["capacity"]

        return (consumed_energy, delta_SOC)
