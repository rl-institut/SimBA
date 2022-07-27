import numpy as np
import csv


class Consumption:

    def __init__(self, vehicle_types, **kwargs) -> None:
        # load temperature of the day, now dummy winter day
        self.temperatures_by_hour = {}
        with open(kwargs.get("outside_temperatures", "data/examples/default_temp_winter.csv")) as f:
            reader = csv.DictReader(f)
            for row in reader:
                self.temperatures_by_hour.update({int(row['hour']): float(row['temperature'])})

        self.consumption_files = {}
        self.vehicle_types = vehicle_types

    def calculate_consumption(self, time, distance, vehicle_type, charging_type):
        """ Calculates consumed amount of energy for a given distance.

        :param time: The date and time at which the trip ends
        :type time: datetime.datetime
        :param distance: Distance travelled [km]
        :type distance: float
        :param vehicle_type: The vehicle type for which to calculate consumption
        :type vehicle_type: str
        :param charging_type: Charging type for the trip. Consumption differs between
                              distinct types.
        :type charging_type: str
        :return: Consumed energy [kWh] and delta SOC as tuple
        :rtype: (float, float)
        """

        vt_ct = f"{vehicle_type}_{charging_type}"

        # in case a constant mileage is provided
        if isinstance(self.vehicle_types[vt_ct]['mileage'], float):
            return self.vehicle_types[vt_ct]['mileage'] * distance / 1000

        temp = np.interp(time.hour,
                         list(self.temperatures_by_hour.keys()),
                         list(self.temperatures_by_hour.values()))

        # load consumption csv
        consumption_file = self.vehicle_types[vt_ct]["mileage"]
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

        mileage = np.interp(temp, xp, fp)  # kWh / m
        consumed_energy = mileage * distance / 1000  # kWh

        return consumed_energy

    def get_delta_soc(self, consumed_energy, vehicle_type, charging_type):
        vt_ct = f"{vehicle_type}_{charging_type}"
        return -1 * (consumed_energy / self.vehicle_types[vt_ct]["capacity"])
