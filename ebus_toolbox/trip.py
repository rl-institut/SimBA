import datetime
import random

from ebus_toolbox.consumption import Consumption


class Trip:
    consumption = Consumption()

    def __init__(self, rotation, departure_time, departure_name,
                 arrival_time, arrival_name, distance, **kwargs):
        self.departure_name = departure_name
        self.departure_time = departure_time
        self.arrival_time = datetime.datetime(2021, 8, 3, random.randint(0, 23), 30)  # testing
        self.arrival_name = arrival_name
        self.distance = float(distance)

        self.rotation = rotation

        self.consumption = 0  # kWh
        self.delta_SOC = 0

    def calculate_consumption(self):
        self.consumption, self.delta_SOC = \
            Trip.consumption.calculate_consumption(self.arrival_time,
                                                   self.distance,
                                                   self.rotation.vehicle_type,
                                                   self.rotation.charging_type)

        return self.consumption
