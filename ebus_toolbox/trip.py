import datetime
import random

from consumption import Consumption


class Trip:
    consumption = Consumption()

    def __init__(self, departure_time, departure_name,
                 arrival_time, arrival_name, distance, **kwargs):
        self.departure_name = departure_name
        self.departure_time = departure_time
        self.arrival_time = datetime.datetime(2021, 8, 3, random.randint(0, 23), 30)
        self.arrival_name = arrival_name
        self.distance = float(distance)

        self.vehicle_type = kwargs.get('vehicle_type') + '_opp'  # dummy for testing purpose
        self.charging_type = None  # maybe make charging type a member?
        self.vehicle_id = kwargs.get('vehicle_id', None)

        self.consumption = 0
        self.delta_SOC = 0

    def calculate_consumption(self):
        self.consumption, self.delta_SOC = \
            Trip.consumption.calculate_consumption(self.arrival_time,
                                                   self.distance,
                                                   self.vehicle_type)

        return self.consumption
