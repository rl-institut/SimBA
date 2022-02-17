import consumption


class Trip:

    def __init__(self, departure_time, departure_name,
                 arrival_time, arrival_name, distance, **kwargs):
        self.departure_name = departure_name
        self.departure_time = departure_time
        self.arrival_time = arrival_time
        self.arrival_name = arrival_name
        self.distance = distance

        self.vehicle_id = kwargs.get('vehicle_id', None)

        self.consumption = getattr(consumption, kwargs['consumption_func'], 'naive')(distance,
                                                                                     departure_time,
                                                                                     arrival_time)

    def calculate_consumption(self, consumption_func=None):
        self.consumption = getattr(consumption, consumption_func, 'naive')(self.distance,
                                                                           self.departure_time,
                                                                           self.arrival_time)
