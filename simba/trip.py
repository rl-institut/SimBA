from datetime import timedelta


class Trip:
    def __init__(self, rotation, departure_time, departure_name,
                 arrival_time, arrival_name, distance, temperature, level_of_loading,
                 height_difference, **kwargs):
        self.departure_name = departure_name
        self.departure_time = departure_time
        self.arrival_time = arrival_time
        self.arrival_name = arrival_name
        self.distance = distance
        self.line = kwargs.get('line', None)
        self.temperature = temperature
        self.height_difference = height_difference
        self.level_of_loading = level_of_loading
        # mean speed in km/h from distance and travel time or from initialization
        # travel time is at least 1 min

        mean_speed = kwargs.get("mean_speed", (self.distance / 1000) /
                                max(1 / 60, ((self.arrival_time - self.departure_time) / timedelta(
                                    hours=1))))
        self.mean_speed = mean_speed

        # Attention: Circular reference!
        # While a rotation carries a references to this trip, this trip
        # itself carries the below reference to that rotation as well.
        # For Python versions < 3.4 this created a problem for the garbage
        # collector.
        self.rotation = rotation

        self.consumption = None  # kWh
        self.delta_soc = None
