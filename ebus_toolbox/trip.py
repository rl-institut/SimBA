from datetime import datetime, timedelta


class Trip:
    def __init__(self, rotation, departure_time, departure_name,
                 arrival_time, arrival_name, distance, **kwargs):
        self.departure_name = departure_name
        self.departure_time = datetime.fromisoformat(departure_time)
        self.arrival_time = datetime.fromisoformat(arrival_time)
        self.arrival_name = arrival_name
        self.distance = float(distance)
        self.line = kwargs.get('line', None)
        self.temperature = kwargs.get('temperature', None)
        # In case of empty temperature column
        try:
            self.temperature = float(self.temperature)
            # In Case of None or empty string
        except ValueError:
            self.temperature = None

        height_diff = kwargs.get("height_difference", None)
        if height_diff is None:
            station_data = kwargs.get("station_data", dict())
            try:
                height_diff = station_data[self.arrival_name]["elevation"] \
                              - station_data[self.departure_name]["elevation"]
            except KeyError or TypeError:
                height_diff = 0
        self.height_diff = height_diff
        self.level_of_loading = kwargs.get('level_of_loading', None)
        try:
            # Clip level of loading to [0,1]
            self.level_of_loading = max(0, min(float(self.level_of_loading), 1))
        # In Case of None or empty string
        except ValueError:
            self.level_of_loading = None
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

    def calculate_consumption(self):
        """ Compute consumption for this trip.

        :return: Consumption of trip [kWh]
        :rtype: float
        """

        try:
            self.consumption, self.delta_soc = \
                Trip.consumption.calculate_consumption(self.arrival_time,
                                                       self.distance,
                                                       self.rotation.vehicle_type,
                                                       self.rotation.charging_type,
                                                       temp=self.temperature,
                                                       height_diff=self.height_diff,
                                                       level_of_loading=self.level_of_loading,
                                                       mean_speed=self.mean_speed)
        except AttributeError:
            print("""To calculate consumption, a consumption object needs to be constructed
                   and linked to Trip class.""")

        return self.consumption
