from datetime import datetime


class Trip:

    def __init__(self, rotation, departure_time, departure_name,
                 arrival_time, arrival_name, distance, **kwargs):
        self.departure_name = departure_name
        self.departure_time = datetime.fromisoformat(departure_time)
        self.arrival_time = datetime.fromisoformat(arrival_time)
        self.arrival_name = arrival_name
        self.distance = float(distance)
        self.line = kwargs.get('line', None)

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
            self.consumption = \
                Trip.consumption.calculate_consumption(self.arrival_time,
                                                       self.distance,
                                                       self.rotation.vehicle_type,
                                                       self.rotation.charging_type)
        except AttributeError:
            print("""To calculate consumption, a consumption object needs to be constructed
                   and linked to Trip class.""")

        return self.consumption

    def get_delta_soc(self):
        """ Compute change in state of charge (SOC) for this trip."""
        if self.consumption is None:
            self.calculate_consumption()

        self.delta_soc = Trip.consumption.get_delta_soc(self.consumption,
                                                        self.rotation.vehicle_type,
                                                        self.rotation.charging_type)
