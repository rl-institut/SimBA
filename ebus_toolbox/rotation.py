from ebus_toolbox.trip import Trip


class Rotation:

    def __init__(self, id, vehicle_type) -> None:
        self.id = id
        self.trips = []

        self.vehicle_type = vehicle_type
        self.vehicle_id = None
        self.charging_type = None

        self.consumption = 0
        self.distance = 0
        self.lines = set()

        self.departure_time = None
        self.departure_name = None
        self.arrival_time = None
        self.arrival_name = None

    def add_trip(self, trip):
        """Create a trip object and append to rotations trip set

        :param trip: Information on trip to be added to rotation
        :type trip: dict
        """
        new_trip = Trip(self, **trip)

        # set charging type if given
        if ('charging_type' in trip and
                any(trip['charging_type'] == t for t in ['depb', 'oppb'])):
            self.charging_type = trip['charging_type']

        self.distance += new_trip.distance
        if new_trip.line:
            self.lines.add(new_trip.line)

        if self.departure_time is None and self.arrival_time is None:
            # first trip added
            self.departure_time = new_trip.departure_time
            self.departure_name = new_trip.departure_name
            self.arrival_time = new_trip.arrival_time
            self.arrival_name = new_trip.arrival_name
        else:
            if self.departure_time > new_trip.departure_time:
                # first of rotation found (for now)
                self.departure_time = new_trip.departure_time
                self.departure_name = new_trip.departure_name
            # '<=' instead of '<' since last trip of rotation has no duration
            # in the sample data. The trips are however chronologically
            # sorted which is why this approach works for sample data.
            # Will also work if one only relies on timestamps!
            elif self.arrival_time <= new_trip.arrival_time:
                # last trip of rotation (for now)
                self.arrival_time = new_trip.arrival_time
                self.arrival_name = new_trip.arrival_name

        self.trips.append(new_trip)

    def calculate_consumption(self):
        """ Calculate consumption of this rotation and all its trips.

        :return: Consumption of rotation [kWh]
        :rtype: float
        """
        rotation_consumption = 0
        for trip in self.trips:
            rotation_consumption += trip.calculate_consumption()

        self.consumption = rotation_consumption

        return rotation_consumption

    def delta_soc_all_trips(self):
        """ Compute change in state of charge (SOC) for every trip
            of this rotation. Stored in the trip objects.
        """
        for trip in self.trips:
            trip.get_delta_soc()
