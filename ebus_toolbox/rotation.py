from ebus_toolbox.trip import Trip


class Rotation:

    def __init__(self, id, vehicle_type, ) -> None:
        self.id = id
        self.trips = []

        self.vehicle_type = vehicle_type
        self.vehicle_id = None
        self.charging_type = 'depot'

        self.consumption = 0

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
            elif self.arrival_time < new_trip.arrival_time:
                # last trip of rotation (for now)
                self.arrival_time = new_trip.arrival_time
                self.arrival_name = new_trip.arrival_name

        self.trips.append(new_trip)

    def calculate_consumption(self):
        rotation_consumption = 0
        for trip in self.trips:
            rotation_consumption += trip.calculate_consumption()

        self.consumption = rotation_consumption

        return rotation_consumption
