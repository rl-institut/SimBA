from ebus_toolbox.trip import Trip


class Rotation:

    def __init__(self, id) -> None:
        self.id = id
        self.trips = []
        self.consumption = 0

    def add_trip(self, trip):
        """Create a trip object and append to rotations trip set

        :param trip: Information on trip to be added to rotation
        :type trip: dict
        """

        self.trips.append(Trip(**trip))

    def calculate_consumption(self):
        rotation_consumption = 0
        for trip in self.trips:
            rotation_consumption += trip.calculate_consumption()

        self.consumption = rotation_consumption

        return rotation_consumption
