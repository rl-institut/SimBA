from trip import Trip


class Rotation:

    def __init__(self, id) -> None:
        self.id = id
        self.trips = []

    def add_trip(self, trip):
        """Create a trip object and append to rotations trip set

        :param trip: Information on trip to be added to rotation
        :type trip: dict
        """

        self.trips.append(Trip(**trip))
