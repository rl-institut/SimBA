import csv

from rotation import Rotation


class Schedule:

    def __init__(self) -> None:
        """Constructs Schedule object from CSV file containing all trips of schedule"""
        self.rotations = {}

    @classmethod
    def from_csv(cls, path_to_csv):
        """Constructs Schedule object from CSV file containing all trips of schedule.

        :param path_to_csv: Path to csv file containing trip data
        :type path_to_csv: str
        """
        schedule = cls()

        with open(path_to_csv, 'r') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                rotation_id = trip['rotation_id']
                if rotation_id not in schedule.rotations.keys():
                    schedule.rotations.update({rotation_id: Rotation(rotation_id)})

                schedule.rotations[rotation_id].add_trip(trip)

    def remove_rotation_by_id(self, id):
        del self.rotations[id]

    def filter_rotations(self, filter_definition):
        """Based on a given filter definition (tbd), rotations will be dropped from schedule."""
        pass

    def add_charging_types(self):
        """Iterate across all rotations/trips and append charging type if not given"""
        pass

    def assign_vehicles():
        """ Assign vehicle IDs to rotations.
            Just randomly? Match consumption with battery sizes?
        """
        pass

    def calculate_consumption(self, consumption_func):
        for rot in self.rotations.values():
            for trip in rot.trips:
                trip.consumption(consumption_func)
