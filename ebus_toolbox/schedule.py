import csv
from ebus_toolbox.rotation import Rotation


class Schedule:

    def __init__(self) -> None:
        """Constructs Schedule object from CSV file containing all trips of schedule"""
        self.rotations = {}
        self.consumption = 0

    @classmethod
    def from_csv(cls, path_to_csv):
        """Constructs Schedule object from CSV file containing all trips of schedule.

        :param path_to_csv: Path to csv file containing trip data
        :type path_to_csv: str
        :return: Returns a new instance of Schedule with all trips from csv loaded.
        :rtype: Schedule
        """
        schedule = cls()

        with open(path_to_csv, 'r') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                rotation_id = trip['rotation_id']
                if rotation_id not in schedule.rotations.keys():
                    schedule.rotations.update({
                        rotation_id: Rotation(id=rotation_id, vehicle_type=trip['vehicle_type'])})
                schedule.rotations[rotation_id].add_trip(trip)

        return schedule

    def remove_rotation_by_id(self, id):
        """ Removes a rotation from schedule

        :param id: Rotation ID to be removed
        :type id: int
        """
        del self.rotations[id]

    def filter_rotations(self):
        """Based on a given filter definition (tbd), rotations will be dropped from schedule."""
        pass

    def set_charging_type(self, ct, rotation_ids=None):
        """Iterate across all rotations/trips and append charging type if not given"""
        if rotation_ids is None:
            rotation_ids = self.rotations.keys()

        for id in rotation_ids:
            self.rotations[id].charging_type = ct

    def _set_vehicle_ids(self):
        pass

    def assign_vehicles(self, preferred_ct):
        """ Depending on preferred charging type and consumption for each rotation
            first assign a charging type to each rotation and then assign specific
            vehicle IDs to each rotation.
        """
        self.set_charging_type(ct=preferred_ct)
        self.calculate_consumption()

        if preferred_ct == 'depot':
            # evaluate for which busses rotation too long (consumed energy > capacity)
            # and change charging types to opp for those
            # self.calculate_consumption()
            pass

        self._set_vehicle_ids()

    def calculate_consumption(self):
        self.consumption = 0
        for rot in self.rotations.values():
            self.consumption += rot.calculate_consumption()

        return self.consumption
