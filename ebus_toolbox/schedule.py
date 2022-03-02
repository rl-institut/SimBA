import csv
from ebus_toolbox.rotation import Rotation


class Schedule:

    def __init__(self, vehicle_types) -> None:
        """Constructs Schedule object from CSV file containing all trips of schedule"""
        self.rotations = {}
        self.consumption = 0
        self.vehicle_types = vehicle_types

    @classmethod
    def from_csv(cls, path_to_csv, vehicle_types):
        """Constructs Schedule object from CSV file containing all trips of schedule.

        :param path_to_csv: Path to csv file containing trip data
        :type path_to_csv: str
        :return: Returns a new instance of Schedule with all trips from csv loaded.
        :rtype: Schedule
        """
        schedule = cls(vehicle_types)

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

    def set_charging_type(self, preferred_ct, rotation_ids=None):
        """Iterate across all rotations/trips and append charging type if not given"""
        assert preferred_ct in ["opp", "depot"], f"Invalid charging type: {preferred_ct}"
        if rotation_ids is None:
            rotation_ids = self.rotations.keys()

        for id in rotation_ids:
            rot = self.rotations[id]
            vehicle_type = self.vehicle_types[f"{rot.vehicle_type}_{rot.charging_type}"]
            if preferred_ct == "opp" or vehicle_type["capacity"] < rot.consumption:
                self.rotations[id].charging_type = "opp"
            else:
                self.rotations[id].charging_type = "depot"

    def assign_vehicles(self):
        """ Depending on preferred charging type and consumption for each rotation
            first assign a charging type to each rotation and then assign specific
            vehicle IDs to each rotation.
        """
        pass

    def calculate_consumption(self):
        self.consumption = 0
        for rot in self.rotations.values():
            self.consumption += rot.calculate_consumption()

        return self.consumption
