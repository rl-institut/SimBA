import csv
import bisect
from datetime import timedelta
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

    def assign_vehicles(self, minimum_standing_time):
        """ Assign vehicle IDs to rotations. A FIFO approach is used.
        For every rotation it is checked whether vehicles with matching type are idle, in which
        case the one with longest standing time since last rotation is used.
        If no vehicle is available a new vehicle ID is generated.

        :param minimum_standing_time: Amount of hours after arrival from previous rotation after
                                      which a vehicle become avaibable for dispatch again.
        :type minimum_standing_time: int
        """
        rotations_in_progress = []
        idle_vehicles = []
        vehicle_type_counts = {vehicle_type: 0 for vehicle_type in self.vehicle_types.keys()}

        rotations = sorted(self.rotations.values(), key=lambda rot: rot.departure_time)

        for rot in rotations:
            # find vehicles that have completed rotation and stood for a minimum staning time
            # mark those vehicle as idle
            for r in rotations_in_progress:
                if rot.departure_time > r.arrival_time + timedelta(hours=minimum_standing_time):
                    idle_vehicles.append(r.vehicle_id)
                    rotations_in_progress.pop(0)
                else:
                    break

            # find idle vehicle for rotation if exists
            # else generate new vehicle id
            vt_ct = f"{rot.vehicle_type}_{rot.charging_type}"
            id = next((id for id in idle_vehicles if vt_ct in id), None)
            if id is None:
                vehicle_type_counts[vt_ct] += 1
                id = f"{vt_ct}_{vehicle_type_counts[vt_ct]}"
            else:
                idle_vehicles.remove(id)

            rot.vehicle_id = id
            arrival_times = [r.arrival_time for r in rotations_in_progress]
            rotations_in_progress.insert(bisect.bisect(arrival_times, rot.arrival_time), rot)

    def calculate_consumption(self):
        self.consumption = 0
        for rot in self.rotations.values():
            self.consumption += rot.calculate_consumption()

        return self.consumption
