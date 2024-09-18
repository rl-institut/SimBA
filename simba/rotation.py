import datetime
import logging

from simba.trip import Trip


class Rotation:

    def __init__(self, id, vehicle_type, schedule) -> None:
        self.id = id
        self.trips = []
        self.schedule = schedule

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

        # Tracks if a warning has been logged if the rotation ends at a non-electrified station.
        # Further warnings will be turned off.
        self.logged_warning = False

    def add_trip(self, trip):
        """ Create a trip object and append to rotations trip set.

        :param trip: Information on trip to be added to rotation
        :type trip: dict
        :raises Exception: if charging type of trip and rotation differ
        """
        new_trip = Trip(self, **trip)

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

        # set charging type if given
        charging_type = trip.get('charging_type')
        self.trips.append(new_trip)
        if charging_type in ['depb', 'oppb']:
            assert self.schedule.vehicle_types.get(
                self.vehicle_type, {}).get(charging_type) is not None, (
                f"The required vehicle type {self.vehicle_type}({charging_type}) "
                "is not given in the vehicle_types.json file.")
            if self.charging_type is None:
                # set CT for whole rotation
                self.set_charging_type(charging_type)
            elif self.charging_type == charging_type:
                # same CT as other trips: just add trip consumption
                self.consumption += self.schedule.calculate_trip_consumption(new_trip)
            else:
                # different CT than rotation: error
                raise Exception(f"Two trips of rotation {self.id} have distinct charging types")

    def set_charging_type(self, ct):
        """ Change charging type of either all or specified rotations.

        This may also change the minimum standing time at depot after completion of rotation.

        :param ct: Choose this charging type if possible. Either 'depb' or 'oppb'.
        :type ct: str
        """
        assert ct in ["oppb", "depb"], f"Invalid charging type: {ct}"

        if ct == self.charging_type:
            return

        assert self.schedule.vehicle_types.get(self.vehicle_type, {}).get(ct), (
            f"Combination of vehicle type {self.vehicle_type} and {ct} not defined.")

        old_consumption = self.consumption
        self.charging_type = ct
        # consumption may have changed with new charging type
        self.consumption = self.schedule.calculate_rotation_consumption(self)

        # recalculate schedule consumption: update for new rotation consumption
        self.schedule.consumption += self.consumption - old_consumption

    @property
    def earliest_departure_next_rot(self):
        """Earliest possible departure of next rotation as datetime."""
        # noqa: DAR201
        return self.arrival_time + datetime.timedelta(hours=self.min_standing_time)

    @property
    def min_standing_time(self):
        """Minimum duration of standing time in minutes.

        No consideration of depot buffer time or charging curve.

        :return: Minimum duration of standing time in minutes.
        """
        # noqa: DAR201
        ct = self.charging_type
        assert ct in ["depb", "oppb"]

        min_recharge_soc = vars(self.schedule)[f"min_recharge_deps_{ct}"]
        stations = self.schedule.stations
        try:
            charge_power = stations[self.arrival_name].get(
                f"cs_power_deps_{ct}", vars(self.schedule)[f"cs_power_deps_{ct}"])
        except KeyError:
            # log a warning once for this. Since min_standing_time is called many times during
            # vehicle assignment, this would clutter the console / log.
            if not self.logged_warning:
                self.logged_warning = True
                logging.warning(f"Rotation {self.id} ends at a non-electrified station.")
                # min_standing_time set to zero, so if another rotation starts here,
                # the vehicle can always be used.
            return 0

        capacity = self.schedule.vehicle_types[self.vehicle_type][ct]["capacity"]

        # minimum time needed to recharge consumed power
        min_standing_time = (self.consumption / charge_power)
        # time to charge battery from 0 to desired SOC
        desired_max_standing_time = ((capacity / charge_power) * min_recharge_soc)
        min_standing_time = min(min_standing_time, desired_max_standing_time)
        return min_standing_time
