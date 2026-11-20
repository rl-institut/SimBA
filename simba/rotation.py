import datetime
import logging

import numpy as np


class Rotation:

    def __init__(self, id, vehicle_type, schedule) -> None:
        self.id = id
        self.trips = []
        self.schedule = schedule

        self.allow_opp_charging_for_oppb: bool = True
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

        :param trip: trip to be added to rotation
        :type trip: Trip
        """
        self.distance += trip.distance
        if trip.line:
            self.lines.add(trip.line)

        if self.departure_time is None and self.arrival_time is None:
            # first trip added
            self.departure_time = trip.departure_time
            self.departure_name = trip.departure_name
            self.arrival_time = trip.arrival_time
            self.arrival_name = trip.arrival_name
        else:
            if self.departure_time > trip.departure_time:
                # first of rotation found (for now)
                self.departure_time = trip.departure_time
                self.departure_name = trip.departure_name
            # '<=' instead of '<' since last trip of rotation has no duration
            # in the sample data. The trips are however chronologically
            # sorted which is why this approach works for sample data.
            # Will also work if one only relies on timestamps!
            elif self.arrival_time <= trip.arrival_time:
                # last trip of rotation (for now)
                self.arrival_time = trip.arrival_time
                self.arrival_name = trip.arrival_name

        # set charging type if given
        self.trips.append(trip)

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
        """Minimum duration of standing time in hours.

        Uses the numeric charging curves of the schedule (charging curve and efficiency, set by
        Schedule.assign_vehicles_w_min_recharge_soc). No consideration of depot buffer time.

        :return: Minimum duration of standing time in hours.
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

        recharge_curves = self.schedule.recharge_curves
        assert recharge_curves is not None, "Schedule.recharge_curves must be set beforehand"
        # numeric charging curve (time in minutes, SoC), ends at desired SoC
        curve = recharge_curves[self.vehicle_type][ct][charge_power]
        times, socs = curve[:, 0], curve[:, 1]
        # rotation started at desired SoC, recharge consumption (at most min_recharge_soc)
        start_soc = socs[-1] - self.consumption / capacity
        target_soc = min(start_soc + min_recharge_soc, socs[-1])
        # below SoC 0: charging power of SoC 0 (as in optimizer_util.get_delta_soc)
        gradient = (socs[1] - socs[0]) / (times[1] - times[0])
        start_min, target_min = [
            np.interp(soc, socs, times) if soc >= 0 else soc / gradient
            for soc in (start_soc, target_soc)]
        return max(target_min - start_min, 0) / 60
