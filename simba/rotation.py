import datetime

from simba.trip import Trip


class Rotation:

    def __init__(self, id, schedule) -> None:
        self.id = id
        self.trips = []
        self.schedule = schedule

        self.vehicle_id = None
        self.vehicle_type = None
        self.charging_type = None

        self.consumption = 0
        self.distance = 0
        self.lines = set()

        self.departure_time = None
        self.departure_name = None
        self.arrival_time = None
        self.arrival_name = None

    def add_trip(self, trip):
        """ Create a trip object and append to rotations trip set.

        :param trip: Information on trip to be added to rotation
        :type trip: simba.trip.Trip
        :raises ValueError: if rotation has two trips of different charging type
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

        if self.vehicle_type is None:
            self.vehicle_type = trip.vehicle_type
            for prior_trip in self.trips:
                prior_trip.vehicle_type = trip.vehicle_type
        else:
            assert self.vehicle_type == trip.vehicle_type, (
                f"Two trips of rotation {self.id} have distinct vehicle types")

        if self.charging_type is not None and self.charging_type != trip.charging_type:
            raise ValueError(f"Two trips of rotation {self.id} have distinct charging types")

        # (re)calculate consumption
        # might change if charging_type changes
        # therefore, add new trip only after charging_type is set
        trip.rotation = self
        if trip.charging_type is None:
            trip.charging_type = self.charging_type
        if trip.charging_type in ["oppb", "depb"]:
            self.set_charging_type(trip.charging_type)
            self.consumption += trip.calculate_consumption()

        self.trips.append(trip)

    def add_trip_from_dict(self, trip_dict):
        """ Create a trip object from dictionary and append to rotations trip set.

        See :py:class:`~simba.trip.Trip` for arguments and defaults.

        :param trip_dict: vars of Trip
        :type trip_dict: dict
        """
        trip = Trip(
            departure_time=trip_dict["departure_time"],
            departure_name=trip_dict["departure_name"],
            arrival_time=trip_dict["arrival_time"],
            arrival_name=trip_dict["arrival_name"],
            distance=trip_dict["distance"],
            line=trip_dict.get("line"),
            temperature=trip_dict.get("temperature"),
            level_of_loading=trip_dict.get("level_of_loading"),
            mean_speed=trip_dict.get("mean_speed"),
            height_diff=trip_dict.get("height_diff"),
            station_data=trip_dict.get("station_data"),
            charging_type=trip_dict.get("charging_type"),
            vehicle_type=trip_dict.get("vehicle_type"),
        )
        self.add_trip(trip)

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

    def set_charging_type(self, ct):
        """ Change charging type of this rotation.

        This may also change the minimum standing time at depot after completion of rotation.

        :param ct: Choose this charging type if possible. Either 'depb' or 'oppb'.
        :type ct: str
        """
        if ct == self.charging_type:
            return

        assert ct in ["oppb", "depb"], f"Invalid charging type '{ct}' for rotation {self.id}"

        assert self.schedule.vehicle_types.get(self.vehicle_type, {}).get(ct), (
            f"Combination of vehicle type {self.vehicle_type} and {ct} not defined.")

        # update associated trips
        for trip in self.trips:
            trip.charging_type = ct

        old_consumption = self.consumption
        self.charging_type = ct
        # consumption may have changed with new charging type
        self.consumption = self.calculate_consumption()

        # recalculate schedule consumption: update for new rotation consumption
        self.schedule.consumption += self.consumption - old_consumption

    @property
    def earliest_departure_next_rot(self):
        """Earliest possible departure of next rotation as datetime."""
        # noqa: DAR201
        return self.arrival_time + datetime.timedelta(hours=self.min_standing_time)

    @property
    def min_standing_time(self):
        """Minimum duration of standing time in minutes."""
        # noqa: DAR201
        assert self.charging_type in ["depb", "oppb"]
        if self.charging_type == "depb":
            capacity_depb = self.schedule.vehicle_types[self.vehicle_type]["depb"]["capacity"]
            # minimum time needed to recharge consumed power from depot charger
            min_standing_time = (self.consumption / self.schedule.cs_power_deps_depb)
            # time to charge battery from 0 to desired SOC
            desired_max_standing_time = ((capacity_depb / self.schedule.cs_power_deps_depb)
                                         * self.schedule.min_recharge_deps_depb)
            if min_standing_time > desired_max_standing_time:
                min_standing_time = desired_max_standing_time
        elif self.charging_type == "oppb":
            capacity_oppb = self.schedule.vehicle_types[self.vehicle_type]["oppb"]["capacity"]
            min_standing_time = ((capacity_oppb / self.schedule.cs_power_deps_oppb)
                                 * self.schedule.min_recharge_deps_oppb)
        return min_standing_time
