import datetime

from ebus_toolbox.trip import Trip


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

    def add_trip(self, trip):
        """Create a trip object and append to rotations trip set

        :param trip: Information on trip to be added to rotation
        :type trip: dict
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
        if ('charging_type' in trip and
                any(trip['charging_type'] == t for t in ['depb', 'oppb'])):
            assert (self.charging_type is None or self.charging_type == trip['charging_type']),\
                f"Two trips of rotation {self.id} have distinct charging types"
            self.set_charging_type(trip['charging_type'])

        self.trips.append(new_trip)

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

    def delta_soc_all_trips(self):
        """ Compute change in state of charge (SOC) for every trip
            of this rotation. Stored in the trip objects.
        """
        for trip in self.trips:
            trip.get_delta_soc()

    def set_charging_type(self, ct):
        """ Change charging type of either all or specified rotations. Adjust minimum standing time
            at depot after completion of rotation.

        :param ct: Choose this charging type wheneever possible. Either 'depb' or 'oppb'.
        :type ct: str
        """
        assert ct in ["oppb", "depb"], f"Invalid charging type: {ct}"

        capacity_depb = self.schedule.vehicle_types[f"{self.vehicle_type}_depb"]["capacity"]
        if ct == "oppb" or capacity_depb < self.consumption:
            self.charging_type = "oppb"
            capacity_oppb = self.schedule.vehicle_types[f"{self.vehicle_type}_oppb"]["capacity"]
            min_standing_time = ((capacity_oppb / self.schedule.cs_power_deps_oppb)
                                 * self.schedule.min_recharge_deps_oppb)
        else:
            self.charging_type = "depb"
            min_standing_time = (self.consumption / self.schedule.cs_power_deps_depb)
            desired_max_standing_time = ((capacity_depb / self.schedule.cs_power_deps_depb)
                                         * self.schedule.min_recharge_deps_depb)
            if min_standing_time > desired_max_standing_time:
                min_standing_time = desired_max_standing_time

        self.earliest_departure_next_rot = \
            self.arrival_time + datetime.timedelta(hours=min_standing_time)
