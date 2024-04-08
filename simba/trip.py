from datetime import datetime


class Trip:
    def __init__(self, departure_time, departure_name,
                 arrival_time, arrival_name, distance,
                 line=None, temperature=None, level_of_loading=None, mean_speed=None,
                 height_diff=None, station_data=None,
                 charging_type=None, vehicle_type=None,
                 ):
        self.departure_time = departure_time
        self.departure_name = departure_name
        self.arrival_time = arrival_time
        self.arrival_name = arrival_name
        self.distance = float(distance)

        # optional / computed
        self.height_diff = height_diff  # m
        self.level_of_loading = level_of_loading  # [0-1]
        self.line = line
        self.mean_speed = mean_speed  # km/h
        self.temperature = temperature  # deg C

        # needed to update rotation properties when adding trip
        self.charging_type = charging_type  # oppb/depb/None
        self.vehicle_type = vehicle_type

        # set when added to rotation
        # needed for consumption calculation and other references, like schedule
        # !Circular reference! problematic for garbage collector in older Python versions
        self.rotation = None

        # set after calculating consumption
        self.consumption = None  # kWh
        self.delta_soc = None

        # ---- cast / compute options ---- #
        # cast times to datetime
        if type(departure_time) is str:
            self.departure_time = datetime.fromisoformat(departure_time)
        if type(arrival_time) is str:
            self.arrival_time = datetime.fromisoformat(arrival_time)

        # cast temperature to float
        try:
            self.temperature = float(self.temperature)
            # In case of empty temperature column or no column at all
        except (TypeError, ValueError):
            self.temperature = None

        # get height difference from station data if not given
        assert height_diff is not None or station_data is not None, (
            "New trip: need either height_diff or station_data")
        if height_diff is None:
            try:
                self.height_diff = (station_data[self.arrival_name]["elevation"]
                                    - station_data[self.departure_name]["elevation"])
            except (KeyError, TypeError):
                self.height_diff = 0

        # clip level of loading to [0,1],
        # but may also be None in case of empty temperature column or no column at all
        try:
            self.level_of_loading = max(0, min(float(level_of_loading), 1))
        except (TypeError, ValueError):
            self.level_of_loading = None

        # mean speed in km/h from distance and travel time if not given
        if mean_speed is None:
            travel_time = (self.arrival_time - self.departure_time).total_seconds()
            # travel time is at least 1 min
            travel_time = max(travel_time, 60)
            mean_speed = self.distance * travel_time  # m/s
            self.mean_speed = mean_speed * 3.6  # km/h

    def calculate_consumption(self):
        """ Compute consumption for this trip.

        :return: Consumption of trip [kWh]
        :rtype: float
        :raises with_traceback: if consumption cannot be constructed
        """

        try:
            self.consumption, self.delta_soc = Trip.consumption.calculate_consumption(
                self.arrival_time, self.distance,
                self.vehicle_type, self.charging_type,
                temp=self.temperature,  height_diff=self.height_diff,
                level_of_loading=self.level_of_loading, mean_speed=self.mean_speed)
        except AttributeError as e:
            raise Exception(
                'To calculate consumption, a consumption object needs to be constructed '
                'and linked to Trip class.').with_traceback(e.__traceback__)

        return self.consumption

    def get_buffer_time(self, default=0):
        """ Get buffer time at arrival station of a trip.

        Buffer time is an abstraction of delays like
        docking procedures and is added to the planned arrival time.

        :param default: Default buffer time if no station specific buffer time is given. [minutes]
        :type default: dict, numeric
        :return: buffer time in minutes
        :rtype: dict or int

        NOTE: Buffer time dictionaries map hours of the day to a buffer time.
        Keys are ranges of hours and corresponding values provide buffer time in
        minutes for that time range.
        An entry with key "else" is a must if not all hours of the day are covered.
        Example: ``buffer_time = {"10-22": 2, "22-6": 3, "else": 1}``
        """

        schedule = self.rotation.schedule
        buffer_time = schedule.stations.get(self.arrival_name, {}).get('buffer_time', default)

        # distinct buffer times depending on time of day can be provided
        # in that case buffer time is of type dict instead of int
        if isinstance(buffer_time, dict):
            # sort dict to make sure 'else' key is last key
            buffer_time = {key: buffer_time[key] for key in sorted(buffer_time)}
            current_hour = self.arrival_time.hour
            for time_range, buffer in buffer_time.items():
                if time_range == 'else':
                    buffer_time = buffer
                    break
                else:
                    start_hour, end_hour = [int(t) for t in time_range.split('-')]
                    if end_hour < start_hour:
                        if current_hour >= start_hour or current_hour < end_hour:
                            buffer_time = buffer
                            break
                    else:
                        if start_hour <= current_hour < end_hour:
                            buffer_time = buffer
                            break
        return buffer_time
