# imports

def check_trip_data(trips):
    """ The provided trips file is checked against the following requirements:
    - Folloing columns are provided:
        Required:
            Rotation,
            Line,
            Leg,
            Departure time,
            Departure day,
            Arrival time,
            Arrival day,
            Break at arrival,
            Vehicle type,
            Departure name,
            Departure ID,
            Arrival name,
            Arrival ID,
            Distance,
        Optional:
            Empty trip,
            Departure (short name),
            Arrival (short name),
    -
    -

    :param trips: Dictionaries for every trip
    :type trips: list

    :return: True if data is valid, else False.
    :rtype: bool
    """
    return True


def add_charging_types():
    """Add charging types to each vehicle if not already given by trips or vehicle type.
    """
    pass
