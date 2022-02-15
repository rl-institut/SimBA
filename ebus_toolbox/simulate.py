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


def consumption():
    """ Calculate the consumptions for each trip and rotation.
    Various consumption schemes possible depending on available data.
    """
    pass


def assign_vehicles():
    """ Assign vehicle IDs to rotations.
    Just randomly? Match consumption with battery sizes?
    """
    pass


def filter_rotations(trips):
    """ Iterate across trips and remove trips that fit criteria given in config.
    Possible filters:
    - specific rotation IDs
    - distance thresholds
    -

    :param trips: Dictionaries for every trip
    :type trips: list
    """
    pass


def simulate(args=None):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    # read CSV input file
    # check for correct data format / standard
    check_trip_data(None)

    # filter trips according to args
    filter_rotations(None)

    # initialize optimizer

    while(True):
        # construct szenario and simulate in spice ev until optimizer is happy
        # if optimizer None, quit after single iteration

        add_charging_types()
        consumption()
        assign_vehicles()

        # write trips to csv in spiceEV format

        # RUN SPICE EV

        # Quit if optimizer is not defined
        # (LATER) Run optimizer, continue from top or quit based on optimizer output

    # create report


if __name__ == 'main':
    # parse args from cmd or file
    simulate()
