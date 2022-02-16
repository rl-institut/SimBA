# imports
import preprocessing
import consumption
import vehicle_assignment
import rotation_filter
import optimizer
import report


def simulate(args=None):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    # read CSV input file
    # check for correct data format / standard
    preprocessing.validate(None)

    # filter trips according to args
    rotation_filter.filter(None)

    # initialize optimizer

    while(True):
        # construct szenario and simulate in spice ev until optimizer is happy
        # if optimizer None, quit after single iteration

        preprocessing.add_charging_types()
        consumption.constant()
        vehicle_assignment.fishgrid()

        # write trips to csv in spiceEV format

        # RUN SPICE EV

        # Quit if optimizer is not defined
        # (LATER) Run optimizer, continue from top or quit based on optimizer output
        optimizer.no_optimization()

    # create report
    report.generate()


if __name__ == 'main':
    # parse args from cmd or file
    simulate()
