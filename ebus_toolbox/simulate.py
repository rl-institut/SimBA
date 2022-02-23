# imports
from ebus_toolbox import optimizer, report
from ebus_toolbox.schedule import Schedule


def simulate(args=None):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    schedule = Schedule.from_csv("./data/private_examples/trips_example-bvg.csv")
    # filter trips according to args
    schedule.filter_rotations(args.filters)

    # initialize optimizer

    while(True):
        # construct szenario and simulate in spice ev until optimizer is happy
        # if optimizer None, quit after single iteration

        schedule.add_charging_types()
        schedule.calculate_consumption(args.consumption_func)
        schedule.assign_vehicles()
        # write trips to csv in spiceEV format

        # RUN SPICE EV

        # Quit if optimizer is not defined
        # (LATER) Run optimizer, continue from top or quit based on optimizer output
        optimizer.no_optimization()

    # create report
    report.generate()


if __name__ == '__main__':
    # parse args from cmd or file
    simulate()
