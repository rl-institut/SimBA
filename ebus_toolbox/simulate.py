# imports
from ebus_toolbox import optimizer, report
from ebus_toolbox.schedule import Schedule


def simulate(args=None):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    schedule = Schedule.from_csv(args.input)
    # filter trips according to args
    schedule.filter_rotations()

    # initialize optimizer

    while(True):
        # construct szenario and simulate in spice ev until optimizer is happy
        # if optimizer None, quit after single iteration
        schedule.calculate_consumption()
        schedule.set_charging_type(preferred_ct=args.preferred_charging_type)
        schedule.assign_vehicles()
        # write trips to csv in spiceEV format

        # RUN SPICE EV

        # Quit if optimizer is not defined
        # (LATER) Run optimizer, continue from top or quit based on optimizer output
        if optimizer.no_optimization() == 'converged':
            break

    # create report
    report.generate()
