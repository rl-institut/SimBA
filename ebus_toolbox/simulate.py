# imports
import json
import warnings
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox import optimizer, report


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """
    try:
        with open(args.vehicle_types) as f:
            vehicle_types = json.load(f)
    except FileNotFoundError:
        warnings.warn("Invalid path for vehicle type JSON. Using default types from EXAMPLE dir.")
        with open("data/examples/vehicle_types.json") as f:
            vehicle_types = json.load(f)

    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(vehicle_types)

    schedule = Schedule.from_csv(args.input, vehicle_types)
    # filter trips according to args
    schedule.filter_rotations()

    # initialize optimizer

    while(True):
        # construct szenario and simulate in spice ev until optimizer is happy
        # if optimizer None, quit after single iteration
        schedule.calculate_consumption()
        schedule.set_charging_type(preferred_ct=args.preferred_charging_type)
        schedule.assign_vehicles(args.min_standing_time)
        # write trips to csv in spiceEV format

        # RUN SPICE EV

        # Quit if optimizer is not defined
        # (LATER) Run optimizer, continue from top or quit based on optimizer output
        if optimizer.no_optimization() == 'converged':
            break

    # create report
    report.generate()
