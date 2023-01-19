# imports
import json
import warnings
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox import report  # , optimizer
from ebus_toolbox.run_sensitivity import prepare_sensitivity
from ebus_toolbox.run_sensitivity import run_sensitivity


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    if args.flag_sensitivity:
        prepare_sensitivity(args)

    try:
        with open(args.vehicle_types) as f:
            vehicle_types = json.load(f)
    except FileNotFoundError:
        warnings.warn("Invalid path for vehicle type JSON. Using default types from EXAMPLE dir.")
        with open("data/uebertragbarkeitsanalyse/vehicle_types.json") as f:
            vehicle_types = json.load(f)

    # parse strategy options for Spice EV
    if args.strategy_option is not None:
        for opt_key, opt_val in args.strategy_option:
            try:
                # option may be number
                opt_val = float(opt_val)
            except ValueError:
                # or not
                pass
            setattr(args, opt_key, opt_val)

    schedule = Schedule.from_csv(args.input_schedule, vehicle_types)
    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(vehicle_types)
    # filter trips according to args
    schedule.filter_rotations()
    schedule.calculate_consumption()
    schedule.set_charging_type(preferred_ct=args.preferred_charging_type, args=args)

    for i in range(args.iterations):
        # (re)calculate the change in SoC for every trip
        # charging types may have changed which may impact battery capacity
        # while mileage is assumed to stay constant
        schedule.delta_soc_all_trips()

        # each rotation is assigned a vehicle ID
        schedule.assign_vehicles()

        scenario = schedule.generate_scenario(args)

        print("Running Spice EV...")
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            # for analyzes
            if args.flag_sensitivity:
                for ix in range(1):
                    run_sensitivity(args, ix)  # args gets changed
                    # Todo create new output folder for every scenario
                    scenario.run('distributed', vars(args).copy())
            else:
                scenario.run('distributed', vars(args).copy())

            # scenario.run('distributed', vars(args).copy())
        print(f"Spice EV simulation complete. (Iteration {i})")

        if i < args.iterations - 1:
            # TODO: replace with optimizer step in the future
            schedule.readjust_charging_type(args, scenario)

    print(f"Rotations {schedule.get_negative_rotations(scenario)} have negative SoC.")

    # create report
    report.generate(schedule, scenario, args)
