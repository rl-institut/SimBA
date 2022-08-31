# imports
import json
import warnings
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox import report, optimization


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """
    try:
        with open(args.vehicle_types) as f:
            vehicle_types = json.load(f)
            del args.vehicle_types
    except FileNotFoundError:
        warnings.warn("Invalid path for vehicle type JSON. Using default types from EXAMPLE dir.")
        with open("data/examples/vehicle_types.json") as f:
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

    schedule = Schedule.from_csv(args.input_schedule,
                                 vehicle_types,
                                 args.electrified_stations,
                                 **vars(args))
    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(vehicle_types)
    schedule.calculate_consumption()
    # set charging type for all rotations without explicitly specified charging type
    for rot in schedule.rotations.values():
        if rot.charging_type is None:
            rot.set_charging_type(ct=args.preferred_charging_type)

    # run the mode specified in config
    if args.mode == 'service_optimization':
        scenario = optimization.service_optimization(schedule, args)
    elif args.mode == "sim":
        # DEFAULT if mode argument is not specified by user

        # calculate the change in SoC for every trip
        # charging types may have changed which may impact battery capacity
        # while mileage is assumed to stay constant
        schedule.delta_soc_all_trips()
        # each rotation is assigned a vehicle ID
        schedule.assign_vehicles()
        scenario = schedule.generate_scenario(args)
        print("Running Spice EV...")
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            scenario.run('distributed', vars(args).copy())
        print("Spice EV simulation complete.")

        print(f"Rotations {schedule.get_negative_rotations(scenario)} have negative SoC.")

    # create report
    report.generate(schedule, scenario, args)
