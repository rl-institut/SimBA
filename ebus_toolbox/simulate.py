# imports
import json
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox import report, optimization


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace

    :raises SystemExit: If an input file does not exist, exit the program.
    """
    # load vehicle types
    try:
        with open(args.vehicle_types) as f:
            vehicle_types = json.load(f)
            del args.vehicle_types
    except FileNotFoundError:
        raise SystemExit(f"Path to vehicle types ({args.vehicle_types}) "
                         "does not exist. Exiting...")

    # load stations file
    try:
        with open(args.electrified_stations) as f:
            stations = json.load(f)
    except FileNotFoundError:
        raise SystemExit(f"Path to electrified stations ({args.electrified_stations}) "
                         "does not exist. Exiting...")

    # load cost parameters
    if args.cost_params is not None:
        try:
            with open(args.cost_params) as f:
                cost_params = json.load(f)
        except FileNotFoundError:
            raise SystemExit(f"Path to cost parameters ({args.cost_params}) "
                             "does not exist. Exiting...")

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
                                 stations,
                                 **vars(args))

    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(vehicle_types)
    # set charging type for all rotations without explicitly specified charging type
    for rot in schedule.rotations.values():
        if rot.charging_type is None:
            rot.set_charging_type(ct=args.preferred_charging_type)

    # run the mode specified in config
    if args.mode == 'service_optimization':
        scenario = optimization.service_optimization(schedule, args)
    elif args.mode == "sim":
        # DEFAULT if mode argument is not specified by user
        # Scenario simulated once
        scenario = schedule.run(args)

    if args.cost_params is not None:
        # Calculate Costs of Iteration
        costs = calculate_costs(cost_params, schedule)
        opex_energy_annual = 0  # ToDo: Import annual energy costs from SpiceEV
        cost_invest = costs["c_invest"]
        cost_annual = costs["c_invest_annual"] + costs["c_maintenance_annual"] + opex_energy_annual
        print(f"Investment cost: {cost_invest} €. Total annual cost: {cost_annual} €.")

    # create report
    report.generate(schedule, scenario, args)
