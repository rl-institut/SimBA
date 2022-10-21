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

    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(vehicle_types)

    schedule = Schedule.from_csv(args.input_schedule,
                                 vehicle_types,
                                 stations,
                                 **vars(args))
    scenario = None

    # run the mode(s) specified in config
    if type(args.mode) != list:
        # backwards compatibility: run single mode
        args.mode = [args.mode]

    for i, mode in enumerate(args.mode):
        if mode == "sim":
            # scenario simulated once
            # default if mode argument is not specified by user
            scenario = schedule.run(args)
        elif mode == "neg_depb_to_oppb":
            # simple optimization
            # change charging type of negative depb rotations from depot to opportunity
            if scenario is None:
                # not simulated yet
                scenario = schedule.run(args)
            neg_rot = schedule.get_negative_rotations(scenario)
            # only depot rotations relevant
            neg_rot = [r for r in neg_rot if schedule.rotations[r].charging_type == "depb"]
            if neg_rot:
                print("Changing charging type to oppb for rotations " + ', '.join(neg_rot))
                schedule.set_charging_type("oppb", neg_rot)
                # simulate again
                scenario = schedule.run(args)
                neg_rot = schedule.get_negative_rotations(scenario)
                if neg_rot:
                    print(f"Rotations {', '.join(neg_rot)} remain negative.")
        elif mode == 'service_optimization':
            # find largest set of rotations that produce no negative SoC
            schedule, scenario = optimization.service_optimization(schedule, args)["optimized"]
            if scenario is None:
                print("*"*49 + "\nNo optimization possible (all rotations negative)")
        elif mode == 'report':
            # create report based on all previous modes
            assert scenario is not None, "Can't report without simulation"
            report_name = '__'.join([m for m in args.mode[:i] if m not in ["report", "cost"]])
            report.generate(schedule, scenario, args, prefix=report_name + '_')
        elif mode == 'cost' and args.cost_params is not None:
            # calculate costs of iteration
            costs = calculate_costs(cost_params, schedule)
            opex_energy_annual = 0  # ToDo: Import annual energy costs from SpiceEV
            cost_invest = costs["c_invest"]
            cost_annual = (opex_energy_annual
                           + costs["c_invest_annual"]
                           + costs["c_maintenance_annual"])
            print(f"Investment cost: {cost_invest} €. Total annual cost: {cost_annual} €.")
