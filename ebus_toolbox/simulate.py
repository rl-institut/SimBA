# imports
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox import report, optimization, util


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace

    :raises SystemExit: If an input file does not exist, exit the program.
    """
    # load vehicle types
    try:
        with open(args.vehicle_types) as f:
            vehicle_types = util.uncomment_json_file(f)
            del args.vehicle_types
    except FileNotFoundError:
        raise SystemExit(f"Path to vehicle types ({args.vehicle_types}) "
                         "does not exist. Exiting...")

    # load stations file
    try:
        with open(args.electrified_stations) as f:
            stations = util.uncomment_json_file(f)
    except FileNotFoundError:
        raise SystemExit(f"Path to electrified stations ({args.electrified_stations}) "
                         "does not exist. Exiting...")

    # load cost parameters
    if args.cost_params is not None:
        try:
            with open(args.cost_params) as f:
                cost_params = util.uncomment_json_file(f)
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
    Trip.consumption = Consumption(vehicle_types,
                                   outside_temperatures=args.outside_temperature_over_day_path,
                                   level_of_loading_over_day=args.level_of_loading_over_day_path)

    schedule = Schedule.from_csv(args.input_schedule,
                                 vehicle_types,
                                 stations,
                                 **vars(args))
    schedule.calculate_consumption()

    # run the mode specified in config
    if args.mode == 'service_optimization':
        schedule, scenario = optimization.service_optimization(schedule, args)["optimized"]
    elif args.mode in ["sim", "neg_depb_to_oppb"]:
        # DEFAULT if mode argument is not specified by user
        # Scenario simulated once
        scenario = schedule.run(args)
    if args.mode == "neg_depb_to_oppb":
        # simple optimization: change charging type from depot to opportunity, simulate again
        neg_rot = schedule.get_negative_rotations(scenario)
        # only depot rotations relevant
        neg_rot = [r for r in neg_rot if schedule.rotations[r].charging_type == "depb"]
        if neg_rot:
            print("Changing charging type from depb to oppb for rotations " + ', '.join(neg_rot))
            schedule.set_charging_type("oppb", neg_rot)
            # simulate again
            scenario = schedule.run(args)
            neg_rot = schedule.get_negative_rotations(scenario)
            if neg_rot:
                print(f"Rotations {', '.join(neg_rot)} remain negative.")

    if args.cost_params is not None:
        # Calculate Costs of Iteration
        costs = calculate_costs(cost_params, schedule)
        opex_energy_annual = 0  # ToDo: Import annual energy costs from SpiceEV
        cost_invest = costs["c_invest"]
        cost_annual = costs["c_invest_annual"] + costs["c_maintenance_annual"] + opex_energy_annual
        print(f"Investment cost: {cost_invest} €. Total annual cost: {cost_annual} €.")

    import pickle
    with open("schedule_buffered_depots.pickle", "wb") as f:
        pickle.dump(schedule, f)
    with open("scenario_buffered_depots.pickle", "wb") as f:
        pickle.dump(scenario, f)
    with open("args_buffered.pickle", "wb") as f:
        pickle.dump(args, f)

    print("pickled")

    # create report
    report.generate(schedule, scenario, args)
