# imports
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox import report, optimization, util
from ebus_toolbox.run_sensitivity import prepare_sensitivity


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace

    :raises SystemExit: If an input file does not exist, exit the program.
    """
    # load vehicle types
    if args.mode == 'sensitivity_analysis':
        args = prepare_sensitivity(args)

    try:
        with open(args.vehicle_types, encoding='utf-8') as f:
            vehicle_types = util.uncomment_json_file(f)
            del args.vehicle_types
    except FileNotFoundError:
        raise SystemExit(f"Path to vehicle types ({args.vehicle_types}) "
                         "does not exist. Exiting...")

    # load stations file
    try:
        with open(args.electrified_stations, encoding='utf-8') as f:
            stations = util.uncomment_json_file(f)
    except FileNotFoundError:
        raise SystemExit(f"Path to electrified stations ({args.electrified_stations}) "
                         "does not exist. Exiting...")

    # load cost parameters
    if args.cost_parameters_file is not None:
        try:
            with open(args.cost_parameters_file, encoding='utf-8') as f:
                cost_parameters_file = util.uncomment_json_file(f)
        except FileNotFoundError:
            raise SystemExit(f"Path to cost parameters ({args.cost_parameters_file}) "
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

    schedule.mode = 0

    # run the mode specified in config
    if args.mode == 'service_optimization':
        schedule, scenario = optimization.service_optimization(schedule, args)["optimized"]
    elif args.mode in ["sim", "neg_depb_to_oppb", "neg_oppb_to_depb"]:
        # DEFAULT if mode argument is not specified by user
        # Scenario simulated once
        scenario = schedule.run(args)
    elif args.mode == 'sensitivity_analysis':
        schedule.mode = 1
        schedule.range_sensitivity = 1
        schedule.run(args)
        return
    if args.mode in ["neg_depb_to_oppb", "neg_oppb_to_depb"]:
        # simple optimization: change charging type, simulate again
        change_from = args.mode[4:8]
        change_to = args.mode[-4:]
        # get negative rotations
        neg_rot = schedule.get_negative_rotations(scenario)
        # check which rotations are relevant and if vehicle with other charging type exists
        neg_rot = [r for r in neg_rot if schedule.rotations[r].charging_type == change_from
                   if change_to in vehicle_types[schedule.rotations[r].vehicle_type]]
        if neg_rot:
            print(f"Changing charging type from {change_from} to {change_to} for rotations "
                  + ', '.join(neg_rot))
            schedule.set_charging_type(change_to, neg_rot)
            # simulate again
            scenario = schedule.run(args)
            neg_rot = schedule.get_negative_rotations(scenario)
            if neg_rot:
                print(f"Rotations {', '.join(neg_rot)} remain negative.")

    if args.cost_calculation:
        # cost calculation following directly after simulation
        calculate_costs(cost_parameters_file, scenario, schedule, args)

    # create report
    if args.generate_report:
        args.save_timeseries = args.output_directory / "simulation_timeseries.csv"
        args.save_results = args.output_directory / "simulation.json"
        args.save_soc = args.output_directory / "vehicle_socs.csv"

        report.generate(schedule, scenario, args)

        args.save_timeseries = None
        args.save_results = None
        args.save_soc = None
