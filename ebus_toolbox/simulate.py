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


    # read filter rotation file
    if args.rotation_filter_variable:
        try:
            with open(args.rotation_filter, encoding='utf-8') as f:
                # TODO: file zu dict umwandeln
                rotation_filter = util.file_wrapper_to_dict(f)
                print(rotation_filter)
        except FileNotFoundError:
            print(f"Path to rotation filter ({args.rotation_filter}) does not exist.")

        # TODO: filter rotations
        if args.rotation_filter_variable == "exclude_rotation":
            excluding_rotations = rotation_filter
            for rotation in rotation_filter:
                if rotation in rotation_filter:
                    try:
                        schedule.rotations.pop(rotation['rotation_id'])
                    except KeyError:
                        print(f"Key of rotation {rotation['rotation_id']} does not exist in schedule.")
        elif args.rotation_filter_variable == "include_rotation":
            included_rotations = rotation_filter
            staying_rotations = dict()
            for rot_filter in rotation_filter:
                for rot_schedule in schedule.rotations:
                    if rot_filter == rot_schedule:
                        staying_rotations.update(rot_schedule)




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
        # only depot rotations relevant and check if oppb-version of vehicle_type exists
        neg_rot = [r for r in neg_rot if schedule.rotations[r].charging_type == "depb"
                   if "oppb" in vehicle_types[schedule.rotations[r].vehicle_type]]
        if neg_rot:
            print("Changing charging type from depb to oppb for rotations " + ', '.join(neg_rot))
            schedule.set_charging_type("oppb", neg_rot)
            # simulate again
            scenario = schedule.run(args)
            neg_rot = schedule.get_negative_rotations(scenario)
            if neg_rot:
                print(f"Rotations {', '.join(neg_rot)} remain negative.")

    if args.cost_calculation:
        # cost calculation following directly after simulation
        calculate_costs(cost_parameters_file, scenario, schedule, args)

    # create report
    report.generate(schedule, scenario, args)
