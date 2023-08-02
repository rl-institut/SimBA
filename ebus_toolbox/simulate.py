import warnings
from warnings import warn

from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox import report, optimization, util
from ebus_toolbox.station_optimization import run_optimization
from ebus_toolbox.optimizer_util import read_config as read_optimizer_config


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

    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(
        vehicle_types,
        outside_temperatures=args.outside_temperature_over_day_path,
        level_of_loading_over_day=args.level_of_loading_over_day_path)

    schedule = Schedule.from_csv(args.input_schedule,
                                 vehicle_types,
                                 stations,
                                 **vars(args))
    schedule.calculate_consumption()
    # scenario simulated once
    scenario = schedule.run(args)

    # run the mode(s) specified in config
    if not isinstance(args.mode, list):
        # backwards compatibility: run single mode
        args.mode = [args.mode]

    for i, mode in enumerate(args.mode):
        # scenario must be set from initial run / prior modes
        assert scenario is not None, f"Scenario became None after mode {args.mode[i-1]} (index {i})"

        if mode == 'service_optimization':
            # find largest set of rotations that produce no negative SoC
            result = optimization.service_optimization(schedule, scenario, args)
            schedule, scenario = result['optimized']
            if scenario is None:
                print('*'*49 + '\nNo optimization possible (all rotations negative), reverting')
                schedule, scenario = result['original']
        elif mode in ['neg_depb_to_oppb', 'neg_oppb_to_depb']:
            # simple optimization: change charging type, simulate again
            change_from = mode[4:8]
            change_to = mode[-4:]
            # get negative rotations
            neg_rot = schedule.get_negative_rotations(scenario)
            # check which rotations are relevant and if vehicle with other charging type exists
            neg_rot = [r for r in neg_rot if schedule.rotations[r].charging_type == change_from
                       if change_to in vehicle_types[schedule.rotations[r].vehicle_type]]
            if neg_rot:
                print(f'Changing charging type from {change_from} to {change_to} for rotations '
                      + ', '.join(neg_rot))
                schedule.set_charging_type(change_to, neg_rot)
                # simulate again
                scenario = schedule.run(args)
                neg_rot = schedule.get_negative_rotations(scenario)
                if neg_rot:
                    print(f'Rotations {", ".join(neg_rot)} remain negative.')
        elif mode == "station_optimization":
            if not args.optimizer_config:
                warnings.warn("Station optimization needs an optimization config file. "
                              "Since no path was given, station optimization is skipped")
                continue
            conf = read_optimizer_config(args.optimizer_config)
            try:
                create_results_directory(args, i+1)
                schedule, scenario = run_optimization(conf, sched=schedule, scen=scenario,
                                                      args=args)
            except Exception as err:
                warnings.warn('During Station optimization an error occurred {0}. '
                              'Optimization was skipped'.format(err))
        elif mode == 'remove_negative':
            neg_rot = schedule.get_negative_rotations(scenario)
            if neg_rot:
                schedule.rotations = {
                    k: v for k, v in schedule.rotations.items() if k not in neg_rot}
                print('Rotations ' + ', '.join(neg_rot) + ' removed')
                # re-run schedule
                scenario = schedule.run(args)
            else:
                print('No negative rotations to remove')
        elif mode == 'report':
            # create report based on all previous modes
            if args.cost_calculation:
                # cost calculation part of report
                calculate_costs(cost_parameters_file, scenario, schedule, args)
            # name: always start with sim, append all prior optimization modes
            create_results_directory(args, i)
            report.generate(schedule, scenario, args)
        elif mode == 'sim':
            if i > 0:
                # ignore anyway, but at least give feedback that this has no effect
                warn('Intermediate sim ignored')
        else:
            warn(f'Unknown mode {mode} ignored')


def create_results_directory(args, i):
    """ Create directory for results.

    :param args: arguments
    :type args: Namespace
    :param i: iteration number of loop
    :type i: int
    """
    prior_modes = ['sim'] + [m for m in args.mode[:i] if m not in ['sim', 'report']]
    report_name = '__'.join(prior_modes)
    args.results_directory = args.output_directory.joinpath(report_name)
    args.results_directory.mkdir(parents=True, exist_ok=True)
