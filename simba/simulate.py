import logging
import traceback

from simba import report, optimization, util
from simba.consumption import Consumption
from simba.costs import calculate_costs
from simba.optimizer_util import read_config as read_optimizer_config
from simba.schedule import Schedule
from simba.station_optimization import run_optimization
from simba.trip import Trip


def simulate(args):
    """ High-level function to create and run a scenario.

    Use the parsed scenario arguments to prepare the schedule,
    run the basic simulation and different modes of SimBA.
    Returns final schedule and scenario after running all modes.

    :param args: configuration
    :type args: Namespace
    :return: final schedule and scenario
    :rtype: tuple
    """
    schedule = pre_simulation(args)
    scenario = schedule.run(args)
    return modes_simulation(schedule, scenario, args)


def pre_simulation(args):
    """ Prepare simulation.

    Read in files, generate consumption info, create and filter schedule.

    :param args: arguments
    :type args: Namespace
    :raises Exception: If an input file does not exist, exit the program.
    :return: schedule
    :rtype: simba.schedule.Schedule
    """
    try:
        with open(args.vehicle_types, encoding='utf-8') as f:
            vehicle_types = util.uncomment_json_file(f)
            del args.vehicle_types
    except FileNotFoundError:
        raise Exception(f"Path to vehicle types ({args.vehicle_types}) "
                        "does not exist. Exiting...")

    # load stations file
    try:
        with open(args.electrified_stations, encoding='utf-8') as f:
            stations = util.uncomment_json_file(f)
    except FileNotFoundError:
        raise Exception(f"Path to electrified stations ({args.electrified_stations}) "
                        "does not exist. Exiting...")

    # load cost parameters
    if args.cost_parameters_file is not None:
        try:
            with open(args.cost_parameters_file, encoding='utf-8') as f:
                args.cost_parameters = util.uncomment_json_file(f)
        except FileNotFoundError:
            raise Exception(f"Path to cost parameters ({args.cost_parameters_file}) "
                            "does not exist. Exiting...")

    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(
        vehicle_types,
        outside_temperatures=args.outside_temperature_over_day_path,
        level_of_loading_over_day=args.level_of_loading_over_day_path)

    # generate schedule from csv
    schedule = Schedule.from_csv(args.input_schedule, vehicle_types, stations, **vars(args))

    # filter rotations
    schedule.rotation_filter(args)

    # calculate consumption of all trips
    schedule.calculate_consumption()

    return schedule


def modes_simulation(schedule, scenario, args):
    """ Run the mode(s) specified in config.

    Ignores unknown and "sim" modes.
    On error, create a report and continue with next mode.

    :param schedule: input schedule
    :type schedule: simba.schedule.Schedule
    :param scenario: corresponding scenario
    :type scenario: spice_ev.scenario.Scenario
    :param args: arguments
    :type args: Namespace
    :return: final schedule and scenario
    :rtype: tuple
    :raises Exception: if args.propagate_mode_errors is set, re-raises error instead of continuing
    """

    if not isinstance(args.mode, list):
        # backwards compatibility: run single mode
        args.mode = [args.mode]

    for i, mode in enumerate(args.mode):
        # scenario must be set from initial run / prior modes
        assert scenario is not None, f"Scenario became None after mode {args.mode[i-1]} (index {i})"

        if mode == "sim":
            if i > 0:
                # ignore anyway, but at least give feedback that this has no effect
                logging.info('Intermediate sim ignored')
            continue

        logging.debug("Starting mode " + mode)
        try:
            func = getattr(Mode, mode)
        except AttributeError:
            logging.error(f'Unknown mode {mode} ignored')
            continue

        try:
            schedule, scenario = func(schedule, scenario, args, i)
            logging.debug("Finished mode " + mode)
        except Exception as e:
            if args.propagate_mode_errors:
                raise
            msg = f"{e.__class__.__name__} during {mode}: {e}"
            logging.error('*'*len(msg))
            logging.error(msg)
            logging.error('*'*len(msg))
            logging.error(traceback.format_exc())
            if scenario is not None and scenario.step_i > 0:
                # generate plot of failed scenario
                args.mode = args.mode[:i] + ["ABORTED"]
                if args.output_directory is None:
                    create_results_directory(args, i+1)
                    report.generate_plots(scenario, args)
                    logging.info(f"Created plot of failed scenario in {args.results_directory}")
            # continue with other modes after error

    # all modes done
    return schedule, scenario


class Mode:
    """ Container for simulation modes.

    Each function takes a schedule, the corresponding scenario and arguments Namespace.
    Optionally, an index of the current mode in the modelist can be given.
    A function must return the updated schedule and scenario objects.
    """

    def service_optimization(schedule, scenario, args, _i):
        # find largest set of rotations that produce no negative SoC
        result = optimization.service_optimization(schedule, scenario, args)
        schedule, scenario = result['optimized']
        if scenario is None:
            logging.warning('No optimization possible (all rotations negative), reverting')
            schedule, scenario = result['original']
        return schedule, scenario

    def neg_depb_to_oppb(schedule, scenario, args, _i):
        return Mode.switch_type(schedule, scenario, args, "depb", "oppb")

    def neg_oppb_to_depb(schedule, scenario, args, _i):
        return Mode.switch_type(schedule, scenario, args, "oppb", "depb")

    def switch_type(schedule, scenario, args, from_type, to_type):
        # simple optimization: change charging type, simulate again
        # get negative rotations
        neg_rot = schedule.get_negative_rotations(scenario)
        # check which rotations are relevant and if vehicle with other charging type exists
        relevant_rotations = []
        for rot_id in neg_rot:
            if schedule.rotations[rot_id].charging_type == from_type:
                new_vehicle_type = f"{schedule.rotations[rot_id].vehicle_type}_{to_type}"
                if scenario.components.vehicle_types.get(new_vehicle_type) is not None:
                    relevant_rotations.append(rot_id)
                else:
                    logging.debug(f"Rotation {rot_id} negative and of type {from_type}, "
                                  f"but vehicle type {new_vehicle_type} does not exist.")
        if relevant_rotations:
            logging.info(
                f'Changing charging type from {from_type} to {to_type} for rotations '
                + ', '.join(sorted(relevant_rotations)))
            schedule.set_charging_type(to_type, relevant_rotations)
            # simulate again
            scenario = schedule.run(args)
            neg_rot = schedule.get_negative_rotations(scenario)
            if neg_rot:
                logging.info(f'Rotations {", ".join(neg_rot)} remain negative.')
        return schedule, scenario

    def station_optimization(schedule, scenario, args, i):
        if not args.optimizer_config:
            logging.warning("Station optimization needs an optimization config file. "
                            "Since no path was given, station optimization is skipped")
            return schedule, scenario
        conf = read_optimizer_config(args.optimizer_config)
        try:
            create_results_directory(args, i+1)
            return run_optimization(conf, sched=schedule, scen=scenario, args=args)
        except Exception as err:
            logging.warning('During Station optimization an error occurred {0}. '
                            'Optimization was skipped'.format(err))
            return schedule, scenario

    def remove_negative(schedule, scenario, args, _i):
        neg_rot = schedule.get_negative_rotations(scenario)
        if neg_rot:
            schedule.rotations = {
                k: v for k, v in schedule.rotations.items() if k not in neg_rot}
            logging.info('Rotations ' + ', '.join(sorted(neg_rot)) + ' removed')
            # re-run schedule
            scenario = schedule.run(args)
        else:
            logging.info('No negative rotations to remove')
        return schedule, scenario

    def recombination(schedule, scenario, args, _i):
        trips, depot_trips = optimization.prepare_trips(schedule)
        recombined_schedule = optimization.recombination(schedule, args, trips, depot_trips)
        # re-run schedule
        scenario = recombined_schedule.run(args)
        return recombined_schedule, scenario

    def report(schedule, scenario, args, i):
        if args.output_directory is None:
            return schedule, scenario

        # create report based on all previous modes
        if args.cost_calculation:
            # cost calculation part of report
            try:
                calculate_costs(args.cost_parameters, scenario, schedule, args)
            except Exception:
                logging.warning(f"Cost calculation failed due to {traceback.print_exc()}")
                if args.propagate_mode_errors:
                    raise
        # name: always start with sim, append all prior optimization modes
        create_results_directory(args, i)
        report.generate(schedule, scenario, args)
        return schedule, scenario


def create_results_directory(args, i):
    """ Create directory for results.

    :param args: arguments
    :type args: Namespace
    :param i: iteration number of loop
    :type i: int
    """

    if args.output_directory is None:
        return

    prior_reports = sum([m.count('report') for m in args.mode[:i]])
    report_name = f"report_{prior_reports+1}"
    args.results_directory = args.output_directory.joinpath(report_name)
    args.results_directory.mkdir(parents=True, exist_ok=True)
    # save used modes in report version
    used_modes = ['sim'] + [m for m in args.mode[:i] if m not in ['sim', 'report']]
    with open(args.results_directory / "used_modes.txt", "w", encoding='utf-8') as f:
        f.write(f"Used modes in this scenario: {', '.join(used_modes)}")
