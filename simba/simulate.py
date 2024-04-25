import logging
import traceback
from copy import deepcopy

from simba import report, optimization, util
from simba.data_container import DataContainer
from simba.costs import calculate_costs
from simba.optimizer_util import read_config as read_optimizer_config
from simba.schedule import Schedule
from simba.station_optimization import run_optimization
from simba.trip import Trip


def simulate(args):
    """
    High-level function to create and run a scenario.

    Use the parsed scenario arguments to prepare the schedule,
    run the basic simulation and different modes of SimBA.
    Returns final schedule and scenario after running all modes.

    :param args: configuration
    :type args: Namespace
    :return: final schedule and scenario
    :rtype: tuple
    """
    # The data_container stores various input data.
    data_container = create_and_fill_data_container(args)

    schedule, args = pre_simulation(args, data_container)
    scenario = schedule.run(args)
    schedule, scenario = modes_simulation(schedule, scenario, args)
    return schedule, scenario


def create_and_fill_data_container(args):
    data_container = DataContainer()
    # Add the vehicle_types from a json file
    data_container.add_vehicle_types_from_json(args.vehicle_types_path)
    # Add consumption data, which is found in the vehicle_type data
    data_container.add_consumption_data_from_vehicle_type_linked_files()
    return data_container


def pre_simulation(args, data_container: DataContainer):
    """
    Prepare simulation.

    Read in files, generate consumption info, create and filter schedule.

    :param args: arguments
    :type args: Namespace
    :param data_container: data needed for simulation
    :type data_container: DataContainer
    :raises Exception: If an input file does not exist, exit the program.
    :return: schedule, args
    :rtype: simba.schedule.Schedule, Namespace
    """
    # Deepcopy args so original args do not get mutated, i.e. deleted
    args = deepcopy(args)

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

    # Add consumption calculator to trip class
    Trip.consumption = data_container.to_consumption()

    # generate schedule from csv
    schedule = Schedule.from_csv(args.input_schedule, data_container.vehicle_types_data, stations,
                                 **vars(args))

    # filter rotations
    schedule.rotation_filter(args)

    # calculate consumption of all trips
    schedule.calculate_consumption()

    # Create soc dispatcher
    schedule.init_soc_dispatcher(args)

    # each rotation is assigned a vehicle ID
    schedule.assign_vehicles()

    return schedule, args


def modes_simulation(schedule, scenario, args):
    """
    Run the mode(s) specified in config.

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
                create_results_directory(args, i+1)
                if not args.skip_plots:
                    report.generate_plots(scenario, args)
                logging.info(f"Created plot of failed scenario in {args.results_directory}")
            # continue with other modes after error

    # all modes done
    return schedule, scenario


class Mode:
    """
    Container for simulation modes.

    Each function takes a schedule, the corresponding scenario and arguments Namespace.
    Optionally, an index of the current mode in the modelist can be given.
    A function must return the updated schedule and scenario objects.
    """

    def sim(schedule, scenario, args, _i):# Noqa
        scenario = schedule.run(args)
        return schedule, scenario

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
        # Work on copies of the original schedule and scenario. In case of an exception the outer
        # schedule and scenario stay intact.
        original_schedule = deepcopy(schedule)
        original_scenario = deepcopy(scenario)
        try:
            create_results_directory(args, i+1)
            return run_optimization(conf, sched=schedule, scen=scenario, args=args)
        except Exception as err:
            logging.warning('During Station optimization an error occurred {0}. '
                            'Optimization was skipped'.format(err))
            return original_schedule, original_scenario

    def station_optimization_single_step(schedule, scenario, args, i):
        """ Electrify only the station with the highest potential"""  # noqa
        if not args.optimizer_config:
            logging.warning("Station optimization needs an optimization config file. "
                            "Since no path was given, station optimization is skipped")
            return schedule, scenario
        conf = read_optimizer_config(args.optimizer_config)
        conf.early_return = True
        # Work on copies of the original schedule and scenario. In case of an exception the outer
        # schedule and scenario stay intact.
        original_schedule = deepcopy(schedule)
        original_scenario = deepcopy(scenario)
        try:
            create_results_directory(args, i+1)
            return run_optimization(conf, sched=schedule, scen=scenario, args=args)
        except Exception as err:
            logging.warning('During Station optimization an error occurred {0}. '
                            'Optimization was skipped'.format(err))
            return original_schedule, original_scenario

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

    def report(schedule, scenario, args, i):
        # create report based on all previous modes
        if args.cost_calculation:
            # cost calculation part of report
            calculate_costs(args.cost_parameters, scenario, schedule, args)
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

    prior_reports = sum([m.count('report') for m in args.mode[:i]])
    report_name = f"report_{prior_reports+1}"
    args.results_directory = args.output_directory.joinpath(report_name)
    args.results_directory.mkdir(parents=True, exist_ok=True)
    # save used modes in report version
    used_modes = ['sim'] + [m for m in args.mode[:i] if m not in ['sim', 'report']]
    with open(args.results_directory / "used_modes.txt", "w", encoding='utf-8') as f:
        f.write(f"Used modes in this scenario: {', '.join(used_modes)}")
