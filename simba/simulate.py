import logging
import dill as pickle
import traceback
from copy import deepcopy
from pathlib import Path

from simba import report, optimization, optimizer_util, util
from simba.data_container import DataContainer
from simba.costs import calculate_costs
from simba.optimizer_util import read_config as read_optimizer_config
from simba.schedule import Schedule
from simba.station_optimization import run_optimization


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
    if vars(args).get("load_pickle"):
        # load pickle file: skip pre_simulation
        if isinstance(args.mode, list):
            first_mode = args.mode[0]
        else:
            first_mode = args.mode
        assert first_mode == "load_pickle", "Load pickle: first mode must be load_pickle"
        # schedule and scenario read out from pickle file in first mode
        # DataContainer is part of schedule
        schedule = None
        scenario = "pickle"  # must not be None
    else:
        # DataContainer stores various input data
        data_container = DataContainer().fill_with_args(args)
        schedule, args = pre_simulation(args, data_container)
        scenario = schedule.run(args)
    return modes_simulation(schedule, scenario, args)


def pre_simulation(args, data_container: DataContainer):
    """
    Prepare simulation.

    Read in files, generate consumption info, create and filter schedule.

    :param args: arguments
    :type args: Namespace
    :param data_container: data needed for simulation
    :type data_container: DataContainer
    :return: schedule, args
    :rtype: simba.schedule.Schedule, Namespace
    """

    # Deepcopy args so original args do not get mutated
    args = deepcopy(args)

    # generate schedule from csv
    schedule = Schedule.from_datacontainer(data_container, args)

    # filter rotations
    schedule.rotation_filter(args)

    # calculate consumption of all trips
    schedule.calculate_consumption()

    # Create soc dispatcher
    schedule.init_soc_dispatcher(args)

    # each rotation is assigned a vehicle ID
    schedule.assign_vehicles(args)

    return schedule, args


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
                if args.output_path is None:
                    create_results_directory(args, i+1)
                    if not args.skip_plots:
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
    @staticmethod
    def sim(schedule, scenario, args, _i):# Noqa
        # Base simulation function for external access.
        # No effect when used directly in SimBA"
        scenario = schedule.run(args, mode="distributed")
        return schedule, scenario

    @staticmethod
    def sim_greedy(schedule, scenario, args, _i):# Noqa
        # Run a basic greedy simulation without depb/oppb distinction
        scenario = schedule.run(args, mode="greedy")
        return schedule, scenario

    @staticmethod
    def service_optimization(schedule, scenario, args, _i):
        """ Find largest set of rotations that produce no negative SoC.
        """ # noqa
        result = optimization.service_optimization(schedule, scenario, args)
        schedule, scenario = result['optimized']
        if scenario is None:
            logging.warning('No optimization possible (all rotations negative), reverting')
            schedule, scenario = result['original']
        return schedule, scenario

    @staticmethod
    def neg_depb_to_oppb(schedule, scenario, args, _i):
        return Mode.switch_type(schedule, scenario, args, "depb", "oppb")

    @staticmethod
    def neg_oppb_to_depb(schedule, scenario, args, _i):
        return Mode.switch_type(schedule, scenario, args, "oppb", "depb")

    @staticmethod
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

    @staticmethod
    def _station_optimization(schedule, scenario, args, i, single_step: bool):
        if not args.optimizer_config_path:
            logging.warning("Station optimization needs an optimization config file. "
                            "Default Config is used.")
            conf = optimizer_util.OptimizerConfig()
        else:
            conf = read_optimizer_config(args.optimizer_config_path)
            util.save_input_file(args.optimizer_config_path, args)
        if single_step:
            conf.early_return = True
        # Work on copies of the original schedule and scenario. In case of an exception the outer
        # schedule and scenario stay intact.
        original_schedule = deepcopy(schedule)
        original_scenario = deepcopy(scenario)
        try:
            create_results_directory(args, i + 1)
            return run_optimization(conf, sched=schedule, scen=scenario, args=args)
        except Exception as err:
            logging.error('During Station optimization, an error occurred {0}. '
                          'Optimization was skipped'.format(err))
            return original_schedule, original_scenario

    @staticmethod
    def station_optimization(schedule, scenario, args, i):
        return Mode._station_optimization(schedule, scenario, args, i, single_step=False)

    @staticmethod
    def station_optimization_single_step(schedule, scenario, args, i):
        """ Electrify only the station with the highest potential

        :param schedule: Schedule
        :type schedule: simba.schedule.Schedule
        :param scenario: Scenario
        :type scenario: spice_ev.scenario.Scenario
        :param args: Arguments
        :type args: argparse.Namespace
        :param i: counter of modes for directory creation
        :return: schedule, scenario

        """  # noqa
        return Mode._station_optimization(schedule, scenario, args, i, single_step=True)

    @staticmethod
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

    @staticmethod
    def split_negative_depb(schedule, scenario, args, _i):
        negative_rotations = schedule.get_negative_rotations(scenario)
        trips, depot_trips = optimization.prepare_trips(schedule, negative_rotations)
        recombined_schedule = optimization.recombination(schedule, args, trips, depot_trips)
        # re-run schedule
        scenario = recombined_schedule.run(args)
        return recombined_schedule, scenario

    @staticmethod
    def load_pickle(_schedule=None, _scenario=None, args=None, _i=None):
        with open(args.load_pickle, 'rb') as f:
            unpickle = pickle.load(f)
        schedule = unpickle["schedule"]
        scenario = unpickle["scenario"]
        # DataContainer is part of schedule
        # However, cost parameters are supposed to be mutable after loading from pickle
        if args.cost_parameters_path:
            schedule.data_container.add_cost_parameters_from_json(args.cost_parameters_path)

        return schedule, scenario

    @staticmethod
    def report(schedule, scenario, args, i):
        if args.output_path is None:
            logging.warning("No output path given, skipping report")
            return schedule, scenario

        # create report based on all previous modes
        if args.cost_calculation:
            # cost calculation part of report
            try:
                cost_parameters = schedule.data_container.cost_parameters_data
                scenario.costs = calculate_costs(cost_parameters, scenario, schedule, args)
            except Exception:
                logging.warning(f"Cost calculation failed due to {traceback.format_exc()}")
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

    if args.output_path is None:
        return

    prior_reports = sum([m.count('report') for m in args.mode[:i]])
    report_name = f"report_{prior_reports+1}"
    args.results_directory = Path(args.output_path).joinpath(report_name)
    args.results_directory.mkdir(parents=True, exist_ok=True)
    # save used modes in report version
    used_modes = ['sim'] + [m for m in args.mode[:i] if m not in ['sim', 'report']]
    file_path = args.results_directory / "used_modes.txt"
    if vars(args).get("scenario_name"):
        file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)
    with open(file_path, "w", encoding='utf-8') as f:
        f.write(f"Used modes in this scenario: {', '.join(used_modes)}")
