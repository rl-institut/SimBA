""" Try to minimize the amount of electrified stations to achieve full electrification."""
from copy import deepcopy
import json
import sys
from pathlib import Path
import logging
import simba.station_optimizer
from simba.station_optimizer import opt_util
from spice_ev.report import generate_soc_timeseries

config = opt_util.OptimizerConfig()


def setup_logger(conf):
    """ Setup file and stream logging by config and args arguments.

    :param conf: configuration object
    :type conf: simba.optimizer_util.OptimizerConfig
    :return: logger
    :rtype: Logger
    """
    this_logger = logging.getLogger(__name__)
    this_logger.setLevel(logging.DEBUG)

    # logging to one file which keeps track of optimization over many runs

    # use the temporary folder in case of testing. In the use case the general log file with logging
    # for all optimizations is located 2 folders above the opt folder, e.g. data/sim_outputs for the
    # directory data/sim_outputs/XXX_SimBA_results/sim__station_optimization/
    if "pytest" in sys.modules:
        general_logger_path = Path(conf.optimizer_output_dir) / "optimizer.log"
    else:
        general_logger_path = Path(conf.optimizer_output_dir) / "../../optimizer.log"
    file_handler_all_opts = logging.FileHandler(general_logger_path)

    file_handler_all_opts.setLevel(conf.debug_level)

    # and logging to a file which is put in the folder with the other optimizer results
    file_handler_this_opt = logging.FileHandler(Path(conf.optimizer_output_dir) /
                                                Path('optimizer.log'))
    file_handler_this_opt.setLevel(conf.debug_level)

    formatter = logging.Formatter('%(asctime)s:%(message)s',
                                  "%m%d %H%M%S")

    file_handler_all_opts.setFormatter(formatter)
    file_handler_this_opt.setFormatter(formatter)

    formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(conf.console_level)

    # Log to an optimization-specific file, a general file and to console
    this_logger.addHandler(file_handler_this_opt)
    this_logger.addHandler(file_handler_all_opts)
    this_logger.addHandler(stream_handler)
    return this_logger


def prepare_filesystem(args, conf):
    """ Prepare files and folders in the optimization results folder.

    :param conf: configuration
    :type conf:  simba.optimizer_util.OptimizerConfig
    :param args: arguments for SimBA toolbox
    :type args: Namespace
    """
    conf.optimizer_output_dir = Path(args.results_directory)
    # create sub folder for specific sim results with timestamp.
    # if folder doesnt exists, create folder.
    conf.optimizer_output_dir.mkdir(parents=True, exist_ok=True)


def run_optimization(conf, sched=None, scen=None, args=None):
    """ Add electrified stations until there are no more negative rotations.

    Configured with arguments from optimizer config file.

    :param conf: Configuration object of optimization
    :type conf: OptimizerConfig
    :param sched: Simulation schedule containing buses, rotations etc.
    :type sched: simba.schedule.Schedule
    :param scen: Simulation scenario containing simulation results
                     including the SoC of all vehicles over time
    :type scen: spice_ev.Scenario
    :param args: Simulation arguments for manipulation of generated outputs
    :type args: Namespace
    :raises Exception: if no rotations can be optimized

    :return: optimized schedule and Scenario
    :rtype: tuple(simba.schedule.Schedule, spice_ev.Scenario)
    """

    # load pickle files if they are given in the optimizer.config
    if conf.schedule:
        # either all optional arguments are given or none are
        error_message = ("To optimize from .pickle files, schedule, scenario and arguments need to "
                         "be provided together")
        assert conf.scenario, error_message
        assert conf.args, error_message
        sched, scen, args = opt_util.toolbox_from_pickle(conf.schedule, conf.scenario, conf.args)

    original_schedule = deepcopy(sched)

    # setup folders, paths and copy config
    prepare_filesystem(args, conf)

    new_ele_stations_path = conf.optimizer_output_dir / "optimized_stations.json"

    logger = setup_logger(conf)

    if args.desired_soc_deps != 1 and conf.solver == "quick":
        logger.error("Fast calculation is not yet optimized for desired socs different to 1")
    optimizer = simba.station_optimizer.StationOptimizer(sched, scen, args, conf, logger)

    # set battery and charging curves through config file
    optimizer.set_battery_and_charging_curves()

    # filter out depot chargers if option is set
    if conf.run_only_oppb:
        optimizer.config.exclusion_rots = optimizer.config.exclusion_rots.union(
            r for r in sched.rotations if "depb" == sched.rotations[r].charging_type)
        sched.rotations = {r: sched.rotations[r] for r in sched.rotations
                           if "oppb" == sched.rotations[r].charging_type}
        if len(sched.rotations) < 1:
            raise Exception("No rotations left after removing depot chargers")

    # rebasing the scenario meaning simulating it again with SpiceEV and the given conditions of
    # included stations, excluded stations, filtered rotations and changed battery sizes
    if conf.rebase_scenario:
        must_include_set, ele_stations = optimizer.rebase_spice_ev()
    else:
        # no new SpiceEV simulation will take place, but some variables need to be adjusted
        must_include_set, ele_stations = optimizer.rebase_simple()
        # make sure vehicle socs exist
        generate_soc_timeseries(scen)

    # create charging dicts which contain soc over time, which is numerically calculated
    optimizer.create_charging_curves()

    # remove none values from socs in the vehicle_socs
    optimizer.replace_socs_from_none_to_value()

    # Remove already electrified stations from possible stations
    optimizer.not_possible_stations = set(optimizer.electrified_stations.keys()).union(
        optimizer.not_possible_stations)

    # all stations electrified: are there still negative rotations?
    if conf.remove_impossible_rotations:
        neg_rots = optimizer.get_negative_rotations_all_electrified()
        optimizer.config.exclusion_rots.update(neg_rots)
        optimizer.schedule.rotations = {
            r: optimizer.schedule.rotations[r] for r in optimizer.schedule.rotations if
            r not in optimizer.config.exclusion_rots}

        logger.warning(f"{len(neg_rots)} negative rotations {neg_rots} were removed from schedule "
                       "because they cannot be electrified")
        assert len(optimizer.schedule.rotations) > 0, (
            "Schedule cannot be optimized, since rotations cannot be electrified.")

    # if the whole network can not be fully electrified if even just a single station is not
    # electrified, this station must be included in a fully electrified network
    # this can make solving networks much simpler. Some information during the path to full
    # electrification gets lost though, e.g. priority of stations.
    if conf.check_for_must_stations:
        must_stations = optimizer.get_critical_stations_and_rebase(relative_soc=False)
        logger.warning("%s must stations %s", len(must_stations), must_stations)

    # Store the rotations to be optimized. optimizer.loop() will mutate the dictionary but not the
    # rotations itself.
    rotations_for_opt = optimizer.schedule.rotations.copy()

    logger.log(msg="Starting greedy station optimization", level=100)

    # start a timer to later check how long the optimization took
    opt_util.get_time()

    # Go into the optimization loop, where stations are subsequently electrified
    ele_stations, ele_station_set = optimizer.loop()

    ele_station_set = ele_station_set.union(must_include_set)
    logger.debug("%s electrified stations : %s", len(ele_station_set), ele_station_set)
    logger.debug("%s total stations", len(ele_stations))
    logger.debug("These rotations could not be electrified: %s", optimizer.could_not_be_electrified)

    # remove none values from socs in the vehicle_socs so timeseries_calc can work
    optimizer.replace_socs_from_none_to_value()

    # Restore the rotations and the scenario which where the goal of optimization.
    optimizer.schedule.rotations = rotations_for_opt
    optimizer.scenario = optimizer.base_scenario

    # Check if the rotations which were part of the optimization are not negative anymore
    vehicle_socs = optimizer.timeseries_calc(optimizer.electrified_station_set)

    new_events = optimizer.get_low_soc_events(soc_data=vehicle_socs)
    if len(new_events) > 0:
        logger.debug("Estimation of network still shows negative rotations")
        for event in new_events:
            logger.debug(event.rotation.id)

    # Change Paths to datatype String, since Paths are not json serializable
    with open(new_ele_stations_path, "w", encoding="utf-8") as file:
        output_dict = {key: value for key, value in ele_stations.items()}
        opt_util.recursive_dict_updater(
            output_dict, lambda key, value: isinstance(value, Path), lambda key, value: str(value))
        json.dump(output_dict, file, ensure_ascii=False, indent=2)

    # Calculation with SpiceEV is more accurate and will show if the optimization is viable or not
    logger.debug("Detailed calculation of an optimized case as a complete scenario")

    # Restore original rotations
    for rotation_id in original_schedule.rotations:
        optimizer.schedule.rotations[rotation_id] = original_schedule.rotations[rotation_id]

    # remove exclusion since internally these would not be simulated
    optimizer.config.exclusion_rots = set()
    _, __ = optimizer.preprocessing_scenario(
        electrified_stations=ele_stations, run_only_neg=False)
    neg_rotations = optimizer.schedule.get_negative_rotations(optimizer.scenario)
    if len(neg_rotations) > 0:
        logger.log(msg=f"Still {len(neg_rotations)} negative rotations: {neg_rotations}", level=100)
    logger.log(msg="Station optimization finished after " + opt_util.get_time(), level=100)

    return optimizer.schedule, optimizer.scenario
