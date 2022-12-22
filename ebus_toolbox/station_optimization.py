""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
from datetime import datetime
import json
from pathlib import Path
from time import time
import logging
import shutil
import matplotlib

import ebus_toolbox.station_optimizer
import ebus_toolbox.optimizer_util as opt_util

config = opt_util.OptimizerConfig()
matplotlib.use("TkAgg")


def time_it(function, timers={}):
    """decorator function to time the duration function calls
    take and count how often they happen
    :param function: function do be decorated
    :type function: function
    :param timers: storage for cumulated time and call number
    :type timers: dict
    :return: decorated function or timer if given function is None
    :rtype function or dict

    """
    if function:
        def decorated_function(*this_args, **kwargs):
            key = function.__name__
            start_time = time()
            return_value = function(*this_args, **kwargs)
            delta_time = time() - start_time
            try:
                timers[key]["time"] += delta_time
                timers[key]["calls"] += 1
            except KeyError:
                timers[key] = dict(time=0, calls=1)
                timers[key]["time"] += delta_time
            return return_value
        return decorated_function

    sorted_timer = dict(sorted(timers.items(), key=lambda x: x[1]["time"] / x[1]["calls"]))
    return sorted_timer


def setup_logger(this_args, conf):
    """ setup file and stream logging by config and args arguments
    :param conf: configuration object
    :param this_args: Namespace object of arguments for ebus toolbox
    :return: logger
    :rtype: Logger
    """
    this_logger = logging.getLogger(__name__)
    this_logger.setLevel(logging.DEBUG)

    # logging to one file which keeps track of optimization over many runs
    file_handler_all_opts = logging.FileHandler('optimizer.log')
    file_handler_all_opts.setLevel(conf.debug_level)

    # and logging to a file which is put in the folder with the other optimizer results
    file_handler_this_opt = logging.FileHandler(Path(this_args.output_directory) /
                                                Path('optimizer.log'))
    file_handler_this_opt.setLevel(conf.debug_level)

    formatter = logging.Formatter('%(asctime)s:%(message)s',
                                  "%m%d %H%M%S")

    file_handler_all_opts.setFormatter(formatter)
    file_handler_this_opt.setFormatter(formatter)

    formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(conf.debug_level)
    this_logger.addHandler(file_handler_this_opt)
    this_logger.addHandler(file_handler_all_opts)
    this_logger.addHandler(stream_handler)
    return this_logger


def main():
    """ main call"""
    config_path = "./data/examples/optimizer.cfg"
    run_optimization(config_path)


def prepare_filesystem(this_args, conf):
    """ Prepare files and folders in the optimization results folder
    :param conf: configuration object
    :param this_args: Namespace object of arguments for ebus toolbox
    """
    now = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    this_args.output_directory = Path(conf.output_path) / str(now + "_optimizer")
    # create sub folder for specific sim results with timestamp.
    # if folder doesnt exists, create folder.
    this_args.output_directory.mkdir(parents=True, exist_ok=True)

    # copy paste the config file
    copy_file = Path(conf.path)
    destination = this_args.output_directory / Path("optimizer_config.cfg")
    shutil.copy(copy_file, destination)


def run_optimization(config_path, sched=None, scen=None, this_args=None):
    """    Optimizes scenario by adding electrified stations sparingly
    until scenario has no below 0 soc events.

    The optimizer config file is used for setting the optimization up.

    :param config_path: path to optimizer.cfg file.
    :type config_path: str

    :param sched: Simulation schedule containing buses, rotations etc.
    :type sched: ebus_toolbox.Schedule

    :param scen: Simulation scenario containing simulation results
                     including the SoC of all vehicles over time
    :type scen: spice_ev.Scenario
    :param this_args: Simulation arguments for manipulation or generated outputs
    :type this_args: object

    :return: (Schedule,Scenario) optimized schedule and Scenario
    :rtype: tuple(ebus_toolbox.Schedule, spice_ev.Scenario)
    """
    conf = opt_util.read_config(config_path)

    # load pickle files
    if sched is None or scen is None or this_args is None:
        # if no schedule was given as argument, make sure no scenario
        # and args input was given as well.
        assert sched == scen == this_args is None
        sched, scen, this_args = opt_util.toolbox_from_pickle(conf.schedule,
                                                              conf.scenario, conf.args)
    args = this_args

    # prepare Filesystem with folders and paths and copy config
    prepare_filesystem(args, conf)

    if args.save_soc:
        args.save_soc = args.output_directory / "simulation_soc_spiceEV.csv"
    new_ele_stations_path = args.output_directory / Path("optimized_stations" + ".json")

    logger = setup_logger(args, conf)

    # remove those args, since they lead to file creation, which is not
    # needed.
    del args.save_timeseries
    del args.save_results
    if args.desired_soc_deps != 1 and conf.opt_type == "quick":
        logger.error("Fast calc is not yet optimized for desired socs unequal to 1")

    optimizer = ebus_toolbox.station_optimizer.StationOptimizer(sched, scen, args, conf, logger)

    # set battery and charging curves through config file if wished for
    optimizer.set_battery_and_charging_curves()

    # create a decision tree or load one from a previous run
    optimizer.set_up_decision_tree()

    # rebasing the scenario meaning simulating it again with the given conditions of
    # included and excluded stations and rotations
    if conf.rebase_scenario:
        must_include_set, ele_stations = optimizer.rebase_spice_ev()
    else:
        must_include_set, ele_stations = optimizer.rebase_simple()

    # create charging dicts which contain soc over time, which is numerically calculated
    optimizer.create_charging_curves()

    # Remove none Values from socs in the vehicle_socs an
    optimizer.remove_none_socs()
    if conf.remove_impossible_rots:
        neg_rots = optimizer.get_negative_rotations_all_electrified()
        optimizer.config.exclusion_rots.update(neg_rots)
        optimizer.schedule.rotations = {r: optimizer.schedule.rotations[r]
                                        for r in optimizer.schedule.rotations if
                                        r not in optimizer.config.exclusion_rots}
        logger.warning("%s negative rotations %s", len(neg_rots), neg_rots)

    if conf.check_for_must_stations:
        must_stations = optimizer.get_must_stations_and_rebase(relative_soc=False)
        logger.warning("%s must stations %s", len(must_stations), must_stations)

    logger.debug("Starting greedy optimization")
    ele_stations, ele_station_set = optimizer.loop()
    ele_station_set = ele_station_set.union(must_include_set)
    logger.debug("%s electrified stations : %s", len(ele_station_set), ele_station_set)
    logger.debug("%s total stations", len(ele_stations))
    logger.debug("These rotations could not be electrified: %s", optimizer.could_not_be_electrified)

    scen.vehicle_socs = optimizer.timeseries_calc()

    new_events = optimizer.get_low_soc_events(soc_data=scen.vehicle_socs)

    if len(new_events) > 0:
        logger.debug("Still not electrified with abs. soc with fast calc")
        for event in new_events:
            logger.debug(event.rotation.id)
    with open(new_ele_stations_path, "w", encoding="utf-8", ) as file:
        json.dump(ele_stations, file, indent=2)

    logger.debug("Spice EV is calculating optimized case as a complete scenario")
    _, __ = optimizer.preprocessing_scenario(
        electrified_stations=ele_stations, run_only_neg=False, cost_calc=True)

    logger.warning("Still negative rotations: %s", optimizer.schedule.
                   get_negative_rotations(optimizer.scenario))
    print("Finished")
    logger.debug(time_it(None))

    return optimizer.schedule, optimizer.scenario


if __name__ == "__main__":
    main()
