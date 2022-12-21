""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
from datetime import datetime, timedelta
import json
import os
import sys
import warnings
import pickle
from copy import copy
from pathlib import Path
from time import time
import logging
import math
import shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


import ebus_toolbox.station_optimizer
from ebus_toolbox.low_soc_events import get_low_soc_events, get_rotation_soc

from ebus_toolbox.consumption import Consumption
from ebus_toolbox.trip import Trip
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox.util import uncomment_json_file, get_buffer_time as get_buffer_time_spice_ev
from optimizer_config import read_config, OptimizerConfig

CONFIG = OptimizerConfig()
matplotlib.use("TkAgg")


def time_it(function, timers={}):
    """decorator function to time the duration function calls
    take and count how often they happen
    :param function: function do be decorated
    :type function: function
    :param timers: storage for cumlulated time and call number
    :type timers: dict
    :return decorated function or timer if given function is None
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


def setup_logger():
    """ setup file and stream logging by config and args arguments
    :return: logger
    :rtype: Logger
    """
    this_logger = logging.getLogger(__name__)
    this_logger.setLevel(logging.DEBUG)

    # logging to one file which keeps track of optimization over many runs
    file_handler_all_opts = logging.FileHandler('optimizer.log')
    file_handler_all_opts.setLevel(CONFIG.debug_level)

    # and logging to a file which is put in the folder with the other optimizer results
    file_handler_this_opt = logging.FileHandler(Path(ARGS.output_directory) / Path('optimizer.log'))
    file_handler_this_opt.setLevel(CONFIG.debug_level)

    formatter = logging.Formatter('%(asctime)s:%(message)s',
                                  "%m%d %H%M%S")

    file_handler_all_opts.setFormatter(formatter)
    file_handler_this_opt.setFormatter(formatter)

    formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(CONFIG.debug_level)
    this_logger.addHandler(file_handler_this_opt)
    this_logger.addHandler(file_handler_all_opts)
    this_logger.addHandler(stream_handler)
    return this_logger


def main():
    """ main call"""
    config_path = "./data/examples/optimizer.cfg"
    run_optimization(config_path)
    LOGGER.debug(time_it(None))


def prepare_filesystem():
    """ Prepare files and folders in the optimization results folder"""
    now = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    ARGS.output_directory = Path(CONFIG.output_path) / str(now + "_optimizer")
    # create subfolder for specific sim results with timestamp.
    # if folder doesnt exists, create folder.
    ARGS.output_directory.mkdir(parents=True, exist_ok=True)

    # copy paste the config file
    copy_file = Path(CONFIG.path)
    destination = ARGS.output_directory / Path("optimizer_config.cfg")
    shutil.copy(copy_file, destination)


@time_it
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
    global CONFIG
    global LOGGER
    global ARGS
    CONFIG = read_config(config_path)

    # load pickle files
    if sched is None or scen is None or this_args is None:
        # if no schedule was given as argument, make sure no scenario
        # and args input was given as well.
        assert sched == scen == this_args is None
        sched, scen, this_args = toolbox_from_pickle(CONFIG.schedule, CONFIG.scenario, CONFIG.args)
    ARGS = this_args

    # prepare Filesystem with folders and paths and copy config
    prepare_filesystem()

    if ARGS.save_soc:
        ARGS.save_soc = ARGS.output_directory / "simulation_soc_spiceEV.csv"
    new_ele_stations_path = ARGS.output_directory / Path("optimized_stations" + ".json")

    LOGGER = setup_logger()

    # remove those args, since they lead to file creation, which is not
    # needed.
    del ARGS.save_timeseries
    del ARGS.save_results
    if ARGS.desired_soc_deps != 1 and CONFIG.opt_type == "quick":
        LOGGER.error("Fast calc is not yet optimized for desired socs unequal to 1")

    optimizer = ebus_toolbox.station_optimizer.StationOptimizer(sched, scen, ARGS, CONFIG, LOGGER)


    # set battery and charging curves through config file if wished for
    optimizer.set_battery_and_charging_curves()

    # create a decision tree or load one from a previous run
    optimizer.set_up_decision_tree()

    # rebasing the scenario meaning simulating it again with the given conditions of
    # included and excluded stations and rotations
    if CONFIG.rebase_scenario:
        must_include_set, ele_stations = optimizer.rebase_spice_ev()
    else:
        must_include_set, ele_stations = optimizer.rebase_simple()

    # create charging dicts which contain soc over time, which is numerically calculated
    optimizer.create_charging_curves()

    # Remove none Values from socs in the vehicle_socs an
    optimizer.remove_none_socs()
    if CONFIG.remove_impossible_rots:
        neg_rots = optimizer.get_negative_rotations_all_electrified()
        optimizer.config.exclusion_rots.update(neg_rots)
        optimizer.schedule.rotations = {r: optimizer.schedule.rotations[r]
                                        for r in optimizer.schedule.rotations if
                                        r not in optimizer.config.exclusion_rots}
        LOGGER.warning("%s negative rotations %s", len(neg_rots), neg_rots)




    if CONFIG.check_for_must_stations:
        must_stations = optimizer.get_must_stations_and_rebase(relative_soc=False)
        LOGGER.warning("%s must stations %s", len(must_stations), must_stations)

    LOGGER.debug("Starting greedy optimization")
    print_time()
    ele_stations, ele_station_set = optimizer.loop()

    print_time()
    ele_station_set = ele_station_set.union(must_include_set)
    LOGGER.debug("Functions took these times in seconds")
    LOGGER.debug("%s electrified stations : %s", len(ele_station_set), ele_station_set)
    LOGGER.debug("%s total stations", len(ele_stations))
    LOGGER.debug(optimizer.could_not_be_electrified)

    scen.vehicle_socs = optimizer.timeseries_calc()

    new_events = get_low_soc_events(optimizer, soc_data=scen.vehicle_socs)

    if len(new_events) > 0:
        LOGGER.debug("Still not electrified with abs. soc with fast calc")
        for event in new_events:
            LOGGER.debug(event.rotation.id)
        LOGGER.debug("#####")

    with open(new_ele_stations_path, "w", encoding="utf-8", ) as file:
        json.dump(ele_stations, file, indent=2)

    LOGGER.debug("Spice EV is calculating optimized case as a complete scenario")
    _, __ = optimizer.preprocessing_scenario(
        electrified_stations=ele_stations, run_only_neg=False, cost_calc=True)

    LOGGER.warning("Still negative rotations: %s", optimizer.schedule.
                   get_negative_rotations(optimizer.scenario))
    print("Finished")
    return optimizer.schedule, optimizer.scenario


@time_it
def print_time(start=[]):
    """ Print the time and automatically set start time to the first time the function getting
    called"""
    if not start:
        start.append(time())
    print(round(time() - start[0], 2), " seconds till start")


@time_it
def plot_(data):
    """ Simple plot of data without having to create subplots"""
    fig, ax = plt.subplots()
    ax.plot(data, linewidth=2.0)
    return ax


def plot_rot(rot_id, this_sched, this_scen, ax=None, rot_only=True):
    """ Simple plot of data without having to create subplots"""
    soc, start, end = get_rotation_soc(rot_id, this_sched, this_scen)
    if not rot_only:
        start = 0
        end = -1
    if ax is None:
        fig, ax = plt.subplots()
        ax.plot(soc[start:end], linewidth=2.0)
        return ax
    ax.plot(soc[start:end], linewidth=2.0)
    return ax



@time_it
def preprocess_schedule(this_sched, this_args, electrified_stations=None):
    Trip.consumption = Consumption(this_sched.vehicle_types,
                                   outside_temperatures=this_args.outside_temperature_over_day_path,
                                   level_of_loading_over_day=this_args.level_of_loading_over_day_path)

    this_sched.stations = electrified_stations
    this_sched.calculate_consumption()
    this_sched.assign_vehicles()

    return this_sched, this_sched.generate_scenario(this_args)


@time_it
def run_schedule(this_sched, this_args, electrified_stations=None, cost_calc=False):
    """Run a given schedule and electrify stations if need be"""
    this_sched2 = copy(this_sched)
    this_sched2.stations = electrified_stations
    this_sched2, new_scen = preprocess_schedule(this_sched2, this_args,
                                                electrified_stations=electrified_stations)
    # Dont print output from spice ev to reduce clutter
    sys.stdout = open(os.devnull, 'w')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        new_scen.run('distributed', vars(this_args).copy())
    sys.stdout = sys.__stdout__
    if this_args.cost_calculation and cost_calc:
        # cost calculation following directly after simulation
        try:
            with open(this_args.cost_parameters_file, encoding='utf-8') as f:
                cost_parameters_file = uncomment_json_file(f)
        except FileNotFoundError:
            raise SystemExit(f"Path to cost parameters ({ARGS.cost_parameters_file}) "
                             "does not exist. Exiting...")
        calculate_costs(cost_parameters_file, new_scen, this_sched2, this_args)
    return this_sched2, new_scen


@time_it
def combs_unordered_no_putting_back(n: int, k: int):
    """ Returns amount of combinations for pulling k elements out of n, without putting elements
    back or looking at the order. this is equal to n over k
    :param n: number of elements in the base group
    :type n: int
    :param k: number of elements in the sub group of picked elements
    :type k: int
    :return: number of combinations
    :rtype: int
    """
    try:
        return math.factorial(n) / ((math.factorial(n - k)) * math.factorial(k))
    except ValueError:
        LOGGER.warning("Value Error")
        return 0





def toolbox_from_pickle(sched_name, scen_name, args_name):
    """ Load the 3 files from pickle"""
    with open(args_name, "rb") as f:
        this_args = pickle.load(f)
    with open(scen_name, "rb") as f:
        scen = pickle.load(f)
    with open(sched_name, "rb") as f:
        sched = pickle.load(f)
    return sched, scen, this_args

if __name__ == "__main__":
    main()
