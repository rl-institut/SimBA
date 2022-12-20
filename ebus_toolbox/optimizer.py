""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
from datetime import datetime, timedelta
import json
import os
import sys
import warnings
import pickle
from copy import copy, deepcopy
from pathlib import Path
from time import time
import logging
import math
import shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import report
from ebus_toolbox.schedule_optimizer import ScheduleOptimizer
from src import scenario
import schedule
import rotation
from optimizer_config import read_config, OptimizerConfig
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
    if ARGS.desired_soc_deps != 1 and CONFIG.opt_type=="quick":
        LOGGER.error("Fast calc is not yet optimized for desired socs unequal to 1")

    optimizer = ScheduleOptimizer(sched, scen, ARGS, CONFIG, LOGGER)

    remove_impossible_rots = CONFIG.remove_impossible_rots

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

    soc_charge_curve_dict = optimizer.create_charging_curves(sched)

    # Remove none Values from socs in the vehicle_socs an
    optimizer.remove_none_socs()

    if CONFIG.check_for_must_stations:
        must_stations = optimizer.get_must_stations_and_rebase(relative_soc=False)
        LOGGER.warning("%s must stations %s", len(must_stations), must_stations)

    LOGGER.debug("Starting greedy optimization")
    ele_station_set = set()
    print_time()
    ele_stations, ele_station_set, could_not_be_electrified, list_greedy_sets =\
        \
        optimization_loop(ele_stations, ele_station_set, scen, sched,
                          not_possible_stations, soc_upper_thresh=ARGS.desired_soc_opps,
                          soc_lower_thresh=CONFIG.min_soc,
                          solver=CONFIG.solver, opt_type=CONFIG.opt_type,
                          node_choice=CONFIG.node_choice,
                          soc_charge_curve_dict=soc_charge_curve_dict,
                          decision_tree=decision_tree)
    print_time()
    ele_station_set = ele_station_set.union(must_include_set)
    LOGGER.debug("Functions took these times in seconds")
    LOGGER.debug("%s Stations : %s", len(ele_station_set), ele_station_set)
    LOGGER.debug(could_not_be_electrified)
    scen.vehicle_socs = timeseries_calc(sched.rotations.values(),
                                        scen.vehicle_socs,
                                        scen, ele_station_set,
                                        soc_charge_curve_dict=soc_charge_curve_dict,
                                        soc_upper_thresh=ARGS.desired_soc_deps)
    new_events = get_low_soc_events(scen, sched)
    if len(new_events) > 0:
        LOGGER.debug("Still not electrified with fast calc")
        for event in new_events:
            LOGGER.debug(event["rotation"].id)
        LOGGER.debug("#####")

    with open(new_ele_stations_path, "w", encoding="utf-8", ) as file:
        json.dump(ele_stations, file, indent=2)

    LOGGER.debug("Spice EV is calculating optimized case as a complete scenario")
    final_sched, final_scen, ele_station_set, ele_stations = preprocessing_scenario(
        sched, scen, ARGS, electrified_stations=ele_stations, run_only_neg=False,
        electrified_station_set=ele_station_set, cost_calc=True)

    LOGGER.warning("Still negative rotations: %s", final_sched.get_negative_rotations(final_scen))
    print("Finished")
    return final_sched, final_scen





@time_it
def get_groups_from_events(events, not_possible_stations=None, could_not_be_electrified=None):
    """ First create simple list of station sets for single events.
    #Electrified and other not possible to electrify stations should not connect groups"""

    # Making sure default arguments are none and not mutable
    if not not_possible_stations:
        not_possible_stations = set()

    if not could_not_be_electrified:
        could_not_be_electrified = set()

    possible_stations = [
        {station for station in event["stations"] if station not in not_possible_stations}
        for event
        in events]
    # If stations overlap join them
    station_subsets = join_all_subsets(possible_stations)
    event_groups = [[] for __ in range(len(station_subsets))]

    # Group the events in the same manner as stations, so that the same amount of event groups are
    # created as station subsets
    for event in events:
        for i, subset in enumerate(station_subsets):
            # Every station from an event must share the same station sub_set. Therefore its enough
            # to check only the first element
            if len(event["stations"]) > 0:
                if next(iter(event["stations"])) in subset:
                    event_groups[i].append(event)
                    break
        else:
            LOGGER.warning('Didnt find rotation %sin any subset'
                           'of possible electrifiable stations', event["rotation"].id)
            # this event will no show up in an event_group.
            # therefore it needs to be put into this set
            could_not_be_electrified.update([event["rotation"].id])

    groups = list(zip(event_groups, station_subsets))
    return sorted(groups, key=lambda x: len(x[1]))


@time_it





def print_time(start=[]):
    """ Print the time and automatically set start time to the first time the function getting
    called"""
    if not start:
        start.append(time())
    print(round(time() - start[0], 2), " seconds till start")


@time_it
def node_to_tree(decision_tree, electrified_station_set, delta_base_energy):
    """ Fill decision tree with the given info
    :return decision tree
    :rtype dict()
    """
    node_name = stations_hash(electrified_station_set)
    try:
        decision_tree[node_name]["missing_energy"] = delta_base_energy
        decision_tree[node_name]["visit_counter"] += 1
        # todo add is_viable
        LOGGER.debug("already visited")
    except KeyError:
        decision_tree[node_name] = {}
        decision_tree[node_name]["missing_energy"] = delta_base_energy
        decision_tree[node_name]["visit_counter"] = 1

    return decision_tree


@time_it



@time_it
def combination_generator(iterable, amount: int):
    """ Generator which yields all possible combinations of choosing
    an amount out of an iterable without putting them back and without caring about the
    order of elements
    :param iterable: Any collection which can be cast to a list
    :param amount: Number of elements which should be drawn from iterable
    :type amount: int
    """
    iterable = list(iterable)

    for i, item in enumerate(iterable):
        # Recursive calling of generator with clock like behavior, e.g right-most item changes until
        # end of list is reached. This leads to a change in the item left to it and so on. Elements
        # on the right can only change to a subset of bigger indicies than their left counter-part.
        # This is due to the ignoring of order, which reduces the amount of possibilities.
        if amount <= 1:
            yield [item]
        else:
            for gen in combination_generator(iterable[i + 1:], amount - 1):
                yield [item] + gen


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


def join_all_subsets(subsets):
    """ join sets for as long as needed until no elements share any intersections"""
    joined_subset = True
    while joined_subset:
        joined_subset, subsets = join_subsets(subsets)
    return subsets




@time_it
def timeseries_calc(rotations, soc_dict, eval_scen, ele_station_set,
                    soc_charge_curve_dict, soc_upper_thresh=1, electrify_stations=None):
    """ A quick estimation of socs by mutating the soc data accordingly
    :return soc dict
    :rtype dict()
    """
    if not electrify_stations:
        electrify_stations=set()
    ele_stations = {*ele_station_set, electrify_stations}
    soc_dict = copy(soc_dict)

    for rot in rotations:
        ch_type = (rot.vehicle_id.find("oppb") > 0) * "oppb" + (
                rot.vehicle_id.find("depb") > 0) * "depb"
        v_type = rot.vehicle_id.split("_" + ch_type)[0]
        soc_over_time_curve = soc_charge_curve_dict[v_type][ch_type]
        soc = soc_dict[rot.vehicle_id]
        for i, trip in enumerate(rot.trips):
            if trip.arrival_name not in ele_stations:
                continue
            idx = get_index_by_time(trip.arrival_time, eval_scen)
            try:
                standing_time_min = get_charging_time(trip, rot.trips[i + 1], ARGS)
            except IndexError:
                standing_time_min = 0

            d_soc = get_delta_soc(soc_over_time_curve, soc[idx], standing_time_min)
            buffer_idx = int((get_buffer_time(trip)) / timedelta(minutes=1))
            delta_idx = int(standing_time_min) + 1
            old_soc = soc[idx + buffer_idx:idx + buffer_idx + delta_idx].copy()
            soc[idx + buffer_idx:] += d_soc
            soc[idx + buffer_idx:idx + buffer_idx + delta_idx] = old_soc
            soc[idx + buffer_idx:idx + buffer_idx + delta_idx] += np.linspace(0, d_soc, delta_idx)

            soc_pre = soc[:idx]
            soc = soc[idx:]
            soc_max = np.max(soc)
            while soc_max > soc_upper_thresh:
                desc = np.arange(len(soc), 0, -1)
                diff = np.hstack((np.diff(soc), -1))
                # masking of socs >1 and negative gradient for local maximum
                idc_loc_max = np.argmax(desc * (soc > 1) * (diff < 0))

                soc_max = soc[idc_loc_max]
                # Reducing everything after local maximum
                soc[idc_loc_max:] = soc[idc_loc_max:] - (soc_max - 1)

                # Capping everything before local maximum
                soc[:idc_loc_max][soc[:idc_loc_max] > 1] = 1
                soc_max = np.max(soc)
            soc = np.hstack((soc_pre, soc))
        soc_dict[rot.vehicle_id] = soc
    return soc_dict


@time_it

@time_it
def get_charging_time(trip1, trip2, args):
    standing_time_min = (trip2.departure_time - trip1.arrival_time) \
                        / timedelta(minutes=1)
    buffer_time = (get_buffer_time(trip1) / timedelta(minutes=1))
    standing_time_min -= buffer_time

    if args.min_charging_time > standing_time_min:
        return 0
    return max(0, standing_time_min)


@time_it
def join_subsets(subsets):
    subsets = [s.copy() for s in subsets]
    for i in range(len(subsets)):
        for ii in range(len(subsets)):
            if i == ii:
                continue
            intersec = subsets[i].intersection(subsets[ii])
            if len(intersec) > 0:
                subsets[i] = subsets[i].union(subsets[ii])
                subsets.remove(subsets[ii])
                return True, subsets
    return False, subsets


@time_it
def get_index_by_time(search_time, this_scen: scenario.Scenario):
    start_time = this_scen.start_time
    delta_time = timedelta(minutes=60 / this_scen.stepsPerHour)
    idx = (search_time - start_time) // delta_time
    return idx



@time_it
def get_rotation_soc(rot_id, this_sched, this_scen, soc_data: dict = None):
    rot = this_sched.rotations[rot_id]
    rot_start_idx = get_index_by_time(rot.departure_time, this_scen)
    rot_end_idx = get_index_by_time(rot.arrival_time, this_scen)
    if soc_data:
        return soc_data[rot.vehicle_id], rot_start_idx, rot_end_idx
    return this_scen.vehicle_socs[rot.vehicle_id], rot_start_idx, rot_end_idx


@time_it
def preprocess_schedule(this_sched, this_args, electrified_stations=None):

    Trip.consumption = Consumption(this_sched.vehicle_types,
                                   outside_temperatures=ARGS.outside_temperature_over_day_path,
                                   level_of_loading_over_day=ARGS.level_of_loading_over_day_path)

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
    if ARGS.cost_calculation and cost_calc:
        # cost calculation following directly after simulation
        try:
            with open(ARGS.cost_parameters_file, encoding='utf-8') as f:
                cost_parameters_file = uncomment_json_file(f)
        except FileNotFoundError:
            raise SystemExit(f"Path to cost parameters ({ARGS.cost_parameters_file}) "
                             "does not exist. Exiting...")
        calculate_costs(cost_parameters_file, new_scen, this_sched2, ARGS)
    return this_sched2, new_scen


@time_it
def charging_curve_to_soc_over_time(charging_curve, capacity, max_charge_from_grid=float('inf'),
                                    time_step=0.1, efficiency=1):
    """create charging curve as nested list of SOC, Power[kW] and capacity in [kWh]"""
    # simple numeric creation of power over time --> to energy over time
    normalized_curve = np.array([[soc, power / capacity] for soc, power in charging_curve])
    soc = 0
    time = 0
    socs = []
    times = []
    while soc < ARGS.desired_soc_opps:
        times.append(time)
        socs.append(soc)
        power1 = min(np.interp(soc, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        soc2 = soc + time_step / 60 * power1
        power2 = min(np.interp(soc2, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        power = (power1 + power2) / 2 * efficiency
        soc += time_step / 60 * power
        time += time_step
    # Fill the soc completely in last timestep
    times.append(time)
    socs.append(ARGS.desired_soc_opps)
    return np.array((times, socs)).T


@time_it



@time_it
def get_buffer_time(trip):
    """  Return the buffer time as timedelta object"""
    return timedelta(minutes=get_buffer_time_spice_ev(trip, ARGS.default_buffer_time_opps))


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


@time_it


def get_vehicle(rot_id, this_sched):
    return this_sched.rotations[rot_id].vehicle_id


def get_vehicle_rots(rot_id, this_sched):
    return [rot for rot in this_sched.rotations if rot_id == get_vehicle(rot, this_sched)]





@time_it
def toolbox_from_pickle(sched_name, scen_name, args_name):
    """ Load the 3 files from pickle"""
    with open(args_name, "rb") as f:
        this_args = pickle.load(f)
    with open(scen_name, "rb") as f:
        scen = pickle.load(f)
    with open(sched_name, "rb") as f:
        sched = pickle.load(f)
    return sched, scen, this_args


@time_it
def toolbox_to_pickle(name, sched, scen, this_args):
    """ Dump the 3 files to pickle files"""
    args_name = "args_" + name + ".pickle"
    with open(args_name, "wb") as f:
        pickle.dump(this_args, f)
    scen_name = "scenario_" + name + ".pickle"
    with open(scen_name, "wb") as f:
        pickle.dump(scen, f)
    sched_name = "schedule_" + name + ".pickle"
    with open(sched_name, "wb") as f:
        pickle.dump(sched, f)
    return sched_name, scen_name, args_name


@time_it
def get_negative_rotations_all_electrified(this_scen, this_sched, soc_upper_thresh,
                                           soc_charge_curve_dict,
                                           not_possible_stations=None, rel_soc=False):
    """Get the rotation ids for the rotations which show negative socs, even when everything
    is electrified"""
    if not not_possible_stations:
        not_possible_stations = set()

    events = get_low_soc_events(this_scen=this_scen, this_sched=this_sched,
                                soc_upper_thresh=soc_upper_thresh,
                                not_possible_stations=not_possible_stations, rel_soc=rel_soc)

    stats = {s for e in events for s in e["stations_list"] if s not in not_possible_stations}
    electrified_station_set = set(stats)
    vehicle_socs = timeseries_calc(this_sched.rotations.values(),
                                   this_scen.vehicle_socs, this_scen, electrified_station_set,
                                   soc_charge_curve_dict=soc_charge_curve_dict,
                                   soc_upper_thresh=ARGS.desired_soc_deps)
    new_events = get_low_soc_events(this_scen, this_sched, soc_data=vehicle_socs)
    return [e["rotation"].id for e in new_events]


def create_charging_curves(this_sched):
    """Cycle through vehicles and create numerically created charging curves with energy supplied
    over time"""
    soc_charge_curve_dict = dict()
    for v_type_name in this_sched.vehicle_types:
        soc_charge_curve_dict[v_type_name] = dict()
    for name, v_type in this_sched.vehicle_types.items():
        for ch_type, data in v_type.items():
            soc_charge_curve_dict[name][ch_type] = charging_curve_to_soc_over_time(
                data["charging_curve"], data["capacity"],
                this_sched.cs_power_opps, efficiency=CONFIG.charge_eff,
                time_step=0.1)
    return soc_charge_curve_dict


@time_it
def get_must_stations(this_scen, this_sched, soc_charge_curve_dict,
                      not_possible_stations=None,
                      relative_soc=False, **kwargs):
    """Electrify everything minus 1 station. If without the stations there are below zero events
     it is a must have station"""

    soc_lower_thresh = kwargs.get("soc_lower_thresh", CONFIG.min_soc)
    events = get_low_soc_events(this_scen=this_scen, this_sched=this_sched,
                                not_possible_stations=not_possible_stations, rel_soc=relative_soc)

    if not not_possible_stations:
        not_possible_stations = set()
    stats = {s for e in events for s in e["stations_list"] if s not in not_possible_stations}
    electrified_station_set_all = set(stats)

    must_stations = set()
    for station in electrified_station_set_all:
        print(".", end="")
        electrified_station_set = electrified_station_set_all.difference([station])
        vehicle_socs = timeseries_calc(this_sched.rotations.values(),
                                       this_scen.vehicle_socs,
                                       this_scen, electrified_station_set,
                                       soc_charge_curve_dict=soc_charge_curve_dict,
                                       soc_upper_thresh=ARGS.desired_soc_deps)

        for rot in this_sched.rotations:
            soc, start, end = get_rotation_soc(rot, this_sched, this_scen, vehicle_socs)
            if np.min(soc[start:end]) < soc_lower_thresh:
                LOGGER.debug("%s , with min soc: %s", station, np.min(soc[start:end]))
                must_stations.add(station)
                break

    return must_stations


if __name__ == "__main__":
    main()
