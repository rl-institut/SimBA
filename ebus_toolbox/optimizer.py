""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
from datetime import datetime, timedelta
import json
import os
import sys
import warnings
from copy import copy, deepcopy
from pathlib import Path
from time import time
import logging
import math
import shutil
import report
# spice_ev scenario
from src import scenario
import schedule
import rotation
import pickle
import matplotlib
import matplotlib.pyplot as plt
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.trip import Trip
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox.util import uncomment_json_file, get_buffer_time as get_buffer_time_spice_ev
from optimizer_config import read_config, OptimizerConfig
import numpy as np

matplotlib.use("TkAgg")

global LOGGER
global args
config = OptimizerConfig()
# global list to be able to keep track of runtimes for optimizing components which are slow
timers = [0] * 10


def setup_logger():
    this_logger = logging.getLogger(__name__)
    this_logger.setLevel(logging.DEBUG)

    # logging to one file which keeps track of optimization over many runs
    file_handler_all_opts = logging.FileHandler('optimizer.log')
    file_handler_all_opts.setLevel(config.debug_level)

    # and logging to a file which is put in the folder with the other optimizer results
    file_handler_this_opt = logging.FileHandler(Path(args.output_directory) / Path('optimizer.log'))
    file_handler_this_opt.setLevel(config.debug_level)

    formatter = logging.Formatter('%(asctime)s:%(message)s',
                                  "%m%d %H%M%S")

    file_handler_all_opts.setFormatter(formatter)
    file_handler_this_opt.setFormatter(formatter)

    formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(config.debug_level)
    this_logger.addHandler(file_handler_this_opt)
    this_logger.addHandler(file_handler_all_opts)
    this_logger.addHandler(stream_handler)
    return this_logger


def main():
    config_path = "./data/examples/optimizer.cfg"
    run_optimization(config_path)
    pass


def run_optimization(config_path, sched=None, scen=None, this_args=None):
    """    Optimizes scenario by adding electrified stations sparingly
    until scenario has no below 0 soc events.

    :param config_path: path to optimizer.cfg file.
    :type config_path: str

    :param sched: Simulation schedule containing buses, rotations etc.
    :type sched: ebus_toolbox.Schedule

    :param scen: Simulation scenario containing simulation results
                     including the SoC of all vehicles over time
    :type scen: spice_ev.Scenario
    :param this_args: Simulation arguments for manipulation or generated outputs
    :type this_args: object

    :return: (Schedule,Scenario) optimizied schedule and Scenario
    :rtype: tuple(ebus_toolbox.Schedule, spice_ev.Scenario)
    """
    global config
    global LOGGER
    config = read_config(config_path)

    # load pickle files
    if sched is None or scen is None or this_args is None:
        # if no schedule was given as argument, make sure no scenario
        # and args input was given as well.
        assert sched == scen == this_args is None
        sched, scen, args = toolbox_from_pickle(config.schedule, config.scenario, config.args)

    global args
    args = this_args

    # prepare Filesystem with folders and paths
    now = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    args.output_directory = Path(config.output_path) / str(now + "_optimizer")
    # create subfolder for specific sim results with timestamp.
    # if folder doesnt exists, create folder.
    # needs to happen after set_options_from_config since
    # args.output_directory can be overwritten by config
    args.output_directory.mkdir(parents=True, exist_ok=True)

    # copy paste the config file
    c = Path(config_path)
    o = args.output_directory / Path("optimizer_config.cfg")
    shutil.copy(c, o)

    if args.save_soc:
        args.save_soc = args.output_directory / "simulation_soc_spiceEV.csv"
    new_ele_stations_path = args.output_directory / Path("optimized_stations" + ".json")
    LOGGER = setup_logger()

    # remove those args, since they lead to file creation, which is not
    # needed.
    del args.save_timeseries
    del args.save_results
    if args.desired_soc_deps != 1:
        LOGGER.error("Fast calc is not yet optimized for desired socs unequal to 1")

    # which rotations should be excluded?
    exclusion_rots = config.exclusion_rots

    # which stations has to be electrified
    inclusion_stations = config.inclusion_stations

    # which station can not be electrified
    exclusion_stations = config.exclusion_stations

    if config.reduce_rots:
        sched.rotations = {rot: sched.rotations[rot] for rot in config.rots}

    opt_type = config.opt_type
    solver = config.solver
    remove_impossible_rots = config.remove_impossible_rots
    rebase_scenario = config.rebase_scenario
    node_choice = config.node_choice

    t = time()
    # set battery and charging curves through config file if wished for
    for name, v_type in sched.vehicle_types.items():
        for charge_type, vehicle in v_type.items():
            if config.battery_capacity is not None:
                vehicle["capacity"] = config.battery_capacity
            if config.charging_curve is not None:
                vehicle["charging_curve"] = config.charging_curve

    # create a decision tree or load one from a previous run
    decision_tree_path = config.decision_tree_path
    if decision_tree_path is not None:
        with open(decision_tree_path, "rb") as file:
            decision_tree = pickle.load(file)
    else:
        decision_tree = dict()

    # stations which are included can not be included again. Therefore they get into
    # the set of not possible stations
    not_possible_stations = inclusion_stations.union(exclusion_stations)
    s = time()

    # if stations have to be included they are stored in this set
    must_include_set = set()

    # rebasing the scenario meaning simulating it again with the given conditions of
    # included and excluded stations and rotations
    if rebase_scenario:
        LOGGER.debug("Spice EV Rebasing Scenario")
        sched, scen, must_include_set, ele_stations = preprocessing_scenario(
            sched, scen, args,
            inclusion_stations,
            exclusion_rots=exclusion_rots, run_only_neg=config.run_only_neg)
        LOGGER.debug(f"Rebasing took {time() - s} sec")
        if config.pickle_rebased:
            toolbox_to_pickle(config.pickle_rebased_name, sched, scen, args)
            LOGGER.debug(f"Rebased scenario pickled  as {config.pickle_rebased_name}")
    else:
        with open(args.electrified_stations, "r", encoding="utf-8", ) as file:
            ele_stations = json.load(file)
        # Electrify inclusion stations
        for stat in inclusion_stations:
            electrify_station(stat, ele_stations, must_include_set)
        rots = {r: sched.rotations[r] for r in sched.rotations if r not in exclusion_rots}
        sched.rotations = rots

    s = time()

    global timer_for_calc
    # create charging dicts which contain energy over time, which is numerically created
    soc_charge_curve_dict = dict()
    for v_type_name in sched.vehicle_types:
        soc_charge_curve_dict[v_type_name] = dict()
    for name, v_type in sched.vehicle_types.items():
        for ch_type, data in v_type.items():
            soc_charge_curve_dict[name][ch_type] = charging_curve_to_soc_over_time(
                data["charging_curve"], data["capacity"],
                sched.cs_power_opps, efficiency=config.charge_eff,
                timestep=0.1)

    ###########
    # Remove none Values from socs in the vehicle_socs an
    remove_none_socs(scen)

    if config.check_for_must_stations:
        must_stations = get_must_stations(scen, sched, args.desired_soc_opps,
                                          soc_charge_curve_dict,
                                          not_possible_stations=not_possible_stations,
                                          soc_lower_thresh=config.min_soc, relative_soc=False)
        LOGGER.warning("%s must stations %s", len(must_stations), must_stations)
        not_possible_stations = not_possible_stations.union(must_stations)
        for stat in must_stations:
            # dont put must stations in electrified set --> therefore create trash set
            electrify_station(stat, ele_stations, must_include_set)
        scen.vehicle_socs = timeseries_calc('best_station_ids[0]', sched.rotations.values(),
                                            scen.vehicle_socs,
                                            scen, must_stations,
                                            soc_charge_curve_dict=soc_charge_curve_dict,
                                            soc_upper_thresh=args.desired_soc_deps)

    LOGGER.debug("Starting greedy optimization")
    ele_station_set = set()
    ele_stations, ele_station_set, could_not_be_electrified, list_greedy_sets = \
        optimization_loop(ele_stations, ele_station_set, scen, sched,
                          not_possible_stations, soc_upper_thresh=args.desired_soc_opps,
                          soc_lower_thresh=config.min_soc,
                          solver=solver, opt_type=opt_type,
                          node_choice=node_choice,
                          soc_charge_curve_dict=soc_charge_curve_dict,
                          decision_tree=decision_tree)

    ele_station_set = ele_station_set.union(must_include_set)
    LOGGER.debug(f"Non Spice Ev methods took {time() - s - timer_for_calc} sec")
    LOGGER.debug(f"Solver: {solver} took {timer_for_calc} s")
    LOGGER.debug(f"Electrified Stations: {len(ele_station_set)}")
    LOGGER.debug(ele_station_set)
    LOGGER.debug(could_not_be_electrified)
    scen.vehicle_socs = timeseries_calc('best_station_ids[0]', sched.rotations.values(),
                                        scen.vehicle_socs,
                                        scen, ele_station_set,
                                        soc_charge_curve_dict=soc_charge_curve_dict,
                                        soc_upper_thresh=args.desired_soc_deps)
    new_events = get_below_zero_soc_events(scen, sched.rotations, sched,
                                           soc_lower_thresh=config.min_soc)

    LOGGER.debug("Still not electrified with fast calc")
    for event in new_events:
        LOGGER.debug(event["rotation"].id)
    LOGGER.debug("###")

    global timers
    LOGGER.debug(timers)
    with open(new_ele_stations_path, "w", encoding="utf-8", ) as file:
        json.dump(ele_stations, file, indent=2)

    LOGGER.debug("Spice EV is calculating optimized case as a complete scenario")
    final_sched, final_scen, ele_station_set, ele_stations = preprocessing_scenario(
        sched, scen, args, electrified_stations=ele_stations, run_only_neg=False,
        electrified_station_set=ele_station_set, cost_calc=True)

    LOGGER.warning(f"Still negative rotations:{final_sched.get_negative_rotations(final_scen)}")
    print("Finished")
    print(f"Opt took {time() - t}")

    return final_sched, final_scen


timer_for_calc = 0


def optimization_loop(electrified_stations, electrified_station_set, new_scen, new_sched,
                      not_possible_stations, soc_charge_curve_dict, soc_upper_thresh=1,
                      soc_lower_thresh=0,
                      decision_tree=None, pre_optimized_set=None, opt_type="greedy",
                      node_choice="step-by-step",
                      **kwargs):
    # Base stations for optimization, so inclusion of stations can be skipped
    base_stations = electrified_stations.copy()
    base_electrified_station_set = electrified_station_set.copy()

    # Connect events with shared stations
    ########
    could_not_be_electrified = set()
    base_scen = copy(new_scen)
    base_sched = copy(new_sched)

    # Get events where soc fell below 0. The events contain info about the problematic
    # timespan, which includes stations which could provide a soc lift
    base_events = get_below_zero_soc_events(this_scen=base_scen,
                                            rotations=list(new_sched.rotations.keys()),
                                            this_sched=base_sched,
                                            soc_upper_thresh=soc_upper_thresh,
                                            filter_standing_time=True,
                                            not_possible_stations=not_possible_stations,
                                            soc_lower_thresh=soc_lower_thresh, relative_soc=False)

    # Check if the events can be divided into subgroups which are independent
    # this makes optimization in smaller groups possible
    groups = get_groups_from_events(base_events, not_possible_stations)

    # Baseline greedy Optimization
    list_greedy_sets = [set()] * len(groups)

    # get_negative_rotations_all_electrified(base_scen, base_sched, soc_upper_thresh, soc_charge_curve_dict,
    #                                        not_possible_stations=not_possible_stations,
    #                                        soc_lower_thresh=soc_lower_thresh, relative_soc=False)

    new_scen.vehicle_socs = timeseries_calc('best_station_ids[0]', new_sched.rotations.values(),
                                            new_scen.vehicle_socs,
                                            new_scen, electrified_station_set,
                                            soc_charge_curve_dict=soc_charge_curve_dict,
                                            soc_upper_thresh=args.desired_soc_deps)

    # Base line is created simply by not having a decision tree and not a pre optimized_set yet
    for group_nr, group in enumerate(groups[:]):
        events, stations = group
        linien = {lne for e in events for lne in e["rotation"].lines}
        electrified_stations = base_stations.copy()
        electrified_station_set = base_electrified_station_set.copy()
        LOGGER.warning(f"Optimizing {group_nr + 1} out of {len(groups)}. This includes these Lines")
        LOGGER.warning(linien)
        LOGGER.warning("%s events", (len(events)))
        solver = kwargs.get("solver", "spiceev")
        if node_choice == "brute":
            choice_func = choose_station_brute
        else:
            choice_func = choose_station_step_by_step

        # first run is always step by step
        group_optimization(group, base_scen, base_sched,
                           electrified_stations, electrified_station_set,
                           could_not_be_electrified, not_possible_stations,
                           choose_station_step_by_step, soc_charge_curve_dict,
                           pre_optimized_set=None,
                           decision_tree=decision_tree, soc_lower_thresh=soc_lower_thresh,
                           soc_upper_thresh=soc_upper_thresh,
                           events_remaining=[len(events)], type=solver,
                           **kwargs)

        LOGGER.warning("Greedy Result ++++++++ %s stations out of %s", len(electrified_station_set),
                       len(stations))
        LOGGER.warning(electrified_station_set)

        if opt_type == "deep":
            sols = []
            i = 0
            cont_loop = True
            t = time()
            combinations = combs_unordered_no_putting_back(len(stations),
                                                           len(electrified_station_set))

            LOGGER.debug(f"There are {combinations} combinations")
            while i < config.max_brute_loop and cont_loop:
                i += 1
                if i % 10 == 0:
                    print(time() - t)
                    t = time()
                    print(len(decision_tree), " Knotenpunkte durchwandert")
                    print(f"Optimal solution has length {len(electrified_station_set)}")
                not_possible_stations = copy(not_possible_stations)
                pre_optimized_set = copy(electrified_station_set)
                could_not_be_electrified_copy = could_not_be_electrified.copy()
                electrified_stations = base_stations.copy()
                # electrified_station_set = base_electrified_station_set.copy()
                new_electrified_set = set()
                new_stations, cont_loop = \
                    group_optimization(group, base_scen, base_sched,
                                       electrified_stations, new_electrified_set,
                                       could_not_be_electrified_copy,
                                       not_possible_stations,
                                       choice_func, soc_charge_curve_dict,
                                       pre_optimized_set=pre_optimized_set,
                                       decision_tree=decision_tree,
                                       events_remaining=[len(events)],
                                       soc_upper_thresh=soc_upper_thresh,
                                       soc_lower_thresh=soc_lower_thresh, type=solver)
                # if a new set was found, print it and save it in sols

                if new_electrified_set != pre_optimized_set and new_stations is not None:
                    LOGGER.warning(
                        f"Optimized with {len(new_electrified_set)}  stations {str('#' * 20)} \
                    {stations_hash(new_electrified_set)}")
                    sols.append(new_electrified_set)
                    if len(new_electrified_set) < len(pre_optimized_set):
                        electrified_station_set = new_electrified_set

            LOGGER.debug(sols)
        list_greedy_sets[group_nr] = electrified_station_set.copy()
        LOGGER.debug("Optimized with {} stations out of {}".format(len(electrified_station_set),
                                                                   len(stations)))
    if config.save_decision_tree:
        with open(args.output_directory / Path("decision_tree.pickle"), "wb") as file:
            pickle.dump(decision_tree, file)
    for single_set in list_greedy_sets:
        for stat in single_set:
            electrify_station(stat, electrified_stations, electrified_station_set)

    return electrified_stations, electrified_station_set, could_not_be_electrified, list_greedy_sets


def get_groups_from_events(events, not_possible_stations=None, could_not_be_electrified=None):
    # First create simple list of station sets for single events.
    # Electrified and other not possible to electrify stations should not connect groups

    # Making sure default arguments are none and not mutable
    if not not_possible_stations:
        not_possible_stations = set()

    if not could_not_be_electrified:
        could_not_be_electrified = set()

    possible_stations = [
        {station for station in e["stations"] if station not in not_possible_stations}
        for e
        in events]
    # If stations overlap join them
    station_subsets = join_all_subsets(possible_stations)
    event_groups = [[] for __ in range(len(station_subsets))]

    # Group the events in the same manner as stations, so that the same amount of event groups are
    # created as station subsets
    for e in events:
        for i, subset in enumerate(station_subsets):
            # Every station from an event must share the same station sub_set. Therefore its enough
            # to check only the first element
            if len(e["stations"]) > 0:
                if next(iter(e["stations"])) in subset:
                    event_groups[i].append(e)
                    break
        else:
            LOGGER.warning(f'Didnt find rotation {e["rotation"].id} in any subset'
                           f'of possible electrifiable stations')
            # this event will no show up in an event_group.
            # therefore it needs to be put into this set
            could_not_be_electrified.update([e["rotation"].id])

    groups = list(zip(event_groups, station_subsets))
    return sorted(groups, key=lambda x: len(x[1]))


def group_optimization(group, base_scen, base_sched,
                       electrified_stations, electrified_station_set,
                       could_not_be_electrified,
                       not_possible_stations, choose_station_function, soc_curve_dict,
                       pre_optimized_set=None, decision_tree=None,
                       soc_lower_thresh=0, soc_upper_thresh=1,
                       tree_position=None, type="spiceev", **kwargs):
    if not tree_position:
        tree_position = []

    # Base socs are the socs without electrification, they get passed through the stacks
    # so the quick calculation can mutate them since quick calculation cant take place iteratly
    if kwargs.get("base_socs"):
        base_socs = kwargs.get("base_socs")
    else:
        base_socs = {id: soc for id, soc in base_scen.vehicle_socs.items()}

    # give tree position
    LOGGER.debug("%s with length of %s", tree_position, len(tree_position))

    # Unpack events and possible stations
    event_group, possible_stations = group

    # Get the events which are still negative
    events_remaining = kwargs.get("events_remaining", [99999])

    # Copy Scen and sched and deepcopy vehicle socs so upper levels are not changed by the following
    # optimization
    new_scen = copy(base_scen)
    new_sched = copy(base_sched)

    # Get rotations from event dict and calculate the missing energy
    rotation_dict = {e["rotation"].id: e["rotation"] for e in event_group}
    missing_energy = get_missing_energy(event_group)

    if missing_energy >= 0:
        LOGGER.debug("Already electrified: Returning set")
        return electrified_stations, True

    station_eval = evaluate(event_group, new_scen, soc_curve_dict,
                            soc_upper_thresh=soc_upper_thresh, soc_lower_thresh=soc_lower_thresh)

    # Debug.Log the missing energy and station evaluation
    station_eval_to_logger(missing_energy, station_eval)

    # Get the best stations, brute or step_by_step
    best_station_ids, recursive = choose_station_function(station_eval, electrified_station_set,
                                                          pre_optimized_set, decision_tree,
                                                          missing_energy=missing_energy)
    stat_eval_dict = {stat_id[0]: stat_id[1]["pot_sum"] for stat_id in station_eval}
    if best_station_ids is None:
        LOGGER.warning(
            f"No useful station found with "
            f"{events_remaining} rotations not electrified yet. "
            f" Stopped after electrifying {len(electrified_station_set)}")

        if pre_optimized_set is not None:
            # Remove electrified stations in this run
            c = electrified_station_set.copy()
            for stat in c:
                electrified_stations.pop(stat)
                electrified_station_set.remove(stat)
            # Overwrite with preoptimized set
            for stat in pre_optimized_set:
                electrify_station(stat, electrified_stations, electrified_station_set)
        else:
            could_not_be_electrified.update(list(rotation_dict.keys()))
        return None, False
    LOGGER.debug("%s, with first pot of %s", best_station_ids, stat_eval_dict[best_station_ids[0]])

    # Electrify station(s)
    for stat_id in best_station_ids:
        electrify_station(stat_id, electrified_stations, electrified_station_set)

    s = time()

    event_rotations = {x["rotation"] for x in event_group}
    if type == "quick":
        # Quick calculation has to electrify everything in one step,
        # or the lifting of socs is not correct
        new_scen.vehicle_socs = deepcopy(base_socs)
        new_scen.vehicle_socs = timeseries_calc(best_station_ids[0], event_rotations,
                                                new_scen.vehicle_socs,
                                                new_scen, electrified_station_set, soc_curve_dict,
                                                soc_upper_thresh=args.desired_soc_deps)
        lifted_socs = deepcopy(new_scen.vehicle_socs)
    else:
        new_sched.rotations = rotation_dict
        new_sched, new_scen = run_schedule(new_sched, args,
                                           electrified_stations=electrified_stations)
        lifted_socs = None

    global timer_for_calc
    global timers
    timer_for_calc += time() - s
    timers[2] += time() - s
    not_possible_stations = set(electrified_stations.keys()).union(not_possible_stations)
    event_rotations_id = {event["rotation"].id for event in event_group}
    new_events = get_below_zero_soc_events(new_scen, event_rotations_id,
                                           new_sched,
                                           soc_upper_thresh=soc_upper_thresh,
                                           filter_standing_time=True,
                                           not_possible_stations=not_possible_stations,
                                           soc_lower_thresh=soc_lower_thresh, relative_soc=True)

    delta_energy = get_missing_energy(new_events)

    events_remaining[0] -= len(event_group) - len(new_events)
    LOGGER.debug("Last electrification electrified %s/%s."
                 " %s remaining events in the base group.",
                 len(event_group) - len(new_events), len(event_group), events_remaining[0])
    delta_base_energy = delta_energy  # get_missing_energy(base_events)

    # Put this node intro the decision tree including the missing energy
    if decision_tree is not None:
        decision_tree = node_to_tree(decision_tree, electrified_station_set, delta_base_energy)

    # Everything electrified
    if delta_energy >= 0:
        return electrified_stations, True
    # Some choice functions might not need a recursive call, they return here. recursive is set
    # by the choose_station_function
    elif not recursive:
        return None, True

    # Check if the events can be divided into subgroups which are independent
    groups = get_groups_from_events(new_events, not_possible_stations, could_not_be_electrified)

    for k, group in enumerate(groups):
        this_tree = tree_position.copy()
        this_tree.append(k)
        new_stations, _ = group_optimization(group, new_scen, new_sched,
                                             electrified_stations,
                                             electrified_station_set,
                                             could_not_be_electrified,
                                             not_possible_stations, choose_station_function,
                                             soc_curve_dict,
                                             pre_optimized_set, decision_tree,
                                             lifted_socs=lifted_socs,
                                             events_remaining=events_remaining,
                                             soc_lower_thresh=soc_lower_thresh,
                                             soc_upper_thresh=soc_upper_thresh,
                                             tree_position=this_tree, type=type,
                                             base_socs=base_socs
                                             )

        if new_stations is not None:
            electrified_stations.update(new_stations)
        else:
            return None, True

        # If there is no pre optimized set function can return
        if pre_optimized_set is None:
            return electrified_stations, True

        # If there is a preoptimized set check if the preoptimized set leads to pruning
        # i.e. if if following  branch makes sense
        if len(pre_optimized_set) - len(electrified_station_set) < 10:
            new_scen.vehicle_socs = timeseries_calc(best_station_ids[0], event_rotations,
                                                    new_scen.vehicle_socs,
                                                    new_scen, electrified_station_set,
                                                    soc_curve_dict,
                                                    soc_upper_thresh=args.desired_soc_deps)

            prune_events = get_below_zero_soc_events(new_scen, event_rotations_id,
                                                     new_sched,
                                                     soc_upper_thresh=soc_upper_thresh,
                                                     filter_standing_time=True,
                                                     not_possible_stations=not_possible_stations,
                                                     soc_lower_thresh=soc_lower_thresh,
                                                     relative_soc=True)

            station_eval = evaluate(prune_events, new_scen, soc_curve_dict,
                                    soc_upper_thresh=soc_upper_thresh,
                                    soc_lower_thresh=soc_lower_thresh)
            prune_missing_energy = get_missing_energy(prune_events)
            if not is_branch_promising(station_eval, electrified_station_set,
                                       pre_optimized_set, prune_missing_energy):
                print("Branch pruned early")
                is_branch_promising(station_eval, electrified_station_set,
                                    pre_optimized_set, missing_energy)
                return None, True


def get_missing_energy(events):
    missing_energy = 0
    for e in events:
        missing_energy += e["min_soc"] * e["capacity"]
    return missing_energy


def node_to_tree(decision_tree, electrified_station_set, delta_base_energy):
    node_name = stations_hash(electrified_station_set)
    try:
        decision_tree[node_name]["missing_energy"] = delta_base_energy
        decision_tree[node_name]["visit_counter"] += 1
        # todo add is_viable
        LOGGER.debug("already visited")
    except KeyError:
        decision_tree[node_name] = dict()
        decision_tree[node_name]["missing_energy"] = delta_base_energy
        decision_tree[node_name]["visit_counter"] = 1

    return decision_tree


def preprocessing_scenario(this_sched, this_scen, this_args,
                           inclusion_stations=None,
                           electrified_stations=None,
                           electrified_station_set=None,
                           exclusion_rots=None,
                           run_only_neg=False,
                           cost_calc=False):
    if not inclusion_stations:
        inclusion_stations = set()

    if not electrified_station_set:
        electrified_station_set = set()

    if not exclusion_rots:
        exclusion_rots = set()

    if electrified_stations is None:
        with open(this_args.electrified_stations, "r", encoding="utf-8", ) as f:
            electrified_stations = json.load(f)

    # Electrify inclusion stations
    for stat in inclusion_stations:
        electrify_station(stat, electrified_stations, electrified_station_set)

    # Calc new but only prev. negative rotations
    if run_only_neg:
        rots = this_sched.get_negative_rotations(this_scen)
        rots = {r: this_sched.rotations[r] for r in rots if r not in exclusion_rots}
        this_sched.rotations = rots
    else:
        rots = {r: this_sched.rotations[r] for r in this_sched.rotations if r not in exclusion_rots}
        this_sched.rotations = rots

    new_sched, new_scen = run_schedule(this_sched, this_args, electrified_stations,
                                       cost_calc=cost_calc)
    report.generate(new_sched, new_scen, args)

    return new_sched, new_scen, electrified_station_set, electrified_stations


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


def plot_rot(id, this_sched, this_scen, ax=None, rot_only=True):
    """ Simple plot of data without having to create subplots"""
    a, b, c = get_rotation_soc(id, this_sched, this_scen)
    if not rot_only:
        b = 0
        c = -1
    if ax is None:
        fig, ax = plt.subplots()
        ax.plot(a[b:c], linewidth=2.0)
        return ax
    else:
        ax.plot(a[b:c], linewidth=2.0)
        return ax


def join_all_subsets(subsets):
    joined_subset = True
    while joined_subset:
        joined_subset, subsets = join_subsets(subsets)
    return subsets


def choose_station_brute(station_eval, electrified_station_set,
                         pre_optimized_set=None, decision_tree=None, missing_energy=0):
    station_ids = [x[0] for x in station_eval]
    a = combination_generator(station_ids, len(pre_optimized_set) - 1)
    station_eval_dict = {stat[0]: stat[1] for stat in station_eval}
    for comb in a:
        node_name = stations_hash(comb)
        if node_name not in decision_tree:
            # Only check the brute force station if they have the remote chance of fullfilling
            # the missing energy
            # Potential>missing energy * 80%
            potential = sum([station_eval_dict[stat]["pot_sum"] for stat in comb])
            if potential > -missing_energy * config.estimation_threshold:
                return comb, False
            else:
                LOGGER.debug("skipped %s since potential is too low %s %%", comb,
                             round(potential / -missing_energy * 100, 0))
    else:
        LOGGER.debug("calculated all viable possibilities")
        return None, False


def is_branch_promising(station_eval, electrified_station_set,
                        pre_optimized_set, missing_energy):
    delta = len(pre_optimized_set) - len(electrified_station_set)
    pot = 0
    for i in range(0, min(delta, len(station_eval))):
        pot += station_eval[i][1]["pot_sum"]
    if pot < -missing_energy * config.estimation_threshold:
        print(f"Not enough potential {round(pot, 0)} / {round(-missing_energy, 0)} after ",
              stations_hash(electrified_station_set))
        return False
    return True


def choose_station_step_by_step(station_eval, electrified_station_set,
                                pre_optimized_set=None, decision_tree=None, missing_energy=0):
    # Filter functions to stop simulating cases which have no hope of being optimal.
    # If in optimization mode, optimization can break if station amount is superseded
    # This filter is done better by the next
    # if pre_optimized_set is not None:
    #     if len(electrified_station_set)>len(pre_optimized_set):
    #         return pre_optimized_set
    # Potentials have to be at least as promising as the pre-optimized case
    if pre_optimized_set is not None:
        if not is_branch_promising(station_eval, electrified_station_set,
                                   pre_optimized_set, missing_energy):
            # Best station id is none and dont go deeper in recursion
            return None, False

    min_count_visited = float('inf')
    for station in station_eval:
        # Create a station combination from already electrified stations and possible new station
        check_stations = electrified_station_set.union([station[0]])

        if decision_tree is not None:
            node_name = stations_hash(check_stations)
            if node_name in decision_tree.keys():
                min_count_visited = min(min_count_visited,
                                        decision_tree[node_name]["visit_counter"])
                # If already checked skip to next one
                continue
        if station[0] not in electrified_station_set:
            best_station_id = station[0]
            return [best_station_id], True
    else:
        # Every possible node from here was evaluated already
        # to do what now?
        # Simply visit the least visited node
        for station in station_eval:
            # Create a station combination from already electrified stations
            # and possible new station
            check_stations = electrified_station_set.union([station[0]])
            if decision_tree[stations_hash(check_stations)]["visit_counter"] == min_count_visited:
                best_station_id = station[0]
                return [best_station_id], True
    return None, True


def stations_hash(stations_set):
    return str(sorted(list(stations_set)))


def timeseries_calc(station, rotations, soc_dict, eval_scen, ele_station_set,
                    soc_charge_curve_dict, soc_upper_thresh=1):
    global timers
    ele_stations = {*ele_station_set, station}
    s2 = time()
    soc_dict = copy(soc_dict)

    for rot in rotations:
        ch_type = (rot.vehicle_id.find("oppb") > 0) * "oppb" + (
                rot.vehicle_id.find("depb") > 0) * "depb"
        v_type = rot.vehicle_id.split("_" + ch_type)[0]
        soc_over_time_curve = soc_charge_curve_dict[v_type][ch_type]
        soc = soc_dict[rot.vehicle_id]
        for i, trip in enumerate(rot.trips):
            s = time()
            if trip.arrival_name not in ele_stations:
                continue
            idx = get_index_by_time(trip.arrival_time, eval_scen)
            try:
                standing_time_min = get_charging_time(trip, rot.trips[i + 1], args)
            except IndexError:
                standing_time_min = 0

            d_soc = get_delta_soc(soc_over_time_curve, soc[idx], standing_time_min)
            buffer_idx = int((get_buffer_time(trip, args)) / timedelta(minutes=1))
            delta_idx = int(standing_time_min) + 1
            old_soc = soc[idx + buffer_idx:idx + buffer_idx + delta_idx].copy()
            soc[idx + buffer_idx:] += d_soc
            soc[idx + buffer_idx:idx + buffer_idx + delta_idx] = old_soc
            soc[idx + buffer_idx:idx + buffer_idx + delta_idx] += np.linspace(0, d_soc, delta_idx)

            soc_pre = soc[:idx]
            soc = soc[idx:]
            soc_max = np.max(soc)
            timers[0] += time() - s
            s = time()
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
            timers[1] += time() - s
        soc_dict[rot.vehicle_id] = soc
    timers[3] += (time() - s2) * 100
    return soc_dict


def evaluate(events, eval_scen, soc_curve_dict, soc_upper_thresh=1, soc_lower_thresh=0,
             electrified_station_set=None, decision_tree=None):
    # Analyse stations for "helpful" energy supply. Energy supply is helpful the minimal soc of an
    # event is raised up to a minimal soc (probably zero). The supplied energy is approximated by
    # loading power, standing time at a station, soc at station, minimal soc of the event
    if not electrified_station_set:
        electrified_station_set = set()

    station_eval = dict()
    for e in events:
        soc_over_time = soc_curve_dict[e["v_type"]][e["ch_type"]]
        for i, trip in enumerate(e["trip"]):
            # Station is only evaluated if station name is part of event stations
            # Only these stations showed potential in electrification, e.g enough standing time
            if trip.arrival_name not in e["stations"]:
                continue
            idx = get_index_by_time(trip.arrival_time, eval_scen)
            soc = eval_scen.vehicle_socs[e["vehicle_id"]][idx]

            max_soc = soc_upper_thresh
            min_soc = soc_lower_thresh
            # Potential is the minimal amount of
            delta_soc_pot = min(max_soc - soc,
                                min_soc - e["min_soc"],
                                soc - e["min_soc"],
                                max_soc - min_soc)

            try:
                standing_time_min = get_charging_time(trip, e["trip"][i + 1], args)
            except IndexError:
                standing_time_min = 0

            # energy_charging_potential = standing_time_min *60 * ch_power
            e_charging_pot = get_delta_soc(soc_over_time, soc, standing_time_min) * e["capacity"]

            # Potential is at max the minimum between the useful delta soc * capacity or the
            # energy provided by charging for the full standing time
            delta_E_pot = min(delta_soc_pot * e["capacity"], e_charging_pot)
            d = dict(E_pot=delta_E_pot,
                     standing_time=timedelta(minutes=standing_time_min))
            try:
                station_eval[trip.arrival_name]["pot_list"].append(d)
            except KeyError:
                station_eval[trip.arrival_name] = dict(pot_list=[], pot_sum=0)
                station_eval[trip.arrival_name]["pot_list"].append(d)

    # time_list = []
    for station_name, stat_dict in station_eval.items():
        if decision_tree is not None:
            check_stations = electrified_station_set.union(station_name)
            if check_stations in decision_tree.keys():
                stat_dict["pot_sum"] = decision_tree[str(check_stations)]["missing_energy"] - \
                                       decision_tree[str(electrified_station_set)]["missing_energy"]
                # decision_tree[str(electrified_station_set)]["children"].append(check_stations)
                continue

        sum = 0
        standing_time = timedelta(minutes=0)
        for pot in stat_dict["pot_list"]:
            sum += pot["E_pot"]
            standing_time += pot["standing_time"]
        stat_dict["pot_sum"] = sum
        # time_list.append(standing_time / timedelta(minutes=1))
    # Sort by pot_sum
    station_eval = list(dict(sorted(station_eval.items(), key=lambda x: x[1]["pot_sum"])).items())
    station_eval.reverse()
    return station_eval


def get_charging_time(trip1, trip2, args):
    standing_time_min = (trip2.departure_time - trip1.arrival_time) \
                        / timedelta(minutes=1)
    buffer_time = (get_buffer_time(trip1, args) / timedelta(minutes=1))
    standing_time_min -= buffer_time

    if args.min_charging_time > standing_time_min:
        return 0
    return max(0, standing_time_min)


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


def get_index_by_time(search_time, this_scen: scenario.Scenario):
    start_time = this_scen.start_time
    delta_time = timedelta(minutes=60 / this_scen.stepsPerHour)
    idx = (search_time - start_time) // delta_time
    return idx


def get_time_by_index(idx, this_scen: scenario.Scenario):
    start_time = this_scen.start_time
    delta_time = timedelta(minutes=60 / this_scen.stepsPerHour)
    searched_time = start_time + delta_time * idx
    return searched_time


def get_trips(rot: rotation.Rotation, start_idx: int, end_idx: int, scen: scenario.Scenario):
    # return trips in a rotation from a start to an end index, if the arrival time is in between
    # the start and end idx
    start_time_event = get_time_by_index(start_idx, scen)
    end_time_event = get_time_by_index(end_idx, scen)

    trips = []
    for i, trip in enumerate(rot.trips):
        if end_time_event > trip.arrival_time > start_time_event:
            # standing_time = 0
            # try:
            #     standing_time = trip.arrival_time - trips[i + 1].departure_time
            # except IndexError:
            #     standing_time = 0
            trips.append(trip)

    return trips


def get_rotation_soc(rot_id, this_sched, this_scen, soc_data: dict = None):
    rot = this_sched.rotations[rot_id]
    rot_start_idx = get_index_by_time(rot.departure_time, this_scen)
    rot_end_idx = get_index_by_time(rot.arrival_time, this_scen)
    if soc_data:
        return soc_data[rot.vehicle_id], rot_start_idx, rot_end_idx
    return this_scen.vehicle_socs[rot.vehicle_id], rot_start_idx, rot_end_idx


def get_below_zero_soc_events(this_scen: scenario.Scenario, rotations,
                              this_sched: schedule.Schedule,
                              soc_upper_thresh=0.9, filter_standing_time=True,
                              not_possible_stations=None, soc_lower_thresh=0, relative_soc=False,
                              soc_data=None):
    if not not_possible_stations:
        not_possible_stations = set()
    # Create list of events which describe trips which end in a soc below zero
    # The event is bound by the lowest soc and an upper soc threshhold which is naturally 1
    # Properties before and after these points have no effect on the event itself, similar to
    # an event horizon
    events = []
    SOC_UPPER_THRESH = soc_upper_thresh
    count_electrified_rot = 0

    for rot_id in rotations:
        rot = this_sched.rotations[rot_id]
        soc, rot_start_idx, rot_end_idx = get_rotation_soc(rot_id, this_sched, this_scen, soc_data)
        soc = [s if s is not None else 999 for s in soc]
        idx = range(0, len(soc))

        comb = list(zip(soc, idx))[rot_start_idx:rot_end_idx]
        min_soc, min_idx = min(comb, key=lambda x: x[0])
        reduced_list = comb.copy()
        soc_lower_thresh_cur = soc_lower_thresh
        # if rotation gets a start soc below 1 this should change below 0 soc events since fixing
        # the rotation before would lead to fixing this rotation

        # if using relative SOC, SOC lookup has to be adjusted
        if relative_soc:
            start_soc = comb[0][0]
            soc_lower_thresh_cur = min(start_soc, soc_upper_thresh) - (
                    soc_upper_thresh - soc_lower_thresh)
            SOC_UPPER_THRESH = soc_lower_thresh_cur + soc_upper_thresh
        if min_soc >= soc_lower_thresh_cur:
            count_electrified_rot += 1
        while min_soc < soc_lower_thresh_cur:
            i = min_idx
            idx = [x[1] for x in reduced_list]
            while soc[i] < SOC_UPPER_THRESH:
                if i == rot_start_idx:
                    break
                i -= 1
            start_comb = idx.index(i)
            start = i
            i = min_idx
            while soc[i] < SOC_UPPER_THRESH:
                if i >= rot_end_idx - 1:
                    break
                i += 1
            end_comb = idx.index(i)
            trips = get_trips(rot=rot, start_idx=start, end_idx=min_idx, scen=this_scen)
            possible_stations = set()
            possible_stations_list = []
            if not filter_standing_time:
                possible_stations = {t.arrival_name for t in trips}
                possible_stations_list = [t.arrival_name for t in trips]
            else:
                for ii, trip in enumerate(trips):
                    try:
                        try:
                            standing_time_min = get_charging_time(trip, trips[ii + 1], args)
                        except IndexError:
                            standing_time_min = 0
                        if standing_time_min > 0:
                            possible_stations.add(trip.arrival_name)
                            possible_stations_list.append(trip.arrival_name)
                    except IndexError:
                        pass

            possible_stations = possible_stations.difference(not_possible_stations)
            cht = rot.vehicle_id.find("depb")
            ch_type = (cht > 0) * "depb" + (cht <= 0) * "oppb"
            type = rot.vehicle_id.split("_" + ch_type)[0]
            event = dict(start_idx=start, end_idx=min_idx,
                         min_soc=min_soc, stations=possible_stations,
                         vehicle_id=rot.vehicle_id, trip=trips,
                         rotation=rot, stations_list=possible_stations_list,
                         capacity=this_sched.vehicle_types[type][ch_type]['capacity'],
                         v_type=type, ch_type=ch_type)

            events.append(event)
            copy_list = reduced_list.copy()
            reduced_list = reduced_list[:start_comb]

            # event_df = pd.DataFrame(soc[start:min_idx])
            # soc_df = pd.DataFrame(soc)

            if end_comb + 1 <= len(copy_list):
                reduced_list.extend(copy_list[end_comb + 1:])
            if len(reduced_list) > 0:
                min_soc, min_idx = min(reduced_list, key=lambda x: x[0])
            else:
                break
    return events


def preprocess_schedule(this_sched, this_args, electrified_stations=None):
    Trip.consumption = Consumption(this_sched.vehicle_types,
                                   outside_temperatures=args.outside_temperature_over_day_path,
                                   level_of_loading_over_day=args.level_of_loading_over_day_path)
    this_sched.stations = electrified_stations
    # filter trips according to args
    this_sched.calculate_consumption()
    # this_sched.set_charging_type(this_args.preferred_charging_type)
    #
    # # (re)calculate the change in SoC for every trip
    # # charging types may have changed which may impact battery capacity
    # # while mileage is assumed to stay constant
    # this_sched.delta_soc_all_trips()
    # each rotation is assigned a vehicle ID
    this_sched.assign_vehicles()

    return this_sched, this_sched.generate_scenario(this_args)


def run_schedule(this_sched, this_args, electrified_stations=None, cost_calc=False):
    this_sched2 = copy(this_sched)
    this_sched2.stations = electrified_stations
    this_sched2.assign_vehicles()

    this_sched2, new_scen = preprocess_schedule(this_sched2, this_args,
                                                electrified_stations=electrified_stations)

    # Dont print output from spice ev to reduce clutter
    print(".", end="")
    print("Running Spice EV...", end="")

    sys.stdout = open(os.devnull, 'w')

    print("Running Spice EV...")
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)

        new_scen.run('distributed', vars(this_args).copy())

    print(" with % s rotations and % s vehicles", this_sched2.rotations,
          len(new_scen.vehicle_socs))
    sys.stdout = sys.__stdout__

    if args.cost_calculation and cost_calc:
        # cost calculation following directly after simulation
        try:
            with open(args.cost_parameters_file, encoding='utf-8') as f:
                cost_parameters_file = uncomment_json_file(f)
        except FileNotFoundError:
            raise SystemExit(f"Path to cost parameters ({args.cost_parameters_file}) "
                             "does not exist. Exiting...")
        calculate_costs(cost_parameters_file, new_scen, this_sched2, args)
    print(".")

    return this_sched2, new_scen


def charging_curve_to_soc_over_time(charging_curve, capacity, max_charge_from_grid=float('inf'),
                                    timestep=0.1, efficiency=1):
    # Charging curve as nested list of SOC, Power[kW] and capacity in [kWh]
    # Simple numeric creation of power over time --> to energy over time
    normalized_curve = np.array([[soc, power / capacity] for soc, power in charging_curve])
    soc = 0
    time = 0
    socs = []
    times = []
    while soc < args.desired_soc_opps:
        times.append(time)
        socs.append(soc)
        power1 = min(np.interp(soc, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        soc2 = soc + timestep / 60 * power1
        power2 = min(np.interp(soc2, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        power = (power1 + power2) / 2 * efficiency
        soc += timestep / 60 * power
        time += timestep
    # Fill the soc completely in last timestep
    times.append(time)
    socs.append(args.desired_soc_opps)
    return np.array((times, socs)).T


def get_delta_soc(soc_over_time_curve, soc, time_delta):
    # Returns expected soc lift for a given start_soc and time_delta.
    # Units for time_delta and time_curve are assumed to be the same, e.g. minutes
    # First element which is bigger than current soc
    if time_delta == 0:
        return 0
    soc = max(min(args.desired_soc_opps, soc), 0)
    first_time, start_soc = soc_over_time_curve[soc_over_time_curve[:, 1] >= soc][0, :]
    second_time = first_time + time_delta
    # Catch out of bounds if time of charging end is bigger than table values

    if second_time >= soc_over_time_curve[-1, 0]:
        end_soc = args.desired_soc_opps
    else:
        end_soc = soc_over_time_curve[soc_over_time_curve[:, 0] >= second_time][0, 1]

    # Make sure to limit delta soc to 1 if negative socs are given. They are possible during
    # the optimization process but will be continuously raised until they are >0.
    return min(args.desired_soc_opps, end_soc - start_soc)


def get_buffer_time_old(search_time, args):
    try:
        return timedelta(minutes=float(args.default_buffer_time_opps))
    except TypeError:
        for window, buffer_time in args.default_buffer_time_opps.items():
            try:
                start, end = window.split("-")
                if float(end) > search_time.hour >= float(start):
                    return timedelta(minutes=buffer_time)
            except ValueError:
                return timedelta(minutes=0)


def electrify_station(stat, stations, electrified_set):
    stations[stat] = config.standard_opp_station
    electrified_set.add(stat)


def get_buffer_time(trip, args):
    # a= get_buffer_time_old(trip.arrival_time, args)
    b = timedelta(minutes=get_buffer_time_spice_ev(trip, args.default_buffer_time_opps))
    return b


def combs_unordered_no_putting_back(n: int, k: int):
    try:
        return math.factorial(n) / ((math.factorial(n - k)) * math.factorial(k))
    except ValueError:
        LOGGER.warning("Value Error")
        return 0


def remove_none_socs(this_scen):
    # Make sure no None values exists in SOCs. Fill later values with last value
    # which was not None, eg, [1,5,4,None,None] becomes [1,5,4,4,4]
    for id, soc in this_scen.vehicle_socs.items():
        soc = np.array(soc)
        last_not_none = "not found"
        if None in soc:
            for ii in range(len(soc) - 1, -1, -1):
                if soc[ii] is not None:
                    last_not_none = soc[ii]
                    break
            soc[soc == np.array(None)] = last_not_none
        this_scen.vehicle_socs[id] = soc


def get_vehicle(rot_id, this_sched):
    return this_sched.rotations[rot_id].vehicle_id


def get_vehicle_rots(rot_id, this_sched):
    return [rot for rot in this_sched.rotations if rot_id == get_vehicle(rot, this_sched)]


def station_eval_to_logger(missing_energy, station_eval):
    LOGGER.debug("Missing energy: %s", missing_energy)
    if LOGGER.getEffectiveLevel() > logging.DEBUG:
        for stat_id in station_eval:
            LOGGER.debug("%s , %s", stat_id[0], stat_id[1]["pot_sum"])


def toolbox_from_pickle(sched_name, scen_name, args_name):
    with open(args_name, "rb") as f:
        args = pickle.load(f)
    with open(scen_name, "rb") as f:
        scen = pickle.load(f)
    with open(sched_name, "rb") as f:
        sched = pickle.load(f)
    return sched, scen, args


def toolbox_to_pickle(name, sched, scen, args):
    args_name = "args_" + name + ".pickle"
    with open(args_name, "wb") as f:
        pickle.dump(args, f)
    scen_name = "scenario_" + name + ".pickle"
    with open(scen_name, "wb") as f:
        pickle.dump(scen, f)
    sched_name = "schedule_" + name + ".pickle"
    with open(sched_name, "wb") as f:
        pickle.dump(sched, f)
    return sched_name, scen_name, args_name


def get_negative_rotations_all_electrified(this_scen, this_sched, soc_upper_thresh,
                                           soc_charge_curve_dict,
                                           not_possible_stations=None,
                                           soc_lower_thresh=0, relative_soc=False):
    if not not_possible_stations:
        not_possible_stations = set()

    events = get_below_zero_soc_events(this_scen=this_scen,
                                       rotations=this_sched.rotations.keys(),
                                       this_sched=this_sched,
                                       soc_upper_thresh=soc_upper_thresh,
                                       filter_standing_time=True,
                                       not_possible_stations=not_possible_stations,
                                       soc_lower_thresh=soc_lower_thresh, relative_soc=relative_soc)

    stats = {s for e in events for s in e["stations_list"] if s not in not_possible_stations}
    electrified_station_set = set(stats)
    vehicle_socs = timeseries_calc('best_station_ids[0]', this_sched.rotations.values(),
                                   this_scen.vehicle_socs,
                                   this_scen, electrified_station_set,
                                   soc_charge_curve_dict=soc_charge_curve_dict,
                                   soc_upper_thresh=args.desired_soc_deps)
    new_events = get_below_zero_soc_events(this_scen, this_sched.rotations,
                                           this_sched, soc_lower_thresh=soc_lower_thresh,
                                           soc_data=vehicle_socs)
    return [e["roation"].id for e in new_events]


# Electrify everything minus 1 station. If without the stations there are below zero events
# it is a must have station
def get_must_stations(this_scen, this_sched, soc_upper_thresh, soc_charge_curve_dict,
                      not_possible_stations=None,
                      soc_lower_thresh=0, relative_soc=False):
    events = get_below_zero_soc_events(this_scen=this_scen,
                                       rotations=this_sched.rotations.keys(),
                                       this_sched=this_sched,
                                       soc_upper_thresh=soc_upper_thresh,
                                       filter_standing_time=True,
                                       not_possible_stations=not_possible_stations,
                                       soc_lower_thresh=soc_lower_thresh, relative_soc=relative_soc)
    if not not_possible_stations:
        not_possible_stations = set()
    stats = {s for e in events for s in e["stations_list"] if s not in not_possible_stations}
    electrified_station_set_all = set(stats)

    must_stations = set()
    for s in electrified_station_set_all:
        print(".", end="")
        electrified_station_set = electrified_station_set_all.difference([s])
        vehicle_socs = timeseries_calc('best_station_ids[0]', this_sched.rotations.values(),
                                       this_scen.vehicle_socs,
                                       this_scen, electrified_station_set,
                                       soc_charge_curve_dict=soc_charge_curve_dict,
                                       soc_upper_thresh=args.desired_soc_deps)

        new_events = get_below_zero_soc_events(this_scen, this_sched.rotations,
                                               this_sched, soc_lower_thresh=soc_lower_thresh,
                                               soc_data=vehicle_socs)
        if len(new_events) > 0:
            must_stations.add(s)
    return must_stations


if __name__ == "__main__":
    main()
