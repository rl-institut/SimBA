""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
from datetime import datetime, timedelta
import json
import os
import sys
import warnings
from copy import copy, deepcopy
import time
import logging

import src.scenario as scenario
import schedule
import rotation
import pickle
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Pool, freeze_support
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.trip import Trip

# Todo this implementation in c ase of changes in ebustoolbox
from ebus_toolbox.util import get_buffer_time as get_buffer_time_spice_ev

matplotlib.use("TkAgg")
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(message)s')
file_handler = logging.FileHandler('optimizer.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)
logger.addHandler(stream_handler)

with open("args_bvg_full_no_ele.pickle", "rb") as f: args = pickle.load(f)

with open("scen_bvg_test.pickle", "rb") as f: scen = pickle.load(f)
with open("sched_bvg_test.pickle", "rb") as f: sched = pickle.load(f)
ROT = '6808730'
# set battery and charging power
BATTERY_CAPACITY = 400
CHARGING_CURVE = [[0, 450], [0.8, 296], [0.9, 210], [1, 20]]
# CHARGING_CURVE = [[0, 450], [0.99, 20],[1,20]]
for name, type in sched.vehicle_types.items():
    for charge_type, vehicle in type.items():
        vehicle["capacity"] = BATTERY_CAPACITY
        vehicle["charging_curve"] = CHARGING_CURVE

CHARGE_EFF = 0.95
CHARGING_POWER = 250

timers = [0] * 10
args.min_charging_time = 0
args.default_buffer_time_opps = {"else": 0}

del args.save_soc
del args.save_timeseries
del args.save_results


def run_optimization(this_sched, this_scen, args,
                     opt_type="greedy", exclusion_rots=set(), inclusion_stations=set(),
                     exclusion_stations=set(), remove_impossible_rots=False, rebase_scenario=False,
                     node_choice="step-by-step",
                     **kwargs):
    """
    Optimizes scenario by adding electrified stations sparingly
    until scenario has no below 0 soc events.

    :param remove_impossible_rots: repeat optimization without impossible to electrify rotations
    default=False
    :type remove_impossible_rots: bool

    :param rebase_scenario: should the given scenario be simulated with the given boundries before
    optimization or can schedule and scenario be used right away. rebasing can be time consuming
    but can speed up optimization afterwards and increase optimization efficiency
    default=False
    :type rebase_scenario: bool

    :param this_sched: Simulation schedule containing buses, rotations etc.
    :type this_sched: ebus_toolbox.Schedule

    :param this_scen: Simulation scenario containing simulation results
                     including the SoC of all vehicles over time
    :type this_scen: spice_ev.Scenario
    :param args: Simulation arguments for manipulation or generated outputs
    :type args: object
    :param opt_type: Type of the following optimizations ["greedy"]
    :type opt_type: str
    :param exclusion_rots: Rotations to be excluded from optimization
    :type exclusion_rots: set(str)
    :param inclusion_stations: Stations which have to be electrified
    :type inclusion_stations: set(str)
    :param exclusion_stations: Stations to be excluded from optimization, eg. cant be electrified
    :type exclusion_stations: set(str)

    :return: (Schedule,Scenario) optimizied schedule and Scenario
    :rtype: tuple(ebus_toolbox.Schedule, spice_ev.Scenario)
    """
    now = datetime.now()
    new_ele_stations_path = "data/examples/optimized_stations_" + \
                            now.strftime("%Y_%m_%d_%H_%M") + ".json"

    not_possible_stations = inclusion_stations.union(exclusion_stations)

    s = time.time()
    # Calculate base scenario
    # this_sched.rotations = {id: rot for id, rot in this_sched.rotations.items() if id in rots}
    if rebase_scenario:
        logger.debug(f"Spice EV Rebasing Scenario")
        new_sched, new_scen, ele_station_set, ele_stations = preprocessing_scenario(
            this_sched, args,
            inclusion_stations,
            exclusion_rots=exclusion_rots, run_only_neg=False)

        logger.debug(f"Rebasing took {time.time() - s} sec")
    else:
        with open(args.electrified_stations, "r") as f:
            ele_stations = json.load(f)
        ele_station_set = set()
        # Electrify inclusion stations
        for stat in inclusion_stations:
            electrify_station(stat, ele_stations, ele_station_set)
        new_scen = this_scen
        new_sched = this_sched

    s = time.time()
    i = 0
    global timer_for_calc
    while True and i < 2:
        if opt_type == "greedy" or opt_type == "deep":
            logger.debug("Starting greedy optimization")
            ele_stations, ele_station_set, could_not_be_electrified, list_greedy_sets = \
                optimization_loop(ele_stations, ele_station_set, new_scen, new_sched,
                                  not_possible_stations, soc_upper_thresh=1, soc_lower_thresh=0,
                                  solver=kwargs.get("solver", "quick"), opt_type=opt_type, node_choice=node_choice)

        i += 1
        if not remove_impossible_rots or len(could_not_be_electrified) == 0:
            break
        else:
            logger.debug(f"Non Spice Ev methods took {time.time() - s - timer_for_calc} sec")
            logger.debug(f"Spice ev took {timer_for_calc} s")
            logger.debug(f"Electrified Stations: {len(ele_station_set)}")
            logger.debug(ele_station_set)
            logger.warning(f"These rotations could not be electrified"
                           f"\n {could_not_be_electrified}")

            logger.warning("Removing impossible rots, rebasing and restarting optimization")
            exclusion_rots.update(could_not_be_electrified)
            new_sched, new_scen, ele_station_set, ele_stations = preprocessing_scenario(
                this_sched, args,
                inclusion_stations,
                exclusion_rots=exclusion_rots, run_only_neg=False)

    logger.debug(f"Non Spice Ev methods took {time.time() - s - timer_for_calc} sec")
    logger.debug(f"Spice ev took {timer_for_calc} s")
    logger.debug(f"Electrified Stations: {len(ele_station_set)}")
    logger.debug(ele_station_set)
    logger.debug(could_not_be_electrified)
    plot_(get_rotation_soc(ROT, new_sched, new_scen)[0])
    new_scen.vehicle_socs = timeseries_calc('best_station_ids[0]', new_sched.rotations.values(),
                                            new_scen.vehicle_socs,
                                            new_scen, ele_station_set)
    plot_(get_rotation_soc(ROT, new_sched, new_scen)[0])
    q = get_rotation_soc(ROT, new_sched, new_scen)[0]
    global timers
    logger.debug(timers)
    with open(new_ele_stations_path, "w") as f:
        json.dump(ele_stations, f, indent=2)

    logger.debug(f"Spice EV is calculating optimized case as a complete scenario")
    new_sched, new_scen, ele_station_set, ele_stations = preprocessing_scenario(
        this_sched, args, electrified_stations=ele_stations, run_only_neg=False,
        electrified_station_set=ele_station_set)

    ax = plot_(get_rotation_soc(ROT, new_sched, new_scen)[0])
    ax.plot(q)
    logger.debug(f"Still negative rotations:{new_sched.get_negative_rotations(new_scen)}")
    print("Finished")
    return new_sched, new_scen


def main():
    # which rotations should be excluded?
    # Teilmenge der 197 die nicht möglich ist, muss ignoriert werden
    exclusion_rots = {'6891222', '6891224', '6891226', '6891228', '7194397',
                      # Linie 901 wird auch ignoriert
                      '6675617', '6675640', '6675626',
                      '6675645', '6675631', '6675648',
                      '6675636', '6675653', '6675661', '6675657'
                      }

    # Ignoieren von umöglichen Umläufen, die existieren ohne Ausschluss von Stationen
    impossible_rots = {'6895194', '6895199', '6895223', '6895226', '6923844', '7049422', '7049423'}

    exclusion_rots.update(impossible_rots)
    exclusion_rots = set()
    exclusion_rots = {'6891228', '7194397', '6891226', '6891222', '7049423', '7049422', '6891224'}
    # which stations have to be electrified?
    inclusion_stations = set()  # {'SWEL01B'} #{'SMZ07B', 'SWEL01B', 'URUD08B'}

    # sched.rotations={ROT:sched.rotations[ROT]}

    exclusion_stations = set()
    t = time.time()
    run_optimization(sched, scen, args, opt_type="deep", exclusion_rots=exclusion_rots,
                     inclusion_stations=inclusion_stations, exclusion_stations=exclusion_stations,
                     remove_impossible_rots=False, rebase_scenario=False, solver="quick",
                     node_choice="brute")
    print(f"Opt took {time.time() - t}")


timer_for_calc = 0


def optimization_loop(electrified_stations, electrified_station_set, new_scen, new_sched,
                      not_possible_stations, soc_upper_thresh=1, soc_lower_thresh=0,
                      decision_tree=dict(), pre_optimzed_set=None, opt_type="greedy", **kwargs):
    # Base stations for optimization, so inclusion of stations can be skipped
    base_stations = electrified_stations.copy()
    base_electrified_station_set = electrified_station_set.copy()

    # Connect events with shared stations
    ########
    could_not_be_electrified = set()
    base_scen = copy(new_scen)
    base_sched = copy(new_sched)

    ###########
    # Make sure no None values exists in SOCs. Fill later values with last value
    # which was not None
    for id, soc in base_scen.vehicle_socs.items():
        soc = np.array(soc)
        if None in soc:
            for ii in range(len(soc) - 1, -1, -1):
                if soc[ii] is not None:
                    last_not_none = soc[ii]
                    break
            soc[soc == np.array(None)] = last_not_none
        base_scen.vehicle_socs[id] = soc

    # Get events where soc fell below 0. The events contain info about the problematic
    # timespan, which includes stations which could provide a soc lift
    events = get_below_zero_soc_events(this_scen=base_scen,
                                       rotations=list(new_sched.rotations.keys()),
                                       this_sched=base_sched,
                                       soc_upper_thresh=0.95, filter_standing_time=True,
                                       not_possible_stations=not_possible_stations,
                                       soc_lower_thresh=0)

    # Check if the events can be divided into subgroups which are independent
    # this makes optimization in smaller groups possible
    groups = get_groups_from_events(events, not_possible_stations)

    # Baseline greedy Optimization
    list_greedy_sets = [set()] * len(groups)

    # Base line is created simply by not having a decision tree and not a pre optimized_set yet
    for group_nr, group in enumerate(groups):
        events, stations = group
        group_rots = {e["rotation"].id: e["rotation"] for e in events}
        group_socs = [e["min_soc"] for e in events]
        linien = {lne for e in events for lne in e["rotation"].lines}
        electrified_stations = base_stations.copy()
        electrified_station_set = base_electrified_station_set.copy()
        print(f"Optimizing {group_nr + 1} out of {len(groups)}. This includes these Lines")
        logger.debug(linien)
        print(len(events), "events")
        solver = kwargs.get("solver", "spiceev")
        node_choice = kwargs.get("node_choice", "step-by-step")
        if node_choice == "brute":
            choice_func = choose_station_brute
        else:
            choice_func = choose_station_step_by_step

        if solver == "quick":
            # first run is always step by step
            group_optimization_quick(group, base_scen, base_sched,
                                     electrified_stations, electrified_station_set,
                                     could_not_be_electrified, not_possible_stations,
                                     choose_station_step_by_step,
                                     pre_optimzed_set=pre_optimzed_set,
                                     decision_tree=decision_tree, base_group=events,
                                     **kwargs)
        if solver == "quick":
            if opt_type == "deep":
                sols = []
                while True:
                    print(len(decision_tree))
                    not_possible_stations = copy(not_possible_stations)
                    pre_optimized_set = copy(electrified_station_set)

                    electrified_stations = base_stations.copy()
                    # electrified_station_set = base_electrified_station_set.copy()
                    new_electrified_set = set()
                    group_optimization_quick(group, base_scen, base_sched,
                                             electrified_stations, new_electrified_set,
                                             could_not_be_electrified
                                             , not_possible_stations,
                                             choice_func,
                                             pre_optimzed_set=pre_optimized_set,
                                             decision_tree=decision_tree,
                                             base_group=events)

                    print("optimized with", len(new_electrified_set), " stations",
                          stations_hash(new_electrified_set))
                    sols.append(new_electrified_set)


        else:
            # use spiceev
            group_optimization(group, base_scen, base_sched,
                               electrified_stations, electrified_station_set,
                               could_not_be_electrified, not_possible_stations,
                               pre_optimzed_set=None,
                               decision_tree=decision_tree, brute=False, **kwargs)

        list_greedy_sets[group_nr] = electrified_station_set.copy()
        logger.debug("Optimized with {} stations out of {}".format(len(electrified_station_set),
                                                                   len(stations)))

    for single_set in list_greedy_sets:
        for stat in single_set:
            electrify_station(stat, electrified_stations, electrified_station_set)

    return electrified_stations, electrified_station_set, could_not_be_electrified, list_greedy_sets




# ########
# # with open("dec_tree.pickle", "rb") as f:
# #     decision_tree=pickle.load(f)
# #
# # with open("pre_opt.pickle", "rb") as f:
# #     electrified_station_set =pickle.load(f)
#
# while True:
#     print(len(decision_tree))
#     not_possible_stations = inclusion_stations.union(exclusion_stations)
#     pre_optimized_set = electrified_station_set
#     group = groups[19]
#     electrified_stations = base_stations.copy()
#     electrified_station_set = base_electrified_station_set.copy()
#
#     group_optimization(group, base_scen, base_sched,
#                        electrified_stations, electrified_station_set, could_not_be_electrified
#                        , not_possible_stations, pre_optimzed_set=pre_optimized_set,
#                        decision_tree=decision_tree, brute=True)
#     print("optimized with", len(electrified_station_set), " stations")
#     sols.append(electrified_station_set)
# ########
#
# # Save Optimization
# for s in list_greedy_sets:
#     for stat in s:
#         electrify_station(stat, electrified_stations, electrified_station_set)
# with open("400kWh_no_inc_new_raus_ohne_197_ohne_imp_ele_station.json", "w") as f:
#     json.dump(electrified_stations, f, indent=2)
#
# electrified_stations = base_stations.copy()
# electrified_station_set = base_electrified_station_set
# #
# # list_opt_sets=[set()]*len(groups)
# #
# # for group_nr, group in enumerate(groups):
#     min_length=float("inf")
#     event_group, _ = group
#     station_eval = list(evaluate(event_group, base_scen).items())
#     station_eval.reverse()
#
#     # Expand greedy Tree from each node of greedy Optimization
#     excl_options = list(list_greedy_sets[group_nr])
#     pre_optimzed_set =list_greedy_sets[group_nr]
#
#     pack = (group, base_scen, base_sched,
#                                 old_stations, old_electrified_station_set,
#                                 could_not_be_electrified, pre_optimzed_set)
#     with Pool(5) as p:
#                      # Outer optimization func,      Repeating input,Stations to be excluded
#         results = p.starmap(outer_group_optimzation,zip(repeat(pack), excl_options))
#
#     for res in results:
#         if len(res)< min_length:
#             min_length=len(res)
#             list_opt_sets[group_nr]=res
#             print(min_length)
#
#
# greedy_stations_count=sum([len(s) for s in list_greedy_sets])
# list_opt_count = sum([len(s) for s in list_opt_sets])
# print(f"Greedy Results in {greedy_stations_count} stations")
# print(f"Opt Results in {list_opt_count} stations")
# print(list_greedy_sets)
# print(list_opt_sets)


def outer_group_optimzation(pack, not_possible_stations):
    group, base_scen, base_sched, old_stations, \
    old_electrified_station_set, could_not_be_electrified, pre_optimzed_set = pack
    electrified_stations = old_stations.copy()
    electrified_station_set = old_electrified_station_set.copy()
    ele_set = group_optimization(group, base_scen, base_sched,
                                 electrified_stations, electrified_station_set,
                                 could_not_be_electrified,
                                 not_possible_stations, pre_optimzed_set=pre_optimzed_set)
    return electrified_station_set.copy()


def get_groups_from_events(events, not_possible_stations=set(), could_not_be_electrified=set()):
    # First create simple list of station sets for single events.
    # Electrified and other not possible to electrify stations should not connect groups
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
            logger.warning(f'Didnt find rotation {e["rotation"].id} in any subset'
                           f'of possible electrifiable stations')
            # this event will no show up in an event_group.
            # Therefor it needs to be put into this set
            could_not_be_electrified.update([e["rotation"].id])

    groups = list(zip(event_groups, station_subsets))
    return sorted(groups, key=lambda x: len(x[1]))


def group_optimization_quick(group, base_scen, base_sched,
                             electrified_stations, electrified_station_set,
                             could_not_be_electrified,
                             not_possible_stations, choose_station_function,
                             pre_optimzed_set=None, decision_tree=None,
                             brute=False, **kwargs):
    event_group, _ = group

    # Loading from pickle faster than deepcopy. Copy enough?
    new_scen = copy(base_scen)
    new_scen.vehicle_socs = deepcopy(base_scen.vehicle_socs)

    lifted_socs = kwargs.get("lifted_socs", None)
    if lifted_socs is not None:
        new_scen.vehicle_socs = lifted_socs
    new_sched = copy(base_sched)

    missing_soc = 0
    rotation_dict = dict()
    for e in event_group:
        rotation_dict[e["rotation"].id] = e["rotation"]
        missing_soc += e["min_soc"]

    if missing_soc >= 0:
        logger.debug("Already electrified: Returning set")
        return electrified_stations
    missing_energy = missing_soc * BATTERY_CAPACITY

    station_eval = evaluate(event_group, new_scen)
    for id in station_eval:
        print(id[0], id[1]["pot_sum"])
    logger.debug(missing_energy)
    best_station_ids, recursive = choose_station_function(station_eval, electrified_station_set,
                                               pre_optimzed_set, decision_tree,
                                               missing_energy=missing_energy)

    logger.debug(best_station_ids)
    if best_station_ids is None:
        logger.warning(
            f"All stations with estimated potential electrified but still missing energy in {len(list(rotation_dict.keys()))} rotations")
        if pre_optimzed_set is not None:
            # Remove electrified stations in this run
            c = electrified_station_set.copy()
            for stat in c:
                electrified_stations.pop(stat)
                electrified_station_set.remove(stat)
            # Overwrite with preoptimized set
            for stat in pre_optimzed_set:
                electrify_station(stat, electrified_stations, electrified_station_set)
        else:
            could_not_be_electrified.update(list(rotation_dict.keys()))
        return electrified_stations

    # Electrify station
    for id in best_station_ids:
        electrify_station(id, electrified_stations, electrified_station_set)

    # pre_opt_sched.rotations = rotation_dict

    s = time.time()

    # Using the base_group for timeseries calculation is a little slower than using the current
    # (smaller) event group which is getting optimized. But it allows for looking at the current
    # base group missing energy see "delta_base_energy" which is put into
    base_group = kwargs.get("base_group", [])
    event_rotations = {x["rotation"] for x in base_group}

    new_scen.vehicle_socs = deepcopy(base_scen.vehicle_socs)
    new_scen.vehicle_socs = timeseries_calc(best_station_ids[0], event_rotations,
                                            new_scen.vehicle_socs,
                                            new_scen, electrified_station_set)
    lifted_socs = deepcopy(new_scen.vehicle_socs)

    global timer_for_calc
    global timers
    timer_for_calc += time.time() - s
    timers[2] += time.time() - s
    not_possible_stations = set(electrified_stations.keys()).union(not_possible_stations)
    event_rotations = {event["rotation"].id for event in event_group}
    new_events = get_below_zero_soc_events(new_scen, event_rotations,
                                           new_sched,
                                           soc_upper_thresh=1,
                                           filter_standing_time=True,
                                           not_possible_stations=not_possible_stations,
                                           soc_lower_thresh=0)

    delta_energy = get_missing_energy(new_events, BATTERY_CAPACITY)

    event_rotations = {event["rotation"].id for event in base_group}
    base_events = get_below_zero_soc_events(new_scen, event_rotations,
                                            new_sched,
                                            soc_upper_thresh=1,
                                            filter_standing_time=True,
                                            not_possible_stations=set(),
                                            soc_lower_thresh=0)

    delta_base_energy = get_missing_energy(base_events, BATTERY_CAPACITY)

    logger.debug(delta_energy)
    if decision_tree is not None:
        node_name = stations_hash(electrified_station_set)
        try:
            decision_tree[node_name]["missing_energy"] = delta_base_energy
            decision_tree[node_name]["visit_counter"] += 1
            logger.debug("already visited")
        except KeyError:
            decision_tree[node_name] = dict()
            decision_tree[node_name]["missing_energy"] = delta_base_energy
            decision_tree[node_name]["visit_counter"] = 1

    # Everything electrified
    if delta_energy >= 0:
        return electrified_stations

    # Some choice functions might not need a recursive call, they return here. recursive is set
    # by the choose_station_function
    if not recursive:
        return



    # Check if the events can be divided into subgroups which are independent
    groups = get_groups_from_events(new_events, not_possible_stations, could_not_be_electrified)

    for group in groups:
        new_stations = group_optimization_quick(group, base_scen, base_sched,
                                                electrified_stations,
                                                electrified_station_set,
                                                could_not_be_electrified,
                                                not_possible_stations, choose_station_function,
                                                pre_optimzed_set, decision_tree,
                                                lifted_socs=lifted_socs, base_group=base_group)
        electrified_stations.update(new_stations)

    return electrified_stations


def get_missing_energy(events, cap):
    missing_soc = 0
    for e in events:
        missing_soc += e["min_soc"]
    return missing_soc * cap


def group_optimization(group, base_scen, base_sched,
                       electrified_stations, electrified_station_set, could_not_be_electrified,
                       not_possible_stations, pre_optimzed_set=None, decision_tree=None,
                       brute=False, **kwargs):
    event_group, _ = group

    # Loading from pickle faster than deepcopy. Copy enough?
    pre_opt_scen = copy(base_scen)
    pre_opt_scen.vehicle_socs = copy(base_scen.vehicle_socs)
    pre_opt_sched = copy(base_sched)
    missing_soc = 0
    rotation_dict = dict()
    for e in event_group:
        rotation_dict[e["rotation"].id] = e["rotation"]
        missing_soc += e["min_soc"]

    if missing_soc >= 0:
        logger.debug("Already electrified: Returning set")
        return electrified_stations

    missing_energy = missing_soc * BATTERY_CAPACITY

    station_eval = evaluate(event_group, pre_opt_scen)
    for id in station_eval:
        logger.debug("%s, %s", id[0], id[1]["pot_sum"])
    logger.debug(missing_energy)

    # best_station_ids = choose_stations_function(station_eval, electrified_station_set,
    #                                             pre_optimzed_set, decision_tree,
    #                                             missing_energy=missing_energy)
    if brute:
        best_station_ids = choose_station_brute(station_eval, electrified_station_set,
                                                pre_optimzed_set, decision_tree,
                                                missing_energy=missing_energy)
    else:
        best_station_ids = choose_station_step_by_step(station_eval, electrified_station_set,
                                                       pre_optimzed_set, decision_tree,
                                                       missing_energy=missing_energy)

    logger.debug(best_station_ids)
    if best_station_ids is None:
        print(
            f"All stations with estimated potential electrified but still missing energy in {len(list(rotation_dict.keys()))} rotations")
        if pre_optimzed_set is not None:
            # Remove electrified stations in this run
            c = electrified_station_set.copy()
            for stat in c:
                electrified_stations.pop(stat)
                electrified_station_set.remove(stat)
            # Overwrite with preoptimized set
            for stat in pre_optimzed_set:
                electrify_station(stat, electrified_stations, electrified_station_set)
        else:
            could_not_be_electrified.update(list(rotation_dict.keys()))
        return electrified_stations

    # Electrify station
    for id in best_station_ids:
        electrify_station(id, electrified_stations, electrified_station_set)
    pre_opt_sched.rotations = rotation_dict

    s = time.time()

    new_sched, new_scen = run_schedule(pre_opt_sched, args,
                                       electrified_stations=electrified_stations)

    global timer_for_calc
    timer_for_calc += time.time() - s

    not_possible_stations = set(electrified_stations.keys()).union(not_possible_stations)
    new_events = get_below_zero_soc_events(new_scen, list(new_sched.rotations.keys()),
                                           new_sched,
                                           soc_upper_thresh=1,
                                           filter_standing_time=True,
                                           not_possible_stations=not_possible_stations,
                                           soc_lower_thresh=0)
    missing_soc = 0
    for e in new_events:
        missing_soc += e["min_soc"]
    delta_energy = missing_soc * BATTERY_CAPACITY

    if decision_tree is not None:
        node_name = str(sorted(list(electrified_station_set)))
        try:
            decision_tree[node_name]["missing_energy"] = delta_energy
            decision_tree[node_name]["visit_counter"] += 1
            print("already visited")
        except KeyError:
            decision_tree[node_name] = dict()
            decision_tree[node_name]["missing_energy"] = delta_energy
            decision_tree[node_name]["visit_counter"] = 1

    # Everything electrified
    if delta_energy >= 0:
        return electrified_stations
    else:
        if brute:
            return

    # Check if the events can be divided into subgroups which are independent
    groups = get_groups_from_events(new_events, not_possible_stations, could_not_be_electrified)

    for group in groups:
        new_stations = group_optimization(group, new_scen, new_sched,
                                          electrified_stations,
                                          electrified_station_set,
                                          could_not_be_electrified,
                                          not_possible_stations,
                                          pre_optimzed_set, decision_tree)
        electrified_stations.update(new_stations)

    return electrified_stations


def preprocessing_scenario(this_sched, this_args,
                           inclusion_stations=set(),
                           electrified_stations=None,
                           electrified_station_set=set(),
                           exclusion_rots=set(),
                           run_only_neg=False):
    if electrified_stations is None:
        with open(this_args.electrified_stations, "r") as f:
            electrified_stations = json.load(f)

    # Electrify inclusion stations
    for stat in inclusion_stations:
        electrify_station(stat, electrified_stations, electrified_station_set)

    # Calc new but only prev. negative rotations
    if run_only_neg:
        rots = this_sched.get_negative_rotations(scen)
        rots = {r: this_sched.rotations[r] for r in rots if r not in exclusion_rots}
        this_sched.rotations = rots
    else:
        rots = {r: this_sched.rotations[r] for r in this_sched.rotations if r not in exclusion_rots}
        this_sched.rotations = rots

    new_sched, new_scen = run_schedule(this_sched, this_args, electrified_stations)
    return new_sched, new_scen, electrified_station_set, electrified_stations


def combination_generator(iterable, amount: int, _pass_list=list(), _counter=0):
    """ Generator which yields all possible combinations of choosing
    an amount out of an iterable without putting them back and without caring about the
    order of elements
    :param iterable: Any collection which can be cast to a list
    :param amount: Number of elements which should be drawn from iterable
    :type amount: int
    """
    iterable = list(iterable)
    # Initialization of list with length of amount
    if len(_pass_list) == 0:
        _pass_list = [0] * amount
    for i, item in enumerate(iterable):
        # Recursive calling of generator with clock like behavior, e.g right-most item changes until
        # end of list is reached. This leads to a change in the item left to it and so on. Elements
        # on the right can only change to a subset of bigger indicies than their left counter-part.
        # This is due to the ignoring of order, which reduces the amount of possibilities.
        _pass_list[_counter] = item
        if amount <= 1:
            yield _pass_list
        else:
            for gen in combination_generator(iterable[i + 1:], amount - 1, _pass_list,
                                             _counter=_counter + 1):
                yield gen


def no_optimization():
    return "converged"


def plot_(data):
    """ Simple plot of data without having to create subplots"""
    fig, ax = plt.subplots()
    ax.plot(data, linewidth=2.0)
    return ax


def join_all_subsets(subsets):
    joined_subset = True
    while joined_subset:
        joined_subset, subsets = join_subsets(subsets)
    return subsets


#  ToDo Further implement
def choose_station_brute(station_eval, electrified_station_set,
                         pre_optimzed_set=None, decision_tree=None, missing_energy=0):
    station_ids = [x[0] for x in station_eval]
    a = combination_generator(station_ids, len(pre_optimzed_set) - 1)
    for comb in a:
        node_name=stations_hash(comb)
        if node_name not in decision_tree:
            return comb, False
    else:
        print("calculated all possibilities")
        return None, False


def choose_station_step_by_step(station_eval, electrified_station_set,
                                pre_optimzed_set=None, decision_tree=None, missing_energy=0):
    # Filter functions to stop simulating cases which have no hope of being optimal.
    # If in optimization mode, optimization can break if station amount is superceded
    # This filter is done better by the next
    # if pre_optimzed_set is not None:
    #     if len(electrified_station_set)>len(pre_optimzed_set):
    #         return pre_optimzed_set
    # Potentials have to be at least as promising as the pre-optimized case
    if pre_optimzed_set is not None:
        delta = len(pre_optimzed_set) - len(electrified_station_set)
        pot = 0
        for i in range(0, min(delta, len(station_eval))):
            pot += station_eval[i][1]["pot_sum"]
        if pot <= -missing_energy:
            print("Not enough potential after ", stations_hash(electrified_station_set))
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
            # Create a station combination from already electrified stations and possible new station
            check_stations = electrified_station_set.union([station[0]])
            if decision_tree[stations_hash(check_stations)]["visit_counter"] == min_count_visited:
                best_station_id = station[0]
                return [best_station_id]
    return None


def stations_hash(stations_set):
    return str(sorted(list(stations_set)))


def timeseries_calc(station, rotations, soc_dict, eval_scen, ele_station_set):
    global timers
    ele_stations = set([*ele_station_set, station])
    s2 = time.time()
    soc_dict = copy(soc_dict)
    CHARGING_CURVE_EFF = [[soc, power] for soc, power in CHARGING_CURVE]

    soc_over_time_curve = charging_curve_to_soc_over_time(CHARGING_CURVE_EFF, BATTERY_CAPACITY,
                                                          CHARGING_POWER, efficiency=CHARGE_EFF,
                                                          timestep=0.1)
    for rot in rotations:
        soc = soc_dict[rot.vehicle_id]
        for i, trip in enumerate(rot.trips):
            s = time.time()
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

            soc_max = np.max(soc)
            timers[0] += time.time() - s
            s = time.time()
            while soc_max > 1:
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
            timers[1] += time.time() - s
        soc_dict[rot.vehicle_id] = soc
    timers[3] += (time.time() - s2) * 100
    return soc_dict


def evaluate(events, eval_scen, soc_upper_thresh=1, soc_lower_thresh=0,
             electrified_station_set=set(), decision_tree=None):
    # Analyse stations for "helpful" energy supply. Energy supply is helpful the minimal soc of an
    # event is raised up to a minimal soc (probably zero). The supplied energy is approximated by
    # loading power, standing time at a station, soc at station, minimal soc of the event
    station_eval = dict()
    soc_over_time_curve = charging_curve_to_soc_over_time(CHARGING_CURVE, BATTERY_CAPACITY,
                                                          CHARGING_POWER)
    for e in events:
        for i, trip in enumerate(e["trip"]):
            # Station is only evaluated if station name is part of event stations
            # Only these stations showed potential in electrification, e.g enough standing time
            if trip.arrival_name not in e["stations"]:
                continue
            idx = get_index_by_time(trip.arrival_time, eval_scen)
            soc = eval_scen.vehicle_socs[e["vehicle_id"]][idx]
            # ToDo define elsewhere
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

            # ToDo get Vehicle Battery Capacity and charging power
            capacity = BATTERY_CAPACITY

            # energy_charging_potential = standing_time_min *60 * ch_power
            energy_charging_potential = get_delta_soc(soc_over_time_curve, soc, standing_time_min) \
                                        * capacity

            # Potential is at max the minimum between the useful delta soc * capacity or the
            # energy provided by charging for the full standing time
            delta_E_pot = min(delta_soc_pot * capacity, energy_charging_potential)
            d = dict(E_pot=delta_E_pot,
                     standing_time=timedelta(minutes=standing_time_min))
            try:
                station_eval[trip.arrival_name]["pot_list"].append(d)
            except:
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
    delay = 0
    standing_time_min = (trip2.departure_time - trip1.arrival_time) \
                        / timedelta(minutes=1)
    if args.min_charging_time > standing_time_min:
        return 0
    # Todo trip1 or trip2
    buffer_time = (get_buffer_time(trip1, args) / timedelta(minutes=1))
    if buffer_time > 0:
        standing_time_min -= buffer_time
    else:
        standing_time_min -= delay
    return max(0, standing_time_min)


def join_subsets(subsets):
    subsets = [s.copy() for s in subsets]
    for i in range(len(subsets)):
        for ii in range(len(subsets)):
            if i == ii: continue
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


def get_rotation_soc(rot_id, this_sched, this_scen):
    rot = this_sched.rotations[rot_id]
    rot_start_idx = get_index_by_time(rot.departure_time, this_scen)
    rot_end_idx = get_index_by_time(rot.arrival_time, this_scen)
    return this_scen.vehicle_socs[rot.vehicle_id], rot_start_idx, rot_end_idx


def get_below_zero_soc_events(this_scen: scenario.Scenario, rotations,
                              this_sched: schedule.Schedule,
                              soc_upper_thresh=0.9, filter_standing_time=True,
                              not_possible_stations=set(), soc_lower_thresh=0):
    # Create list of events which describe trips which end in a soc below zero
    # The event is bound by the lowest soc and an upper soc threshhold which is naturally 1
    # Properties before and after these points have no effect on the event itself, similar to
    # an event horizon
    events = []
    SOC_UPPER_THRESH = soc_upper_thresh
    count_electrified_rot = 0
    # pos_rots=pd.DataFrame(index=range(0,1549))
    # neg_rots=pd.DataFrame(index=range(0,1549))

    for rot_id in rotations:
        rot = this_sched.rotations[rot_id]
        soc, rot_start_idx, rot_end_idx = get_rotation_soc(rot_id, this_sched, this_scen)
        soc = [s if s is not None else 999 for s in soc]
        idx = range(0, len(soc))

        comb = list(zip(soc, idx))[rot_start_idx:rot_end_idx]
        min_soc, min_idx = min(comb, key=lambda x: x[0])
        reduced_list = comb.copy()
        start_soc = comb[0][0]
        # if rotation gets a start soc below 1 this should change below 0 soc events since fixing
        # the rotation before would lead to fixing this rotation
        # soc_lower_thresh_cur = min(start_soc, soc_upper_thresh) - (
        #             soc_upper_thresh - soc_lower_thresh)
        soc_lower_thresh_cur = soc_lower_thresh
        if min_soc >= soc_lower_thresh_cur:
            count_electrified_rot += 1
        #
        # if min_soc<soc_lower_thresh:
        #     neg_rots[rot_id]=None
        #     neg_rots.loc[rot_start_idx:rot_end_idx-1, rot_id]=soc[rot_start_idx:rot_end_idx]
        # else:
        #     pos_rots[rot_id]=None
        #     pos_rots.loc[rot_start_idx:rot_end_idx-1, rot_id]=soc[rot_start_idx:rot_end_idx]

        while min_soc < soc_lower_thresh_cur:
            i = min_idx
            idx = [x[1] for x in reduced_list]
            while soc[i] < SOC_UPPER_THRESH:
                if i == rot_start_idx: break
                i -= 1
            start_comb = idx.index(i)
            start = i
            i = min_idx
            while soc[i] < SOC_UPPER_THRESH:
                if i >= rot_end_idx - 1: break
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
            event = dict(start_idx=start, end_idx=min_idx,
                         min_soc=min_soc, stations=possible_stations,
                         vehicle_id=rot.vehicle_id, trip=trips,
                         rotation=rot, stations_list=possible_stations_list)
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

            # red_df = pd.DataFrame(reduced_list)
            # event_df.plot()
            # soc_df.plot()
            # soc_df[rot_start_idx: rot_end_idx].plot()
            # red_df.loc[:,0].plot()
            # plt.show()

    print(f"Last electrification electrified {count_electrified_rot} rotations")
    # neg_rots.plot()
    # pos_rots.plot()
    # plt.show()

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


def run_schedule(this_sched, this_args, electrified_stations=None):
    this_sched2 = copy(this_sched)
    this_sched2.stations = electrified_stations
    this_sched2.assign_vehicles()

    this_sched2, new_scen = preprocess_schedule(this_sched2, this_args,
                                                electrified_stations=electrified_stations)

    # Dont print output from spice ev to reduce clutter
    print(".", end="")
    sys.stdout = open(os.devnull, 'w')

    print("Running Spice EV...")
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)

        new_scen.run('distributed', vars(this_args).copy())

    sys.stdout = sys.__stdout__
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
    while soc < 1:
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
    socs.append(1)
    return np.array((times, socs)).T


def get_delta_soc(soc_over_time_curve, soc, time_delta):
    # Returns expected soc lift for a given start_soc and time_delta.
    # Units for time_delta and time_curve are assumed to be the same, e.g. minutes
    # First element which is bigger than current soc
    if time_delta == 0:
        return 0
    soc = max(soc, 0)
    first_time, start_soc = soc_over_time_curve[soc_over_time_curve[:, 1] >= soc][0, :]
    second_time = first_time + time_delta
    # Catch out of bounds if time of charging end is bigger than table values

    if second_time >= soc_over_time_curve[-1, 0]:
        end_soc = 1
    else:
        end_soc = soc_over_time_curve[soc_over_time_curve[:, 0] >= second_time][0, 1]

    # Make sure to limit delta soc to 1 if negative socs are given. They are possible during
    # the optimization process but will be continuously raised until they are >0.
    return min(1, end_soc - start_soc)


def get_buffer_time_old(search_time, args):
    for window, buffer_time in args.default_buffer_time_opps.items():
        try:
            start, end = window.split("-")

            if float(end) > search_time.hour >= float(start):
                return timedelta(minutes=buffer_time)
        except ValueError:
            return timedelta(minutes=0)


def electrify_station(stat, stations, electrified_set):
    stations[stat] = {'type': 'opps', 'n_charging_stations': 200}
    electrified_set.add(stat)


def get_buffer_time(trip, args):
    return get_buffer_time_old(trip.arrival_time, args)
    # return timedelta(minutes=get_buffer_time_spice_ev(trip,))


if __name__ == "__main__":
    freeze_support()
    main()
