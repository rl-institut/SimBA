""" Optimizer class which implements the optimizer object and methods needed"""
import json
from copy import copy
from datetime import timedelta
import pickle
from pathlib import Path

import numpy as np
import schedule

from ebus_toolbox import report, rotation
from src import scenario
from optimizer_config import read_config, OptimizerConfig

import logging
import optimizer


class ScheduleOptimizer:
    def __init__(self, sched: schedule.Schedule, scen: scenario.Scenario, args,
                 config: OptimizerConfig(), logger: logging.Logger, **kwargs):
        self.could_not_be_electrified = set()
        self.must_include_set = set()
        self.soc_charge_curve_dict = {}
        self.decision_tree = {}
        self.schedule = sched
        self.scenario = scen
        self.base_schedule = self.schedule
        self.base_scenario = self.scenario
        self.args = args
        self.logger = logger
        self.config = config
        self.electrified_station_set = set()
        with open(self.args.electrified_stations, "r", encoding="utf-8", ) as file:
            self.electrified_stations = json.load(file)
        self.base_stations = self.electrified_stations.copy()
        self.base_electrified_station_set = self.electrified_station_set.copy()

        # stations which are included can not be included again. Therefore they get into
        # the set of not possible stations together with excluded stations
        self.not_possible_stations = config.inclusion_stations.union(config.exclusion_stations)

        if config.reduce_rots:
            self.schedule.rotations = {rot: sched.rotations[rot] for rot in config.rots}

    def loop(self, opt_type="greedy", node_choice="step-by-step", **kwargs):
        # Loops over every base group with the core of group optimization and handles deep or
        # greedy arguments

        # base stations for optimization, so inclusion of stations can be skipped
        base_stations = self.base_stations.copy()
        base_electrified_station_set = self.base_electrified_station_set.copy()

        self.base_scenario = copy(self.scenario)
        self.base_schedule = copy(self.schedule)

        # get events where soc fell below 0. The events contain info about the problematic
        # time span, which includes stations which could provide a soc lift
        base_events = self.get_low_soc_events(rel_soc=False,
                                              soc_data=self.base_scenario.vehicle_socs)

        # check if the events can be divided into subgroups which are independent
        # this makes optimization in smaller groups possible
        groups = optimizer.get_groups_from_events(base_events, self.not_possible_stations)

        # Baseline greedy Optimization
        list_greedy_sets = [set()] * len(groups)

        # get_negative_rotations_all_electrified(base_scen, base_sched, soc_upper_thresh, soc_charge_curve_dict,
        #                                        not_possible_stations=not_possible_stations,
        #                                        rel_soc=False)

        # new_scen.vehicle_socs = timeseries_calc(new_sched.rotations.values(),
        #                                         new_scen.vehicle_socs,
        #                                         new_scen, electrified_station_set,
        #                                         soc_charge_curve_dict=soc_charge_curve_dict,
        #                                         soc_upper_thresh=ARGS.desired_soc_deps)

        # Base line is created simply by not having a decision tree and not a pre optimized_set yet
        for group_nr, group in enumerate(groups[:]):
            events, stations = group
            linien = {lne for e in events for lne in e["rotation"].lines}
            self.electrified_stations = self.base_stations.copy()
            self.electrified_station_set = self.base_electrified_station_set.copy()
            self.logger.warning("Optimizing %s out of %s. This includes these Lines", group_nr + 1,
                                len(groups))
            self.logger.warning(linien)
            self.logger.warning("%s events", (len(events)))

            # first run is always step by step
            self.group_optimization(group, self.choose_station_step_by_step, **kwargs)

            self.logger.warning("Greedy Result ++++++++ %s stations out of %s",
                                len(self.electrified_station_set),
                                len(stations))
            self.logger.warning(self.electrified_station_set)

            list_greedy_sets[group_nr] = self.electrified_station_set.copy()

        # if no deep analysis is needed, loop can continue and part further down is skipped
        if opt_type != "deep":
            for single_set in list_greedy_sets:
                for stat in single_set:
                    self.electrify_station(stat, self.electrified_station_set)
            return self.electrified_stations, self.electrified_station_set

        # from here on only for deep analysis
        for group_nr, group in enumerate(groups[:]):
            events, stations = group
            sols = []
            i = 0
            cont_loop = True
            combinations = optimizer.combs_unordered_no_putting_back(len(stations),
                                                                     len(list_greedy_sets[
                                                                             group_nr]))
            if node_choice == "brute":
                choice_func = self.choose_station_brute
            else:
                choice_func = self.choose_station_step_by_step
            self.logger.debug("There are %s combinations", combinations)
            pre_optimized_set = list_greedy_sets[group_nr]

            while i < self.config.max_brute_loop and cont_loop:
                i += 1
                if i % 10 == 0:
                    print(len(self.decision_tree), " Knotenpunkte durchwandert")
                    print(f"Optimal solution has length {len(self.electrified_station_set)}")
                # not_possible_stations = copy(not_possible_stations)
                self.electrified_stations = self.base_stations.copy()
                new_electrified_set = set()

                new_stations, cont_loop = self.group_optimization(group, choice_func,
                                                                  track_not_possible_rots=False,
                                                                  **kwargs)
                # if a new set was found, print it and save it in sols
                if new_electrified_set != pre_optimized_set and new_stations is not None:
                    self.logger.warning("Optimized with %s  stations %s %s",
                                        len(new_electrified_set), str('#' * 20),
                                        stations_hash(new_electrified_set))
                    sols.append(new_electrified_set)
                    if len(new_electrified_set) < len(pre_optimized_set):
                        electrified_station_set = new_electrified_set

            self.logger.debug("All solutions for this group: %s", sols)
            list_greedy_sets[group_nr] = electrified_station_set.copy()
            self.logger.debug("Optimized with %s stations out of %s", len(electrified_station_set),
                              len(stations))

        # Saving decision tree only in case of deep analysis
        if self.config.save_decision_tree:
            with open(self.args.output_directory / Path("decision_tree.pickle"), "wb") as file:
                pickle.dump(self.decision_tree, file)

        for single_set in list_greedy_sets:
            for stat in single_set:
                self.electrify_station(stat, self.electrified_station_set)

        return self.electrified_stations, self.electrified_station_set

    def group_optimization(self, group, choose_station_function, track_not_possible_rots=True,
                           tree_position=None, **kwargs):
        """ Optimize a single group events
        :return (electrified_stations if optimization continues, None if no further recursive calls
        should happen, bool if further deep analysis should take place)
        :rtype (dict() or None, bool)
        """
        could_not_be_electrified = self.could_not_be_electrified

        # in deep analysis not all nodes are visited if they are suboptimal.
        # therefore could not be electrified is only important in the first run, where every station
        # is possible
        if not track_not_possible_rots:
            could_not_be_electrified = set()

        solver = kwargs.get("solver", "spiceev")

        if not tree_position:
            tree_position = []

        # give tree position
        self.logger.debug("%s with length of %s", tree_position, len(tree_position))

        # Unpack events and possible stations
        event_group, _ = group

        # Get the events which are still negative
        events_remaining = kwargs.get("events_remaining", [99999])

        # Copy Scen and sched
        new_scen = copy(self.base_scenario)
        new_sched = copy(self.base_schedule)

        # Get rotations from event dict and calculate the missing energy
        rotation_dict = {e["rotation"].id: e["rotation"] for e in event_group}
        missing_energy = get_missing_energy(event_group)
        if missing_energy >= 0:
            self.logger.debug("Already electrified: Returning")
            return self.electrified_stations, True

        station_eval = evaluate(event_group, new_scen, soc_curve_dict,
                                kwargs)


        # Get the best stations, brute or step_by_step
        best_station_ids, recursive = choose_station_function(station_eval, electrified_station_set,
                                                              pre_optimized_set, decision_tree,
                                                              missing_energy=missing_energy)
        stat_eval_dict = {stat_id[0]: stat_id[1]["pot_sum"] for stat_id in station_eval}
        if best_station_ids is None:
            LOGGER.warning("No useful station found with %s rotations not electrified yet. "
                           "Stopped after electrifying %s", events_remaining,
                           len(electrified_station_set))

            if pre_optimized_set is None:
                could_not_be_electrified.update(list(rotation_dict.keys()))
                return None, False

            # Remove electrified stations in this run
            copied_set = electrified_station_set.copy()
            for stat in copied_set:
                electrified_stations.pop(stat)
                electrified_station_set.remove(stat)
            # Overwrite with preoptimized set
            for stat in pre_optimized_set:
                self.electrify_station(stat, self.electrified_station_set)
            return None, False

        LOGGER.debug("%s, with first pot of %s", best_station_ids,
                     stat_eval_dict[best_station_ids[0]])

        # Electrify station(s)
        for stat_id in best_station_ids:
            self.electrify_station(stat_id, self.electrified_station_set)

        event_rotations = {x["rotation"] for x in event_group}
        if solver == "quick":
            # Quick calculation has to electrify everything in one step,
            # or the lifting of socs is not correct
            new_scen.vehicle_socs = deepcopy(base_socs)
            new_scen.vehicle_socs = timeseries_calc(event_rotations,
                                                    new_scen.vehicle_socs,
                                                    new_scen, electrified_station_set,
                                                    soc_curve_dict,
                                                    soc_upper_thresh=ARGS.desired_soc_deps,
                                                    electrify_stations=best_station_ids[0])
        else:
            new_sched.rotations = rotation_dict
            new_sched, new_scen = run_schedule(new_sched, ARGS,
                                               electrified_stations=electrified_stations)

        not_possible_stations = set(electrified_stations.keys()).union(not_possible_stations)
        event_rotations_ids = {event["rotation"].id for event in event_group}
        new_events = get_low_soc_events(new_scen,
                                        new_sched, rotations=event_rotations_ids,
                                        not_possible_stations=not_possible_stations, rel_soc=True,
                                        **kwargs)

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
        if not recursive:
            return None, True

        # Check if the events can be divided into subgroups which are independent
        groups = get_groups_from_events(new_events, not_possible_stations, could_not_be_electrified)

        for k, this_group in enumerate(groups):
            this_tree = tree_position.copy()
            this_tree.append(k)
            new_stations, _ = group_optimization(this_group, new_scen, new_sched,
                                                 electrified_stations,
                                                 electrified_station_set,
                                                 could_not_be_electrified,
                                                 not_possible_stations, choose_station_function,
                                                 soc_curve_dict,
                                                 pre_optimized_set, decision_tree,
                                                 tree_position=this_tree,
                                                 **kwargs)

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
                new_scen.vehicle_socs = timeseries_calc(event_rotations,
                                                        new_scen.vehicle_socs,
                                                        new_scen, electrified_station_set,
                                                        soc_curve_dict,
                                                        soc_upper_thresh=ARGS.desired_soc_deps,
                                                        electrify_stations=best_station_ids[0])

                prune_events = get_low_soc_events(new_scen, new_sched,
                                                  rotations=event_rotations_ids,
                                                  not_possible_stations=not_possible_stations,
                                                  rel_soc=True, **kwargs)

                station_eval = evaluate(prune_events, new_scen, soc_curve_dict)
                prune_missing_energy = get_missing_energy(prune_events)
                if not is_branch_promising(station_eval, electrified_station_set,
                                           pre_optimized_set, prune_missing_energy):
                    print("Branch pruned early")
                    is_branch_promising(station_eval, electrified_station_set,
                                        pre_optimized_set, missing_energy)
                    return None, True


    def is_branch_promising(self, station_eval, electrified_station_set,
                            pre_optimized_set, missing_energy):
        """ quickly evaluates if following a branch is promising by summing up estimated potentials
        :return is the branch promising
        :rtype bool"""

        delta = len(pre_optimized_set) - len(electrified_station_set)
        pot = 0
        for i in range(0, min(delta, len(station_eval))):
            pot += station_eval[i][1]["pot_sum"]
        if pot < -missing_energy * self.config.estimation_threshold:
            print(f"Not enough potential {round(pot, 0)} / {round(-missing_energy, 0)} after ",
                  stations_hash(electrified_station_set))
            return False
        return True

    def choose_station_brute(self, station_eval, electrified_station_set,
                             pre_optimized_set=None, decision_tree=None, missing_energy=0,
                             gens=dict()):
        """ Gives back a possible set of Stations to electrify which shows potential and hasnt been
        tried yet. The set of stations is smaller than the best optimized set so far."""

        station_ids = [x[0] for x in station_eval]
        try:
            generator = gens[str(station_ids) + str(len(pre_optimized_set) - 1)]
        except KeyError:
            generator = combination_generator(station_ids, len(pre_optimized_set) - 1)
            gens[str(station_ids) + str(len(pre_optimized_set) - 1)] = generator
        station_eval_dict = {stat[0]: stat[1] for stat in station_eval}
        for comb in generator:
            node_name = stations_hash(comb)
            if node_name not in decision_tree:
                # only check the brute force station if they have the remote chance of fulfilling
                # the missing energy
                # potential>missing energy * 80%
                potential = sum([station_eval_dict[stat]["pot_sum"] for stat in comb])
                if potential > -missing_energy * CONFIG.estimation_threshold:
                    return comb, False
                else:
                    LOGGER.debug("skipped %s since potential is too low %s %%", comb,
                                 round(potential / -missing_energy * 100, 0))
        LOGGER.debug("calculated all viable possibilities")
        return None, False

    def choose_station_step_by_step(self, station_eval,
                                    pre_optimized_set=None, missing_energy=0):
        """ Gives back a station of possible stations to electrify which shows the biggest potential
        and hasnt been picked before."""

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
                if decision_tree[stations_hash(check_stations)][
                    "visit_counter"] == min_count_visited:
                    best_station_id = station[0]
                    return [best_station_id], True
        return None, True

    def set_battery_and_charging_curves(self):
        for v_type in self.schedule.vehicle_types.values():
            for vehicle in v_type.values():
                if self.config.battery_capacity is not None:
                    vehicle["capacity"] = self.config.battery_capacity
                if self.config.charging_curve is not None:
                    vehicle["charging_curve"] = self.config.charging_curve

    def set_up_decision_tree(self):
        if self.config.decision_tree_path is not None:
            with open(self.config.decision_tree_path, "rb") as file:
                self.decision_tree = self.config.load(file)

    def rebase_spice_ev(self):
        self.logger.debug("Spice EV Rebasing Scenario")
        must_include_set, ele_stations = self.preprocessing_scenario(
            run_only_neg=self.config.run_only_neg)
        self.logger.debug("Rebasing finished")
        if self.config.pickle_rebased:
            toolbox_to_pickle(self.config.pickle_rebased_name, self.schedule, self.scenario,
                              self.args)
            self.logger.debug("Rebased scenario pickled  as %s", self.config.pickle_rebased_name)

        self.must_include_set = must_include_set
        self.base_stations = self.electrified_stations.copy()
        return must_include_set, ele_stations

    def rebase_simple(self):
        must_include_set = set()
        # Electrify inclusion stations
        for stat in self.config.inclusion_stations:
            self.electrify_station(stat, must_include_set)
        self.schedule.rotations = {r: self.schedule.rotations[r] for r in self.schedule.rotations if
                                   r not in self.config.exclusion_rots}
        self.must_include_set = must_include_set
        self.base_stations = self.electrified_stations.copy()
        return must_include_set, self.electrified_stations

    def preprocessing_scenario(self, electrified_stations=None, run_only_neg=False,
                               cost_calc=False):
        """ Prepare scenario and run schedule
        :return schedule, scenario, electrified_station_set, electrified_stations
        :rtype (schedule.Schedule(), scenario.Scneario(), set(), dict())
        """
        must_include_set = set()
        # Electrify inclusion stations
        for stat in self.config.inclusion_stations:
            self.electrify_station(stat, must_include_set)

        # Calc new but only prev. negative rotations
        if run_only_neg:
            rots = self.schedule.get_negative_rotations(self.scenario)
            rots = {r: self.schedule.rotations[r] for r in rots if r not in
                    self.config.exclusion_rots}
            self.schedule.rotations = rots
        else:
            rots = {r: self.schedule.rotations[r] for r in self.schedule.rotations if
                    r not in self.config.exclusion_rots}
            self.schedule.rotations = rots

        new_sched, new_scen = optimizer.run_schedule(self.schedule, self.args,
                                                     self.electrified_stations,
                                                     cost_calc=cost_calc)
        report.generate(new_sched, new_scen, self.args)
        self.schedule = new_sched
        self.scenario = new_scen
        self.must_include_set = must_include_set
        return must_include_set, electrified_stations

    def electrify_station(self, stat, electrified_set):
        """electrify a station and keep track of it in the electrified set file"""
        self.electrified_stations[stat] = self.config.standard_opp_station
        electrified_set.add(stat)

    def create_charging_curves(self):
        """Cycle through vehicles and create numerically created charging curves with energy supplied
        over time"""
        soc_charge_curve_dict = {}
        for v_type_name in self.schedule.vehicle_types:
            soc_charge_curve_dict[v_type_name] = {}
        for name, v_type in self.schedule.vehicle_types.items():
            for ch_type, data in v_type.items():
                soc_charge_curve_dict[name][ch_type] = charging_curve_to_soc_over_time(
                    data["charging_curve"], data["capacity"], self.args,
                    self.schedule.cs_power_opps, efficiency=self.config.charge_eff,
                    time_step=0.1)
        self.soc_charge_curve_dict = soc_charge_curve_dict

    def get_must_stations_and_rebase(self, relative_soc=False):
        """Electrify everything minus 1 station. If without the stations there are below zero events
         it is a must have station"""

        events = self.get_low_soc_events(rel_soc=relative_soc)

        stats = {station for event in events for station in event["stations_list"]
                 if station not in self.not_possible_stations}

        electrified_station_set_all = set(stats)
        must_stations = set()
        for station in electrified_station_set_all:
            self.logger.debug(".", end="")
            electrified_station_set = electrified_station_set_all.difference([station])
            vehicle_socs = optimizer. \
                timeseries_calc(self.schedule.rotations.values(),
                                self.scenario.vehicle_socs,
                                self.scenario, electrified_station_set,
                                soc_charge_curve_dict=self.soc_charge_curve_dict,
                                soc_upper_thresh=self.args.desired_soc_deps)

            for rot in self.schedule.rotations:
                soc, start, end = self.get_rotation_soc(rot, vehicle_socs)
                if np.min(soc[start:end]) < self.config.min_soc:
                    self.logger.debug("%s , with min soc: %s", station, np.min(soc[start:end]))
                    must_stations.add(station)
                    break
        self.not_possible_stations = self.not_possible_stations.union(must_stations)

        # rebasing
        for stat in must_stations:
            # dont put must stations in electrified set, but in extra set must_include_set
            self.electrify_station(stat, self.must_include_set)
        self.scenario.vehicle_socs = optimizer. \
            timeseries_calc(self.schedule.rotations.values(),
                            self.scenario.vehicle_socs,
                            self.scenario, must_stations,
                            soc_charge_curve_dict=self.soc_charge_curve_dict,
                            soc_upper_thresh=self.args.desired_soc_deps)
        return must_stations

    def remove_none_socs(self):
        """ Removes soc values of None by filling them with the last value which is not None
        :param this_scen: scenario from which the None socs should be removed
        :type this_scen: scenario.Scenario()
        """
        # make sure no None values exists in SOCs. Fill later values with last value
        # which was not None, eg, [1,5,4,None,None] becomes [1,5,4,4,4]
        for v_id, soc in self.scenario.vehicle_socs.items():
            soc = np.array(soc)
            last_not_none = "not found"
            if None in soc:
                for ii in range(len(soc) - 1, -1, -1):
                    if soc[ii] is not None:
                        last_not_none = soc[ii]
                        break
                soc[soc == np.array(None)] = last_not_none
            self.scenario.vehicle_socs[v_id] = soc


    def get_rotation_soc(self, rot_id, soc_data: dict = None):
        rot = self.schedule.rotations[rot_id]
        rot_start_idx = self.get_index_by_time(rot.departure_time)
        rot_end_idx = self.get_index_by_time(rot.arrival_time)
        if soc_data:
            return soc_data[rot.vehicle_id], rot_start_idx, rot_end_idx
        return self.scenario.vehicle_socs[rot.vehicle_id], rot_start_idx, rot_end_idx

    def get_index_by_time(self, search_time):
        start_time = self.scenario.start_time
        delta_time = timedelta(minutes=60 / self.scenarioc.stepsPerHour)
        idx = (search_time - start_time) // delta_time
        return idx

    def get_trips(self, rot: rotation.Rotation, start_idx: int, end_idx: int):
        # return trips in a rotation from a start to an end index, if the arrival time is in between
        # the start and end idx
        start_time_event = self.get_time_by_index(start_idx, self.scenario)
        end_time_event = self.get_time_by_index(end_idx, self.scenario)

        trips = []
        for i, trip in enumerate(rot.trips):
            if end_time_event > trip.arrival_time > start_time_event:
                trips.append(trip)

        return trips

    def get_time_by_index(self, idx):
        start_time = self.scenario.start_time
        delta_time = timedelta(minutes=60 / self.scenario.stepsPerHour)
        searched_time = start_time + delta_time * idx
        return searched_time








def toolbox_from_pickle(sched_name, scen_name, args_name):
    """ Load the 3 files from pickle"""
    with open(args_name, "rb") as f:
        this_args = pickle.load(f)
    with open(scen_name, "rb") as f:
        scen = pickle.load(f)
    with open(sched_name, "rb") as f:
        sched = pickle.load(f)
    return sched, scen, this_args


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


def charging_curve_to_soc_over_time(charging_curve, capacity, args,
                                    max_charge_from_grid=float('inf'),
                                    time_step=0.1, efficiency=1):
    """create charging curve as nested list of SOC, Power[kW] and capacity in [kWh]"""
    # simple numeric creation of power over time --> to energy over time
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
        soc2 = soc + time_step / 60 * power1
        power2 = min(np.interp(soc2, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        power = (power1 + power2) / 2 * efficiency
        soc += time_step / 60 * power
        time += time_step
    # Fill the soc completely in last timestep
    times.append(time)
    socs.append(args.desired_soc_opps)
    return np.array((times, socs)).T


def get_missing_energy(events):
    """ Sum up all the missing energies of the given events"""
    missing_energy = 0
    for event in events:
        missing_energy += event["min_soc"] * event["capacity"]
    return missing_energy


def stations_hash(stations_set):
    """ Create a simple str as hash for a set of stations"""
    return str(sorted(list(stations_set)))
