""" Optimizer class which implements the optimizer object and methods needed"""

from copy import deepcopy, copy
from datetime import datetime, timedelta
import logging
import numpy as np
from pathlib import Path
import pickle
from typing import Iterable

import simba.optimizer_util as opt_util
from spice_ev import scenario
from simba import rotation, schedule, util


class StationOptimizer:
    """Class for station optimization"""

    def __init__(self, sched: schedule.Schedule, scen: scenario.Scenario, args,
                 config: 'opt_util.OptimizerConfig', logger: logging.Logger):
        self.base_not_possible_stations = set()
        self.could_not_be_electrified = set()
        self.must_include_set = set()
        self.soc_charge_curve_dict = {}
        self.decision_trees = []
        self.current_tree = {}
        self.schedule = sched
        self.scenario = scen
        self.base_schedule = self.schedule
        self.base_scenario = self.scenario
        self.args = args
        self.logger = logger
        self.config = config
        self.electrified_station_set = set()
        self.electrified_stations = deepcopy(self.schedule.stations)
        self.base_stations = self.electrified_stations.copy()
        self.base_electrified_station_set = set()

        # Generator for producing combinations for the brute force algorithm
        self.brute_generator = dict()

        # stations which are included can not be included again. Therefore they go into
        # the set of not possible stations together with excluded stations
        self.not_possible_stations = config.inclusion_stations.union(config.exclusion_stations)

        # if the config defines special rotations, only use them in the schedule
        if config.reduce_rotations:
            self.schedule.rotations = {rot: sched.rotations[rot] for rot in config.rotations}

    def loop(self, **kwargs):
        """ Loops recursively over every base group and handles deep or greedy arguments.

        :param kwargs: arguments handled by the function and config.
        :type kwargs: dict
        :return: tuple of the dict electrified_stations and the set electrified_station_set
            in the optimal found case
        :rtype: (dict, set)
        """

        node_choice = kwargs.get("node_choice", self.config.node_choice)
        opt_type = kwargs.get("opt_type", self.config.opt_type)
        electrified_station = self.must_include_set.union(self.electrified_station_set)
        self.scenario.vehicle_socs = self.timeseries_calc(electrified_station)
        self.base_scenario = deepcopy(self.scenario)
        self.base_schedule = copy(self.schedule)

        self.base_stations = self.electrified_stations.copy()
        self.base_not_possible_stations = self.not_possible_stations.copy()

        # get events where SoC fell below the minimal
        # SoC. The events contain info about the problematic
        # time span, which includes stations which could provide a SoC lift
        base_events = self.get_low_soc_events(
            rel_soc=False, soc_data=self.base_scenario.vehicle_socs)

        # check if the events can be divided into subgroups which are independent
        # this makes optimization in smaller groups possible
        groups = opt_util.get_groups_from_events(
            base_events, self.not_possible_stations,
            could_not_be_electrified=self.could_not_be_electrified, optimizer=self)

        if len(groups) == 0:
            self.logger.info("The scenario is already optimized, as it has no low SoC.")
            return self.electrified_stations, self.electrified_station_set

        # sort groups by highest potential of a single electrification. Used if partial
        # electrification is relevant
        groups = sorted(groups, key=lambda group: opt_util.evaluate(
            group[0], self, soc_data=self.scenario.vehicle_socs)[0][1])

        # storage for the sets of stations which will be generated
        list_greedy_sets = [set()] * len(groups)

        # create a decision tree or load one from a previous run.
        # there are as many trees as groups.
        self.set_up_decision_tree(len(groups))

        # baseline greedy optimization
        # base line is created simply by not having a decision tree and not a pre optimized_set yet
        for group_nr, group in enumerate(groups[:]):
            # unpack the group
            events, stations = group
            # define the current tree from the decision trees
            self.current_tree = self.decision_trees[group_nr]
            lines = {lne for e in events for lne in e.rotation.lines}
            self.logger.warning(
                "Optimizing %s out of %s groups. This includes these Lines: %s",
                group_nr + 1, len(groups), lines)
            self.logger.warning("%s low soc events with %s potential stations",
                                (len(events)), len(stations))

            # the electrified_station_set get mutated by the group optimization
            # the results are stored in the list of greedy sets, but before each group optimization
            # the base stations are restored, so the resulting sets do not stack on top of each
            # other
            self.electrified_station_set = self.base_electrified_station_set.copy()

            if self.config.solver == "spiceev":
                # running SpiceEV changes the scenario and vehicle ids. Since the low SoC events
                # point to the base vehicle ids, these have to be used in the group optimization
                self.scenario = deepcopy(self.base_scenario)

            # first run is always step by step
            self.group_optimization(
                group, self.choose_station_step_by_step, events_remaining=[len(events)], **kwargs)

            self.logger.warning("Greedy Result ++++++++ %s stations out of %s",
                                len(self.electrified_station_set), len(stations))
            self.logger.warning(self.electrified_station_set)

            list_greedy_sets[group_nr] = self.electrified_station_set.copy()

        # if no deep analysis is needed, return here and part further down is skipped
        # if deep analysis is needed run can continue
        if opt_type != "deep":
            for single_set in list_greedy_sets:
                for stat in single_set:
                    self.electrify_station(stat, self.electrified_station_set)
            return self.electrified_stations, self.electrified_station_set
        self.logger.warning("Starting deep analysis with mode: %s", self.config.node_choice)
        self.logger.debug(opt_util.get_time())
        # from here on only for deep analysis
        for group_nr, group in enumerate(groups[:]):
            events, stations = group
            solutions = []
            self.current_tree = self.decision_trees[group_nr]
            if node_choice == "brute":
                choice_func = self.choose_station_brute
            else:
                choice_func = self.choose_station_step_by_step

            combinations = opt_util.combs_unordered_no_putting_back(
                len(stations), len(list_greedy_sets[group_nr]) - 1)
            self.logger.debug(
                "There are %s combinations with 1 station less than the current solution with %s "
                "stations out of %s", combinations, len(list_greedy_sets[group_nr]), len(stations))

            pre_optimized_set = list_greedy_sets[group_nr]
            for i in range(self.config.max_brute_loop):
                if i % 10 == 0:
                    self.logger.log(msg=f"{len(self.current_tree)} nodes checked", level=100)
                    self.logger.log(msg=f"{i}/{self.config.max_brute_loop} simulations for group "
                                        f"{group_nr}/{len(groups)}", level=100)
                    self.logger.log(
                        msg=f"Optimal solution has length {len(pre_optimized_set)}", level=100)

                self.electrified_stations = self.base_stations.copy()
                self.electrified_station_set = self.base_electrified_station_set.copy()
                self.scenario = copy(self.base_scenario)
                self.schedule = copy(self.base_schedule)

                self.scenario.vehicle_socs = deepcopy(self.base_scenario.vehicle_socs)
                # deep copy of schedule.rotations is very slow. Not needed for quick calculation.
                if self.config.solver == "spiceev":
                    self.schedule.rotations = deepcopy(self.base_schedule.rotations)
                self.not_possible_stations = self.base_not_possible_stations.copy()

                # create a new electrified set
                try:
                    new_stats = self.group_optimization(
                        group, choice_func, track_not_possible_rots=False,
                        pre_optimized_set=pre_optimized_set, events_remaining=[len(events)],
                        **kwargs)

                    # The terminal node is checked therefore not viable anymore
                    node_name = opt_util.stations_hash(self.electrified_station_set)
                    self.current_tree[node_name]["viable"] = False
                    self.current_tree[node_name]["completed"] = True
                except opt_util.SuboptimalSimulationException:
                    # exception if the optimization encountered a dead end, where it is clear
                    # that this case will not lead to an optimal solution
                    continue
                except opt_util.AllCombinationsCheckedException:
                    # exception if all combinations in the brute force method are checked or at
                    # least estimated or all nodes show no potential in the step-by-step method
                    self.logger.warning("The brute force method checked all possible combinations")
                    break

                new_electrified_set = self.electrified_station_set
                # if a new set was found, print it and save it in solutions
                if new_electrified_set != pre_optimized_set and new_stats is not None:
                    self.logger.warning(
                        "Optimized with %s  stations %s %s", len(new_electrified_set),
                        str('#' * 20), opt_util.stations_hash(new_electrified_set))
                    solutions.append(new_electrified_set)
                    if len(new_electrified_set) < len(pre_optimized_set):
                        list_greedy_sets[group_nr] = new_electrified_set.copy()
                        pre_optimized_set = new_electrified_set.copy()
            else:
                self.logger.warning(f"Ran all {self.config.max_brute_loop} loops")

            # print all solutions
            self.logger.debug("All solutions for this group: %s", solutions)
            # how many stations do these solutions need in comparison to the available stations?
            self.logger.debug(
                "Optimized with %s stations out of %s", len(pre_optimized_set), len(stations))

        # make sure it gets printed in the console through high level
        self.logger.log(msg=opt_util.get_time(), level=100)
        # saving decision tree only in case of deep analysis
        if self.config.save_decision_tree:
            with open(self.config.optimizer_output_dir / Path("decision_tree.pickle"), "wb") as f:
                pickle.dump(self.current_tree, f)

        self.electrified_stations = self.base_stations.copy()
        self.electrified_station_set = self.base_electrified_station_set.copy()
        for single_set in list_greedy_sets:
            for stat in single_set:
                self.electrify_station(stat, self.electrified_station_set)
        # dump the measured running times of the functions
        self.logger.debug("Durations of Optimization: ", opt_util.time_it(None))
        return self.electrified_stations, self.electrified_station_set

    def get_negative_rotations_all_electrified(self, rel_soc=False):
        """ Get the ids for the rotations which show negative SoCs when everything is electrified.

        :param rel_soc: if true, the start SoC is handled like it has the desired deps SoC
        :type rel_soc: bool
        :return: rotation ids which are negative even with all stations electrified
        :rtype: set
        """
        events = self.get_low_soc_events(rel_soc=rel_soc)

        stats = {stat for event in events for stat in event.stations_list
                 if stat not in self.not_possible_stations}
        electrified_station_set = set(stats)
        vehicle_socs = self.timeseries_calc(electrified_station_set)
        new_events = self.get_low_soc_events(soc_data=vehicle_socs)
        return {event.rotation.id for event in new_events}

    @opt_util.time_it
    def group_optimization(self, group, choose_station_function, track_not_possible_rots=True,
                           pre_optimized_set=None, **kwargs):
        """ Optimize a single group events and returns the electrified stations and a flag.

        :param group: tuple of events and station_ids to be optimized together
        :type group: (list(simba.optimizer_util.LowSocEvent), list(str))
        :param choose_station_function: function to be used for choosing the next function
        :type choose_station_function: function
        :param track_not_possible_rots: not possible tracks need to be tracked in the first run
        :type track_not_possible_rots: bool
        :param pre_optimized_set: stations which were optimized before hand used
            as upper threshold for optimization
        :type pre_optimized_set: set
        :param kwargs: internally used functionality
        :type kwargs: dict
        :return: electrified_stations
        :rtype dict
        :raises AllCombinationsCheckedException: If all combinations have been checked
        :raises SuboptimalSimulationException: if a suboptimal Simulation has been identified
        """
        could_not_be_electrified = self.could_not_be_electrified

        # could_not_be_electrified only important in first loop, where every station is possible"
        if not track_not_possible_rots:
            # could_not_be_electrified gets decoupled with optimizer attribute
            could_not_be_electrified = set()

        solver = self.config.solver

        tree_position = kwargs.get("tree_position", [1])

        # give tree position
        self.logger.debug(
            "%s with length of %s and a total of %s stations", tree_position, len(tree_position),
            len(self.electrified_station_set))

        # unpack events and possible stations
        event_group, _ = group

        # get the events which are still negative
        events_remaining = kwargs.get("events_remaining", [99999])

        kwargs["lifted_socs"] = kwargs.get("lifted_socs", self.scenario.vehicle_socs)

        # evaluate with current scenario and schedule and lifted values from the
        # call above. Important not to use the scenario in optimizer, since it gets overwritten
        # every call and only relevant trips are adjusted. This leads to problems
        # when the recursive call goes on a higher level eg. level 0 - level 1 - level 2 - level 1
        station_eval = opt_util.evaluate(event_group, self, soc_data=kwargs["lifted_socs"])

        missing_energy = opt_util.get_missing_energy(event_group, min_soc=self.config.min_soc)
        if missing_energy >= 0:
            self.logger.debug("Already electrified: Returning")
            return self.electrified_stations

        # get rotations from event dict and calculate the missing energy
        rotation_dict = {e.rotation.id: e.rotation for e in event_group}

        self.expand_tree(station_eval)

        # check if the children are viable or not. If no children are viable this node is not
        # viable itself for further inspection
        # if the node is the root of the tree, i.e. no electrified stations yet, than all
        # combinations have been checked
        if not self.is_node_viable():
            node_name = opt_util.stations_hash(self.electrified_station_set)
            self.current_tree[node_name]["viable"] = False
            if len(self.electrified_station_set) == 0:
                raise opt_util.AllCombinationsCheckedException
            raise opt_util.SuboptimalSimulationException

        # get the best stations, brute or step_by_step
        best_station_ids, recursive = choose_station_function(
            station_eval, pre_optimized_set=pre_optimized_set, missing_energy=missing_energy)

        stat_eval_dict = {stat_id[0]: stat_id[1] for stat_id in station_eval}

        self.logger.debug("%s, with first pot of %s", best_station_ids,
                          round(stat_eval_dict[best_station_ids[0]], 1))

        # electrify station(s)
        for stat_id in best_station_ids:
            self.electrify_station(stat_id, self.electrified_station_set)

        event_rotations = {event.rotation for event in event_group}

        # copy scen and sched from base scenario and base schedule. This way we can be sure
        # to have the correct original basis for simulation. If we would not copy, there can be
        # unintended, changes in scenario start times for example, since some values do no stay the
        # same if we mutate the rotation dictionary and run the scenario (e.g changing vehicle ids)
        # This led to problems in the past
        self.copy_scen_sched()

        if solver == "quick":
            # quick calculation has to electrify everything in one step, that is chronologically.
            # or the lifting of SoCs is not correct.
            # Since its only adding charge on top of SoC time series, electrification that took
            # place before in the optimization but later in the time series in regards to the
            # current electrification would show to much charge, since charging curves decrease over
            # SoC.
            self.deepcopy_socs()
            electrified_stations = self.electrified_station_set.union(best_station_ids)
            self.scenario.vehicle_socs = self.timeseries_calc(electrified_stations, event_rotations)
        else:
            self.schedule.rotations = rotation_dict
            self.schedule, self.scenario = opt_util.run_schedule(
                self.schedule, self.args,  electrified_stations=self.electrified_stations)

        kwargs["lifted_socs"] = self.scenario.vehicle_socs.copy()

        self.not_possible_stations = set(self.electrified_stations.keys()).union(
            self.not_possible_stations)
        event_rotations_ids = {event.rotation.id for event in event_group}
        new_events = self.get_low_soc_events(rotations=event_rotations_ids,
                                             rel_soc=True, **kwargs)

        delta_energy = opt_util.get_missing_energy(new_events, self.config.min_soc)
        events_remaining[0] -= len(event_group) - len(new_events)

        new_rotation = {event.rotation for event in new_events}
        r_electrified = len(event_rotations_ids) - len(new_rotation)
        self.logger.debug(
            f"Last electrification electrified {r_electrified}/{len(event_rotations_ids)} "
            f"Rotations. {events_remaining[0]} remaining events in the base group.")
        delta_base_energy = delta_energy

        # put this node into the decision tree including the missing energy
        self.node_to_tree(delta_base_energy)

        # everything electrified
        if delta_energy >= 0:
            return self.electrified_stations
        # some choice functions might not need a recursive call, they return here. recursive is set
        # by the choose_station_function
        if not recursive:
            raise opt_util.SuboptimalSimulationException

        # check if the events can be divided into subgroups which are independent
        groups = opt_util.get_groups_from_events(
            new_events, self.not_possible_stations, could_not_be_electrified, optimizer=self)

        for k, this_group in enumerate(groups):
            this_tree = tree_position + [k + 1]
            kwargs["tree_position"] = this_tree
            new_stations = self.group_optimization(
                this_group, choose_station_function,
                track_not_possible_rots=track_not_possible_rots,
                pre_optimized_set=pre_optimized_set, **kwargs)

            if new_stations is not None:
                self.electrified_stations.update(new_stations)

            # if there is no pre optimized set function can return
            if pre_optimized_set is None:
                continue
            # if there is a pre optimized set check if the pre optimized set leads to pruning
            # i.e. if following branch makes sense. Since the prognosis is not very accurate
            # over electrifying many stations, pruning is very unlikely early on. The pruning
            # threshold should be set not to high, since checking for pruning is likely to be an
            # unnecessary calculation
            thresh = self.config.pruning_threshold
            if len(pre_optimized_set) - len(self.electrified_station_set) < thresh:
                self.copy_scen_sched()
                self.deepcopy_socs()
                electrified_stations = self.electrified_station_set.union(best_station_ids)
                self.scenario.vehicle_socs = self.timeseries_calc(
                    electrified_stations, event_rotations)
                prune_events = self.get_low_soc_events(
                    rotations=event_rotations_ids, rel_soc=True, **kwargs)
                station_eval = opt_util.evaluate(prune_events, self)
                prune_missing_energy = \
                    opt_util.get_missing_energy(prune_events, self.config.min_soc)
                if not self.is_branch_promising(station_eval, self.electrified_station_set,
                                                pre_optimized_set, prune_missing_energy):
                    self.logger.debug("Branch pruned early")
                    node_name = opt_util.stations_hash(self.electrified_station_set)
                    self.current_tree[node_name]["viable"] = False
                    raise opt_util.SuboptimalSimulationException

        return self.electrified_stations

    @opt_util.time_it
    def deepcopy_socs(self):
        """ Deepcopy of the SoCs in the scenario."""
        self.scenario.vehicle_socs = deepcopy(self.base_scenario.vehicle_socs)

    @opt_util.time_it
    def copy_scen_sched(self):
        """Copy of the base scenario and base schedule."""
        self.scenario = copy(self.base_scenario)
        self.schedule = copy(self.base_schedule)

    def sort_station_events(self, charge_events_single_station):
        """ Sort events by arrival times.

        :param charge_events_single_station: charge events
        :type charge_events_single_station: list
        :return: sorted charging events
        :rtype: list
        """
        return sorted(charge_events_single_station, key=lambda x: x.arrival_time)

    @opt_util.time_it
    def timeseries_calc(self, electrified_stations: set, rotations=None) -> dict:
        """ A quick estimation of SoCs.

        Iterates through rotations and calculates the SoC.
        The start value is assumed to be desired_soc_deps.
        The function subtracts the consumption of each trip
        and numerically estimates the charge at the electrified station.
        Charging powers depend on the soc_charge_curve_dict, which were earlier created using the
        vehicle charge curve and the schedule.cs_power_opps.

        :param electrified_stations: stations which are calculated as electrified with cs_power_opps
        :type electrified_stations: set
        :param rotations: Rotations to be calculated. Defaults to optimizer.schedule.rotations
        :type rotations: iterable[Rotation]
        :param electrified_stations: Stations which are electrified. Default None leads to using the
            so far optimized optimizer.electrified_station_set
        :return: Returns SoC dict with lifted SoCs
        :rtype dict()
        """
        if rotations is None:
            rotations = self.schedule.rotations.values()

        vehicle_socs = deepcopy(self.scenario.vehicle_socs)
        for rot in rotations:
            ch_type = (rot.vehicle_id.find("oppb") > 0) * "oppb" + (
                    rot.vehicle_id.find("depb") > 0) * "depb"
            v_type = rot.vehicle_id.split("_" + ch_type)[0]
            capacity = self.schedule.vehicle_types[v_type][ch_type]["capacity"]
            soc_over_time_curve = self.soc_charge_curve_dict[v_type][ch_type]
            soc = vehicle_socs[rot.vehicle_id]
            last_soc = self.args.desired_soc_deps
            for i, trip in enumerate(rot.trips):
                # Handle consumption during trip
                # Find start, end and delta_soc of the current trip
                idx_start = opt_util.get_index_by_time(self.scenario, trip.departure_time)
                idx_end = opt_util.get_index_by_time(self.scenario, trip.arrival_time)
                delta_idx = idx_end + 1 - idx_start
                d_soc = trip.consumption / capacity

                # Linear interpolation of SoC during trip
                soc[idx_start:idx_end + 1] = np.linspace(last_soc, last_soc - d_soc, delta_idx)

                # Update last known SoC with current value
                last_soc = last_soc - d_soc

                # Fill the values while the vehicle is standing waiting for the next trip
                idx = opt_util.get_index_by_time(self.scenario, trip.arrival_time)
                try:
                    idx_end = opt_util.get_index_by_time(
                        self.scenario, rot.trips[i + 1].departure_time)
                except IndexError:
                    # No next trip. Rotation finished.
                    break
                # Set SoC at arrival time for whole standing time
                soc[idx:idx_end+1] = last_soc

                if trip.arrival_name not in electrified_stations:
                    # Arrival station is not electrified: skip station
                    continue

                # Arrival station is not the last station of the rotation and electrified.
                # Calculate the charge at this station

                # Get the standing time in minutes. Buffer time is already subtracted
                standing_time_min = opt_util.get_charging_time(
                    trip, rot.trips[i + 1], self.args)

                # Assume that the start SoC for charging is at least 0, since this results in the
                # largest possible charge for a monotonous decreasing charge curve in the SoC range
                # of [0%,100%].
                search_soc = max(0, soc[idx])

                # Get the SoC lift
                d_soc = opt_util.get_delta_soc(soc_over_time_curve, search_soc, standing_time_min)

                # Clip the SoC lift to the desired_soc_deps,
                # which is the maximum that can reached when the rotation stays positive
                d_soc = min(d_soc, self.args.desired_soc_opps)

                # Add the charge as linear interpolation during the charge time, but only start
                # after the buffer time
                buffer_minutes = util.get_buffer_time(trip, self.args.default_buffer_time_opps)
                buffer_idx = int(timedelta(minutes=buffer_minutes) / self.scenario.interval)
                delta_idx = int(timedelta(minutes=standing_time_min) / self.scenario.interval + 1)
                soc[idx + buffer_idx:idx + buffer_idx + delta_idx] += np.linspace(0, d_soc,
                                                                                  delta_idx)
                # Keep track of the last SoC as starting point for the next trip
                last_soc = soc[idx + buffer_idx + delta_idx-1]

            start_idx = opt_util.get_index_by_time(self.scenario, rot.trips[0].departure_time)
            end_idx = opt_util.get_index_by_time(self.scenario, rot.trips[-1].arrival_time)
            vehicle_socs[rot.vehicle_id][start_idx:end_idx+1] = soc[start_idx:end_idx+1]
        return vehicle_socs

    @opt_util.time_it
    def expand_tree(self, station_eval):
        """ Check which nodes can be reached from the current node.

        :param station_eval: station evaluation with station names and potential
        :type station_eval: list(str, float)
        """
        try:
            parent_name = opt_util.stations_hash(self.electrified_station_set)
            self.current_tree[parent_name]
        except KeyError:
            self.current_tree[parent_name] = get_init_node()
        for stat, ev_score in station_eval:
            node_name = opt_util.stations_hash(self.electrified_station_set.union([stat]))
            self.current_tree[parent_name]["children"].add(node_name)
            try:
                self.current_tree[node_name]
            except KeyError:
                self.current_tree[node_name] = get_init_node()

    @opt_util.time_it
    def is_node_viable(self):
        """ Check if a node has viable children.

        Viable children are not terminal and did not show
        lack of potential yet

        :return: True if the node still has potential
        """
        parent_name = opt_util.stations_hash(self.electrified_station_set)
        for child in self.current_tree[parent_name]["children"]:
            if self.current_tree[child]["viable"]:
                # at least one viable node is found --> return True
                return True
        # not a single viable node was found --> return False
        return False

    @opt_util.time_it
    def is_branch_promising(self, station_eval, electrified_station_set,
                            pre_optimized_set, missing_energy):
        """ Return if following a branch is promising by summing up estimated potentials.

        :param station_eval: sorted station evaluation with name and potential of station
        :type station_eval: list(str, float)
        :param electrified_station_set: electrified stations
        :type electrified_station_set: set(str)
        :param pre_optimized_set: already optimized stations
        :type pre_optimized_set: set(str)
        :param missing_energy: the amount of energy missing in this branch
        :type missing_energy: float
        :return: is the branch promising
        :rtype bool
        """

        delta = len(pre_optimized_set) - len(electrified_station_set)
        if delta < 0:
            return False
        pot = 0
        # sum up highest potentials until you would reach the number of electrified stations
        # in the pre optimized case
        for i in range(0, min(delta, len(station_eval))):
            pot += station_eval[i][1]
        if pot < -missing_energy * self.config.estimation_threshold:
            self.logger.debug(
                "Not enough potential: %s  after: %s ", round(pot / -missing_energy, 0),
                opt_util.stations_hash(electrified_station_set))
            return False
        return True

    def node_to_tree(self, delta_base_energy):
        """ Fill decision tree with the given info.

        :param delta_base_energy: missing energy
        :type delta_base_energy: float
        """
        node_name = opt_util.stations_hash(self.electrified_station_set)
        self.current_tree[node_name]["missing_energy"] = delta_base_energy
        self.current_tree[node_name]["visit_counter"] += 1
        if self.current_tree[node_name]["visit_counter"] > 1:
            self.logger.debug("already visited this node")

    @opt_util.time_it
    def choose_station_brute(self, station_eval,
                             pre_optimized_set=None, missing_energy=0):
        """ Return a possible set of stations to electrify which has not been tried yet.

        Gives back a possible set of stations to electrify which shows potential and has not been
        tried yet. The set of stations is smaller than the best optimized set so far.

        :param station_eval: sorted station evaluation with names and potential
        :type station_eval list(str, float)
        :param pre_optimized_set: optimized stations thus far
        :type pre_optimized_set: set
        :param missing_energy: missing energy in this branch before electrification
        :type missing_energy: float
        :return: combination of stations to electrify and
            false since this function does not support recursive calling
        :raises AllCombinationsCheckedException: If all combinations have been checked
        """

        station_ids = [x[0] for x in station_eval]
        # get combination generator (reuse if possible)
        generator_key = f"{station_ids} {len(pre_optimized_set) - 1}"
        try:
            generator = self.brute_generator[generator_key]
        except KeyError:
            # create new generator
            generator = opt_util.combination_generator(station_ids, len(pre_optimized_set) - 1)
            self.brute_generator[generator_key] = generator
        station_eval_dict = {stat[0]: stat[1] for stat in station_eval}
        for comb in generator:
            node_name = opt_util.stations_hash(comb)
            if node_name not in self.current_tree:
                # only check the combination of stations if they have the remote chance of
                # fulfilling the missing energy. The default for the estimation_threshold is 80%
                potential = sum([station_eval_dict[stat] for stat in comb])
                if potential > -missing_energy * self.config.estimation_threshold:
                    self.current_tree[node_name] = get_init_node()
                    return comb, False
                else:
                    self.logger.debug(
                        "skipped %s since potential is too low %s %%", comb,
                        round(potential / -missing_energy * 100, 0))
                    try:
                        self.current_tree[node_name]["viable"] = False
                    except KeyError:
                        self.current_tree[node_name] = get_init_node()
                        self.current_tree[node_name]["viable"] = False
        self.logger.debug("Calculated all viable possibilities")
        raise opt_util.AllCombinationsCheckedException

    @opt_util.time_it
    def choose_station_step_by_step(self, station_eval,
                                    pre_optimized_set=None, missing_energy=0):
        """ Returns a station of possible stations to electrify which shows the biggest potential.

        It checks that this station has not been picked before and the path of electrification is
        viable. If all stations have been picked, the minimum amount of picks is determined. From
        the group of stations with this number of picks, the one with the highest amount of
        potential is chosen.

        :param station_eval: sorted station evaluation with names and potential
        :type station_eval list(str, float)
        :param pre_optimized_set: optimized stations thus far
        :type pre_optimized_set: set
        :param missing_energy: missing energy in this branch before electrification
        :type missing_energy: float
        :return: a station to electrify and True since this function needs recursive calling
        :rtype: str, bool
        :raises SuboptimalSimulationException: if a suboptimal Simulation has been identified
        """

        # filter function to stop simulating cases which have no hope of being optimal.
        # potentials have to be at least as promising as the pre-optimized case.
        if pre_optimized_set is not None:
            if not self.is_branch_promising(
                    station_eval, self.electrified_station_set, pre_optimized_set, missing_energy):
                # best station id is none and do not go deeper in recursion
                node_name = opt_util.stations_hash(self.electrified_station_set)
                self.current_tree[node_name]["viable"] = False
                raise opt_util.SuboptimalSimulationException

        min_nr_visited = float('inf')
        # find the lowest amount of visits in the possible stations /  children
        # not viable stations are skipped
        for station, _ in station_eval:
            # create a station combination from already electrified
            # stations and possible new stations
            stats = self.electrified_station_set.union([station])
            node_name = opt_util.stations_hash(stats)
            if not self.current_tree[node_name]["viable"]:
                continue
            min_nr_visited = min(min_nr_visited, self.current_tree[node_name]["visit_counter"])
        # simply visit the least visited node, starting from the one with the highest potential
        # not viable stations are skipped
        for stat_eval in station_eval:
            station, _ = stat_eval
            # create a station combination from already electrified stations
            # and possible new station
            stats = self.electrified_station_set.union([station])
            node_name = opt_util.stations_hash(stats)
            if not self.current_tree[node_name]["viable"]:
                continue

            if self.current_tree[opt_util.stations_hash(stats)]["visit_counter"] == min_nr_visited:
                best_station_id = station
                return [best_station_id], True

        # not a single children node was deemed viable. Therefore the current node itself can not be
        # viable
        node_name = opt_util.stations_hash(self.electrified_station_set)
        self.current_tree[node_name]["viable"] = False
        raise opt_util.SuboptimalSimulationException

    def set_battery_and_charging_curves(self):
        """ Set battery and charging curves from config. """
        for v_type in self.schedule.vehicle_types.values():
            for vehicle in v_type.values():
                if self.config.battery_capacity is not None:
                    vehicle["capacity"] = self.config.battery_capacity
                if self.config.charging_curve is not None:
                    vehicle["charging_curve"] = self.config.charging_curve

    def set_up_decision_tree(self, group_amount):
        """ Load decision tree if given in the config.

        :param group_amount: amount of groups for decision tree
        :type group_amount: int
        """
        if self.config.decision_tree_path is not None:
            with open(self.config.decision_tree_path, "rb") as file:
                self.decision_trees = pickle.load(file)
        else:
            self.decision_trees = [{} for _ in range(group_amount)]

    def rebase_spice_ev(self):
        """ Run SpiceEV simulation with new schedule and parameters.

        Configure various variables according to the input data and run a SpiceEV simulation.

        :return: must_include_set and electrified_stations
        :rtype: (set,dict)
        """
        self.logger.debug("Rebasing Scenario")
        must_include_set, ele_stations = self.preprocessing_scenario(
            run_only_neg=self.config.run_only_neg)
        self.logger.debug("Rebasing finished")
        if self.config.pickle_rebased:
            opt_util.toolbox_to_pickle(
                self.config.pickle_rebased_name, self.schedule, self.scenario, self.args)

            self.logger.debug("Rebased scenario pickled  as %s", self.config.pickle_rebased_name)

        self.must_include_set = must_include_set
        self.base_stations = self.electrified_stations.copy()
        return must_include_set, ele_stations

    def rebase_simple(self):
        """ Configure various variables according to the input data.

        Rebuild scenario with new schedule parameters.

        :return: must_include_set and electrified_stations
        :rtype: (set,dict)
        """
        must_include_set = set()
        # electrify inclusion stations
        for stat in self.config.inclusion_stations:
            self.electrify_station(stat, must_include_set)
        self.schedule.rotations = {r: self.schedule.rotations[r] for r in self.schedule.rotations
                                   if r not in self.config.exclusion_rots}
        self.must_include_set = must_include_set
        self.base_stations = self.electrified_stations.copy()
        return must_include_set, self.electrified_stations

    def preprocessing_scenario(self, electrified_stations=None, run_only_neg=False):
        """ Prepare scenario and run schedule.

        :param electrified_stations: optional electrified stations to use in simulation.
            Default None leads to using optimizer.electrified_stations
        :type electrified_stations: dict
        :param run_only_neg: should only negative rotations be simulated
        :type run_only_neg: bool
        :return: stations that must be included and stations which are electrified
        :rtype: (set, dict)
        """
        if electrified_stations is None:
            electrified_stations = self.electrified_stations
        must_include_set = set()
        # electrify inclusion stations
        for stat in self.config.inclusion_stations:
            self.electrify_station(stat, must_include_set)

        # calc new but only prev. negative rotations
        if run_only_neg:
            rots = self.schedule.get_negative_rotations(self.scenario)
            rots = {r: self.schedule.rotations[r] for r in rots
                    if r not in self.config.exclusion_rots}
            self.schedule.rotations = rots
        else:
            rots = {r: self.schedule.rotations[r] for r in self.schedule.rotations
                    if r not in self.config.exclusion_rots}
            self.schedule.rotations = rots

        new_sched, new_scen = opt_util.run_schedule(
            self.schedule, self.args, electrified_stations)
        self.schedule = new_sched
        self.scenario = new_scen
        self.must_include_set = must_include_set
        return must_include_set, electrified_stations

    def electrify_station(self, stat, electrified_set):
        """ Electrify a station and keep track of it in the electrified set file.

        :param stat: station id to be electrified
        :param stat: str
        :param electrified_set: the group of stations which is mutated along the optimization
        :type electrified_set: set(str)
        """
        self.electrified_stations[stat] = self.config.standard_opp_station
        electrified_set.add(stat)

    def create_charging_curves(self):
        """ Create charging curves with energy supplied over time for all vehicles.

        Cycle through vehicles and create numerically created charging curves with energy
        supplied over time, taking efficiencies into consideration """
        soc_charge_curve_dict = {}
        for v_type_name in self.schedule.vehicle_types:
            soc_charge_curve_dict[v_type_name] = {}
        for name, v_type in self.schedule.vehicle_types.items():
            for ch_type, data in v_type.items():
                eff = self.scenario.components.vehicle_types[f"{name}_{ch_type}"].battery_efficiency
                soc_charge_curve_dict[name][ch_type] = opt_util.charging_curve_to_soc_over_time(
                    data["charging_curve"], data["capacity"], self.args.desired_soc_opps,
                    self.schedule.cs_power_opps, efficiency=eff,
                    time_step=0.1, eps=self.config.eps)
        self.soc_charge_curve_dict = soc_charge_curve_dict

    @opt_util.time_it
    def get_critical_stations_and_rebase(self, relative_soc=False):
        """ Get the stations that must be electrified for full electrification of the system.

        Get the stations that must be electrified for full electrification of the system and put
        them into the electrified stations and an extra set of critical stations.
        Electrify every station but one. If without this single station there are below zero SoC
        events it is a critical station.

        :param relative_soc: should the evaluation use the relative or absolute SoC
        :param relative_soc: bool
        :return: Group of stations which have to be part of a fully electrified system
        :rtype: set(str)
        """

        events = self.get_low_soc_events(rel_soc=relative_soc)

        stats = {station for event in events for station in event.stations_list
                 if station not in self.not_possible_stations}

        electrified_station_set_all = set(stats)
        critical_stations = set()
        nr_of_all_stations = len(electrified_station_set_all)
        self.logger.log(msg=f"Electrifying {nr_of_all_stations-1}/{nr_of_all_stations} stations",
                        level=100)
        for station in sorted(electrified_station_set_all):
            electrified_station_set = electrified_station_set_all.difference([station])
            electrified_stations = self.electrified_station_set.union(electrified_station_set)
            vehicle_socs = self.timeseries_calc(electrified_stations)
            soc_min = 1
            for rot in self.schedule.rotations:
                soc, start, end = self.get_rotation_soc(rot, vehicle_socs)

                soc_min = min(np.min(soc[start:end]), soc_min)
                if soc_min < self.config.min_soc:
                    critical_stations.add(station)
            self.logger.debug("%s , with min soc: %s", station, soc_min)

        self.not_possible_stations = self.not_possible_stations.union(critical_stations)

        # rebasing
        for stat in critical_stations:
            # do not put must stations in electrified set, but in extra set must_include_set
            self.electrify_station(stat, self.must_include_set)
        return critical_stations

    def replace_socs_from_none_to_value(self):
        """ Removes SoC values of None by filling them with the last value which is not None.

        The function only changes None values which are at the end of SoC time series. Data type
        is switched to np.array for easier handling at later stages.
        """
        # make sure no None values exist in SOCs. Fill later values with last value
        # which was not None, eg, [1,5,4,None,None] becomes [1,5,4,4,4]
        for v_id, soc in self.scenario.vehicle_socs.items():
            soc = np.array(soc)
            if soc[-1] is None:
                for i in range(len(soc) - 1, -1, -1):
                    if soc[i] is not None:
                        # value is the last value which is not None. Overwrite the None values until
                        # the end of the list with this value
                        soc[i+1:] = soc[i]
                        break
            self.scenario.vehicle_socs[v_id] = soc

    def get_rotation_soc(self, rot_id, soc_data: dict = None):
        """ Gets the SoC object with start and end index for a given rotation id.

        :param rot_id: rotation_id
        :param soc_data: optional soc_data if not the scenario data should be used
        :return: tuple with SoC array, start index and end index
        """
        return opt_util.get_rotation_soc(
            rot_id, self.schedule, self.scenario, soc_data=soc_data)

    def get_index_by_time(self, search_time: datetime):
        """Get the index for a given time.

        :param search_time: The time for which to return the index as datetime object
        :type search_time: datetime
        :return: interval index. In case of no exact match, the index before is given.
        :rtype: int
        """
        return opt_util.get_index_by_time(self.scenario, search_time)

    @opt_util.time_it
    def get_trips(self, rot: rotation.Rotation, start_idx: int, end_idx: int):
        """ Return trips from rotation with start to end index.

        Get trips in a rotation from a start to an end index, if the arrival time is in between
        the start and end idx.

        :param rot: The rotation object containing the trips.
        :type rot: simba.rotation.Rotation
        :param start_idx: start index, representing the start time.
        :type start_idx: int
        :param end_idx: end index, representing the end time.
        :type end_idx: int
        :return: trips that arrive between the start and end time.
        :rtype: list(simba.trip.Trip)
        """

        start_time_event = self.get_time_by_index(start_idx)
        end_time_event = self.get_time_by_index(end_idx)

        trips = []
        for i, trip in enumerate(rot.trips):
            if end_time_event >= trip.arrival_time > start_time_event:
                trips.append(trip)

        return trips

    @opt_util.time_it
    def get_time_by_index(self, idx):
        """ Get the time for a given index.

        :param idx: The index for which to return the time.
        :type idx: int
        :return: the time corresponding to the given index.
        """

        start_time = self.scenario.start_time

        delta_time = self.scenario.interval
        searched_time = start_time + delta_time * idx
        return searched_time

    @opt_util.time_it
    def get_low_soc_events(self, rotations: Iterable = None, filter_standing_time=True,
                           rel_soc=False, soc_data=None, **kwargs):
        """ Return low SoC events below the config threshold.

        :param rotations: rotation_ids to be searched for low SoC events. Default None means whole
            schedule is searched
        :type rotations: Iterable
        :param filter_standing_time: Should the stations be filtered by standing time. True leads to
            an output with only stations with charging potential
        :type filter_standing_time: bool
        :param rel_soc: Defines if the start SoC should be seen as full even when not.
            i.e. a drop from a rotation start SoC from 0.1 to -0.3 is not critical since the start
            SoC from the rotation will be raised
            If the rel_soc
            is false, it means coupled rotations might be prone to errors due to impossible lifts.
        :param soc_data: SoC data to be used. Default None means the SoC data from the optimizer
            scenario is used
        :type soc_data: dict
        :param kwargs: optional soc_lower_thresh or soc_upper_thresh if from optimizer differing
            values should be used
        :raises InfiniteLoopException: If while loop does not change during iterations
        :return: low SoC events
        :rtype: list(simba.optimizer_util.LowSocEvent)
        """
        if not rotations:
            rotations = self.schedule.rotations

        soc_lower_thresh = kwargs.get("soc_lower_thresh", self.config.min_soc)
        soc_upper_thresh = kwargs.get("soc_upper_thresh", self.args.desired_soc_opps)
        # create list of events which describe trips which end in a SoC below zero
        # the event is bound by the lowest SoC and an upper SoC threshold which is naturally 1
        # properties before and after these points have no effect on the event itself, similar to
        # an event horizon
        events = []

        for rot_id in rotations:
            events_per_rotation = 0
            rot = self.schedule.rotations[rot_id]
            soc, rot_start_idx, rot_end_idx = self.get_rotation_soc(rot_id, soc_data)
            rot_end_idx += 1
            idx = range(0, len(soc))

            # if rotation gets a start SoC below the args.desired_soc_deps, this should change
            # below 0 SoC events, since fixing the rotation before leads to fixing this rotation.
            # If using relative SOC, SOC lookup has to be adjusted
            if rel_soc:
                soc = np.array(soc)
                # if the rotation starts with lower than desired soc, lift the SoC
                if soc[rot_start_idx] < self.args.desired_soc_deps:
                    soc = self.lift_and_clip_positive_gradient(rot_start_idx, soc, soc_upper_thresh)

            # combined data of the SoC data of the rotation and the original index
            # the array will get masked later, so the index of the array might be different
            # to the index data column
            soc_idx = np.array((soc, idx))[:, rot_start_idx:rot_end_idx]

            # Mask for soc_idx which is used to know if an index has been checked or not
            mask = np.ones(len(soc_idx[0])).astype(bool)

            # get the minimum SoC and index of this value. The mask does nothing yet
            min_soc, min_idx = get_min_soc_and_index(soc_idx, mask)

            # Used to check if an infinite loop is happening
            old_idx = -1

            # while the minimal SoC of the SoC time series is below the threshold find the events
            # which are connected with this event. Every iteration the data of the found events is
            # removed. At some point the reduced time series is either empty or does not have a
            # SoC below the lower threshold.
            while min_soc < soc_lower_thresh:
                # soc_idx is of type float, so the index needs to be cast to int
                min_idx = int(min_idx)

                if min_idx == old_idx:
                    raise opt_util.InfiniteLoopException
                else:
                    old_idx = min_idx
                i = min_idx

                # find the first index by going back from the minimal SoC index, where the SoC
                # was above the upper threshold OR the index where the rotation started.
                while soc[i] < soc_upper_thresh:
                    if i <= rot_start_idx:
                        break
                    i -= 1

                # which index in the original data is the found index?
                start = i

                i = min_idx
                # do the same as before but find the index after the minimal SoC.
                while soc[i] < soc_upper_thresh:
                    if i >= rot_end_idx - 1:
                        break
                    i += 1
                end = i

                # with the start index and the minimal index find trips, in this time span
                trips = self.get_trips(rot=rot, start_idx=start, end_idx=min_idx)
                possible_stations = set()
                possible_stations_list = []

                # check which stations the bus arrives at with these trips. Depending on
                # filter_standing_time, every station is returned or only the stations where the bus
                # stops long enough to be able to charge.
                if not filter_standing_time:
                    possible_stations = {t.arrival_name for t in trips}
                    possible_stations_list = [t.arrival_name for t in trips]
                else:
                    for ii, trip in enumerate(trips):
                        try:
                            standing_time_min = opt_util.get_charging_time(
                                trip, trips[ii + 1], self.args)
                        except IndexError:
                            standing_time_min = 0
                        if standing_time_min > 0:
                            possible_stations.add(trip.arrival_name)
                            possible_stations_list.append(trip.arrival_name)

                # the stations which have potential should be filtered by the not_possible_stations,
                # which are stations which for some reason should not be electrified.
                possible_stations = possible_stations.difference(self.not_possible_stations)
                ch_type = "depb" if "depb" in rot.vehicle_id else "oppb"
                v_type = rot.vehicle_id.split("_" + ch_type)[0]

                # with the gathered data create the event object
                event = opt_util.LowSocEvent(
                    start_idx=start, end_idx=min_idx, min_soc=min_soc, stations=possible_stations,
                    vehicle_id=rot.vehicle_id, trip=trips, rot=rot,
                    stations_list=possible_stations_list,
                    capacity=self.schedule.vehicle_types[v_type][ch_type]['capacity'],
                    v_type=v_type, ch_type=ch_type)

                events.append(event)
                events_per_rotation += 1

                if events_per_rotation > len(rot.trips):
                    self.logger.error("More low-soc-events than trips found for rotation "
                                      f"{rot.id} with vehicle {rot.vehicle_id}")

                # the mask is expanded to the just checked low_soc_event
                mask[start - rot_start_idx:end - rot_start_idx + 1] = False

                # if the mask does not leave any value, the loop is finished
                if not np.any(mask):
                    break
                # Check the remaining unmasked SoCs for the minimal SoC
                min_soc, min_idx = get_min_soc_and_index(soc_idx, mask)

        return events

    def lift_and_clip_positive_gradient(self, start_idx: int, soc: np.array,
                                        soc_upper_thresh: float):
        """ Lift a SoC array to a given soc_upper_thresh and handle clipping

        Lifts a given SoC array to the value of soc_upper_thresh from the start_index to the end.
        Values which exceed soc_upper_thresh and have a positive gradient are clipped to
        soc_upper_thresh. If values exist which exceed soc_upper_thresh and have negative gradients,
        they and the following SoCs have this delta subtracted instead. This is the behavior of
        approximate behaviour of charging a vehicle above 100%. Values can not exceed 100%. As soon
        as discharge happens, the value drops.

        :param start_idx: index from where the array is lifted
        :type start_idx: int
        :param soc: array which contains SoCs to be lifted
        :type soc: np.array
        :param soc_upper_thresh: max value of SoCs used for clipping
        :type soc_upper_thresh: float
        :return: None

        """

        delta = self.args.desired_soc_deps - soc[start_idx]
        # lift all beyond the start index (including the start_idx)
        soc[start_idx:] += delta
        # values before the start index stay unaffected
        soc_pre = soc[:start_idx]

        # to speed up computation, only look at the affected part
        soc = soc[start_idx:]
        soc_max = np.max(soc)
        # Clip SoC values to soc_upper_thresh.
        # Multiple local maxima above this value can exist which,
        # if they are rising over time, need multiple clippings.
        while soc_max > soc_upper_thresh:
            # descending array from len(soc) to 0
            desc = np.arange(len(soc), 0, -1)
            # gradient of SoC i.e. positive if charging, negative if discharging
            diff = np.hstack((np.diff(soc), -1))
            # masking of SoCs >1 and negative gradient for local maximum
            # i.e. after lifting the SoC, it finds the first spot where the SoC is bigger
            # than the upper threshold and descending.
            idc_loc_max = np.argmax(desc * (soc > 1) * (diff < 0))

            # find the SoC value of this local maximum
            soc_max = soc[idc_loc_max]
            # reducing everything after local maximum
            soc[idc_loc_max:] = soc[idc_loc_max:] - (soc_max - soc_upper_thresh)

            # clip everything before local maximum
            soc[:idc_loc_max][soc[:idc_loc_max] > 1] = soc_upper_thresh
            soc_max = np.max(soc)
        # restore SoC with unaffected part
        soc = np.hstack((soc_pre, soc))
        return soc


def get_min_soc_and_index(soc_idx, mask):
    """ Returns the minimal SoC and the corresponding index of a masked soc_idx.

    :param soc_idx: SoC values and their corresponding index of shape (2,nr_values).
    :type soc_idx: np.array
    :param mask: boolean mask. It is false where the soc_idx has been checked already
    :type mask: np.array
    :return: minimal SoC and corresponding index
    :rtype: tuple(Float, Float)
    """
    return soc_idx[:, mask][:, np.argmin(soc_idx[0, mask])]


def get_init_node():
    """ Returns the initialization dictionary / node of the decision tree.

    :return: init node
    :rtype: dict
    """
    return {"completed": False, "viable": True,
            "missing_energy": None, "visit_counter": 0, "children": set()}
