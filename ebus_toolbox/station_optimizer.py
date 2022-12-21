""" Optimizer class which implements the optimizer object and methods needed"""
import json
import typing
from copy import copy, deepcopy
from datetime import timedelta
import pickle
from pathlib import Path
import numpy as np
import logging

import schedule
from ebus_toolbox import report, rotation
from src import scenario
from optimizer_config import OptimizerConfig
from low_soc_events import evaluate, get_low_soc_events, get_groups_from_events, get_delta_soc, \
    get_buffer_time, get_charging_time, get_index_by_time,\
    get_rotation_soc as get_rotation_soc_low_soc
import station_optimization as opt


class StationOptimizer:
    """ Class for station optimization"""

    def __init__(self, sched: schedule.Schedule, scen: scenario.Scenario, args,
                 config: OptimizerConfig(), logger: logging.Logger, **kwargs):
        self.base_not_possible_stations = set()
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

    def loop(self, **kwargs):
        """ Loops over every base group with the core of group optimization and handles deep or
        greedy arguments

        :param opt_type: type of optimization. Default from optimizer object
        :param kwargs: arguments handled by the function and config.
        :return: tuple of the dict electrified_stations and the set electrified_station_set
        in the optimal found case
        """

        node_choice=kwargs.get("node_choice", self.config.node_choice)
        opt_type=kwargs.get("opt_type", self.config.opt_type)
        self.base_scenario = copy(self.scenario)
        self.base_schedule = copy(self.schedule)

        self.base_stations =  self.electrified_stations.copy()
        self.base_electrified_station_set = self.electrified_station_set.copy()
        self.base_not_possible_stations = self.not_possible_stations.copy()

        # get events where soc fell below 0. The events contain info about the problematic
        # time span, which includes stations which could provide a soc lift
        base_events = get_low_soc_events(self, rel_soc=False,
                                         soc_data=self.base_scenario.vehicle_socs)

        # check if the events can be divided into subgroups which are independent
        # this makes optimization in smaller groups possible
        groups = get_groups_from_events(base_events, self.not_possible_stations,
                                        could_not_be_electrified=self.could_not_be_electrified,
                                        optimizer=self)

        # storage for the sets of stations which will be generated
        list_greedy_sets = [set()] * len(groups)

        # baseline greedy optimization
        # base line is created simply by not having a decision tree and not a pre optimized_set yet
        for group_nr, group in enumerate(groups[:]):
            events, stations = group
            linien = {lne for e in events for lne in e.rotation.lines}
            self.logger.warning("Optimizing %s out of %s. This includes these Lines", group_nr + 1,
                                len(groups))
            self.logger.warning(linien)
            self.logger.warning("%s events", (len(events)))

            # first run is always step by step
            self.group_optimization(group, self.choose_station_step_by_step,
                                    events_remaining=[len(events)], **kwargs)

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
            if node_choice == "brute":
                choice_func = self.choose_station_brute
            else:
                choice_func = self.choose_station_step_by_step

            combinations = opt.combs_unordered_no_putting_back(len(stations),
                                                               len(list_greedy_sets[group_nr]) - 1)
            self.logger.debug("There are %s combinations with 1 station less than "
                              "the current solution", combinations)

            pre_optimized_set = list_greedy_sets[group_nr]

            while i < self.config.max_brute_loop and cont_loop:
                i += 1
                if i % 10 == 0:
                    print(len(self.decision_tree), " nodes checked")
                    print(f"Optimal solution has length {len(pre_optimized_set)}")
                self.electrified_stations = self.base_stations.copy()
                self.electrified_station_set= self.base_electrified_station_set.copy()
                self.scenario = copy(self.base_scenario)
                self.schedule = copy(self.base_schedule)

                self.scenario.vehicle_socs = deepcopy(self.base_scenario.vehicle_socs)
                self.schedule.rotations = deepcopy(self.base_schedule.rotations)
                self.not_possible_stations = self.base_not_possible_stations.copy()

                # create a new electrified set
                new_stats, cont_loop = self.group_optimization(group, choice_func,
                                                               track_not_possible_rots=False,
                                                               pre_optimized_set=pre_optimized_set,
                                                                events_remaining = [len(events)],
                                                               **kwargs)

                new_electrified_set = self.electrified_station_set
                # if a new set was found, print it and save it in sols
                if new_electrified_set != pre_optimized_set and new_stats is not None:
                    self.logger.warning("Optimized with %s  stations %s %s",
                                        len(new_electrified_set), str('#' * 20),
                                        stations_hash(new_electrified_set))
                    sols.append(new_electrified_set)
                    if len(new_electrified_set) < len(pre_optimized_set):
                        list_greedy_sets[group_nr] = new_electrified_set.copy()
                        pre_optimized_set = new_electrified_set.copy()

            self.logger.debug("All solutions for this group: %s", sols)
            self.logger.debug("Optimized with %s stations out of %s", len(pre_optimized_set),
                              len(stations))

        # saving decision tree only in case of deep analysis
        if self.config.save_decision_tree:
            with open(self.args.output_directory / Path("decision_tree.pickle"), "wb") as file:
                pickle.dump(self.decision_tree, file)

        self.electrified_stations = self.base_stations.copy()
        self.electrified_station_set = self.base_electrified_station_set.copy()
        for single_set in list_greedy_sets:
            for stat in single_set:
                self.electrify_station(stat, self.electrified_station_set)

        return self.electrified_stations, self.electrified_station_set

    def get_negative_rotations_all_electrified(self, rel_soc=False):
        """Get the rotation ids for the rotations which show negative socs, even when everything
        is electrified

        :param rel_soc: if true, the start soc is handled like it has the desired deps soc
        :return: set of rotation ids which are negative even with all stations electrified
        """
        events = get_low_soc_events(self, rel_soc=rel_soc)

        stats = {stat for event in events for stat in event.stations_list
                 if stat not in self.not_possible_stations}
        electrified_station_set = set(stats)
        vehicle_socs = self.timeseries_calc(ele_station_set=electrified_station_set)
        new_events = get_low_soc_events(self, soc_data=vehicle_socs)
        return {event.rotation.id for event in new_events}

    def group_optimization(self, group, choose_station_function, track_not_possible_rots=True,
                           pre_optimized_set=None, **kwargs):
        """Optimize a single group events

        :param group: tuple of events and station_ids to be optimized together
        :param choose_station_function: function to be used for choosing the next function
        :param track_not_possible_rots: not possible tracks need to be tracked in the first run
        :param pre_optimized_set: set which was optimized before hand used as upper threshold for
        optimization
        :param kwargs: internally used functionality
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

        solver = self.config.solver

        tree_position = kwargs.get("tree_position", [1])

        # give tree position
        self.logger.debug("%s with length of %s and a total of %s "
                          "stations", tree_position, len(tree_position),
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
        station_eval = evaluate(event_group, self, soc_data=kwargs["lifted_socs"])

        # get rotations from event dict and calculate the missing energy
        rotation_dict = {e.rotation.id: e.rotation for e in event_group}
        missing_energy = get_missing_energy(event_group)
        if missing_energy >= 0:
            self.logger.debug("Already electrified: Returning")
            return self.electrified_stations, True

        # get the best stations, brute or step_by_step
        best_station_ids, recursive = choose_station_function(station_eval,
                                                              pre_optimized_set=pre_optimized_set,
                                                              missing_energy=missing_energy)

        stat_eval_dict = {stat_id[0]: stat_id[1] for stat_id in station_eval}
        if best_station_ids is None:
            self.logger.warning("No useful station found with %s rotations not electrified yet. "
                                "Stopped after electrifying %s", events_remaining,
                                len(self.electrified_station_set))

            if pre_optimized_set is None:
                could_not_be_electrified.update(list(rotation_dict.keys()))
                return None, False

            # remove electrified stations in this run
            copied_set = self.electrified_station_set.copy()
            for stat in copied_set:
                self.electrified_stations.pop(stat)
                self.electrified_station_set.remove(stat)
            # overwrite with preo ptimized set
            for stat in pre_optimized_set:
                self.electrify_station(stat, self.electrified_station_set)
            return None, False

        self.logger.debug("%s, with first pot of %s", best_station_ids,
                          round(stat_eval_dict[best_station_ids[0]],1))

        # electrify station(s)
        for stat_id in best_station_ids:
            self.electrify_station(stat_id, self.electrified_station_set)

        event_rotations = {event.rotation for event in event_group}

        # copy scen and sched
        self.scenario = copy(self.base_scenario)
        self.schedule = copy(self.base_schedule)

        if solver == "quick":
            # quick calculation has to electrify everything in one step,
            # or the lifting of socs is not correct
            self.scenario.vehicle_socs = deepcopy(self.base_scenario.vehicle_socs)
            self.scenario.vehicle_socs = self.timeseries_calc(event_rotations,
                                                              electrify_stations=best_station_ids)
        else:
            self.schedule.rotations = rotation_dict
            self.schedule, self.scenario = \
                opt.run_schedule(self.schedule, self.args,
                                 electrified_stations=self.electrified_stations)

        kwargs["lifted_socs"] = self.scenario.vehicle_socs.copy()

        self.not_possible_stations = set(self.electrified_stations.keys()).union(
            self.not_possible_stations)
        event_rotations_ids = {event.rotation.id for event in event_group}
        new_events = get_low_soc_events(self, rotations=event_rotations_ids,
                                        rel_soc=True, **kwargs)

        delta_energy = get_missing_energy(new_events)
        events_remaining[0] -= len(event_group) - len(new_events)
        self.logger.debug("Last electrification electrified %s/%s."
                          " %s remaining events in the base group.",
                          len(event_group) - len(new_events), len(event_group), events_remaining[0])
        delta_base_energy = delta_energy  # get_missing_energy(base_events)

        # put this node intro the decision tree including the missing energy
        self.node_to_tree(delta_base_energy)

        # everything electrified
        if delta_energy >= 0:
            return self.electrified_stations, True
        # some choice functions might not need a recursive call, they return here. recursive is set
        # by the choose_station_function
        if not recursive:
            return None, True

        # check if the events can be divided into subgroups which are independent
        groups = get_groups_from_events(new_events, self.not_possible_stations,
                                        could_not_be_electrified,
                                        optimizer=self)

        for k, this_group in enumerate(groups):
            this_tree = tree_position.copy()
            this_tree.append(k+1)
            kwargs["tree_position"] = this_tree
            new_stations, _ = \
                self.group_optimization(this_group, choose_station_function,
                                        track_not_possible_rots=track_not_possible_rots,
                                        pre_optimized_set=pre_optimized_set,
                                        **kwargs)

            if new_stations is not None:
                self.electrified_stations.update(new_stations)
            else:
                return None, True

            # if there is no pre optimized set function can return
            if pre_optimized_set is None:
                continue

            # if there is a pre optimized set check if the preoptimized set leads to pruning
            # i.e. if if following  branch makes sense
            if len(pre_optimized_set) - len(self.electrified_station_set) < 10:
                self.scenario.vehicle_socs = \
                    self.timeseries_calc(event_rotations, electrify_stations=best_station_ids)

                prune_events = \
                    get_low_soc_events(self, rotations=event_rotations_ids, rel_soc=True, **kwargs)

                station_eval = evaluate(prune_events, self)
                prune_missing_energy = get_missing_energy(prune_events)
                if not self.is_branch_promising(station_eval, self.electrified_station_set,
                                                pre_optimized_set, prune_missing_energy):
                    self.logger.debug("Branch pruned early")
                    return None, True

        return self.electrified_stations, True

    # ToDo implement a fast calculation of timeseries_calc but with a limited amount of charging
    #  points. Implementation could look like this
    # For each possible station find the charging times of each vehicle, e.g
    # PotsdamerPlatz: Vehicle1: 12:30-12:40, 13:20-14:00
    #                 Vehicle5: 9:30-9:40, 13:10-13:25
    # For each station sort the events accordingly, e.g.
    # Step 1
    # PotsdamerPlatz: 9:30-9:40 Vehicle5
    #                 12:30-12:40 Vehicle1
    #                 13:10-13:25 Vehicle5
    #                 13:20-14:00 Vehicle1
    # For a number of charging points mutate the events, e.g. 2 charging points
    #
    # Step2
    # PotsdamerPlatz: 9:30-9:40 Vehicle5
    #                 12:30-12:40 Vehicle1
    #                 13:10-13:25 Vehicle5
    #                 13:20-14:00 Vehicle1
    # or 1 charging piont
    # PotsdamerPlatz: 9:30-9:40 Vehicle5
    #                 12:30-12:40 Vehicle1
    #                 13:10-13:25 Vehicle5
    #                 13:25-14:00 Vehicle1
    # do this for every station
    # Step 3
    # Sort all events of all stations by time
    # PotsdamerPlatz: 9:30-9:40 Vehicle5
    # Hauptbahnhof: 9:35-9:45
    # ...
    # go through events and lift the socs according to timeseries_calc
    # check if socs are clipped. if so mutate the event time and repeat step2 and 3.
    # go further

    def timeseries_calc(self, rotations=None, soc_dict=None, ele_station_set=None,
                        soc_upper_thresh=None, electrify_stations=None) -> object:
        """ A quick estimation of socs by mutating the soc data according to full electrification.
        This means a electrified station has unlimited charging points

        :param rotations: Optional if not optimizer.schedule.rotations should be used
        :param soc_dict: Optional if not optimizer.scenario.vehicle_socs should be used
        :param ele_station_set: set to be electrified. Default None leads to using the
        so far optimized optimizer.electrified_station_set
        :param soc_upper_thresh: optional if not optimizer.config.desired_soc_deps is used
        :param electrify_stations: dict of stations to be electrified
        :return soc dict with lifted socs
        :rtype dict()
        """
        if not rotations:
            rotations = self.schedule.rotations.values()

        if not soc_dict:
            soc_dict = self.scenario.vehicle_socs
        if not ele_station_set:
            ele_station_set = self.electrified_station_set
        if not soc_upper_thresh:
            soc_upper_thresh = self.args.desired_soc_deps
        if not electrify_stations:
            electrify_stations = set()

        if not electrify_stations:
            electrify_stations = set()
        ele_stations = {*ele_station_set, *electrify_stations}
        soc_dict = copy(soc_dict)

        for rot in rotations:
            ch_type = (rot.vehicle_id.find("oppb") > 0) * "oppb" + (
                    rot.vehicle_id.find("depb") > 0) * "depb"
            v_type = rot.vehicle_id.split("_" + ch_type)[0]
            soc_over_time_curve = self.soc_charge_curve_dict[v_type][ch_type]
            soc = soc_dict[rot.vehicle_id]
            for i, trip in enumerate(rot.trips):
                if trip.arrival_name not in ele_stations:
                    continue
                idx = get_index_by_time(self.scenario, trip.arrival_time)
                try:
                    standing_time_min = get_charging_time(trip, rot.trips[i + 1], self.args)
                except IndexError:
                    standing_time_min = 0

                d_soc = get_delta_soc(soc_over_time_curve, soc[idx], standing_time_min, self)
                buffer_idx = int((get_buffer_time(trip, self.args.default_buffer_time_opps))
                                 / timedelta(minutes=1))
                delta_idx = int(standing_time_min) + 1
                old_soc = soc[idx + buffer_idx:idx + buffer_idx + delta_idx].copy()
                soc[idx + buffer_idx:] += d_soc
                soc[idx + buffer_idx:idx + buffer_idx + delta_idx] = old_soc
                soc[idx + buffer_idx:idx + buffer_idx + delta_idx] += np.linspace(0, d_soc,
                                                                                  delta_idx)

                soc_pre = soc[:idx]
                soc = soc[idx:]
                soc_max = np.max(soc)
                while soc_max > soc_upper_thresh:
                    # descending array
                    desc = np.arange(len(soc), 0, -1)
                    # gradient of soc i.e. positive if charging negative if discharging
                    diff = np.hstack((np.diff(soc), -1))
                    # masking of socs >1 and negative gradient for local maximum
                    # i.e. after lifiting the soc, it finds the first spot where the soc is bigger
                    # than the upper threshold and descending.
                    idc_loc_max = np.argmax(desc * (soc > 1) * (diff < 0))

                    # find the soc value of this local maximum
                    soc_max = soc[idc_loc_max]
                    # reducing everything after local maximum
                    soc[idc_loc_max:] = soc[idc_loc_max:] - (soc_max - 1)

                    # capping everything before local maximum
                    soc[:idc_loc_max][soc[:idc_loc_max] > 1] = soc_upper_thresh
                    soc_max = np.max(soc)
                soc = np.hstack((soc_pre, soc))
            soc_dict[rot.vehicle_id] = soc
        return soc_dict

    def is_branch_promising(self, station_eval, electrified_station_set,
                            pre_optimized_set, missing_energy):
        """Quickly evaluates if following a branch is promising by summing up estimated potentials
        :param station_eval: list of  sorted station evaluation with list(station_id, pot)
        :param electrified_station_set: set of electrified stations
        :param pre_optimized_set: set of already optimized stations
        :param missing_energy:
        :return is the branch promising
        :rtype bool
        """

        delta = len(pre_optimized_set) - len(electrified_station_set)
        pot = 0
        for i in range(0, min(delta, len(station_eval))):
            pot += station_eval[i][1]
        if pot < -missing_energy * self.config.estimation_threshold:
            print(f"Not enough potential {round(pot, 0)} / {round(-missing_energy, 0)} after ",
                  stations_hash(electrified_station_set))
            return False
        return True

    def node_to_tree(self, delta_base_energy):
        """ Fill decision tree with the given info
        :param delta_base_energy: missing energy
        :return decision tree
        :rtype dict()
        """
        node_name = stations_hash(self.electrified_station_set)
        try:
            self.decision_tree[node_name]["missing_energy"] = delta_base_energy
            self.decision_tree[node_name]["visit_counter"] += 1
            # todo add is_viable
            self.logger.debug("already visited")
        except KeyError:
            self.decision_tree[node_name] = {}
            self.decision_tree[node_name]["missing_energy"] = delta_base_energy
            self.decision_tree[node_name]["visit_counter"] = 1

    def choose_station_brute(self, station_eval,
                             pre_optimized_set=None, missing_energy=0,
                             gens=dict()):
        """Gives back a possible set of Stations to electrify which shows potential and hasnt been
        tried yet. The set of stations is smaller than the best optimized set so far.

        :param station_eval: list of  sorted station evaluation with list(station_id, pot)
        :param pre_optimized_set: set of optimized stations thus far
        :param missing_energy: missing energy in this branch before electrification
        :param gens: dict of generators for brute force generation
        :return: combination of stations or None if no viable stations exists and
        false since this function does not support recursive calling
        """

        station_ids = [x[0] for x in station_eval]
        try:
            generator = gens[str(station_ids) + str(len(pre_optimized_set) - 1)]
        except KeyError:
            generator = combination_generator(station_ids, len(pre_optimized_set) - 1)
            gens[str(station_ids) + str(len(pre_optimized_set) - 1)] = generator
        station_eval_dict = {stat[0]: stat[1] for stat in station_eval}
        for comb in generator:
            node_name = stations_hash(comb)
            if node_name not in self.decision_tree:
                # only check the brute force station if they have the remote chance of fulfilling
                # the missing energy
                # potential>missing energy * 80%
                potential = sum([station_eval_dict[stat] for stat in comb])
                if potential > -missing_energy * self.config.estimation_threshold:
                    return comb, False
                else:
                    self.logger.debug("skipped %s since potential is too low %s %%", comb,
                                      round(potential / -missing_energy * 100, 0))
        self.logger.debug("calculated all viable possibilities")
        return None, False

    def choose_station_step_by_step(self, station_eval,
                                    pre_optimized_set=None, missing_energy=0):
        """ Gives back a station of possible stations to electrify which shows the biggest potential
        and hasnt been picked before.


        :param station_eval: list of  sorted station evaluation with list(station_id, pot)
        :param pre_optimized_set: set of optimized stations thus far
        :param  missing energy in this branch before electrification
        :return:a station to electrify or None if no viable stations exists and
        false since this function does not support recursive calling
        """

        # filter functions to stop simulating cases which have no hope of being optimal.
        # if in optimization mode, optimization can break if station amount is superseded
        # this filter is done better by the next
        # if pre_optimized_set is not None:
        #     if len(electrified_station_set)>len(pre_optimized_set):
        #         return pre_optimized_set
        # potentials have to be at least as promising as the pre-optimized case
        if pre_optimized_set is not None:
            if not self.is_branch_promising(station_eval, self.electrified_station_set,
                                            pre_optimized_set, missing_energy):
                # best station id is none and dont go deeper in recursion
                return None, False

        min_nr_visited = float('inf')
        for stat_tuple in station_eval:
            station, pot = stat_tuple
            # create a station combination from already electrified
            # stations and possible new station
            check_stations = self.electrified_station_set.union([station])

            if self.decision_tree is not None:
                node_name = stations_hash(check_stations)
                if node_name in self.decision_tree.keys():
                    min_nr_visited = min(min_nr_visited,
                                         self.decision_tree[node_name]["visit_counter"])
                    # if already checked skip to next one
                    continue
            if station not in self.electrified_station_set:
                best_station_id = station
                return [best_station_id], True
        # every possible node from here was evaluated already
        # to do what now?
        # simply visit the least visited node
        for stat_eval in station_eval:
            station, eval = stat_eval
            # create a station combination from already electrified stations
            # and possible new station
            check_stations = self.electrified_station_set.union([station])
            if self.decision_tree[stations_hash(check_stations)]["visit_counter"] == min_nr_visited:
                best_station_id = station
                return [best_station_id], True
        return None, True

    def set_battery_and_charging_curves(self):
        """ Create battery and charging curves for fast lookup
        :return: Void
        """
        for v_type in self.schedule.vehicle_types.values():
            for vehicle in v_type.values():
                if self.config.battery_capacity is not None:
                    vehicle["capacity"] = self.config.battery_capacity
                if self.config.charging_curve is not None:
                    vehicle["charging_curve"] = self.config.charging_curve

    def set_up_decision_tree(self):
        """ Load decision tree if given in the config
        :return: Void
        """
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
        # electrify inclusion stations
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
        # electrify inclusion stations
        for stat in self.config.inclusion_stations:
            self.electrify_station(stat, must_include_set)

        # calc new but only prev. negative rotations
        if run_only_neg:
            rots = self.schedule.get_negative_rotations(self.scenario)
            rots = {r: self.schedule.rotations[r] for r in rots if r not in
                    self.config.exclusion_rots}
            self.schedule.rotations = rots
        else:
            rots = {r: self.schedule.rotations[r] for r in self.schedule.rotations if
                    r not in self.config.exclusion_rots}
            self.schedule.rotations = rots

        new_sched, new_scen = opt.run_schedule(self.schedule, self.args,
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
        """Cycle through vehicles and create numerically created charging curves with energy
        supplied over time """
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

        events = get_low_soc_events(self, rel_soc=relative_soc)

        stats = {station for event in events for station in event.stations_list
                 if station not in self.not_possible_stations}

        electrified_station_set_all = set(stats)
        must_stations = set()
        print(f"Electrifying {len(electrified_station_set_all)} stations minus one", )
        for station in electrified_station_set_all:
            electrified_station_set = electrified_station_set_all.difference([station])
            vehicle_socs = self.timeseries_calc(electrify_stations=electrified_station_set)
            soc_min = 1
            for rot in self.schedule.rotations:
                soc, start, end = self.get_rotation_soc(rot, vehicle_socs)

                soc_min = min(np.min(soc[start:end]), soc_min)
                if soc_min < self.config.min_soc:
                    must_stations.add(station)
            self.logger.debug("%s , with min soc: %s", station, soc_min)

        self.not_possible_stations = self.not_possible_stations.union(must_stations)

        # rebasing
        for stat in must_stations:
            # dont put must stations in electrified set, but in extra set must_include_set
            self.electrify_station(stat, self.must_include_set)

        self.scenario.vehicle_socs = self.timeseries_calc(ele_station_set=must_stations)
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
        """ gets you the soc object with start and end index for a given rotation id
        :param rot_id: rotation_id
        :param soc_data: optional soc_data if not the scenario data should be used
        :return: tuple with soc array, start index and end index
        """
        return get_rotation_soc_low_soc(rot_id, self.schedule, self.scenario, soc_data=soc_data)


    def get_index_by_time(self, search_time):
        start_time = self.scenario.start_time
        delta_time = timedelta(minutes=60 / self.scenario.stepsPerHour)
        idx = (search_time - start_time) // delta_time
        return idx

    def get_trips(self, rot: rotation.Rotation, start_idx: int, end_idx: int):
        # return trips in a rotation from a start to an end index, if the arrival time is in between
        # the start and end idx
        start_time_event = self.get_time_by_index(start_idx)
        end_time_event = self.get_time_by_index(end_idx)

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
    # fill the soc completely in last time step
    times.append(time)
    socs.append(args.desired_soc_opps)
    return np.array((times, socs)).T


def get_missing_energy(events):
    """ Sum up all the missing energies of the given events"""
    missing_energy = 0
    for event in events:
        missing_energy += event.min_soc * event.capacity
    return missing_energy


def stations_hash(stations_set):
    """ Create a simple str as hash for a set of stations"""
    return str(sorted(list(stations_set)))


def combination_generator(iterable: typing.Iterable, amount: int):
    """ Generator which yields all possible combinations of choosing
    an amount out of an iterable without putting them back and without caring about the
    order of elements
    :param iterable: Any collection which can be cast to a list
    :param amount: Number of elements which should be drawn from iterable
    :type amount: int
    """
    iterable = list(iterable)

    for i, item in enumerate(iterable):
        # recursive calling of generator with clock like behavior, e.g right-most item changes until
        # end of list is reached. This leads to a change in the item left to it and so on. Elements
        # on the right can only change to a subset of bigger indices than their left counter-part.
        # this is due to the ignoring of order, which reduces the amount of possibilities.
        if amount <= 1:
            yield [item]
        else:
            for gen in combination_generator(iterable[i + 1:], amount - 1):
                yield [item] + gen
