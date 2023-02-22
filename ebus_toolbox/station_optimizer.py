""" Optimizer class which implements the optimizer object and methods needed"""
import traceback
import logging
import pickle
import warnings
from copy import deepcopy, copy
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np

import ebus_toolbox.optimizer_util as util
from spice_ev import scenario
from ebus_toolbox import report, rotation, schedule
from ebus_toolbox.util import uncomment_json_file


class StationOptimizer:
    """ Class for station optimization"""

    def __init__(self, sched: schedule.Schedule, scen: scenario.Scenario, args,
                 config: 'util.OptimizerConfig', logger: logging.Logger):
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
        with open(self.args.electrified_stations, "r", encoding="utf-8", ) as file:
            self.electrified_stations = uncomment_json_file(file)
        self.base_stations = self.electrified_stations.copy()
        self.base_electrified_station_set = set()

        # stations which are included can not be included again. Therefore they go into
        # the set of not possible stations together with excluded stations
        self.not_possible_stations = config.inclusion_stations.union(config.exclusion_stations)

        # if the config defines special rotations, only use them in the schedule
        if config.reduce_rots:
            self.schedule.rotations = {rot: sched.rotations[rot] for rot in config.rots}

    def loop(self, **kwargs):
        """ Loops over every base group with the core of group optimization and handles deep or
        greedy arguments

        :param kwargs: arguments handled by the function and config.
        :return: tuple of the dict electrified_stations and the set electrified_station_set
            in the optimal found case
        """

        node_choice = kwargs.get("node_choice", self.config.node_choice)
        opt_type = kwargs.get("opt_type", self.config.opt_type)

        self.scenario.vehicle_socs = self.timeseries_calc(electrify_stations=self.must_include_set)
        self.base_scenario = copy(self.scenario)
        self.base_schedule = copy(self.schedule)

        self.base_stations = self.electrified_stations.copy()
        self.base_not_possible_stations = self.not_possible_stations.copy()

        # get events where soc fell below the minimal
        # soc. The events contain info about the problematic
        # time span, which includes stations which could provide a soc lift
        base_events = self.get_low_soc_events(rel_soc=False,
                                              soc_data=self.base_scenario.vehicle_socs)

        # check if the events can be divided into subgroups which are independent
        # this makes optimization in smaller groups possible
        groups = \
            util.get_groups_from_events(base_events, self.not_possible_stations,
                                        could_not_be_electrified=self.could_not_be_electrified,
                                        optimizer=self)

        # storage for the sets of stations which will be generated
        list_greedy_sets = [set()] * len(groups)

        # create a decision tree or load one from a previous run.
        # there are as many trees as groups.
        self.set_up_decision_tree(len(groups))

        # baseline greedy optimization
        # base line is created simply by not having a decision tree and not a pre optimized_set yet
        util.print_time()
        for group_nr, group in enumerate(groups[:]):
            # unpack the group
            events, stations = group
            # define the current tree from the decision trees
            self.current_tree = self.decision_trees[group_nr]
            lines = {lne for e in events for lne in e.rotation.lines}
            self.logger.warning("Optimizing %s out of %s. This includes these Lines", group_nr + 1,
                                len(groups))
            self.logger.warning(lines)
            self.logger.warning("%s events with %s stations", (len(events)), len(stations))

            # first run is always step by step
            self.group_optimization(group, self.choose_station_step_by_step,
                                    events_remaining=[len(events)], **kwargs)

            self.logger.warning("Greedy Result ++++++++ %s stations out of %s",
                                len(self.electrified_station_set),
                                len(stations))
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
        util.print_time()
        # from here on only for deep analysis
        for group_nr, group in enumerate(groups[:]):
            events, stations = group
            sols = []
            self.current_tree = self.decision_trees[group_nr]
            if node_choice == "brute":
                choice_func = self.choose_station_brute
            else:
                choice_func = self.choose_station_step_by_step

            combinations = util.combs_unordered_no_putting_back(len(stations),
                                                                len(list_greedy_sets[group_nr]) - 1)
            self.logger.debug("There are %s combinations with 1 station less than "
                              "the current solution with %s stations out of %s", combinations,
                              len(list_greedy_sets[group_nr]), len(stations))

            pre_optimized_set = list_greedy_sets[group_nr]
            for i in range(self.config.max_brute_loop):
                i += 1
                if i % 10 == 0:
                    print(len(self.current_tree), " nodes checked")
                    print(f"Optimal solution has length {len(pre_optimized_set)}")

                self.electrified_stations = self.base_stations.copy()
                self.electrified_station_set = self.base_electrified_station_set.copy()
                self.scenario = copy(self.base_scenario)
                self.schedule = copy(self.base_schedule)

                self.scenario.vehicle_socs = deepcopy(self.base_scenario.vehicle_socs)
                # deep copy of schedule.rotations is very slow. Not needed for
                # quick calculation.
                if self.config.solver == "spiceev":
                    self.schedule.rotations = deepcopy(self.base_schedule.rotations)
                self.not_possible_stations = self.base_not_possible_stations.copy()

                # create a new electrified set
                try:
                    new_stats = \
                        self.group_optimization(group, choice_func,
                                                track_not_possible_rots=False,
                                                pre_optimized_set=pre_optimized_set,
                                                events_remaining=[len(events)],
                                                **kwargs)

                    # The terminal node is checked therefore not viable anymore
                    node_name = util.stations_hash(self.electrified_station_set)
                    self.current_tree[node_name]["viable"] = False
                    self.current_tree[node_name]["completed"] = True
                except util.SuboptimalSimulationException:
                    # exception if the optimization encountered a situation, where it is clear
                    # that going down the tree further, will not lead to an optimal solution
                    continue
                except util.AllCombinationsCheckedException:
                    # exception if all combinations in the brute force method are checked or at
                    # least estimated or all nodes show no potential in the step-by-step method
                    print("all checked")
                    break

                new_electrified_set = self.electrified_station_set
                # if a new set was found, print it and save it in sols
                if new_electrified_set != pre_optimized_set and new_stats is not None:
                    self.logger.warning("Optimized with %s  stations %s %s",
                                        len(new_electrified_set), str('#' * 20),
                                        util.stations_hash(new_electrified_set))
                    sols.append(new_electrified_set)
                    if len(new_electrified_set) < len(pre_optimized_set):
                        list_greedy_sets[group_nr] = new_electrified_set.copy()
                        pre_optimized_set = new_electrified_set.copy()
            else:
                print(f"Ran all {self.config.max_brute_loop} loops")

            # print all solutions with the length of the
            self.logger.debug("All solutions for this group: %s", sols)
            self.logger.debug("Optimized with %s stations out of %s", len(pre_optimized_set),
                              len(stations))

        util.print_time()
        # saving decision tree only in case of deep analysis
        if self.config.save_decision_tree:
            with open(self.args.output_directory / Path("decision_tree.pickle"), "wb") as file:
                pickle.dump(self.current_tree, file)

        self.electrified_stations = self.base_stations.copy()
        self.electrified_station_set = self.base_electrified_station_set.copy()
        for single_set in list_greedy_sets:
            for stat in single_set:
                self.electrify_station(stat, self.electrified_station_set)
        # dump the measured running times of the functions
        print(util.time_it(None))
        return self.electrified_stations, self.electrified_station_set

    def get_negative_rotations_all_electrified(self, rel_soc=False):
        """Get the rotation ids for the rotations which show negative socs, even when everything
        is electrified

        :param rel_soc: if true, the start soc is handled like it has the desired deps soc
        :return: set of rotation ids which are negative even with all stations electrified
        """
        events = self.get_low_soc_events(rel_soc=rel_soc)

        stats = {stat for event in events for stat in event.stations_list
                 if stat not in self.not_possible_stations}
        electrified_station_set = set(stats)
        vehicle_socs = self.timeseries_calc(ele_station_set=electrified_station_set)
        new_events = self.get_low_soc_events(soc_data=vehicle_socs)
        return {event.rotation.id for event in new_events}

    @util.time_it
    def group_optimization(self, group, choose_station_function, track_not_possible_rots=True,
                           pre_optimized_set=None, **kwargs):
        """Optimize a single group events and returns the electrified stations and a flag

        :param group: tuple of events and station_ids to be optimized together
        :param choose_station_function: function to be used for choosing the next function
        :param track_not_possible_rots: not possible tracks need to be tracked in the first run
        :param pre_optimized_set: set which was optimized before hand used as upper threshold for
            optimization
        :param kwargs: internally used functionality
        :return: Returns (electrified_stations if optimization continues, None if no further
            recursive calls should happen, bool if further deep analysis should take place)
        :rtype (dict() or None, bool)
        :raises AllCombinationsCheckedException: If all combinations have been checked
        :raises SuboptimalSimulationException: if a suboptimal Simulation has been identified
        """
        could_not_be_electrified = self.could_not_be_electrified

        # in deep analysis not all nodes are visited if they are suboptimal.
        # therefore could not be electrified is only important in the first run, where every station
        # is possible
        if track_not_possible_rots is None:
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
        station_eval = util.evaluate(event_group, self, soc_data=kwargs["lifted_socs"])

        missing_energy = util.get_missing_energy(event_group)
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
            node_name = util.stations_hash(self.electrified_station_set)
            self.current_tree[node_name]["viable"] = False
            if len(self.electrified_station_set) == 0:
                raise util.AllCombinationsCheckedException
            raise util.SuboptimalSimulationException

        # get the best stations, brute or step_by_step
        best_station_ids, recursive = choose_station_function(station_eval,
                                                              pre_optimized_set=pre_optimized_set,
                                                              missing_energy=missing_energy)

        stat_eval_dict = {stat_id[0]: stat_id[1] for stat_id in station_eval}

        self.logger.debug("%s, with first pot of %s", best_station_ids,
                          round(stat_eval_dict[best_station_ids[0]], 1))

        # electrify station(s)
        for stat_id in best_station_ids:
            self.electrify_station(stat_id, self.electrified_station_set)

        event_rotations = {event.rotation for event in event_group}

        # copy scen and sched
        self.copy_scen_sched()

        if solver == "quick":
            # quick calculation has to electrify everything in one step,
            # or the lifting of socs is not correct
            self.deepcopy_socs()
            self.scenario.vehicle_socs = self.timeseries_calc(event_rotations,
                                                              electrify_stations=best_station_ids)
        else:
            self.schedule.rotations = rotation_dict
            self.schedule, self.scenario = \
                util.run_schedule(self.schedule, self.args,
                                  electrified_stations=self.electrified_stations)

        kwargs["lifted_socs"] = self.scenario.vehicle_socs.copy()

        self.not_possible_stations = set(self.electrified_stations.keys()).union(
            self.not_possible_stations)
        event_rotations_ids = {event.rotation.id for event in event_group}
        new_events = self.get_low_soc_events(rotations=event_rotations_ids,
                                             rel_soc=True, **kwargs)

        delta_energy = util.get_missing_energy(new_events)
        events_remaining[0] -= len(event_group) - len(new_events)
        self.logger.debug("Last electrification electrified %s/%s."
                          " %s remaining events in the base group.",
                          len(event_group) - len(new_events), len(event_group), events_remaining[0])
        delta_base_energy = delta_energy

        # put this node into the decision tree including the missing energy
        self.node_to_tree(delta_base_energy)

        # everything electrified
        if delta_energy >= 0:
            return self.electrified_stations
        # some choice functions might not need a recursive call, they return here. recursive is set
        # by the choose_station_function
        if not recursive:
            raise util.SuboptimalSimulationException

        # check if the events can be divided into subgroups which are independent
        groups = util.get_groups_from_events(new_events, self.not_possible_stations,
                                             could_not_be_electrified,
                                             optimizer=self)

        for k, this_group in enumerate(groups):
            this_tree = tree_position + [k + 1]
            kwargs["tree_position"] = this_tree
            new_stations = \
                self.group_optimization(this_group, choose_station_function,
                                        track_not_possible_rots=track_not_possible_rots,
                                        pre_optimized_set=pre_optimized_set,
                                        **kwargs)

            if new_stations is not None:
                self.electrified_stations.update(new_stations)

            # if there is no pre optimized set function can return
            if pre_optimized_set is None:
                continue
            # if there is a pre optimized set check if the pre optimized set leads to pruning
            # i.e. if if following branch makes sense. Since the prognosis is not very accurate
            # over many station electrification pruning is very unlikely early on. The pruning
            # threshold should be set not to high, since checking for pruning is likely to be
            # unnecessary calculation
            thresh = self.config.pruning_threshold
            if len(pre_optimized_set) - len(self.electrified_station_set) < thresh:
                self.scenario.vehicle_socs = \
                    self.timeseries_calc(event_rotations, electrify_stations=best_station_ids)
                prune_events = self.get_low_soc_events(rotations=event_rotations_ids,
                                                       rel_soc=True, **kwargs)
                station_eval = util.evaluate(prune_events, self)
                prune_missing_energy = util.get_missing_energy(prune_events)
                if not self.is_branch_promising(station_eval, self.electrified_station_set,
                                                pre_optimized_set, prune_missing_energy):
                    self.logger.debug("Branch pruned early")
                    node_name = util.stations_hash(self.electrified_station_set)
                    self.current_tree[node_name]["viable"] = False
                    raise util.SuboptimalSimulationException

        return self.electrified_stations

    @util.time_it
    def deepcopy_socs(self):
        self.scenario.vehicle_socs = deepcopy(self.base_scenario.vehicle_socs)

    @util.time_it
    def copy_scen_sched(self):
        self.scenario = copy(self.base_scenario)
        self.schedule = copy(self.base_schedule)

# ToDo implement a fast calculation of timeseries_calc but with a limited amount of charging
#  points. Implementation could look like this

# def get_charge_events_per_station(self, station_name, rotations=None):
#     """ Gather low soc events below the config threshold.
#
#     :param rotations: rotations to be searched for low soc events. Default None means whole
#         schedule is searched
#     :param station_name: name of the station to be checked
#     :return: list(ChargingEvents)
#     """
#     if not rotations:
#         rotations = self.schedule.rotations
#
#     # Find relevant charging events for this station
#     charging_events = []
#
#     for rot_id in rotations:
#         rot = self.schedule.rotations[rot_id]
#         for i, trip in enumerate(rot.trips):
#             # if the trip station name is not the searched station continue.
#             if trip.arrival_name != station_name:
#                 continue
#             try:
#                 charging_event_start_time = util.get_charging_start(trip, self.args)
#                 end_time = rot.trips[i + 1].departure_time
#
#
#                 # Do not add the event if there is charging time of at least the defined
#                 # min charging time
#                 if charging_event_start_time + self.args.min_charging_time >= end_time:
#                     continue
#             except IndexError:
#                 warnings.warn("Station to be checked has no following trip. Final destinations"
#                               "can not offer lift to the soc and are therefore discarded")
#                 continue
#             arrival_time= trip.arrival_time
#             start_idx = self.get_index_by_time(charging_event_start_time)
#             end_time = rot.trips[i+1].departure_time
#             buffer_time = util.get_buffer_time(trip, self.args.default_buffer_time_opps)
#             end_idx = self.get_index_by_time(end_time)
#             cht = rot.vehicle_id.find("depb")
#             ch_type = (cht > 0) * "depb" + (cht <= 0) * "oppb"
#             v_type = rot.vehicle_id.split("_" + ch_type)[0]
#             event = util.ChargingEvent(start_idx=start_idx,end_idx=end_idx,
#                                        arrival_time=arrival_time,
#                                        start_time=charging_event_start_time, end_time=end_time,
#                                        buffer_time=buffer_time, vehicle_id=rot.vehicle_id,
#                                        capacity=self.schedule.vehicle_types[v_type][ch_type][
#                                          'capacity'],
#                                        station_name=station_name, rotation=rot)
#             charging_events.append(event)
#     return charging_events

    def sort_station_events(self, charge_events_single_station):
        return sorted(charge_events_single_station, key=lambda x: x.arrival_time)

    def mutate_events_for_n_charging_points(self, charge_events_single_station, nr_charge_points):
        pass

    # mutate_events_for_n_charging_points()
    #
    # sort_all_station_events()
    #
    # calculate_events()
    #
    # get_below_0_soc_events()
    #
    # get_missing_energy()
    #
    # evaluate_soc_lift_vs_objective_function()

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
    # or 1 charging point
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

    @util.time_it
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
        :return: Returns soc dict with lifted socs
        :rtype dict()
        """
        if rotations is None:
            rotations = self.schedule.rotations.values()

        if soc_dict is None:
            soc_dict = self.scenario.vehicle_socs
        if ele_station_set is None:
            ele_station_set = self.electrified_station_set
        if not soc_upper_thresh:
            soc_upper_thresh = self.args.desired_soc_deps
        if electrify_stations is None:
            electrify_stations = set()

        if electrify_stations is None:
            electrify_stations = set()
        ele_stations = {*ele_station_set, *electrify_stations}
        soc_dict = deepcopy(soc_dict)

        for rot in rotations:
            ch_type = (rot.vehicle_id.find("oppb") > 0) * "oppb" + (
                    rot.vehicle_id.find("depb") > 0) * "depb"
            v_type = rot.vehicle_id.split("_" + ch_type)[0]
            soc_over_time_curve = self.soc_charge_curve_dict[v_type][ch_type]
            soc = soc_dict[rot.vehicle_id]
            for i, trip in enumerate(rot.trips):
                if trip.arrival_name not in ele_stations:
                    continue
                idx = util.get_index_by_time(self.scenario, trip.arrival_time)
                try:
                    standing_time_min = util.get_charging_time(trip, rot.trips[i + 1], self.args)
                except IndexError:
                    standing_time_min = 0

                d_soc = util.get_delta_soc(soc_over_time_curve, soc[idx], standing_time_min, self)
                buffer_idx = int((util.get_buffer_time(trip, self.args.default_buffer_time_opps))
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
                    # i.e. after lifting the soc, it finds the first spot where the soc is bigger
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

    @util.time_it
    def expand_tree(self, station_eval):
        try:
            parent_name = util.stations_hash(self.electrified_station_set)
            self.current_tree[parent_name]
        except KeyError:
            self.current_tree[parent_name] = get_init_node()
        for stat, ev_score in station_eval:
            node_name = util.stations_hash(self.electrified_station_set.union([stat]))
            self.current_tree[parent_name]["children"].add(node_name)
            try:
                self.current_tree[node_name]
            except KeyError:
                self.current_tree[node_name] = get_init_node()

    @util.time_it
    def is_node_viable(self):
        """Check if a node has viable children. Viable children are not terminal and did not show
        lack of potential yet
        :return: True if the node still has potential
        """
        parent_name = util.stations_hash(self.electrified_station_set)
        for child in self.current_tree[parent_name]["children"]:
            if not self.current_tree[child]["viable"]:
                continue
            # if one viable node is found  return True
            return True
        # not a single viable node was found, return False
        return False

    @util.time_it
    def is_branch_promising(self, station_eval, electrified_station_set,
                            pre_optimized_set, missing_energy):
        """Quickly evaluates if following a branch is promising by summing up estimated potentials
        :param station_eval: list of  sorted station evaluation with list(station_id, pot)
        :param electrified_station_set: set of electrified stations
        :param pre_optimized_set: set of already optimized stations
        :param missing_energy: the amount of energy missing in this branch
        :return: is the branch promising
        :rtype bool
        """

        delta = len(pre_optimized_set) - len(electrified_station_set)
        pot = 0
        # sum up highest potentials until you would reach the number of electrified stations
        # in the pre optimized case
        for i in range(0, min(delta, len(station_eval))):
            pot += station_eval[i][1]
        if pot < -missing_energy * self.config.estimation_threshold:
            self.logger.debug("Not enough potential: %s  after: %s ",
                              round(pot / -missing_energy, 0),
                              util.stations_hash(electrified_station_set))
            return False
        return True

    def node_to_tree(self, delta_base_energy):
        """ Fill decision tree with the given info
        :param delta_base_energy: missing energy
        :return decision tree
        :rtype dict()
        """
        node_name = util.stations_hash(self.electrified_station_set)
        self.current_tree[node_name]["missing_energy"] = delta_base_energy
        self.current_tree[node_name]["visit_counter"] += 1
        if self.current_tree[node_name]["visit_counter"] > 1:
            self.logger.debug("already visited this node")

    @util.time_it
    def choose_station_brute(self, station_eval,
                             pre_optimized_set=None, missing_energy=0,
                             gens=dict()):
        """Gives back a possible set of Stations to electrify which shows potential and has not been
        tried yet. The set of stations is smaller than the best optimized set so far.

        :param station_eval: list of  sorted station evaluation with list(station_id, pot)
        :param pre_optimized_set: set of optimized stations thus far
        :param missing_energy: missing energy in this branch before electrification
        :param gens: dict of generators for brute force generation
        :return: combination of stations or None if no viable stations exists and
            false since this function does not support recursive calling
        :raises AllCombinationsCheckedException: If all combinations have been checked
        """

        station_ids = [x[0] for x in station_eval]
        try:
            generator = gens[str(station_ids) + str(len(pre_optimized_set) - 1)]
        except KeyError:
            generator = util.combination_generator(station_ids, len(pre_optimized_set) - 1)
            gens[str(station_ids) + str(len(pre_optimized_set) - 1)] = generator
        station_eval_dict = {stat[0]: stat[1] for stat in station_eval}
        for comb in generator:
            node_name = util.stations_hash(comb)
            if node_name not in self.current_tree:
                # only check the brute force station if they have the remote chance of fulfilling
                # the missing energy
                # potential>missing energy * 80%
                potential = sum([station_eval_dict[stat] for stat in comb])
                if potential > -missing_energy * self.config.estimation_threshold:
                    self.current_tree[node_name] = get_init_node()
                    return comb, False
                else:
                    self.logger.debug("skipped %s since potential is too low %s %%", comb,
                                      round(potential / -missing_energy * 100, 0))
                    try:
                        self.current_tree[node_name]["viable"] = False
                    except KeyError:
                        self.current_tree[node_name] = get_init_node()
                        self.current_tree[node_name]["viable"] = False
        self.logger.debug("calculated all viable possibilities")
        raise util.AllCombinationsCheckedException

    @util.time_it
    def choose_station_step_by_step(self, station_eval,
                                    pre_optimized_set=None, missing_energy=0):
        """ Gives back a station of possible stations to electrify which shows the biggest potential
        and has not been picked before.


        :param station_eval: list of  sorted station evaluation with list(station_id, pot)
        :param pre_optimized_set: set of optimized stations thus far
        :param missing_energy: missing energy in this branch before electrification
        :return:a station to electrify or None if no viable stations exists and
            false since this function does not support recursive calling
        :raises SuboptimalSimulationException: if a suboptimal Simulation has been identified
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
                # best station id is none and do not go deeper in recursion
                node_name = util.stations_hash(self.electrified_station_set)
                self.current_tree[node_name]["viable"] = False
                raise util.SuboptimalSimulationException

        min_nr_visited = float('inf')
        # find the lowes amount of visits in the possible stations /  children
        # not viable stations are skipped
        for stat_tuple in station_eval:
            station, _ = stat_tuple
            # create a station combination from already electrified
            # stations and possible new stations
            stats = self.electrified_station_set.union([station])
            node_name = util.stations_hash(stats)
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
            node_name = util.stations_hash(stats)
            if not self.current_tree[node_name]["viable"]:
                continue

            if self.current_tree[util.stations_hash(stats)]["visit_counter"] == min_nr_visited:
                best_station_id = station
                return [best_station_id], True

        node_name = util.stations_hash(self.electrified_station_set)
        self.current_tree[node_name]["viable"] = False
        raise util.SuboptimalSimulationException

    def set_battery_and_charging_curves(self):
        """ Set battery and charging curves from config
        """
        for v_type in self.schedule.vehicle_types.values():
            for vehicle in v_type.values():
                if self.config.battery_capacity is not None:
                    vehicle["capacity"] = self.config.battery_capacity
                if self.config.charging_curve is not None:
                    vehicle["charging_curve"] = self.config.charging_curve

    def set_up_decision_tree(self, group_amount):
        """ Load decision tree if given in the config
        :param group_amount: amount of groups for decision tree
        """
        if self.config.decision_tree_path is not None:
            with open(self.config.decision_tree_path, "rb") as file:
                self.decision_trees = pickle.load(file)
        else:
            self.decision_trees = [{} for _ in range(group_amount)]

    def rebase_spice_ev(self):
        """ Rebase the scenario meaning configuring various variables according to the input data
        and running a spice ev simulation
        :return: must_include_set and dict of electrified_stations
        """
        self.logger.debug("Spice EV Rebasing Scenario")
        must_include_set, ele_stations = self.preprocessing_scenario(
            run_only_neg=self.config.run_only_neg)
        self.logger.debug("Rebasing finished")
        if self.config.pickle_rebased:
            util.toolbox_to_pickle(self.config.pickle_rebased_name, self.schedule, self.scenario,
                                   self.args)
            self.logger.debug("Rebased scenario pickled  as %s", self.config.pickle_rebased_name)

        self.must_include_set = must_include_set
        self.base_stations = self.electrified_stations.copy()
        return must_include_set, ele_stations

    def rebase_simple(self):
        """ Rebase the scenario meaning configuring various variables according to the input data
        :return: must_include_set and dict of electrified_stations
        """
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
        """Prepare scenario and run schedule
        :param electrified_stations: optional dict of electrified stations to use in simulation.
            Default None leads to using optimizer.electrified_stations
        :param run_only_neg: should only negative rotations be simulated
        :param cost_calc: should costs be calculated
        :return: schedule, scenario, electrified_station_set, electrified_stations
        :rtype (schedule.Schedule, scenario.Scenario, set(), dict())
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
            rots = {r: self.schedule.rotations[r] for r in rots if r not in
                    self.config.exclusion_rots}
            self.schedule.rotations = rots
        else:
            rots = {r: self.schedule.rotations[r] for r in self.schedule.rotations if
                    r not in self.config.exclusion_rots}
            self.schedule.rotations = rots

        new_sched, new_scen = util.run_schedule(self.schedule, self.args,
                                                electrified_stations,
                                                cost_calc=cost_calc)
        try:
            report.generate(new_sched, new_scen, self.args)
        except Exception:
            warnings.warn('Report generation failed: " {0}'.format(traceback.format_exc()))
        self.schedule = new_sched
        self.scenario = new_scen
        self.must_include_set = must_include_set
        return must_include_set, electrified_stations

    def electrify_station(self, stat, electrified_set):
        """Electrify a station and keep track of it in the electrified set file
        :param stat: station id to be electrified
        :param electrified_set:the set that is mutated along the electrification
        """
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
                soc_charge_curve_dict[name][ch_type] = util.charging_curve_to_soc_over_time(
                    data["charging_curve"], data["capacity"], self.args,
                    self.schedule.cs_power_opps, efficiency=self.config.charge_eff,
                    time_step=0.1)
        self.soc_charge_curve_dict = soc_charge_curve_dict

    @util.time_it
    def get_must_stations_and_rebase(self, relative_soc=False):
        """ Get the stations that must be electrified and put them into the electrified stations
        and an extra set of must stations
        Electrify everything minus 1 station. If without the stations there are below zero events
        it is a must have station

        :param relative_soc: should the evaluation use the relative or absolute soc
        :return: Set(Station_ids)
        """

        events = self.get_low_soc_events(rel_soc=relative_soc)

        stats = {station for event in events for station in event.stations_list
                 if station not in self.not_possible_stations}

        electrified_station_set_all = set(stats)
        must_stations = set()
        print(f"Electrifying {len(electrified_station_set_all)} stations minus one", )
        for station in sorted(electrified_station_set_all):
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
            # do not put must stations in electrified set, but in extra set must_include_set
            self.electrify_station(stat, self.must_include_set)
        return must_stations

    def remove_none_socs(self):
        """ Removes soc values of None by filling them with the last value which is not None
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
        """ Gets you the soc object with start and end index for a given rotation id
        :param rot_id: rotation_id
        :param soc_data: optional soc_data if not the scenario data should be used
        :return: tuple with soc array, start index and end index
        """
        return util.get_rotation_soc_util(rot_id, self.schedule, self.scenario, soc_data=soc_data)

    def get_index_by_time(self, search_time: datetime):
        """Get the index for a given time.
        :param search_time: The time for which to return the index as datetime object
        :return: the time corresponding to the given index.
        """
        start_time = self.scenario.start_time
        delta_time = timedelta(minutes=60 / self.scenario.stepsPerHour)
        idx = (search_time - start_time) // delta_time
        return idx

    @util.time_it
    def get_trips(self, rot: rotation.Rotation, start_idx: int, end_idx: int):
        """ Get trips in a rotation from a start to an end index, if the arrival time is in between
         the start and end idx

        :param rot: The rotation object containing the trips.
        :param start_idx: The start index, representing the start time.
        :param end_idx:self The end index, representing the end time.
        :return: A list of trip objects that arrive between the start and end time.
        """

        start_time_event = self.get_time_by_index(start_idx)
        end_time_event = self.get_time_by_index(end_idx)

        trips = []
        for i, trip in enumerate(rot.trips):
            if end_time_event >= trip.arrival_time > start_time_event:
                trips.append(trip)

        return trips

    @util.time_it
    def get_time_by_index(self, idx):
        """Get the time for a given index.

        :param idx: The index for which to return the time.
        :return: the time corresponding to the given index.
        """

        start_time = self.scenario.start_time
        delta_time = timedelta(minutes=60 / self.scenario.stepsPerHour)
        searched_time = start_time + delta_time * idx
        return searched_time

    @util.time_it
    def get_low_soc_events(self, rotations=None, filter_standing_time=True,
                           rel_soc=False, soc_data=None, **kwargs):
        """ Gather low soc events below the config threshold.

        :param rotations: rotations to be searched for low soc events. Default None means whole
            schedule is searched
        :param filter_standing_time: Should the stations be filtered by standing time. True leads to
            an output only stations with charging potential
        :param rel_soc: Defines if the start soc should be seen as full even when not.
            i.e. a drop from a rotation start soc from 0.1 to -0.3 is not critical since the start
            soc from the rotation will be raised
            If the rel_soc
            is false, it means coupled rotations might be prone to errors due to impossible lifts.
        :param soc_data: soc data to be used. Default None means the soc data from the optimizer
            scenario is used
        :param kwargs: optional soc_lower_thresh or soc_upper_thresh if from optimizer differing
            values should be used
        :return: list(LowSocEvents)
        """
        if not rotations:
            rotations = self.schedule.rotations

        soc_lower_thresh = kwargs.get("soc_lower_thresh", self.config.min_soc)
        soc_upper_thresh = kwargs.get("soc_upper_thresh", self.args.desired_soc_deps)
        # Create list of events which describe trips which end in a soc below zero
        # The event is bound by the lowest soc and an upper soc threshold which is naturally 1
        # Properties before and after these points have no effect on the event itself, similar to
        # an event horizon
        events = []
        count_electrified_rot = 0

        for rot_id in rotations:
            rot = self.schedule.rotations[rot_id]
            soc, rot_start_idx, rot_end_idx = self.get_rotation_soc(rot_id, soc_data)
            rot_end_idx += 1
            idx = range(0, len(soc))
            comb = list(zip(soc, idx))[rot_start_idx:rot_end_idx]
            min_soc, min_idx = min(comb, key=lambda x: x[0])
            reduced_list = comb.copy()
            soc_lower_thresh_cur = soc_lower_thresh
            # if rotation gets a start soc below 1 this should change below 0 soc
            # events since fixing the rotation before would lead to fixing this rotation

            # if using relative SOC, SOC lookup has to be adjusted
            if rel_soc:
                start_soc = comb[0][0]
                soc_lower_thresh_cur = min(start_soc, soc_upper_thresh) - (
                        soc_upper_thresh - soc_lower_thresh)
                soc_upper_thresh = soc_lower_thresh_cur + soc_upper_thresh
            if min_soc >= soc_lower_thresh_cur:
                count_electrified_rot += 1
            while min_soc < soc_lower_thresh_cur:
                i = min_idx
                idx = [x[1] for x in reduced_list]
                while soc[i] < soc_upper_thresh:
                    if i == rot_start_idx:
                        break
                    i -= 1
                start_comb = idx.index(i)
                start = i
                i = min_idx
                while soc[i] < soc_upper_thresh:
                    if i >= rot_end_idx-1:
                        break
                    i += 1
                end_comb = idx.index(i)
                trips = self.get_trips(rot=rot, start_idx=start, end_idx=min_idx)
                possible_stations = set()
                possible_stations_list = []
                if not filter_standing_time:
                    possible_stations = {t.arrival_name for t in trips}
                    possible_stations_list = [t.arrival_name for t in trips]
                else:
                    for ii, trip in enumerate(trips):
                        try:
                            standing_time_min = util.get_charging_time(trip, trips[ii + 1],
                                                                       self.args)
                        except IndexError:
                            standing_time_min = 0
                        if standing_time_min > 0:
                            possible_stations.add(trip.arrival_name)
                            possible_stations_list.append(trip.arrival_name)

                possible_stations = possible_stations.difference(self.not_possible_stations)
                cht = rot.vehicle_id.find("depb")
                ch_type = (cht > 0) * "depb" + (cht <= 0) * "oppb"
                v_type = rot.vehicle_id.split("_" + ch_type)[0]
                event = util.LowSocEvent(start_idx=start, end_idx=min_idx,
                                         min_soc=min_soc, stations=possible_stations,
                                         vehicle_id=rot.vehicle_id, trip=trips,
                                         rot=rot, stations_list=possible_stations_list,
                                         capacity=self.schedule.vehicle_types[v_type][ch_type][
                                             'capacity'],
                                         v_type=v_type, ch_type=ch_type)

                events.append(event)
                copy_list = reduced_list.copy()
                reduced_list = reduced_list[:start_comb]
                if end_comb + 1 <= len(copy_list):
                    reduced_list.extend(copy_list[end_comb + 1:])
                if len(reduced_list) > 0:
                    min_soc, min_idx = min(reduced_list, key=lambda x: x[0])
                else:
                    break
        return events


def get_init_node():
    """ Returns the initialization dictionary / node of the decision tree
    :return: init node
    :rtype: dict
    """
    return {"completed": False, "viable": True,
            "missing_energy": None, "visit_counter": 0, "children": set()}