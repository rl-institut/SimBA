""" Module for LowSocEvent, ChargingEvent, OptimizerConfig and utility functionality """
import logging
import math
import os
from pathlib import Path
import pickle
import sys
import typing
import warnings
import configparser
import json
from copy import copy
from datetime import timedelta, datetime
from time import time
import numpy as np
from matplotlib import pyplot as plt

if typing.TYPE_CHECKING:
    from simba.station_optimizer import StationOptimizer

from simba.util import get_buffer_time as get_buffer_time_util
from spice_ev.report import generate_soc_timeseries


class ChargingEvent:
    """ Class to gather information about a charging event """

    def __init__(self, start_idx, end_idx, arrival_time, start_time, end_time, buffer_time,
                 vehicle_id, capacity,
                 station_name, rotation):
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.arrival_time = arrival_time
        self.start_time = start_time
        self.end_time = end_time
        self.buffer_time = buffer_time
        self.vehicle_id = vehicle_id
        self.capacity = capacity
        self.station_name = station_name
        self.rotation = rotation


class LowSocEvent:
    """ Class to gather information about a low soc event """
    event_counter = 0

    def __init__(self, start_idx, end_idx, min_soc, stations, vehicle_id, trip, rot,
                 stations_list, capacity, v_type, ch_type):
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.min_soc = min_soc
        self.stations = stations
        self.vehicle_id = vehicle_id
        self.trip = trip
        self.rotation = rot
        self.stations_list = stations_list
        self.capacity = capacity
        self.v_type = v_type
        self.ch_type = ch_type
        self.event_counter = LowSocEvent.event_counter
        LowSocEvent.event_counter += 1


class OptimizerConfig:
    """ Class for the configuration file """

    def __init__(self):
        self.logger_name = ""
        self.debug_level = 0
        self.console_level = 99

        self.exclusion_rots = []
        self.exclusion_stations = set()
        self.inclusion_stations = set()
        self.standard_opp_station = {"type": "opps", "n_charging_stations": None}

        self.schedule = None
        self.scenario = None
        self.args = None
        self.charge_eff = None
        self.battery_capacity = None
        self.charging_curve = None
        self.charging_power = None
        self.min_soc = 0

        self.solver = "spiceev"
        self.rebase_scenario = False
        self.pickle_rebased = False
        # used for gradual scenario analysis in django-simba
        self.early_return = False

        self.pickle_rebased_name = None
        self.opt_type = "greedy"
        self.eps = 0.0001

        self.remove_impossible_rotations = False
        self.node_choice = "step-by-step"
        self.max_brute_loop = 20
        self.run_only_neg = True
        self.run_only_oppb = True
        self.estimation_threshold = 0.8
        self.check_for_must_stations = False
        self.pruning_threshold = 3

        self.decision_tree_path = None
        self.save_decision_tree = False
        self.optimizer_output_dir = None
        self.reduce_rotations = False
        self.rotations = None
        self.path = None
        self.save_all_results = None


def time_it(function, timers={}):
    """ Decorator function to time the duration and number of function calls.

    :param function: function do be decorated
    :type function: function
    :param timers: storage for cumulated time and call number
    :type timers: dict
    :return: decorated function or timer if given function is None
    :rtype: function or dict

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


def read_config(config_path):
    """ Read the config path to a config object.

    :param config_path: path to file
    :type config_path: str
    :return: configuration
    :rtype: simba.optimizer_util.OptimizerConfig
    """
    config_parser = configparser.ConfigParser()

    assert Path(config_path).is_file(), (f"Path to optimizer_config: {config_path} "
                                         "does not lead to file")
    try:
        config_parser.read(config_path, encoding="utf-8")
    except configparser.MissingSectionHeaderError:
        # make sure there is always a DEFAULT section.
        with open(config_path, 'r', encoding='utf-8') as f:
            config_string = '[DEFAULT]\n' + f.read()
        config_parser.read_string(config_string)
    conf = OptimizerConfig()
    conf.path = config_path
    assert len(config_parser.sections()) == 0, "Sections are not allowed for the optimizer config."

    default = config_parser["DEFAULT"]
    conf.debug_level = default.getint("debug_level", 0)
    conf.console_level = default.getint("console_level", 99)

    conf.exclusion_rots = set(json.loads(default.get("exclusion_rots", "[]")))
    conf.exclusion_stations = set(json.loads(default.get("exclusion_stations", "[]")))
    conf.inclusion_stations = set(json.loads(default.get("inclusion_stations", "[]")))
    conf.standard_opp_station = dict(json.loads(default.get("standard_opp_station", "{}")))

    conf.schedule = default.get("schedule", "")
    conf.scenario = default.get("scenario", "")
    conf.args = default.get("args", "")

    conf.battery_capacity = default.getfloat("battery_capacity", 0)
    if conf.battery_capacity == 0:
        conf.battery_capacity = None
    conf.charging_curve = json.loads(default.get("charging_curve", "[]"))
    if not conf.charging_curve:
        conf.charging_curve = None
    conf.charging_power = default.getfloat("charging_power", 0)
    if conf.charging_power == 0:
        conf.charging_power = None
    conf.min_soc = default.getfloat("min_soc", 0)

    conf.solver = default.get("solver", "spiceev")
    conf.rebase_scenario = default.getboolean("rebase_scenario", True)
    conf.pickle_rebased = default.getboolean("pickle_rebased", False)
    conf.pickle_rebased_name = default.get(
        "pickle_rebased_name", "rebased_" + datetime.now().isoformat(sep='-', timespec='seconds'))
    conf.opt_type = default.get("opt_type", "greedy")
    conf.eps = default.getfloat("eps", 0.0001)
    conf.remove_impossible_rotations = default.getboolean("remove_impossible_rotations", False)
    conf.node_choice = default.get("node_choice", "step-by-step")
    conf.max_brute_loop = default.getint("max_brute_loop", 20)
    conf.run_only_neg = default.getboolean("run_only_neg", False)
    conf.run_only_oppb = default.getboolean("run_only_oppb", False)
    conf.estimation_threshold = default.getfloat("estimation_threshold", 0.8)
    conf.check_for_must_stations = default.getboolean("check_for_must_stations", True)
    conf.pruning_threshold = default.getint("pruning_threshold", 3)

    conf.decision_tree_path = default.get("decision_tree_path", None)
    if conf.decision_tree_path in ["", '""', "''"]:
        # path should not be empty or "empty like"
        conf.decision_tree_path = None
    conf.save_decision_tree = default.getboolean("save_decision_tree", False)
    conf.reduce_rotations = default.getboolean("reduce_rotations", False)
    conf.rotations = json.loads(default.get("rotations", "[]"))

    return conf


def get_charging_time(trip1, trip2, args):
    """ Returns the charging time in minutes between trips as numeric value/float.

    :param trip1: first trip
    :type trip1: simba.trip.Trip
    :param trip2: following trip
    :type trip2: simba.trip.Trip
    :param args:  arguments including default buffer time
    :type args: Namespace
    :return: maximum possible charging time in minutes between trips
    :rtype: float
    """
    standing_time_min = (trip2.departure_time - trip1.arrival_time) / timedelta(minutes=1)
    buffer_time = (get_buffer_time(trip1, args.default_buffer_time_opps) / timedelta(minutes=1))
    standing_time_min -= buffer_time

    if args.min_charging_time > standing_time_min:
        return 0
    return max(0, standing_time_min)


def get_charging_start(trip1, args):
    """ Return the possible start time of charging.

    This function considers the buffer times before charging can take place

    :param trip1: trip to be checked
    :type trip1: simba.trip.Trip
    :param args: arguments including default buffer time
    :type args: Namespace
    :return: First possible charging time as datetime object
    """
    buffer_time = get_buffer_time(trip1, args.default_buffer_time_opps)
    return trip1.arrival_time+buffer_time


def get_buffer_time(trip, default_buffer_time_opps):
    """ Return the buffer time as timedelta object.

    :param trip: trip to be checked
    :type trip: simba.trip.Trip
    :param default_buffer_time_opps: the default buffer time at opps charging stations
    :return: buffer time
    :rtype: datetime.timedelta
    """
    return timedelta(minutes=get_buffer_time_util(trip, default_buffer_time_opps))


def get_index_by_time(scenario, search_time):
    """ Get the index for a given time.

    In case the time does not coincide with a simulation index the lower index is returned.

    :param scenario: scenario to be checked
    :type scenario: spice_ev.Scenario
    :param search_time: search time
    :type search_time: datetime.datetime
    :return: index
    :rtype: int
    """
    return (search_time - scenario.start_time) // scenario.interval


def get_rotation_soc(rot_id, schedule, scenario, soc_data: dict = None):
    """ Return the SoC time series with start and end index for a given rotation ID.

    :param rot_id: rotation_id
    :type rot_id: str
    :param schedule: schedule containing rotation information
    :type schedule: simba.schedule.Schedule
    :param scenario: scenario containing the soc data
    :type scenario: spice_ev.Scenario
    :param soc_data: optional soc_data if the scenario data should not be used
    :type soc_data: dict
    :return: soc time series, start index and end index
    :rtype: (list, int, int)
    """
    rot = schedule.rotations[rot_id]
    rot_start_idx = get_index_by_time(scenario, rot.departure_time)
    rot_end_idx = get_index_by_time(scenario, rot.arrival_time)
    if soc_data:
        return soc_data[rot.vehicle_id], rot_start_idx, rot_end_idx
    return scenario.vehicle_socs[rot.vehicle_id], rot_start_idx, rot_end_idx


def get_delta_soc(soc_over_time_curve, soc, duration_min: float):
    """ Return expected SoC lift for a given SoC charging time series, start_soc and time_delta.

    :param soc_over_time_curve: Data with columns: time, soc and n rows
    :type soc_over_time_curve: np.array() with shape(n, 2)
    :param soc: start socs
    :type soc: float
    :param duration_min: duration of charging in minutes
    :type duration_min: float
    :return: positive delta of the soc
    :rtype: float
    """
    # units for time_delta and time_curve are assumed to be the same, e.g. minutes
    # first element which is bigger than current soc

    if duration_min <= 0:
        return 0
    negative_soc = 0

    # find time until vehicle soc is no longer negative, since the charging curve only defines
    # socs between 0 and 1.
    # when charging a negative soc, the charging power in the negative region is considered to
    # be constant, with the power value of the soc at 0.
    if soc < 0:
        # keep track of the negative soc to add it later
        negative_soc = -soc
        try:
            gradient = soc_over_time_curve[1, 1] - soc_over_time_curve[0, 1]
        except IndexError:
            # If no charge curve exists, i.e. less than 2 elements
            return 0
        charge_duration_till_0 = -soc / gradient
        if charge_duration_till_0 > duration_min:
            return duration_min * gradient
        duration_min = duration_min - charge_duration_till_0
        soc = 0

    idx = np.searchsorted(soc_over_time_curve[:, 1], soc, side='left')
    first_time, start_soc = soc_over_time_curve[idx, :]
    second_time = first_time + duration_min
    # catch out of bounds if time of charging end is bigger than table values

    if second_time >= soc_over_time_curve[-1, 0]:
        end_soc = soc_over_time_curve[-1, 1]
        # this makes sure the battery is actually at 100%
        start_soc = soc
    else:
        end_soc = soc_over_time_curve[soc_over_time_curve[:, 0] >= second_time][0, 1]

    return end_soc - start_soc + negative_soc


class InfiniteLoopException(Exception):
    pass


class SuboptimalSimulationException(Exception):
    pass


class AllCombinationsCheckedException(Exception):
    pass


@time_it
def evaluate(events: typing.Iterable[LowSocEvent],
             optimizer: 'StationOptimizer', **kwargs):
    """ Analyse stations for useful energy supply.

    Energy supply is helpful if the minimal soc of an event is raised (up to a minimal soc
    (probably zero)). The supplied energy is approximated by  charging power, standing time at a
    station, soc at station and minimal soc of the event.

    :param events: events to be evaluated
    :type events: list(simba.optimizer_util.LowSocEvent)
    :param optimizer: optimizer with scenario and schedule data
    :type optimizer: simba.station_optimizer.StationOptimizer
    :param kwargs: optional overwriting of soc_lower_thresh, soc_upper_thresh or soc_data
    :return: sorted stations and potentials
    :rtype: list(str(station_id), float(potential))

    """
    soc_lower_thresh = kwargs.get("soc_lower_thresh", optimizer.config.min_soc)
    soc_upper_thresh = kwargs.get("soc_upper_thresh", optimizer.args.desired_soc_deps)
    soc_data = kwargs.get("soc_data", optimizer.scenario.vehicle_socs)

    station_eval = {}
    # Note: Lift describes the positive delta in the soc time series through electrification.
    # cycle through events and determine how much lift can be provided by electrifying a station
    # the lift is determined by the soc position, standing time, power supply and charging curve
    for e in events:
        soc_over_time = optimizer.soc_charge_curve_dict[e.v_type][e.ch_type]
        for i, trip in enumerate(e.trip):
            # station is only evaluated if station name is part of event stations
            # only these stations showed potential in electrification, e.g. enough standing time
            if trip.arrival_name not in e.stations:
                continue
            idx = get_index_by_time(optimizer.scenario, trip.arrival_time, )
            soc = soc_data[e.vehicle_id][idx]

            # potential is the minimal amount of the following boundaries
            # - soc can only be lifted to the upper threshold
            # - useful lift is in between the minimal soc of the event and the
            #      lower threshold of soc
            # - useful lift can only occur in between the current soc and the
            # - the highest useful lift is the amount between a "full" and "empty" battery, where
            # "full" and "empty" are described by the upper and lower thresholds
            delta_soc_pot = min(
                soc_upper_thresh - soc,
                soc_lower_thresh - e.min_soc,
                soc - e.min_soc,
                soc_upper_thresh - soc_lower_thresh)

            try:
                standing_time_min = get_charging_time(trip, e.trip[i + 1], optimizer.args)
            except IndexError:
                standing_time_min = 0

            desired_soc = optimizer.args.desired_soc_opps
            soc = max(soc, 0)
            d_soc = get_delta_soc(soc_over_time, soc, standing_time_min)
            pot_kwh = min(d_soc, desired_soc) * e.capacity

            # potential is at max the minimum between the useful delta soc * capacity or the
            # energy provided by charging for the full standing time
            delta_e_pot = min(delta_soc_pot * e.capacity, pot_kwh)
            try:
                station_eval[trip.arrival_name] += delta_e_pot
            except KeyError:
                station_eval[trip.arrival_name] = delta_e_pot

    # sort by pot_sum
    station_eval_list = sorted(station_eval.items(), key=lambda x: x[1], reverse=True)
    return station_eval_list


def get_groups_from_events(events, impossible_stations=None, could_not_be_electrified=None,
                           optimizer=None):
    """ Create groups from events which need to be optimized together.

    First it creates a simple list of station sets for single events.
    They are connected if they share possible stations.
    Electrified and non-electrifiable stations are ignored,
    i.e. will not show up as a possible station.

    :param events: events for a given state of a scenario
    :type events: list(simba.optimizer_util.LowSocEvent)
    :param impossible_stations: stations to be discarded
    :type impossible_stations: set(str)
    :param could_not_be_electrified: rotations to be discarded
    :type could_not_be_electrified: set(str)
    :param optimizer: Optimizer
    :type optimizer: simba.station_optimizer.StationOptimizer
    :return: Groups of events and stations which have to be optimized together
    :rtype: list((simba.optimizer_util.LowSocEvent, str))
    """

    # making sure default arguments are none and not mutable
    if impossible_stations is None:
        impossible_stations = set()

    if could_not_be_electrified is None:
        could_not_be_electrified = set()

    possible_stations = [
        {station for station in event.stations if station not in impossible_stations}
        for event in events]
    # if station sets have intersections they are returned as single set
    station_subsets = join_all_subsets(possible_stations)
    event_groups = [[] for __ in range(len(station_subsets))]

    # group the events in the same manner as stations, so that the same amount of event groups are
    # created as station subsets
    for event in events:
        for i, subset in enumerate(station_subsets):
            # every station from an event must share the same station subset. Therefore its enough
            # to check only the first element
            if len(event.stations) > 0:
                if next(iter(event.stations)) in subset:
                    event_groups[i].append(event)
                    break
        else:
            if optimizer:
                optimizer.logger.warning(f"Rotation {event.rotation.id} has no possible "
                                         "electrifiable stations and will be removed.")
                # this event will not show up in an event_group.
                # therefore it needs to be put into this set
            could_not_be_electrified.update([event.rotation.id])

    groups = list(zip(event_groups, station_subsets))
    # each event group should have events and stations. If not something went wrong.
    filtered_groups = list(filter(lambda x: len(x[0]) != 0 and len(x[1]) != 0, groups))
    if len(filtered_groups) != len(groups):
        if optimizer:
            optimizer.logger.error("An event group has no possible electrifiable stations and "
                                   "will not be optimized.")
    return sorted(filtered_groups, key=lambda x: len(x[1]))


def join_all_subsets(subsets):
    """ Return sets which are joined together if they have any intersections.

    :param subsets: sets to be joined
    :type subsets: iterable
    :return: joined subsets if they connect with other subsets in some way
    :rtype: list(set)
    """

    # create set of all stations in given subsets
    all_stations = {station for subset in subsets for station in subset}
    # make list from station set to have fixed order
    all_stations_list = list(all_stations)

    # create look-up-table from station name to index in all_stations_list
    all_stations_index = dict()
    for i, station in enumerate(all_stations_list):
        all_stations_index[station] = i

    # station_array: boolean matrix of dimension n*m with n being the number of unique stations and
    # m the amount of subsets. If a station is part of the subset it will be set to True
    # this will turn the subsets
    # Sub_s_0={Station4}
    # Sub_s_1={Station1}
    # Sub_s_2={Station2}
    # Sub_s_3={Station0,Station2}
    # Sub_s_4={Station3}
    # into
    #          Sub_s_0  Sub_s_1  Sub_s_2  Sub_s_3  Sub_s_4
    # Station0 False    True     False    True     False
    # Station1 False    False    True     False    False
    # Station2 False    False    False    True     False
    # Station3 False    False    False    False    True
    # Station4 True     False    False    False    False
    station_array = np.zeros((len(all_stations), len(subsets))).astype(bool)
    for i, subset in enumerate(subsets):
        for station in subset:
            station_array[all_stations_index[station], i] = True

    # Traverse the matrix row wise. Columns with True values in this row are merged to a single one.
    # All rows of these columns are merged. This translates to subsets sharing the same station.
    # The result cannot contain any overlapping columns / sets
    #          Sub_s_0  Sub_s_1  Sub_s_2  Sub_s_3  Sub_s_4
    # Station0 False    True     False    True     False

    # returns the indicies 1 and 3. the columns are merged.
    #          Sub_s_0  Sub_s_1/Sub_s_3  Sub_s_2  Sub_s_4
    # Station0 False    True             False    False
    # Station1 False    False            True     False
    # Station2 False    True             False     False
    # Station3 False    False            False     True
    # Station4 True     False            False    False

    # The other rows are compared as well but share no True values and are not changed
    rows = station_array.shape[0]
    for row in range(rows):
        indicies = np.where(station_array[row, :])[0]
        if len(indicies) > 1:
            station_array[:, indicies[0]] = np.sum(station_array[:, indicies], axis=1).astype(bool)
            station_array = np.delete(station_array, indicies[1:], axis=1)

    # Translate the matrix back to subsets
    columns = station_array.shape[1]
    subsets = []
    for column in range(columns):
        subset = set()
        indicies = np.where(station_array[:, column])[0]
        for ind in indicies:
            subset.add(all_stations_list[ind])
        subsets.append(subset)
    return subsets


def toolbox_from_pickle(sched_name, scen_name, args_name):
    """ Load the 3 files from pickle.

    :param sched_name: name of schedule file
    :type sched_name: str
    :param scen_name: name of scenario file
    :type scen_name: str
    :param args_name: name of args file
    :type args_name: str
    :return: schedule, scenario and arguments
    :rtype: (simba.schedule.Schedule, spice_ev.Scenario, Namespace)
    """
    with open(args_name, "rb") as file:
        args = pickle.load(file)
    with open(scen_name, "rb") as file:
        scen = pickle.load(file)
    with open(sched_name, "rb") as file:
        sched = pickle.load(file)
    return sched, scen, args


def toolbox_to_pickle(name, sched, scen, args):
    """ Dump the 3 files to pickle files.

    :param name: base name of the files
    :param sched: schedule
    :type sched: simba.schedule.Schedule
    :param scen: scenario
    :type scen: spice_ev.Scenario
    :param args: arguments
    :type args: Namespace
    :return: names of the created files
    :rtype: str
    """
    args_name = "args_" + name + ".pickle"
    with open(args_name, "wb") as file:
        pickle.dump(args, file)
    scen_name = "scenario_" + name + ".pickle"
    with open(scen_name, "wb") as file:
        pickle.dump(scen, file)
    sched_name = "schedule_" + name + ".pickle"
    with open(sched_name, "wb") as file:
        pickle.dump(sched, file)
    return sched_name, scen_name, args_name


def charging_curve_to_soc_over_time(
        charging_curve, capacity, final_value: float, max_charge_from_grid=float('inf'),
        time_step=0.1, efficiency=1, eps=0.001, logger: logging.Logger = None):
    """ Create charging curve as np.array with soc and time as two columns of an np.array.

    :param logger: logger
    :type logger: logging.Logger
    :param eps: smallest normalized power, where charging curve will stop
    :type eps: float
    :param charging_curve: the charging curve with power in kw over soc
    :type charging_curve: list(float,float)
    :param capacity: capacity of the vehicle in kwh
    :type capacity: float
    :param final_value: soc value to which the curve is simulated
    :type final_value: float
    :param max_charge_from_grid: maximum amount of charge from grid / connector in kw
    :type max_charge_from_grid: float
    :param time_step: time step in minutes for simulation
    :type time_step: float
    :param efficiency: efficiency of charging
    :type efficiency: float
    :return: soc and time with n values each
    :rtype: np.array of shape (n, 2)
    """
    # simple numeric creation of power over time --> to energy over time
    normalized_curve = np.array([[soc, power / capacity] for soc, power in charging_curve])
    soc = 0
    charge_time = 0
    socs = []
    times = []

    starting_power = min(
        np.interp(soc, normalized_curve[:, 0], normalized_curve[:, 1]),
        max_charge_from_grid / capacity)
    if starting_power <= 0:
        times.append(charge_time)
        socs.append(soc)
        return np.array((times, socs)).T

    while soc < final_value:
        times.append(charge_time)
        socs.append(soc)
        power1 = min(
            np.interp(soc, normalized_curve[:, 0], normalized_curve[:, 1]),
            max_charge_from_grid / capacity)
        soc2 = soc + time_step / 60 * power1
        power2 = min(
            np.interp(soc2, normalized_curve[:, 0], normalized_curve[:, 1]),
            max_charge_from_grid / capacity)
        power = (power1 + power2) / 2 * efficiency
        delta_soc = time_step / 60 * power
        soc += delta_soc
        charge_time += time_step
        if power/starting_power < eps:
            if logger:
                warnings.warn("charging_curve_to_soc_over_time stopped early")
                logger.warning(
                    "charging_curve_to_soc_over_time stopped early, because the charging power of "
                    "%s was to low for eps: %s at an soc of %s and a desired soc of %s", power, eps,
                    soc, final_value)
            final_value = soc
            break
    # fill the soc completely in last time step
    times.append(charge_time)
    socs.append(final_value)
    return np.array((times, socs)).T


def get_missing_energy(events, min_soc=0):
    """ Sum up all the missing energies of the given events.

    :param events: events to be checked
    :type events: list(simba.optimizer_util.LowSocEvent)
    :param min_soc: minimal soc as desired value
    :type min_soc: float
    :return: missing energy
    :rtype: float
    """
    missing_energy = 0
    for event in events:
        # events should be by definition below the min_soc defined by the config
        assert event.min_soc < min_soc
        missing_energy += (event.min_soc-min_soc) * event.capacity
    return missing_energy


def stations_hash(stations_set):
    """ Create a simple string as hash for a set of stations.

    :param stations_set: stations to be hashed
    :type stations_set: set
    :return: hash
    :rtype: str
    """
    return str(sorted(stations_set))


def recursive_dict_updater(dict_to_change, filter_function, modify_function):
    """ Change nested dictionary in place given a filter and modify function.

    Goes through all values of a dictionary and modifies the value when filter criteria are met.
    The filter criteria are checked by the filter_function which gets the arguments key and value.
    The values are updated by the modify_function which gets the arguments key and value and returns
    the updated value.

    :param dict_to_change: nested dictionary that needs to be updated
    :type dict_to_change: dict
    :param filter_function: function that returns True if the value should be changed with key and
        value as arguments
    :type filter_function: function
    :param modify_function: function that returns the dictionary value with key and value as
        arguments
    :type modify_function: function

    """
    # iterate over all items. For every item, try iterating over it as well until an AttributeError
    for key, value in dict_to_change.items():
        try:
            recursive_dict_updater(value, filter_function, modify_function)
        except AttributeError:
            # could not iterate, therefore its a "final" value. check for filter condition
            # and apply the modify function
            if filter_function(key, value):
                dict_to_change[key] = modify_function(key, value)


def combination_generator(iterable: typing.Iterable, amount: int):
    """ Yields all combinations of choosing an amount, without putting back and without order.

    Generator which yields all possible combinations of choosing an amount out of an iterable
    without putting them back and without caring about the order of elements

    :param iterable: Any collection which can be cast to a list
    :rtype iterable: iterable
    :param amount: Number of elements which should be drawn from iterable
    :type amount: int
    :yields: list of items

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


def combs_unordered_no_putting_back(n: int, k: int):
    """ Return number of combinations of choosing an amount, without putting back and without order.

    Returns amount of combinations for pulling k elements out of n,
    without putting elements back or looking at the order. This is equal to n over k.

    :param n: number of elements to chose from
    :type n: int
    :param k: number of elements in the sub group of picked elements
    :type k: int
    :return: number of combinations
    :rtype: int
    """
    try:
        return math.factorial(n) / ((math.factorial(n - k)) * math.factorial(k))
    except ValueError:
        warnings.warn(f"Value Error. n={n}, k={k}")
        return 0


def run_schedule(sched, args, electrified_stations=None):
    """ Run a given schedule and electrify stations if needed.

    :param sched: schedule object
    :type sched: simba.schedule.Schedule
    :param args: arguments
    :type args: Namespace
    :param electrified_stations: electrified stations. Default value None means no further
        stations are electrified
    :type electrified_stations: dict or None
    :return: schedule and scenario after SpiceEV simulation
    :rtype: simba.schedule.Schedule, spice_ev.Scenario
    """
    sched_copy = copy(sched)
    sched_copy.stations = electrified_stations
    new_scen = sched_copy.generate_scenario(args)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        if "pytest" not in sys.modules:
            # do not print output from SpiceEV to reduce clutter. Do not do it in testing
            # since it produces errors
            sys.stdout = open(os.devnull, 'w')
        new_scen.run('distributed', vars(args).copy())
        if "pytest" not in sys.modules:
            sys.stdout = sys.__stdout__
    generate_soc_timeseries(new_scen)
    return sched_copy, new_scen


def get_time(start=[]):
    """ Prints the time which passed since the first function call.

    :param start: start time
    :type start: list(float)
    :return: String with seconds which passed since the first call.
    :rtype: str
    """
    if not start:
        start.append(time())
    delta = round(time() - start[0], 2)
    if delta > 0:
        return str(delta) + " seconds since start"


def plot_(data, axis=None):
    """ Simple plot of data without having to create subplots.

    Checks if data is nested, e.g. dict(key, list), or list(list()), and plots the inner values.
    If non-nested data is given, this data is plotted.

    :param data: data to be plotted
    :type data: iterable or iterable of iterables
    :param axis: axis to be plotted to. Default None will create new figure and axis
    :type axis: matplotlib.axes.Axes
    :return: axis of the plot """
    if axis is None:
        fig, axis = plt.subplots()
    try:
        # if the data is iterable, use the lower dimension
        _ = iter(next(iter(data)))
        for dat in data:
            try:
                name, plot_data = dat, data[dat]
            except Exception:
                plot_data = dat
                name = None
            axis.plot(plot_data, label=name, linewidth=2.0)
            axis.legend()
    except Exception:
        axis.plot(data, linewidth=2.0)
    return axis


def plot_rot(rot_id, sched, scen, axis=None, rot_only=True):
    """ Simple plot of data without having to create subplots.

    :param rot_id: id of the rotation
    :type rot_id: str
    :param sched: schedule
    :type sched: simba.schedule.Schedule
    :param scen: scenario
    :type scen: spice_ev.Scenario
    :param axis: axis to be plotted on
    :type axis: int
    :param rot_only: show only the rot or the whole vehicle socs
    :type rot_only: bool
    :return: axis of the plot
    :rtype: matplotlib.axes
    """
    soc, start, end = get_rotation_soc(rot_id, sched, scen)
    if not rot_only:
        start = 0
        end = -1
    if axis is None:
        fig, axis = plt.subplots()
        axis.plot(soc[start:end], linewidth=2.0)
        return axis
    axis.plot(soc[start:end], linewidth=2.0)
    return axis
