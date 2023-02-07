""" Module for the class LowSocEvent which gathers functionality around
LowSocEvents, like evaluating them, gathering them and so on"""
import math
import os
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
    from ebus_toolbox.station_optimizer import StationOptimizer

from ebus_toolbox.consumption import Consumption
from ebus_toolbox.costs import calculate_costs
from ebus_toolbox.trip import Trip
from ebus_toolbox.util import get_buffer_time as get_buffer_time_spice_ev, uncomment_json_file


class ChargingEvent:
    """ Class to gather information about a charging event"""
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
    """Class to gather information about a low soc event"""
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
    """Class for the configuration file"""

    def __init__(self):
        self.debug_level = None
        self.exclusion_rots = None
        self.exclusion_stations = None
        self.inclusion_stations = None
        self.standard_opp_station = None
        self.schedule = None
        self.scenario = None
        self.args = None
        self.charge_eff = None
        self.battery_capacity = None
        self.charging_curve = None
        self.charging_power = None
        self.min_soc = None
        self.solver = None
        self.rebase_scenario = None
        self.pickle_rebased = None
        self.pickle_rebased_name = None
        self.opt_type = None
        self.remove_impossible_rots = None
        self.node_choice = None
        self.max_brute_loop = None
        self.run_only_neg = None
        self.estimation_threshold = None
        self.output_path = None
        self.check_for_must_stations = None
        self.decision_tree_path = None
        self.save_decision_tree = None
        self.reduce_rots = None
        self.rots = None
        self.path = None
        self.pruning_threshold = None
        self.save_all_results = None


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


def read_config(config_path):
    """ Read the config path to a config object
    :param config_path: path to file
    :return: config object
    """
    config_parser = configparser.ConfigParser()
    config_parser.sections()
    config_parser.read(config_path, encoding="utf-8")

    conf = OptimizerConfig()
    conf.path = config_path

    default = config_parser["DEFAULT"]
    conf.debug_level = int(default.get("debug_level", "0"))
    sce = config_parser["SCENARIO"]
    conf.exclusion_rots = set(json.loads(sce.get("exclusion_rots", "[]")))
    conf.exclusion_stations = set(json.loads(sce.get("exclusion_stations", "[]")))
    conf.inclusion_stations = set(json.loads(sce.get("inclusion_stations", "[]")))
    conf.standard_opp_station = dict(json.loads(sce.get("standard_opp_station", "{}")))

    pick = config_parser["PICKLE"]
    conf.schedule = pick.get("schedule", "")
    conf.scenario = pick.get("scenario", "")
    conf.args = pick.get("args", "")

    vehicle = config_parser["VEHICLE"]
    conf.charge_eff = float(vehicle.get("charge_eff", "0.95"))
    conf.battery_capacity = float(vehicle.get("battery_capacity", "0"))
    if conf.battery_capacity == 0:
        conf.battery_capacity = None
    conf.charging_curve = json.loads(vehicle.get("charging_curve", "[]"))
    if not conf.charging_curve:
        conf.charging_curve = None
    conf.charging_power = float(vehicle.get("charging_power", "0"))
    if conf.charging_power == 0:
        conf.charging_power = None
    conf.min_soc = float(vehicle.get("min_soc", "0.0"))

    optimizer = config_parser["OPTIMIZER"]
    conf.solver = optimizer.get("solver", "spiceev")
    conf.rebase_scenario = optimizer.getboolean("rebase_scenario", True)
    conf.pickle_rebased = optimizer.getboolean("pickle_rebased", False)
    conf.pickle_rebased_name = optimizer.get("pickle_rebased_name",
                                             "rebased_" + str(
                                                 datetime.now().strftime("%Y-%m-%d-%H-%M-%S")))
    conf.opt_type = optimizer.get("opt_type", "greedy")
    conf.remove_impossible_rots = optimizer.getboolean("remove_impossible_rots", False)
    conf.node_choice = optimizer.get("node_choice", "step-by-step")
    conf.max_brute_loop = int(optimizer.get("max_brute_loop", "20"))
    conf.run_only_neg = optimizer.getboolean("run_only_neg", False)
    conf.estimation_threshold = float(optimizer.get("estimation_threshold", "0.8"))
    conf.output_path = optimizer.get("output_path")
    conf.check_for_must_stations = optimizer.getboolean("check_for_must_stations", True)
    conf.pruning_threshold = int(optimizer.get("pruning_threshold", "3"))
    conf.save_all_results = optimizer.getboolean("save_all_results", False)

    special = config_parser["SPECIAL"]
    conf.decision_tree_path = special.get("decision_tree_path", None)
    if conf.decision_tree_path in ["", '""', "''"]:
        conf.decision_tree_path = None
    conf.save_decision_tree = special.getboolean("save_decision_tree", False)
    conf.reduce_rots = special.getboolean("reduce_rots", False)
    conf.rots = json.loads(special.get("rots", []))

    return conf


def get_charging_time(trip1, trip2, args):
    """ Returns the charging time between trips as numeric value

    :param trip1: First trip
    :param trip2: Following trip
    :param args:  arguments Namespace with default buffer time
    :return: maximum possible charging time in minutes between trips
    """
    standing_time_min = (trip2.departure_time - trip1.arrival_time) / timedelta(minutes=1)
    buffer_time = (get_buffer_time(trip1, args.default_buffer_time_opps) / timedelta(minutes=1))
    standing_time_min -= buffer_time

    if args.min_charging_time > standing_time_min:
        return 0
    return max(0, standing_time_min)


def get_charging_start(trip1, args):
    """ Returns the possible start of charging consindering buffer times
    :param trip1: First trip
    :param args:  arguments Namespace with default buffer time
    :return: First possible charging time as datetime object
    """
    buffer_time = get_buffer_time(trip1, args.default_buffer_time_opps)
    return trip1.arrival_time+buffer_time


def get_buffer_time(trip, default_buffer_time_opps):
    """  Return the buffer time as timedelta object
    :param trip: trip object
    :param default_buffer_time_opps: the default buffer time at opps charging stations
    :return: timedelta object for the buffer time
    """
    return timedelta(minutes=get_buffer_time_spice_ev(trip, default_buffer_time_opps))


def get_index_by_time(scenario, search_time):
    """ Get the index for a given time
    :param scenario: scenario object
    :param search_time: search time as datetime object
    :return: index as int
    """
    start_time = scenario.start_time
    delta_time = timedelta(minutes=60 / scenario.stepsPerHour)
    idx = (search_time - start_time) // delta_time
    return idx


def get_rotation_soc_util(rot_id, this_sched, this_scen, soc_data: dict = None):
    """Gets you the soc object with start and end index for a given rotation id
    :param rot_id: rotation_id
    :param this_sched: schedule object contain rotation information
    :param this_scen: scenario object containing the soc data
    :param soc_data: optional soc_data if not the scenario data should be used
    :return: tuple with soc array, start index and end index
    """
    rot = this_sched.rotations[rot_id]
    rot_start_idx = get_index_by_time(this_scen, rot.departure_time)
    rot_end_idx = get_index_by_time(this_scen, rot.arrival_time)
    if soc_data:
        return soc_data[rot.vehicle_id], rot_start_idx, rot_end_idx
    return this_scen.vehicle_socs[rot.vehicle_id], rot_start_idx, rot_end_idx


def get_delta_soc(soc_over_time_curve, soc, time_delta, optimizer: 'StationOptimizer'):
    """get expected soc lift for a given start_soc and time_delta.

    :param soc_over_time_curve: array with socs over time
    :param soc: start socs
    :param time_delta: time of charging
    :param optimizer: optimizer object
    :return: (float) the lift of the soc
    """
    # units for time_delta and time_curve are assumed to be the same, e.g. minutes
    # first element which is bigger than current soc
    if time_delta == 0:
        return 0
    soc = max(min(optimizer.args.desired_soc_opps, soc), 0)
    idx = np.searchsorted(soc_over_time_curve[:, 1], soc, side='left')
    first_time, start_soc = soc_over_time_curve[idx, :]
    second_time = first_time + time_delta
    # catch out of bounds if time of charging end is bigger than table values

    if second_time >= soc_over_time_curve[-1, 0]:
        end_soc = soc_over_time_curve[-1, 1]
    else:
        end_soc = soc_over_time_curve[soc_over_time_curve[:, 0] >= second_time][0, 1]

    # make sure to limit delta soc to 1 if negative socs are given. They are possible during
    # the optimization process but will be continuously raised until they are >0.
    return min(optimizer.args.desired_soc_opps, optimizer.args.desired_soc_opps - start_soc,
               end_soc - start_soc)


class SuboptimalSimulationException(Exception):
    pass


class AllCombinationsCheckedException(Exception):
    pass


@time_it
def evaluate(events: typing.Iterable[LowSocEvent],
             optimizer: 'StationOptimizer', **kwargs):
    """Analyse stations for "helpful" energy supply. Energy supply is helpful if the minimal soc of
    an event is raised (up to a minimal soc (probably zero)). The supplied energy is approximated
    by  loading power, standing time at a station, soc at station and minimal soc of the event

    :param events: events to be evaluated
    :param optimizer: StationOptimizer object with scenario and schedule data
    :param kwargs: optional overwriting of soc_lower_thresh, soc_upper_thresh or soc_data
    :return: sorted list with the best station and its potential on index 0
    :rtype: list(str(station_id), float(potential))
    """
    soc_lower_thresh = kwargs.get("soc_lower_thresh", optimizer.config.min_soc)
    soc_upper_thresh = kwargs.get("soc_upper_thresh", optimizer.args.desired_soc_deps)
    soc_data = kwargs.get("soc_data", optimizer.scenario.vehicle_socs)

    station_eval = {}
    # cycle through events and determine how much lift can be provided by electrifying a station
    # the lift is determined by the soc position, standing time, power supply and charging curve
    for e in events:
        soc_over_time = optimizer.soc_charge_curve_dict[e.v_type][e.ch_type]
        for i, trip in enumerate(e.trip):
            # station is only evaluated if station name is part of event stations
            # only these stations showed potential in electrification, e.g enough standing time
            if trip.arrival_name not in e.stations:
                continue
            idx = get_index_by_time(optimizer.scenario, trip.arrival_time, )
            soc = soc_data[e.vehicle_id][idx]

            # potential is the minimal amount of
            delta_soc_pot = min(soc_upper_thresh - soc,
                                soc_lower_thresh - e.min_soc,
                                soc - e.min_soc,
                                soc_upper_thresh - soc_lower_thresh)

            try:
                standing_time_min = get_charging_time(trip, e.trip[i + 1], optimizer.args)
            except IndexError:
                standing_time_min = 0

            pot_kwh = get_delta_soc(soc_over_time, soc, standing_time_min, optimizer) * e.capacity

            # potential is at max the minimum between the useful delta soc * capacity or the
            # energy provided by charging for the full standing time
            delta_e_pot = min(delta_soc_pot * e.capacity, pot_kwh)
            try:
                station_eval[trip.arrival_name] += delta_e_pot
            except KeyError:
                station_eval[trip.arrival_name] = delta_e_pot

    # sort by pot_sum
    station_eval = list(dict(sorted(station_eval.items(), key=lambda x: x[1])).items())
    station_eval.reverse()
    return station_eval


def get_groups_from_events(events, not_possible_stations=None, could_not_be_electrified=None,
                           optimizer=None):
    """ Create groups from events which need to be optimized together

    First it creates a simple list of station sets for single events. They are connected if they
    share possible stations.
    Electrified and other not possible to electrify stations should not connect groups

    :param events: events for a given state of a scenario
    :param not_possible_stations: stations to be discarded
    :param could_not_be_electrified: rotations to be discarded
    :param optimizer: StationOptimizer object
    :return: list((events, stations)) which give you events to optimize together with the stations
        which might help them
    """

    # making sure default arguments are none and not mutable
    if not_possible_stations is None:
        not_possible_stations = set()

    if could_not_be_electrified is None:
        could_not_be_electrified = set()

    possible_stations = [
        {station for station in event.stations if station not in not_possible_stations}
        for event
        in events]
    # if stations overlap join them
    station_subsets = join_all_subsets(possible_stations)
    event_groups = [[] for __ in range(len(station_subsets))]

    # group the events in the same manner as stations, so that the same amount of event groups are
    # created as station subsets
    for event in events:
        for i, subset in enumerate(station_subsets):
            # every station from an event must share the same station sub_set. Therefore its enough
            # to check only the first element
            if len(event.stations) > 0:
                if next(iter(event.stations)) in subset:
                    event_groups[i].append(event)
                    break
        else:
            if optimizer:
                optimizer.logger.warning('Did not find rotation %s in any subset'
                                         'of possible electrifiable stations', event.rotation.id)
                # this event will no show up in an event_group.
                # therefore it needs to be put into this set
            could_not_be_electrified.update([event.rotation.id])

    groups = list(zip(event_groups, station_subsets))
    return sorted(groups, key=lambda x: len(x[1]))


def join_all_subsets(subsets):
    """ join sets for as long as needed until no elements share any intersections
    :param subsets: iterable of sets
    :return: joined subsets if they connect with other subsets in some way
    """
    joined_subset = True
    while joined_subset:
        joined_subset, subsets = join_subsets(subsets)
    return subsets


def join_subsets(subsets: typing.Iterable[set]):
    """ Run through every subset and check with every other subset if there is an intersection
    If an intersection is found. The subsets are joined and returned with a boolean of True.
    If not intersection is found over all subsets False is returned which will cancel the outer
    call in join_all_subsets

    :param subsets: iterable of sets
    :return: boolean if joining subsets is finished (i.e. False if all subsets are connected
        and thus far connected subsets
    """
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


def toolbox_to_pickle(name, sched, scen, this_args):
    """ Dump the 3 files to pickle files

    :param name: base name of the files
    :param sched: schedule object
    :param scen: scenario object
    :param this_args: args Namespace object
    :return: Names of the created files
    """
    args_name = "args_" + name + ".pickle"
    with open(args_name, "wb") as file:
        pickle.dump(this_args, file)
    scen_name = "scenario_" + name + ".pickle"
    with open(scen_name, "wb") as file:
        pickle.dump(scen, file)
    sched_name = "schedule_" + name + ".pickle"
    with open(sched_name, "wb") as file:
        pickle.dump(sched, file)
    return sched_name, scen_name, args_name


def charging_curve_to_soc_over_time(charging_curve, capacity, args,
                                    max_charge_from_grid=float('inf'),
                                    time_step=0.1, efficiency=1):
    """ create charging curve as nested list of SOC, Power[kW] and capacity in [kWh]

    :param charging_curve: the charging curve with power over soc
    :param capacity: capacity of the vehicle
    :param args: args namespace object for simulation arguments
    :param max_charge_from_grid: maximum amount of charge from grid / connector
    :param time_step: time step for simulation
    :param efficiency: efficiency of charging
    :return: np.array with soc over time
    """
    # simple numeric creation of power over time --> to energy over time
    normalized_curve = np.array([[soc, power / capacity] for soc, power in charging_curve])
    soc = 0
    charge_time = 0
    socs = []
    times = []
    final_value = args.desired_soc_opps
    while soc < args.desired_soc_opps:
        times.append(charge_time)
        socs.append(soc)
        power1 = min(np.interp(soc, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        soc2 = soc + time_step / 60 * power1
        power2 = min(np.interp(soc2, normalized_curve[:, 0], normalized_curve[:, 1]),
                     max_charge_from_grid / capacity)
        power = (power1 + power2) / 2 * efficiency
        delta_soc = time_step / 60 * power
        soc += delta_soc
        charge_time += time_step
        if args.desired_soc_opps-soc < 0.0001 and delta_soc < 1e-10:
            final_value = soc
            break
    # fill the soc completely in last time step
    times.append(charge_time)
    socs.append(final_value)
    return np.array((times, socs)).T


def get_missing_energy(events):
    """ Sum up all the missing energies of the given events

    :param events: events to be checked
    :return: (float) missing energy
    """
    missing_energy = 0
    for event in events:
        missing_energy += event.min_soc * event.capacity
    return missing_energy


def stations_hash(stations_set):
    """ Create a simple str as hash for a set of stations

    :param stations_set: set of stations to be hashed
    :return: hash string
    """
    return str(sorted(list(stations_set)))


def combination_generator(iterable: typing.Iterable, amount: int):
    """ Generator which yields all possible combinations of choosing
    an amount out of an iterable without putting them back and without caring about the
    order of elements
    :param iterable: Any collection which can be cast to a list
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


def toolbox_from_pickle(sched_name, scen_name, args_name):
    """ Load the 3 files from pickle

    :param sched_name: name of schedule file
    :param scen_name: name of scenario file
    :param args_name: name of args file
    :return: schedule, scenario and args object
    """
    with open(args_name, "rb") as file:
        this_args = pickle.load(file)
    with open(scen_name, "rb") as file:
        scen = pickle.load(file)
    with open(sched_name, "rb") as file:
        sched = pickle.load(file)
    return sched, scen, this_args


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
        warnings.warn("Value Error")
        return 0


def run_schedule(this_sched, this_args, electrified_stations=None, cost_calc=False):
    """Run a given schedule and electrify stations if need be
    :param this_sched: schedule object
    :param this_args: args namespace object
    :param electrified_stations: dict of electrified stations. Default value None means no further
        stations are electrified
    :param cost_calc: should the cost be calculated
    :raises SystemExit: in case of wrong cost calculation file
    :return: schedule and scenario objects after spiceev simulation
    """
    this_sched2 = copy(this_sched)
    this_sched2.stations = electrified_stations
    this_sched2, new_scen = preprocess_schedule(this_sched2, this_args,
                                                electrified_stations=electrified_stations)
    # do not print output from spice ev to reduce clutter
    sys.stdout = open(os.devnull, 'w')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        new_scen.run('distributed', vars(this_args).copy())
    sys.stdout = sys.__stdout__
    if this_args.cost_calculation and cost_calc:
        # cost calculation following directly after simulation
        try:
            with open(this_args.cost_parameters_file, encoding='utf-8') as file:
                cost_parameters_file = uncomment_json_file(file)
        except FileNotFoundError:
            raise SystemExit(f"Path to cost parameters ({this_args.cost_parameters_file}) "
                             "does not exist. Exiting...")
        calculate_costs(cost_parameters_file, new_scen, this_sched2, this_args)
    return this_sched2, new_scen


def preprocess_schedule(this_sched, this_args, electrified_stations=None):
    """ Prepare the schedule by calculating consumption, setting elctrified stations and assigning
    vehicles

    :param this_sched: schedule containing the rotations
    :param this_args: arguments for simulation
    :param electrified_stations: dict of stations to be electrified
    :return: schedule and scenario to be simulated
    """
    Trip.consumption = \
        Consumption(this_sched.vehicle_types,
                    outside_temperatures=this_args.outside_temperature_over_day_path,
                    level_of_loading_over_day=this_args.level_of_loading_over_day_path)

    this_sched.stations = electrified_stations
    this_sched.calculate_consumption()
    this_sched.assign_vehicles()

    return this_sched, this_sched.generate_scenario(this_args)


def print_time(start=[]):
    """Print the time and automatically set start time to the first time the function getting
    called
    :param start: list for storing times
    """
    if not start:
        start.append(time())
    delta = round(time() - start[0], 2)
    if delta > 0:
        print(delta, " seconds till start")


def plot_(data):
    """ Simple plot of data without having to create subplots
    :param data: data to be plotted
    :return: axis of the plot """
    fig, axis = plt.subplots()
    axis.plot(data, linewidth=2.0)
    return axis


def plot_rot(rot_id, this_sched, this_scen, axis=None, rot_only=True):
    """Simple plot of data without having to create subplots
    :param rot_id: id of the rotation
    :param this_sched: schedule object
    :param this_scen: scenario object
    :param axis: axis to be plotted on
    :param rot_only: show only the rot or the whole vehicle socs
    :return: axis of the plot
    """
    soc, start, end = get_rotation_soc_util(rot_id, this_sched, this_scen)
    if not rot_only:
        start = 0
        end = -1
    if axis is None:
        fig, axis = plt.subplots()
        axis.plot(soc[start:end], linewidth=2.0)
        return axis
    axis.plot(soc[start:end], linewidth=2.0)
    return axis
