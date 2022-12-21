""" Module for the class LowSocEvent which gathers functionality around
LowSocEvents, like evaluating them, gathering them and so on"""

import typing
from datetime import timedelta

from ebus_toolbox.util import get_buffer_time as get_buffer_time_spice_ev


class LowSocEvent:
    event_counter = 0

    def __init__(self, start_idx, end_idx, min_soc, stations, vehicle_id, trip, rotation,
                 stations_list, capacity, v_type, ch_type, soc_curve):
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.min_soc = min_soc
        self.stations = stations
        self.vehicle_id = vehicle_id
        self.trip = trip
        self.rotation = rotation
        self.stations_list = stations_list
        self.capacity = capacity
        self.v_type = v_type
        self.ch_type = ch_type
        self.event_counter = LowSocEvent.event_counter
        LowSocEvent.event_counter += 1


def get_low_soc_events(optimizer, rotations=None, filter_standing_time=True,
                       rel_soc=False, soc_data=None, **kwargs):
    """ Gather low soc events below the config threshold.

    :param optimizer: StationOptimizer which defines the search conditions as well as
    the scenario and schedule
    :param rotations: rotations to be searched for low soc events. Default None means whole
    schedule is searched
    :param filter_standing_time: Should the stations be filtered by standing time. True leads to
    an output only stations with charging potential
    :param rel_soc: Defines if the start soc should be seen as full even when not. If the rel_soc
    is false, it means coupled rotations might be prone to errors due to impossible lifts.
    :param soc_data: soc data to be used. Default None means the soc data from the optimizer
    scenario is used
    :param kwargs: optional soc_lower_tresh or soc_upper_thresh if from optimizer differing values
    should be used
    :return: list(LowSocEvents)
    """
    if not rotations:
        rotations = optimizer.schedule.rotations

    soc_lower_thresh = kwargs.get("soc_lower_thresh", optimizer.config.min_soc)
    soc_upper_thresh = kwargs.get("soc_upper_thresh", optimizer.args.desired_soc_deps)
    # Create list of events which describe trips which end in a soc below zero
    # The event is bound by the lowest soc and an upper soc threshhold which is naturally 1
    # Properties before and after these points have no effect on the event itself, similar to
    # an event horizon
    events = []
    count_electrified_rot = 0

    for rot_id in rotations:
        rot = optimizer.schedule.rotations[rot_id]
        soc, rot_start_idx, rot_end_idx = optimizer.get_rotation_soc(rot_id, soc_data)
        idx = range(0, len(soc))
        comb = list(zip(soc, idx))[rot_start_idx:rot_end_idx]
        min_soc, min_idx = min(comb, key=lambda x: x[0])
        reduced_list = comb.copy()
        soc_lower_thresh_cur = soc_lower_thresh
        # if rotation gets a start soc below 1 this should change below 0 soc events since fixing
        # the rotation before would lead to fixing this rotation

        # if using relative SOC, SOC lookup has to be adjusted
        if rel_soc:
            start_soc = comb[0][0]
            soc_lower_thresh_cur = min(start_soc, soc_upper_thresh) - (
                    soc_upper_thresh - soc_lower_thresh)
            soc_upper_thresh = soc_lower_thresh_cur + soc_upper_thresh
        if min_soc >= soc_lower_thresh_cur:
            count_electrified_rot += 1
        while min_soc < soc_lower_thresh_cur:
            if LowSocEvent.event_counter == 235:
                print(optimizer.scenario)
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
                if i >= rot_end_idx - 1:
                    break
                i += 1
            end_comb = idx.index(i)
            trips = optimizer.get_trips(rot=rot, start_idx=start, end_idx=min_idx)
            possible_stations = set()
            possible_stations_list = []
            if not filter_standing_time:
                possible_stations = {t.arrival_name for t in trips}
                possible_stations_list = [t.arrival_name for t in trips]
            else:
                for ii, trip in enumerate(trips):
                    try:
                        standing_time_min = get_charging_time(trip, trips[ii + 1], optimizer.args)
                    except IndexError:
                        standing_time_min = 0
                    if standing_time_min > 0:
                        possible_stations.add(trip.arrival_name)
                        possible_stations_list.append(trip.arrival_name)

            possible_stations = possible_stations.difference(optimizer.not_possible_stations)
            cht = rot.vehicle_id.find("depb")
            ch_type = (cht > 0) * "depb" + (cht <= 0) * "oppb"
            v_type = rot.vehicle_id.split("_" + ch_type)[0]
            event = LowSocEvent(start_idx=start, end_idx=min_idx,
                                min_soc=min_soc, stations=possible_stations,
                                vehicle_id=rot.vehicle_id, trip=trips,
                                rotation=rot, stations_list=possible_stations_list,
                                capacity=optimizer.schedule.vehicle_types[v_type][ch_type][
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


def get_buffer_time(trip, default_buffer_time_opps):
    """  Return the buffer time as timedelta object
    :param trip: trip object
    :param default_buffer_time_opps:
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


def get_delta_soc(soc_over_time_curve, soc, time_delta, optimizer):
    """get expected soc lift for a given start_soc and time_delta.

    :param soc_over_time_curve:
    :param soc:
    :param time_delta:
    :param optimizer:
    :return:
    """
    # units for time_delta and time_curve are assumed to be the same, e.g. minutes
    # first element which is bigger than current soc
    if time_delta == 0:
        return 0
    soc = max(min(optimizer.args.desired_soc_opps, soc), 0)
    first_time, start_soc = soc_over_time_curve[soc_over_time_curve[:, 1] >= soc][0, :]
    second_time = first_time + time_delta
    # catch out of bounds if time of charging end is bigger than table values

    if second_time >= soc_over_time_curve[-1, 0]:
        end_soc = soc_over_time_curve[-1, 1]
    else:
        end_soc = soc_over_time_curve[soc_over_time_curve[:, 0] >= second_time][0, 1]

    # make sure to limit delta soc to 1 if negative socs are given. They are possible during
    # the optimization process but will be continuously raised until they are >0.
    return min(optimizer.args.desired_soc_opps, end_soc - start_soc)


def evaluate(events: typing.Iterable[LowSocEvent], optimizer, **kwargs):
    """Analyse stations for "helpful" energy supply. Energy supply is helpful if the minimal soc of
    an event is raised (up to a minimal soc (probably zero)). The supplied energy is approximated
    by  loading power, standing time at a station, soc at station and minimal soc of the event

    :param events: events to be evaulated
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
    if not not_possible_stations:
        not_possible_stations = set()

    if not could_not_be_electrified:
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
