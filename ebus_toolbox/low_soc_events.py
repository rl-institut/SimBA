import typing
from datetime import timedelta

from ebus_toolbox.util import uncomment_json_file, get_buffer_time as get_buffer_time_spice_ev

from ebus_toolbox.schedule_optimizer import ScheduleOptimizer


class LowSocEvent:
    def __init__(self, start_idx, end_idx, min_soc, stations, vehicle_id, trip, rotation,
                 stations_list,
                 capacity, v_type, ch_type):
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


def get_low_soc_events(optimizer: ScheduleOptimizer, rotations=None, filter_standing_time=True,
                       rel_soc=False, soc_data=None, **kwargs):
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
    standing_time_min = (trip2.departure_time - trip1.arrival_time) \
                        / timedelta(minutes=1)
    buffer_time = (get_buffer_time(trip1, args.default_buffer_time_opps) / timedelta(minutes=1))
    standing_time_min -= buffer_time

    if args.min_charging_time > standing_time_min:
        return 0
    return max(0, standing_time_min)


def get_buffer_time(trip, default_buffer_time_opps):
    """  Return the buffer time as timedelta object"""
    return timedelta(minutes=get_buffer_time_spice_ev(trip, default_buffer_time_opps))


def get_index_by_time(scenario, search_time):
    start_time = scenario.start_time
    delta_time = timedelta(minutes=60 / scenario.stepsPerHour)
    idx = (search_time - start_time) // delta_time
    return idx


def get_delta_soc(soc_over_time_curve, soc, time_delta, optimizer:ScheduleOptimizer):
    """ get expected soc lift for a given start_soc and time_delta."""
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


def evaluate(events: typing.Iterable[LowSocEvent], optimizer: ScheduleOptimizer, **kwargs):
    """Analyse stations for "helpful" energy supply. Energy supply is helpful the minimal soc of an
     event is raised up to a minimal soc (probably zero). The supplied energy is approximated by
     loading power, standing time at a station, soc at station, minimal soc of the event"""

    soc_lower_thresh = kwargs.get("soc_lower_thresh", optimizer.config.min_soc)
    soc_upper_thresh = kwargs.get("soc_upper_thresh", optimizer.args.desired_soc_deps)

    # electrified_station_set = kwargs.get("electrified_station_set", set())

    station_eval = {}
    for event in events:
        soc_over_time = optimizer.soc_charge_curve_dict[event.v_type][event.ch_type]
        for i, trip in enumerate(event.trip):
            # Station is only evaluated if station name is part of event stations
            # Only these stations showed potential in electrification, e.g enough standing time
            if trip.arrival_name not in event["stations"]:
                continue
            idx = get_index_by_time(optimizer.scenario, trip.arrival_time, )
            soc = optimizer.scenario.vehicle_socs[event.vehicle_id][idx]

            # Potential is the minimal amount of
            delta_soc_pot = min(soc_upper_thresh - soc,
                                soc_lower_thresh - event.min_soc,
                                soc - event.min_soc,
                                soc_upper_thresh - soc_lower_thresh)

            try:
                standing_time_min = get_charging_time(trip, event.trip[i + 1], optimizer.args)
            except IndexError:
                standing_time_min = 0

            # energy_charging_potential = standing_time_min *60 * ch_power
            e_charging_pot = get_delta_soc(soc_over_time, soc, standing_time_min, optimizer)\
                             * event.capacity

            # Potential is at max the minimum between the useful delta soc * capacity or the
            # energy provided by charging for the full standing time
            delta_e_pot = min(delta_soc_pot * event.capacity, e_charging_pot)
            try:
                station_eval[trip.arrival_name] += delta_e_pot
            except KeyError:
                station_eval[trip.arrival_name]
                station_eval[trip.arrival_name] = delta_e_pot

    # ToDo can get true potential from decision tree but is it necessary?
    # # time_list = []
    # for station_name, stat_dict in station_eval.items():
    #     if decision_tree is not None:
    #         check_stations = electrified_station_set.union(station_name)
    #         if check_stations in decision_tree.keys():
    #             stat_dict["pot_sum"] = decision_tree[str(check_stations)]["missing_energy"] - \
    #                                    decision_tree[str(electrified_station_set)]["missing_energy"]
    #             # decision_tree[str(electrified_station_set)]["children"].append(check_stations)
    #             continue

        # time_list.append(standing_time / timedelta(minutes=1))
    # Sort by pot_sum
    station_eval = list(dict(sorted(station_eval.items(), key=lambda x: x[1])).items())
    station_eval.reverse()
    return station_eval
