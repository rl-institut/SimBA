""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
import datetime

import src.scenario as scenario


import schedule
import rotation
import pickle
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
import numpy as np


with open("args_only_depot.pickle", "rb") as f: args = pickle.load(f)
with open("scenario_only_depot.pickle", "rb") as f: scen = pickle.load(f)
with open("schedule_only_depot.pickle", "rb") as f: sched = pickle.load(f)
BATTERY_CAPACITY=400
CHARGING_POWER=250

exclusions = ['BF B A',
              'BF C A',
              'BF I A',
              'BF L A',
              'BF M A',
              'BF SO A']

def get_buffer_time(search_time, args):
    for window, buffer_time in args.default_buffer_time_opps.items():
        try:
            start, end = window.split("-")

            if float(end)>search_time.hour>=float(start):
                return datetime.timedelta(minutes=buffer_time)
        except ValueError:
            return datetime.timedelta(minutes=0)



def main():
    # Find below zero soc events and stations which have charging capabilties
    events = get_below_zero_soc_events(scen=scen,
                                       rotations=sched.get_negative_rotations(scen),
                                       sched=sched,
                                       soc_upper_thresh=1, filter_standing_time=True)

    # Connect events with shared stations

    # First create simple list of station sets for single events
    possible_stations = [{station for station in e["stations"]} for e
                         in events]
    # If stations overlap join them
    station_subsets = join_all_subsets(possible_stations)

    event_groups = [[] for __ in range(len(station_subsets))]

    # Group the events in the same manner as stations, so that the same amount of event group are
    # created as station subsets
    for e in events:
        for i, subset in enumerate(station_subsets):
            # Every station from an event must share the same station sub_set. Therefore its enough
            # to check only the first element
            if len(e["stations"])>0:
                if next(iter(e["stations"])) in subset:
                    event_groups[i].append(e)
                    break

    groups = list(zip(event_groups, station_subsets))

    # Analyse stations for "helpful" energy supply. Energy supply is helpful the minimal soc of an
    # event is raised up to a minimal soc (probably zero). The supplied energy is approximated by
    # loading power, standing time at a station, soc at station, minimal soc of the event

    station_eval = dict()
    for i, group in enumerate(groups):
        event_group, station_subset = group
        for e in event_group:
            for ii, trip in enumerate(e["trip"]):
                idx = get_index_by_time(trip.arrival_time)
                soc = scen.vehicle_socs[e["vehicle_id"]][idx]
                # ToDo define elsewhere
                max_soc=1
                min_soc=0
                delta_E_pot = min(max_soc-soc, min_soc-e["min_soc"],max_soc-min_soc)
                try:
                    standing_time = e["trip"][ii+1].departure_time - trip.arrival_time
                except IndexError:
                    standing_time = datetime.timedelta(minutes=0)
                standing_time -= get_buffer_time(trip.arrival_time, args)
                standing_time_h = max(0,standing_time/datetime.timedelta(minutes=60))
                # ToDo get Vehicle Battery Capacity and charging power
                capacity= BATTERY_CAPACITY
                ch_power=CHARGING_POWER
                delta_E_pot= min(delta_E_pot*capacity, standing_time_h*ch_power)
                d = dict(E_pot=delta_E_pot,
                         standing_time=standing_time)
                try:
                    station_eval[trip.arrival_name]["pot_list"].append(d)
                except:
                    station_eval[trip.arrival_name]=dict(pot_list=[],pot_sum=0)
                    station_eval[trip.arrival_name]["pot_list"].append(d)

    time_list=[]
    for station_name, stat_dict in station_eval.items():
        sum=0
        standing_time=datetime.timedelta(minutes=0)
        for pot in stat_dict["pot_list"]:
            sum +=pot["E_pot"]
            standing_time += pot["standing_time"]
        stat_dict["pot_sum"]=sum
        time_list.append(standing_time/datetime.timedelta(minutes=1))


    # Sort by pot_sum
    station_eval=dict(sorted(station_eval.items(), key=lambda x: x[1]["pot_sum"]))

    pot=[]
    ids=[]
    for stat_id, stat in station_eval.items():
        pot.append(stat["pot_sum"])
        ids.append(stat_id)

    table= np.array((pot,ids))

    # Reverse order
    table=table[:,::-1]

    fig, ax = plt.subplots()
    ax.plot(table[0,:].astype(float), linewidth=2.0)
    ax.set_title("sortierte pot. Ladeenergie in kWh nach Haltepunkt")
    ax.set_ylabel("Pot. Ladeenergie in kWh")
    ax.set_xlabel("Haltestellen")

    # fig, ax = plt.subplots()
    # time_list=sorted(time_list)
    # time_list.reverse()
    # ax.plot(time_list, linewidth=2.0)
    # ax.set_title("sortierte Standzeiten nach Haltestellen")
    # ax.set_ylabel("kumulierte Standzeiten in min")
    # ax.set_xlabel("Haltestellen")

    # event = dict(start_idx=start, end_idx=min_idx,
    #              min_soc=min_soc, stations=possible_stations,
    #              vehicle_id=rot.vehicle_id, trip=trips,
    #              rotation=rot)

    possible_stations_ = {station for e in events for station in e["stations"] if
                          station not in exclusions}
    possible_stations = [{station for station in e["stations"] if station not in exclusions} for e
                         in events]
    station_subsets = join_all_subsets(possible_stations)

    a = set()
    for s in possible_stations:
        a.update(s)

    a = sorted(events, key=lambda x: len(x["stations"]))
    neg_rots = [sched.rotations[rot_id] for rot_id in sched.get_negative_rotations(scen)]
    neg_arrivals = {t.arrival_name for rot in neg_rots for t in rot.trips}
    neg_departure = {t.departure_name for rot in neg_rots for t in rot.trips}
    neg_stops = neg_arrivals.union(neg_departure)


def no_optimization():
    return "converged"


def plot_(data):
    fig, ax = plt.subplots()
    ax.plot(data, linewidth=2.0)

def join_all_subsets(subsets):
    joined_subset = True
    while joined_subset:
        joined_subset, subsets = join_subsets(subsets)
    return subsets


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


def get_subsets(sched: schedule.Schedule, rotations):
    # Betriebshöfe are exluced since they are depot and not opportunity. therefore they connect
    # basically all stops
    subsets = [{t.arrival_name for t in rot.trips if t.arrival_name not in exclusions} for rot in
               rotations]
    return join_all_subsets(subsets)


def get_index_by_time(search_time):
    start_time_sched = datetime.datetime.fromisoformat(sched.scenario["scenario"]["start_time"])
    delta_time = datetime.timedelta(minutes=sched.scenario["scenario"]["interval"])
    idx = (search_time - start_time_sched) // delta_time
    return idx


def get_time_by_index(idx):
    start_time_sched = datetime.datetime.fromisoformat(sched.scenario["scenario"]["start_time"])
    delta_time = datetime.timedelta(minutes=sched.scenario["scenario"]["interval"])
    searched_time = start_time_sched + delta_time * idx
    return searched_time


def get_trips(rot: rotation.Rotation, start_idx: int, end_idx: int, sched: schedule.Schedule):
    # return trips in a rotation from a start to an end index, if the arrival time is in between
    # the start and end idx
    start_time_event = get_time_by_index(start_idx)
    end_time_event = get_time_by_index(end_idx)

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


def get_below_zero_soc_events(scen: scenario.Scenario, rotations, sched: schedule.Schedule,
                              soc_upper_thresh=0.9, filter_standing_time=True):
    # Create list of events which describe trips which end in a soc below zero
    # The event is bound by the lowest soc and an upper soc threshhold which is naturally 1
    # Properties before and after these points have no effect on the event itself, similar to
    # an event horizon
    events = []
    SOC_UPPER_THRESH = soc_upper_thresh
    for rot_id in rotations:
        rot = sched.rotations[rot_id]
        rot_start_idx = get_index_by_time(rot.departure_time)
        rot_end_idx = get_index_by_time(rot.arrival_time)

        soc = scen.vehicle_socs[rot.vehicle_id]
        soc = [s if s is not None else 999 for s in soc]
        idx = range(0, len(soc))

        comb = list(zip(soc, idx))[rot_start_idx:rot_end_idx]
        min_soc, min_idx = min(comb, key=lambda x: x[0])
        reduced_list = comb.copy()
        if min_soc>=0:
            print(rot_id,min_soc)
            a=1
        while min_soc < 0:
            i = min_idx
            idx = [x[1] for x in reduced_list]
            while soc[i] < SOC_UPPER_THRESH:
                if i == rot_start_idx: break
                i -= 1
            start = idx.index(i)

            i = min_idx
            while soc[i] < SOC_UPPER_THRESH:
                if i >= rot_end_idx - 1: break
                i += 1
            end = idx.index(i)

            trips = get_trips(rot=rot, start_idx=start, end_idx=min_idx, sched=sched)
            possible_stations=set()
            if not filter_standing_time:
                possible_stations = {t.arrival_name for t in trips}
            else:
                for ii,t in enumerate(trips):
                    try:
                        standing_time = trips[ii + 1].departure_time - t.arrival_time
                        if standing_time>get_buffer_time(t.arrival_time,args):
                            possible_stations.add(t.arrival_name)
                    except IndexError:
                        pass


            all_stations = {t.arrival_name for t in rot.trips}
            # print(all_stations.difference(possible_stations) , rot_id)

            event = dict(start_idx=start, end_idx=min_idx,
                         min_soc=min_soc, stations=possible_stations,
                         vehicle_id=rot.vehicle_id, trip=trips,
                         rotation=rot)
            events.append(event)
            copy_list = reduced_list.copy()
            reduced_list = reduced_list[:start]
            if end + 1 >= len(copy_list):
                reduced_list.extend(copy_list[end + 1:])
            if len(reduced_list) > 0:
                min_soc, min_idx = min(reduced_list, key=lambda x: x[0])
            else:
                break
    return events


if __name__ == "__main__":
    main()
