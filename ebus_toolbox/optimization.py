""" Collection of procedures optimizing arbitrary parameters of a bus schedule or infrastructure.
"""

import datetime
from copy import deepcopy


def service_optimization(schedule, scenario, args):
    """ Optimize rotations based on feasability.
        Try to find sets of rotations that produce no negative SoC

    :param schedule: Schedule to be optimized
    :type schedule: ebus_toolbox.Schedule
    :param scenario: Simulated schedule
    :type scenario: spice_ev.Scenario
    :param args: Command line arguments
    :type args: argparse.Namespace
    :return: original and most optimized scenario (highest electrification rate)
    :rtype: dict of tuple of schedule and spice_ev.Scenario
    :raises Exception: if the charging type of a rotation with negative SoC is not oppb
    """
    common_stations = schedule.get_common_stations(only_opps=True)

    original = (deepcopy(schedule), deepcopy(scenario))
    optimal = (None, None)

    # single out negative rotations. Try to run these with common non-negative rotations
    negative_rotations = schedule.get_negative_rotations(scenario)
    print(f"Initially, rotations {sorted(negative_rotations)} have neg. SoC.")

    negative_sets = {}
    for rot_key in negative_rotations:
        rotation = schedule.rotations[rot_key]
        if rotation.charging_type != "oppb":
            raise Exception(f"Rotation {rot_key} should be optimized, "
                            f"but is of type {rotation.charging_type}.")
        # oppb: build non-interfering sets of negative rotations
        # (these include the dependent non-negative rotations)
        s = {rot_key}
        vid = rotation.vehicle_id
        last_neg_soc_time = scenario.negative_soc_tracker[vid][-1]
        last_neg_soc_time = datetime.datetime.fromisoformat(last_neg_soc_time)
        dependent_station = {r: t for r, t in common_stations[rot_key].items()
                             if t <= last_neg_soc_time}
        while dependent_station:
            r, t = dependent_station.popitem()
            if r not in negative_rotations:
                s.add(r)
                # add dependencies of r
                dependent_station.update({r2: t2 for r2, t2
                                          in common_stations[r].items() if t2 <= t})
        negative_sets[rot_key] = s
        # remove negative rotation from initial schedule
        del schedule.rotations[rot_key]

    # run scenario with non-negative rotations only
    if schedule.rotations:
        scenario = schedule.run(args)
        optimal = (deepcopy(schedule), deepcopy(scenario))

    # run singled-out negative rotations
    ignored = []
    for i, (rot, s) in enumerate(negative_sets.items()):
        schedule.rotations = {r: original[0].rotations[r] for r in s}
        print(f"{i+1} / {len(negative_sets)} negative schedules: {rot}")
        scenario = schedule.run(args)
        if scenario.negative_soc_tracker:
            # still fail: try just the negative rotation
            schedule.rotations = {rot: original[0].rotations[rot]}
            scenario = schedule.run(args)
            if scenario.negative_soc_tracker:
                # no hope, this just won't work
                print(f"Rotation {rot} will stay negative")
            else:
                # works alone with other non-negative rotations
                print(f"Rotation {rot} works alone")
            ignored.append(rot)

    negative_sets = {k: v for k, v in negative_sets.items() if k not in ignored}

    # now, in negative_sets are just rotations that work with others still
    # try to combine them
    possible = [(a, b) for a in negative_sets for b in negative_sets if a != b]
    while possible:
        print(f"{len(possible)} combinations remain")
        r1, r2 = possible.pop()
        combined = negative_sets[r1].union(negative_sets[r2])
        schedule.rotations = {r: original[0].rotations[r] for r in combined}
        scenario = schedule.run(args)
        if not scenario.negative_soc_tracker:
            # compatible (don't interfere): keep union, remove r2
            negative_sets[r1] = combined
            # simplified. What about triangle? (r1+r2, r1+r3, r2!+r3)?
            possible = [t for t in possible if r2 not in t]

            # save scenario with highest electrification rate
            if optimal is None or len(scenario[0].rotations) > len(optimal[0].rotations):
                optimal = (deepcopy(schedule), deepcopy(scenario))

    print(negative_sets)

    return {
        "original": original,
        "optimized": optimal,
    }
