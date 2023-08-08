""" Collection of procedures optimizing arbitrary parameters of a bus schedule or infrastructure.
"""
from copy import deepcopy
import datetime
import logging


def service_optimization(schedule, scenario, args):
    """ Optimize rotations based on feasability.
        Try to find sets of rotations that produce no negative SoC

    :param schedule: Schedule to be optimized
    :type schedule: simba.Schedule
    :param scenario: Simulated schedule
    :type scenario: spice_ev.Scenario
    :param args: Command line arguments
    :type args: argparse.Namespace
    :return: original and most optimized scenario (highest electrification rate)
    :rtype: dict of tuple of schedule and spice_ev.Scenario
    """
    common_stations = schedule.get_common_stations(only_opps=True)

    original = (deepcopy(schedule), deepcopy(scenario))
    optimal = (None, None)

    # single out negative rotations. Try to run these with common non-negative rotations
    negative_rotations = schedule.get_negative_rotations(scenario)
    logging.info(f"Initially, rotations {sorted(negative_rotations)} have neg. SoC.")

    if not negative_rotations:
        return {
            "original": original,
            "optimized": original,
        }

    negative_sets = {}
    for rot_key in negative_rotations:
        logging.debug(f"Looking at negative rotation {rot_key}")
        # remove negative rotation from initial schedule
        rotation = schedule.rotations.pop(rot_key)
        if rotation.charging_type != "oppb":
            # only oppb rotations are optimized -> skip others
            logging.warning(f"Rotation {rot_key} should be optimized, "
                            f"but is of type {rotation.charging_type}.")
            continue
        # oppb: build non-interfering sets of negative rotations
        # (these include the dependent non-negative rotations)
        s = {rot_key}
        vid = rotation.vehicle_id
        last_neg_soc_time = scenario.negative_soc_tracker[vid][-1]
        last_neg_soc_time = datetime.datetime.fromisoformat(last_neg_soc_time)
        dependent_station = {r: t for r, t in common_stations[rot_key].items()
                             if t <= last_neg_soc_time}
        logging.debug(f"Dependent stations: {dependent_station}")
        while dependent_station:
            # get next rotation ID and last common time
            r, t = dependent_station.popitem()
            # rotation ID -> rotation object
            dependent_rotation = schedule.rotations.get(r)
            logging.debug(f"\tRotation {r} ({dependent_rotation}) @ {t}:")
            if r not in negative_rotations:
                # r not negative: add to set of dependent rotations
                s.add(r)
                # add dependencies of r, avoid circular dependency
                dependencies = {r2: t2 for r2, t2 in common_stations[r].items()
                                if t2 <= t and r2 not in s}
                dependent_station.update(dependencies)
                logging.debug(f"\t{dependencies}")
            elif dependent_rotation is not None and dependent_rotation.charging_type != "obbp":
                logging.warning(f"Rotation {rot_key} depends on negative non-oppb rotation")

        negative_sets[rot_key] = s
        logging.debug(f"Negative sets: {negative_sets}")

    # run scenario with non-negative rotations only
    if schedule.rotations:
        scenario = schedule.run(args)
        optimal = (deepcopy(schedule), deepcopy(scenario))

    # run singled-out negative rotations
    ignored = []
    for i, (rot, s) in enumerate(negative_sets.items()):
        schedule.rotations = {r: original[0].rotations[r] for r in s}
        logging.debug(f"{i+1} / {len(negative_sets)} negative schedules: {rot}")
        scenario = schedule.run(args)
        if scenario.negative_soc_tracker:
            # still fail: try just the negative rotation
            schedule.rotations = {rot: original[0].rotations[rot]}
            scenario = schedule.run(args)
            if scenario.negative_soc_tracker:
                # no hope, this just won't work
                logging.info(f"Rotation {rot} will stay negative")
            else:
                # works alone with other non-negative rotations
                logging.info(f"Rotation {rot} works alone")
            ignored.append(rot)

    negative_sets = {k: v for k, v in negative_sets.items() if k not in ignored}

    # now, in negative_sets are just rotations that work with others still
    # try to combine them
    possible = [(a, b) for a in negative_sets for b in negative_sets if a != b]
    while possible:
        logging.debug(f"{len(possible)} combinations remain")
        r1, r2 = possible.pop()
        combined = negative_sets[r1].union(negative_sets[r2])
        schedule.rotations = {r: original[0].rotations[r] for r in combined}
        scenario = schedule.run(args)
        logging.debug(f"{r1} + {r2}: {not scenario.negative_soc_tracker}")
        if not scenario.negative_soc_tracker:
            # compatible (don't interfere): keep union, remove r2
            negative_sets[r1] = combined
            # simplified. What about triangle? (r1+r2, r1+r3, r2!+r3)?
            possible = [t for t in possible if r2 not in t]
            logging.debug(possible)

            # save scenario with highest electrification rate
            if optimal is None or len(scenario[0].rotations) > len(optimal[0].rotations):
                optimal = (deepcopy(schedule), deepcopy(scenario))

    logging.info("Service optimization finished")

    return {
        "original": original,
        "optimized": optimal,
    }
