# imports
from copy import deepcopy
import datetime
import json
import warnings
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from ebus_toolbox.trip import Trip
from ebus_toolbox import report  # , optimizer


def simulate(args):
    """Simulate the given scenario and eventually optimize for given metric(s).

    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """
    try:
        with open(args.vehicle_types) as f:
            vehicle_types = json.load(f)
    except FileNotFoundError:
        warnings.warn("Invalid path for vehicle type JSON. Using default types from EXAMPLE dir.")
        with open("data/examples/vehicle_types.json") as f:
            vehicle_types = json.load(f)

    # parse strategy options for Spice EV
    if args.strategy_option is not None:
        for opt_key, opt_val in args.strategy_option:
            try:
                # option may be number
                opt_val = float(opt_val)
            except ValueError:
                # or not
                pass
            setattr(args, opt_key, opt_val)

    schedule = Schedule.from_csv(args.input_schedule, vehicle_types, args.electrified_stations)
    # setup consumption calculator that can be accessed by all trips
    Trip.consumption = Consumption(vehicle_types)
    schedule.calculate_consumption()
    schedule.set_charging_type(preferred_ct=args.preferred_charging_type, args=args)

    common_stations = schedule.get_common_stations(only_opps=True)

    # initial run
    scenario = schedule.run(args)

    # single out negative rotations. Try to run these with common non-negative rotations
    original_rotations = deepcopy(schedule.rotations)
    negative_rotations = schedule.get_negative_rotations(scenario)
    print(f"Initially, rotations {sorted(negative_rotations)} have neg. SoC.")

    negative_sets = {}
    for rot_key in negative_rotations:
        rotation = schedule.rotations[rot_key]
        if rotation.charging_type == "depb":
            schedule.set_charging_type("oppb", args, [rot_key])
            # todo: actually this case should not happen, but it still does happen.. why?
        else:
            # oppb: build non-interfering sets of negative rotations
            s = {rot_key}
            vid = rotation.vehicle_id
            last_neg_soc_time = scenario.negative_soc_tracker[vid][-1]
            last_neg_soc_time = datetime.datetime.fromisoformat(last_neg_soc_time)
            dep_stations = {r: t for r, t in common_stations[rot_key].items()
                            if t <= last_neg_soc_time}
            while dep_stations:
                r, t = dep_stations.popitem()
                if r not in negative_rotations:
                    s.add(r)
                    # add dependencies of r
                    dep_stations.update({r2: t2 for r2, t2
                                        in common_stations[r].items() if t2 <= t})
            negative_sets[rot_key] = s

    # run singled-out rotations
    ignored = []
    for rot, s in negative_sets.items():
        schedule.rotations = {r: original_rotations[r] for r in s}
        scenario = schedule.run(args)
        if scenario.negative_soc_tracker:
            # still fail: try just the negative rotation
            schedule.rotations = {rot: original_rotations[rot]}
            scenario = schedule.run(args)
            if scenario.negative_soc_tracker:
                # no hope, this just won't work
                print(f"Rotation {rot} will stay negative")
            else:
                # works alone with other non-negative rotations
                print(f"Rotation {rot} works with alone")
            ignored.append(rot)

    negative_sets = {k: v for k, v in negative_sets.items() if k not in ignored}

    # now, in negative_sets are just rotations that work with others still
    # try to combine them
    possible = [(a, b) for a in negative_sets for b in negative_sets if a != b]
    while possible:
        r1, r2 = possible.pop()
        combined = negative_sets[r1].union(negative_sets[r2])
        schedule.rotations = {r: original_rotations[r] for r in s}
        scenario = schedule.run(args)
        if not scenario.negative_soc_tracker:
            # compatible (don't interfere): keep union, remove r2
            negative_sets[r1] = combined
            # simplified. What about triangle? (r1+r2, r1+r3, r2!+r3)?
            possible = [t for t in possible if t[1] != r2]

    print(negative_sets)

    # create report
    report.generate(schedule, scenario, args)
