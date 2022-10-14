""" Module to generate meaningful output files and/or figures to describe simulation process
"""
import csv
import datetime
import warnings


def generate(schedule, scenario, args):
    rotation_infos = []

    negative_rotations = schedule.get_negative_rotations(scenario)

    interval = datetime.timedelta(minutes=args.interval)
    sim_start_time = \
        schedule.get_departure_of_first_trip() - datetime.timedelta(minutes=args.signal_time_dif)

    incomplete_rotations = []
    rotation_socs = {}
    for id, rotation in schedule.rotations.items():
        # get SOC timeseries for this rotation
        vehicle_id = rotation.vehicle_id

        # get soc timeseries for current rotation
        vehicle_soc = scenario.vehicle_socs[vehicle_id]
        start_idx = (rotation.departure_time - sim_start_time) // interval
        end_idx = start_idx + ((rotation.arrival_time-rotation.departure_time) // interval)
        if end_idx > scenario.n_intervals:
            # SpiceEV stopped before rotation was fully simulated
            incomplete_rotations.append(id)
            continue
        rotation_soc_ts = vehicle_soc[start_idx:end_idx]

        # bus does not return before simulation end
        # replace trailing None values with last numeric value
        for i, soc in enumerate(reversed(rotation_soc_ts)):
            if soc is not None:
                break
        last_known_idx = len(rotation_soc_ts) - 1 - i
        rotation_soc_ts[last_known_idx+1:] = i * [rotation_soc_ts[last_known_idx]]

        rotation_info = {
            "rotation_id": id,
            "start_time": rotation.departure_time.isoformat(),
            "end_time": rotation.arrival_time.isoformat(),
            "vehicle_type": rotation.vehicle_type,
            "vehicle_id": rotation.vehicle_id,
            "depot_name": rotation.departure_name,
            "lines": ':'.join(rotation.lines),
            "total_consumption_[kWh]": rotation.consumption,
            "distance": rotation.distance,
            "charging_type": rotation.charging_type,
            "SOC_at_arrival": rotation_soc_ts[-1],
            "Minimum_SOC": min(rotation_soc_ts),
            "Negative_SOC": 1 if id in negative_rotations else 0
        }
        rotation_infos.append(rotation_info)

        # save SOCs for each rotation
        rotation_socs[id] = [None]*scenario.n_intervals
        rotation_socs[id][start_idx:end_idx] = rotation_soc_ts

    if incomplete_rotations:
        warnings.warn("SpiceEV stopped before simulation of the these rotations were completed:\n"
                      f"{', '.join(incomplete_rotations)}\n"
                      "Omit parameter <days> to simulate entire schedule.",
                      stacklevel=100)

    with open(args.output_directory / "rotation_socs.csv", "w+", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(("time",) + tuple(rotation_socs.keys()))
        for i, row in enumerate(zip(*rotation_socs.values())):
            t = sim_start_time + i * scenario.interval
            csv_writer.writerow((t,) + row)

    with open(args.output_directory / "rotation_summary.csv", "w+", newline='') as f:
        csv_writer = csv.DictWriter(f, list(rotation_infos[0].keys()))
        csv_writer.writeheader()
        csv_writer.writerows(rotation_infos)
