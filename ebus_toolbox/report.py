""" Module to generate meaningful output files and/or figures to describe simulation process
"""
import csv
import datetime
import warnings

from src.report import aggregate_timeseries


def generate_vehicle_socs(scenario, args):

    sim_start_time = scenario.start_time
    v_list = []
    for v in scenario.vehicle_socs:
        v_list.append(v)

    with open(args.output_directory / "vehicle_socs.csv", "w+", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["timestep", "time",] + v_list)
        for i, row in enumerate(zip(*scenario.vehicle_socs.values())):
            t = sim_start_time + i * scenario.interval
            csv_writer.writerow((i, t,) + row)


def generate_station_name(scenario, args):
    for gc in scenario.gcPowerSchedule:
        gc_info = aggregate_timeseries(scenario, gc)
        with open(args.output_directory / f"simulation_spiceEV_{gc}.csv", "w+", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(gc_info["header"])
            for elem in gc_info["timeseries"]:
                csv_writer.writerow(elem)


def generate(schedule, scenario, args):

    # create csv out of vehicle's socs
    generate_vehicle_socs(scenario, args)

    # create csv for all stations
    generate_station_name(scenario, args)

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
        end_idx = start_idx + ((rotation.arrival_time - rotation.departure_time) // interval)
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
        rotation_soc_ts[last_known_idx + 1:] = i * [rotation_soc_ts[last_known_idx]]

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
        rotation_socs[id] = [None] * scenario.n_intervals
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

    # summary of used vehicle types and all costs
    if args.cost_calculation:
        with open(args.output_directory / "summary_vehicles_costs.csv", "w", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(["parameter", "value", "unit"])
            for key, value in schedule.vehicle_type_counts.items():
                if value > 0:
                    csv_writer.writerow([key, value, "vehicles"])
            for key, value in scenario.costs.items():
                if "annual" in key:
                    csv_writer.writerow([key, round(value, 2), "€/year"])
                else:
                    csv_writer.writerow([key, round(value, 2), "€"])
