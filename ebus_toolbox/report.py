""" Module to generate meaningful output files and/or figures to describe simulation process
"""
import csv
import datetime
import warnings
import json
import matplotlib.pyplot as plt

from src.report import aggregate_timeseries, aggregate_local_results, aggregate_global_results, plot


def sanitize(s, chars=''):
    """
    Removes special characters from string.

    Used to make strings safe for file paths.
    :param s: input to be sanitized
    :type s: string
    :param chars: characters to replace
    :type chars: string
    :return: input without special characters in chars
    :rtype: string
    """
    if not chars:
        chars = '</|\\>:"?*'
    return s.translate({ord(c): "" for c in chars})


def generate_vehicle_socs(scenario, args):
    """Generates a csv file from the vehicle's socs in the specified simulation time.

    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    sim_start_time = scenario.start_time
    v_list = list(scenario.vehicle_socs.keys())
    with open(args.output_directory / "vehicle_socs.csv", "w", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["timestep", "time", ] + v_list)
        for i, row in enumerate(zip(*scenario.vehicle_socs.values())):
            t = (sim_start_time + i * scenario.interval).isoformat()
            csv_writer.writerow([i, t] + list(row))


def generate_station_name_csv(scenario, args):
    """Generates a csv file from the grid connectors
    and their header information in the specified simulation time.

    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    for gc in scenario.constants.grid_connectors.keys():
        gc_info = aggregate_timeseries(scenario, gc)
        file_name = f"simulation_{sanitize(gc)}.csv"
        with open(args.output_directory / file_name, "w", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(gc_info["header"])
            for elem in gc_info["timeseries"]:
                csv_writer.writerow(elem)


def generate_station_name_json(scenario, args):
    """Generates a json file from the grid connectors
    and their header information in the specified simulation time.

    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    file_name_prefix = "simulation"
    for gc in scenario.constants.grid_connectors.keys():
        gc_info = aggregate_local_results(scenario, gc)
        file_name = f"{file_name_prefix}_{sanitize(gc)}.json"
        with open(args.output_directory / file_name, 'w') as f:
            json.dump(gc_info, f, indent=2)


def generate_gc_power_overview(scenario, args):
    """Generates a csv file from each grid connectors summed up
    charging station power in the specified simulation time.

    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    gc_list = list(scenario.constants.grid_connectors.keys())

    with open(args.output_directory / "gc_power_overview_timeseries.csv", "w", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["time", ] + gc_list)
        stations = []
        time_col = getattr(scenario, f"{gc_list[0]}_timeseries")["time"]
        for i in range(len(time_col)):
            time_col[i] = time_col[i].isoformat()
        stations.append(time_col)
        for gc in gc_list:
            stations.append([-x for x in getattr(scenario, f"{gc}_timeseries")["grid power [kW]"]])
        gc_power_overview = list(map(list, zip(*stations)))
        csv_writer.writerows(gc_power_overview)


def generate(schedule, scenario, args):
    """Generates all output files/ plots and saves them in the output directory.

    :param schedule: Driving schedule for the simulation.
    :type schedule: eBus-Toolbox.Schedule
    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    # generate csv out of vehicle's socs
    generate_vehicle_socs(scenario, args)

    # generate csv and json for all stations
    generate_station_name_csv(scenario, args)
    generate_station_name_json(scenario, args)

    # generate cs power overview
    generate_gc_power_overview(scenario, args)

    # save plots as png and pdf
    aggregate_global_results(scenario)
    with plt.ion():     # make plotting temporarily interactive, so plt.show does not block
        plot(scenario)
        plt.gcf().set_size_inches(10, 10)
        plt.savefig(args.output_directory / "run_overview.png")
        plt.savefig(args.output_directory / "run_overview.pdf")
        plt.close()

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
