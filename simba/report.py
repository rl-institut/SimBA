""" Module to generate meaningful output files and/or figures to describe simulation process
"""
import csv
import datetime
import logging
from typing import Iterable

import matplotlib.pyplot as plt
from spice_ev.report import aggregate_global_results, plot, generate_reports


def generate_gc_power_overview_timeseries(scenario, args):
    """Generate a csv timeseries with each grid connector's summed up charging station power

    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    gc_list = list(scenario.components.grid_connectors.keys())

    if not gc_list:
        return

    with open(args.results_directory / "gc_power_overview_timeseries.csv", "w", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["time"] + gc_list)
        stations = []
        time_col = getattr(scenario, f"{gc_list[0]}_timeseries")["time"]
        for i in range(len(time_col)):
            time_col[i] = time_col[i].isoformat()
        stations.append(time_col)
        for gc in gc_list:
            stations.append([-x for x in getattr(scenario, f"{gc}_timeseries")["grid supply [kW]"]])
        gc_power_overview = list(map(list, zip(*stations)))
        csv_writer.writerows(gc_power_overview)


def generate_gc_overview(schedule, scenario, args):
    """Generate a csv file with information regarding electrified stations.

    For each electrified station, the name, type, max. power, max. number of occupied
    charging stations, sum of charged energy and use factors of least used stations is saved.

    :param schedule: Driving schedule for the simulation.
    :type schedule: simba.Schedule
    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    all_gc_list = list(schedule.stations.keys())
    used_gc_list = list(scenario.components.grid_connectors.keys())
    stations = getattr(schedule, "stations")

    with open(args.results_directory / "gc_overview.csv", "w", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["station_name",
                             "station_type",
                             "maximum_power",
                             "maximum Nr charging stations",
                             "sum of CS energy",
                             "use factor least CS",
                             "use factor 2nd least CS",
                             "use factor 3rd least CS"])
        for gc in all_gc_list:
            if gc in used_gc_list:
                ts = getattr(scenario, f"{gc}_timeseries")
                max_gc_power = -min(ts["grid supply [kW]"])
                max_nr_cs = max(ts["# CS in use [-]"])
                sum_of_cs_energy = sum(ts["sum CS power [kW]"]) * args.interval/60

                # use factors: to which percentage of time are the three least used CS in use
                num_ts = scenario.n_intervals  # number of timesteps
                # three least used CS. Less if number of CS is lower.
                least_used_num = min(3, max_nr_cs)
                # count number of timesteps with this exact number of occupied CS
                count_nr_cs = [ts["# CS in use [-]"].count(max_nr_cs - i) for i in range(
                    least_used_num)]
                use_factors = [sum(count_nr_cs[:i + 1]) / num_ts for i in range(
                    least_used_num)]  # relative occupancy with at least this number of occupied CS
                use_factors = use_factors + [None] * (3 - least_used_num)  # fill up line with None

            else:
                max_gc_power = 0
                max_nr_cs = 0
                sum_of_cs_energy = 0
                use_factors = [None, None, None]
            station_type = stations[gc]["type"]
            csv_writer.writerow([gc,
                                 station_type,
                                 max_gc_power,
                                 max_nr_cs,
                                 sum_of_cs_energy,
                                 *use_factors])


def generate_plots(scenario, args):
    """Save plots as png and pdf.

    :param scenario: Scenario to plot.
    :type scenario: spice_ev.Scenario
    :param args: Configuration. Uses results_directory and show_plots.
    :type args: argparse.Namespace
    """
    aggregate_global_results(scenario)
    # disable DEBUG logging from matplotlib
    logging.disable(logging.INFO)
    with plt.ion():  # make plotting temporarily interactive, so plt.show does not block
        plt.clf()
        plot(scenario)
        plt.gcf().set_size_inches(10, 10)
        plt.savefig(args.results_directory / "run_overview.png")
        plt.savefig(args.results_directory / "run_overview.pdf")
    if args.show_plots:
        plt.show()
    # revert logging override
    logging.disable(logging.NOTSET)


def generate(schedule, scenario, args):
    """Generates all output files/ plots and saves them in the output directory.

    :param schedule: Driving schedule for the simulation.
    :type schedule: simba.Schedule
    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    :raises Exception: if cost calculation cannot be written to file and
        args.propagate_mode_errors=True
    """

    # generate simulation_timeseries.csv, simulation.json and vehicle_socs.csv in SpiceEV
    # re-route output paths
    args.save_soc = args.results_directory / "vehicle_socs.csv"
    args.save_results = args.results_directory / "info.json"
    args.save_timeseries = args.results_directory / "ts.csv"
    generate_reports(scenario, vars(args).copy())
    args.save_timeseries = None
    args.save_results = None
    args.save_soc = None

    # generate gc power overview
    generate_gc_power_overview_timeseries(scenario, args)

    # generate gc overview
    generate_gc_overview(schedule, scenario, args)

    # save plots as png and pdf
    generate_plots(scenario, args)

    # calculate SOCs for each rotation
    rotation_infos = []

    negative_rotations = schedule.get_negative_rotations(scenario)

    interval = datetime.timedelta(minutes=args.interval)
    sim_start_time = schedule.get_departure_of_first_trip()
    # rotations might be empty, start_time is None
    if sim_start_time:
        sim_start_time -= datetime.timedelta(minutes=args.signal_time_dif)
    else:
        sim_start_time = datetime.datetime.fromtimestamp(0)

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
        logging.warning(
            "SpiceEV stopped before simulation of the these rotations were completed:\n"
            f"{', '.join(incomplete_rotations)}\n"
            "Omit parameter <days> to simulate entire schedule.")

    if rotation_infos:
        try:
            # order rotations naturally
            rotations = sorted(rotation_socs.keys(), key=lambda k: int(k))
        except ValueError:
            # Some strings cannot be cast to int and throw a value error
            rotations = sorted(rotation_socs.keys(), key=lambda k: k)
        data = [["time"] + rotations]
        for i in range(scenario.n_intervals):
            t = sim_start_time + i * scenario.interval
            socs = [str(rotation_socs[k][i]) for k in rotations]
            data.append([str(t)] + socs)
        write_csv(data, args.results_directory / "rotation_socs.csv",
                  propagate_errors=args.propagate_mode_errors)

        with open(args.results_directory / "rotation_summary.csv", "w", newline='') as f:
            csv_writer = csv.DictWriter(f, list(rotation_infos[0].keys()))
            csv_writer.writeheader()
            csv_writer.writerows(rotation_infos)

    # summary of used vehicle types and all costs
    if args.cost_calculation:
        file_path = args.results_directory / "summary_vehicles_costs.csv"
        csv_report = None
        try:
            csv_report = scenario.costs.to_csv_lists()
        except Exception as e:
            logging.warning(f"Generating the cost calculation to {file_path} calculation failed "
                            f"due to {str(e)}")
            if args.propagate_mode_errors:
                raise
        write_csv(csv_report, file_path, propagate_errors=args.propagate_mode_errors)

    logging.info(f"Plots and output files saved in {args.results_directory}")


def write_csv(data: Iterable, file_path, propagate_errors=False):
    """Write iterable data to .csv file.

    Errors are caught if propagate_errors=False. If data is None, no file is created.
    :param data: data which is written to a file
    :type data: Iterable
    :param file_path: file_path for file creation
    :type file_path: str or Path
    :param propagate_errors: should errors be propagated?
    :type propagate_errors: bool
    :raises Exception: if file cannot be written propagate_errors=True
    """

    if data is None:
        logging.debug(f"Writing to {file_path} aborted since no data is passed.")
        return
    try:
        with open(file_path, "w", newline='', encoding='utf-8') as f:
            csv_writer = csv.writer(f)
            for row in data:
                csv_writer.writerow(row)
    except Exception as e:
        logging.warning(f"Writing to {file_path} failed due to {str(e)}")
        if propagate_errors:
            raise
