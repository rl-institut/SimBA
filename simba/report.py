""" Module to generate meaningful output files and/or figures to describe simulation process. """
import csv
import datetime
import logging
from math import ceil
import re
from typing import Iterable

import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from spice_ev.report import aggregate_global_results, plot, generate_reports, aggregate_timeseries
from spice_ev.util import sanitize

from simba import util

DPI = 300


def open_for_csv(file_path):
    """ Create a file handle to write to.

    :param file_path: Path to new file, overwritten if existing
    :type file_path: string or pathlib.Path
    :return: Function to open file
    :rtype: function
    """
    return open(file_path, "w", newline='', encoding='utf-8')


def generate_gc_power_overview_timeseries(scenario, args):
    """ Generate a csv timeseries with each grid connector's summed up charging station power.

    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    gc_list = list(scenario.components.grid_connectors.keys())

    if not gc_list:
        return

    file_path = args.results_directory / "gc_power_overview_timeseries.csv"
    if vars(args).get("scenario_name"):
        file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)

    with open_for_csv(file_path) as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["time"] + gc_list)
        stations = []
        time_col = getattr(scenario, f"{gc_list[0]}_timeseries")["time"]
        for i in range(len(time_col)):
            time_col[i] = str(time_col[i])
        stations.append(time_col)
        for gc in gc_list:
            stations.append([-x for x in getattr(scenario, f"{gc}_timeseries")["grid supply [kW]"]])
        gc_power_overview = list(map(list, zip(*stations)))
        csv_writer.writerows(gc_power_overview)


def generate_gc_overview(schedule, scenario, args):
    """ Generate a csv file with information regarding electrified stations.

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
    file_path = args.results_directory / "gc_overview.csv"
    if vars(args).get("scenario_name"):
        file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)

    with open_for_csv(file_path) as f:
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


def generate_trips_timeseries_data(schedule):
    """
    Build a trip data structure that can be saved to CSV.

    This should be in the form of a valid trips.csv input file.
    :param schedule: current schedule for the simulation
    :type schedule: simba.Schedule
    :return: trip data
    :rtype: Iterable
    """
    header = [
        # identifier
        "rotation_id", "line",
        # stations
        "departure_name", "departure_time", "arrival_name", "arrival_time",
        # types
        "vehicle_type", "charging_type",
        # consumption (minimal)
        "distance", "consumption"
        # consumption (extended). Not strictly needed.
        # "distance", "temperature", "height_diff", "level_of_loading", "mean_speed",
    ]
    data = [header]
    rotations = schedule.rotations.values()
    # sort rotations naturally by ID
    rotations = sorted(rotations, key=lambda r: [
        # natural sort: split by numbers, then sort numbers by value and chars by lowercase
        int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', r.id)])
    for rotation in rotations:
        for trip in rotation.trips:
            # get trip info from trip or trip.rotation (same name in Trip/Rotation as in CSV)
            row = [vars(trip).get(k, vars(trip.rotation).get(k)) for k in header]
            # special case rotation_id
            row[0] = rotation.id
            data.append(row)
    return data


def generate_plots(schedule, scenario, args):
    """ Save plots as png and pdf.

    Optionally create extended output plots as well.

    :param schedule: Driving schedule for the simulation, used in extended plots.
    :type schedule: simba.Schedule
    :param scenario: Scenario to plot.
    :type scenario: spice_ev.Scenario
    :param args: Configuration. Uses results_directory, show_plots and extended_output_plots.
    :type args: argparse.Namespace
    """
    aggregate_global_results(scenario)
    # disable DEBUG logging from matplotlib
    logging.disable(logging.INFO)
    with plt.ion():  # make plotting temporarily interactive, so plt.show does not block
        plt.clf()
        plot(scenario)
        plt.gcf().set_size_inches(10, 10)
        file_path_png = args.results_directory / "run_overview.png"
        if vars(args).get("scenario_name"):
            file_path_png = file_path_png.with_stem(file_path_png.stem + '_' + args.scenario_name)
        plt.savefig(file_path_png, dpi=300)
        file_path_pdf = args.results_directory / "run_overview.pdf"
        if vars(args).get("scenario_name"):
            file_path_pdf = file_path_pdf.with_stem(file_path_pdf.stem + '_' + args.scenario_name)
        plt.savefig(file_path_pdf, dpi=300)
    if args.show_plots:
        plt.show()

    if args.extended_output_plots:
        plt.clf()
        # create directory for extended plots
        extended_plots_path = args.results_directory.joinpath("extended_plots")
        extended_plots_path.mkdir(parents=True, exist_ok=True)

        plot_blocks_dense(schedule, extended_plots_path)
        plot_vehicle_services(schedule, extended_plots_path)

        plot_consumption_per_rotation_distribution(extended_plots_path, schedule)
        plot_distance_per_rotation_distribution(extended_plots_path, schedule)
        plot_charge_type_distribution(extended_plots_path, scenario, schedule)
        plot_gc_power_timeseries(extended_plots_path, scenario)
        plot_active_rotations(extended_plots_path, scenario, schedule)

    # revert logging override
    logging.disable(logging.NOTSET)


def generate(schedule, scenario, args):
    """ Generates all output files and saves them in the args.results_directory.

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
    file_path = args.results_directory / "vehicle_socs.csv"
    if vars(args).get("scenario_name"):
        file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)
    args.save_soc = file_path
    # bundle station-specific output files in subdirectory
    gc_dir = args.results_directory / "gcs"
    gc_dir.mkdir(exist_ok=True)
    args.save_results = gc_dir / "info.json"
    args.save_timeseries = gc_dir / "ts.csv"
    generate_reports(scenario, vars(args).copy())
    args.save_timeseries = None
    args.save_results = None
    args.save_soc = None

    # generate gc power overview
    generate_gc_power_overview_timeseries(scenario, args)

    # generate gc overview
    generate_gc_overview(schedule, scenario, args)

    # save plots as png and pdf, includes optional extended plots
    generate_plots(schedule, scenario, args)

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

        start_idx = (rotation.departure_time - sim_start_time) // interval
        end_idx = start_idx + ((rotation.arrival_time - rotation.departure_time) // interval)
        if end_idx > scenario.n_intervals:
            # SpiceEV stopped before rotation was fully simulated
            incomplete_rotations.append(id)
            continue

        # get soc timeseries for current rotation
        vehicle_soc = scenario.vehicle_socs[vehicle_id]

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
            "start_time": str(rotation.departure_time),
            "end_time": str(rotation.arrival_time),
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
        # count active rotations for each timestep
        num_active_rotations = count_active_rotations(scenario, schedule)
        data = [["time"] + rotations + ['# active rotations']]
        for i in range(scenario.n_intervals):
            t = sim_start_time + i * scenario.interval
            socs = [str(rotation_socs[k][i]) for k in rotations]
            data.append([str(t)] + socs + [num_active_rotations[i]])

        file_path = args.results_directory / "rotation_socs.csv"
        if vars(args).get("scenario_name"):
            file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)
        write_csv(data, file_path, propagate_errors=args.propagate_mode_errors)

        file_path = args.results_directory / "rotation_summary.csv"
        if vars(args).get("scenario_name"):
            file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)
        with open_for_csv(file_path) as f:
            csv_writer = csv.DictWriter(f, list(rotation_infos[0].keys()))
            csv_writer.writeheader()
            csv_writer.writerows(rotation_infos)

    if vars(args).get('create_trips_in_report', False):
        file_path = args.results_directory / "trips.csv"
        write_csv(generate_trips_timeseries_data(schedule), file_path)

    # summary of used vehicle types and all costs
    if args.cost_calculation:
        file_path = args.results_directory / "summary_vehicles_costs.csv"
        if vars(args).get("scenario_name"):
            file_path = file_path.with_stem(file_path.stem + '_' + args.scenario_name)
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
    """ Write iterable data to CSV file.

    Errors are caught if propagate_errors=False. If data is None, no file is created.

    :param data: data which is written to a file
    :type data: Iterable
    :param file_path: file_path for file creation
    :type file_path: str or Path
    :param propagate_errors: should errors be propagated?
    :type propagate_errors: bool
    :raises Exception: if file cannot be written and propagate_errors=True
    """

    if data is None:
        logging.debug(f"Writing to {file_path} aborted since no data is passed.")
        return
    try:
        with open_for_csv(file_path) as f:
            csv_writer = csv.writer(f)
            for row in data:
                csv_writer.writerow(row)
    except Exception as e:
        logging.warning(f"Writing to {file_path} failed due to {str(e)}")
        if propagate_errors:
            raise


# ##### EXTENDED PLOTTING ##### #

def prepare_histogram(rotations, schedule):
    """ Find suitable number of histogram bins for given rotation values.

    :param rotations: Rotation values to create histogram for. ID -> value
    :type rotations: dict
    :param schedule: Driving schedule
    :type schedule: simba.schedule.Schedule
    :return: histogram bins, labels
    """
    # find suitable step size / number of bins
    min_value = int(min(rotations.values()))
    max_value = int(max(rotations.values()))
    # maximum number of bins (some may be empty): between 1 and 20, optimally half of rotations
    max_num_bins = min(max(len(rotations) / 2, 1), 20)
    steps = [1, 2.5, 5]  # extended to 10, 25, 50, 100, 250, ...
    idx = 0
    mult = 1
    while True:
        step = steps[idx] * mult
        min_bin = (min_value // step) * step
        num_bins = int((max_value - min_bin) // step + 1)
        if num_bins <= max_num_bins:
            # first step with large enough step size / small enough number of bins:
            # use this step size
            if step != 2.5:
                # step size is integer: cast to int for better labels
                step = int(step)
                min_bin = int(min_bin)
            break
        # too many bins: increase step size
        idx = (idx + 1) % len(steps)
        # all steps iterated: append a zero, try again
        mult = mult if idx else mult * 10

    # suitable step size found
    labels = [f"{min_bin + i*step} - {min_bin + (i+1)*step}" for i in range(num_bins)]
    # init bins: track bins for each vehicle type individually
    bins = {v_types: [0]*num_bins for v_types in schedule.vehicle_types}

    # fill bins with rotations
    for rot, value in rotations.items():
        position = int((value - min_bin) // step)
        bins[schedule.rotations[rot].vehicle_type][position] += 1

    return bins, labels


def plot_distance_per_rotation_distribution(extended_plots_path, schedule):
    """Plots the distribution of bus types in distance brackets as a stacked bar chart.

    :param extended_plots_path: directory to save plot to
    :type extended_plots_path: Path
    :param schedule: Driving schedule for the simulation, schedule.rotations are used
    :type schedule: simba.schedule.Schedule
    """
    distances = {rot: schedule.rotations[rot].distance / 1000 for rot in schedule.rotations}
    bins, labels = prepare_histogram(distances, schedule)

    # plot
    fig, ax = plt.subplots()
    bar_bottom = [0] * len(labels)
    for v_type in schedule.vehicle_types:
        ax.bar(labels, bins[v_type], width=0.9, label=v_type, bottom=bar_bottom)
        for i in range(len(labels)):
            bar_bottom[i] += bins[v_type][i]
    ax.set_xlabel('Distance [km]')
    plt.xticks(rotation=30)  # slant labels for better readability
    ax.set_ylabel('Number of rotations')
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.yaxis.grid(True)
    ax.set_title('Distribution of rotation length per vehicle type')
    ax.legend()
    plt.tight_layout()
    plt.savefig(extended_plots_path / "distribution_distance.png", dpi=300)
    plt.close()


def plot_consumption_per_rotation_distribution(extended_plots_path, schedule):
    """Plots the distribution of bus types in consumption brackets as a stacked bar chart.

    :param extended_plots_path: directory to save plot to
    :type extended_plots_path: Path
    :param schedule: Driving schedule for the simulation, schedule.rotations are used
    :type schedule: simba.schedule.Schedule
    """
    consumption = {rot: schedule.rotations[rot].consumption for rot in schedule.rotations}
    bins, labels = prepare_histogram(consumption, schedule)

    # plot
    fig, ax = plt.subplots()
    bar_bottom = [0] * len(labels)
    for v_type in schedule.vehicle_types:
        ax.bar(labels, bins[v_type], width=0.9, label=v_type, bottom=bar_bottom)
        for i in range(len(labels)):
            bar_bottom[i] += bins[v_type][i]
    ax.set_xlabel('Energy consumption [kWh]')
    plt.xticks(rotation=30)
    ax.set_ylabel('Number of rotations')
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.yaxis.grid(True)
    ax.set_title('Distribution of energy consumption of rotations per vehicle type')
    ax.legend()
    plt.tight_layout()
    plt.savefig(extended_plots_path / "distribution_consumption", dpi=300)
    plt.close()


def plot_charge_type_distribution(extended_plots_path, scenario, schedule):
    """Plots the number of rotations of each charging type in a bar chart.

    :param extended_plots_path: directory to save plot to
    :type extended_plots_path: Path
    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param schedule: Driving schedule for the simulation. schedule.rotations are used
    :type schedule: simba.schedule.Schedule
    """
    # count charging types (also with regard to negative rotations)
    charging_types = {'oppb': 0, 'oppb_neg': 0, 'depb': 0, 'depb_neg': 0}
    negative_rotations = schedule.get_negative_rotations(scenario)
    for rot in schedule.rotations:
        ct = schedule.rotations[rot].charging_type
        if rot in negative_rotations:
            ct += '_neg'
        try:
            charging_types[ct] += 1
        except KeyError:
            logging.error(f"Rotation {rot}: unknown charging type: '{ct}'")

    # plot
    fig, ax = plt.subplots()
    bars1 = ax.bar(
        ["Opportunity", "Depot"],
        [charging_types["oppb"], charging_types["depb"]],
    )
    bars2 = ax.bar(
        ["Opportunity", "Depot"],
        [charging_types["oppb_neg"], charging_types["depb_neg"]],
        bottom=[charging_types["oppb"], charging_types["depb"]],
    )
    # create labels with counts
    # create empty labels for empty bins
    labels = [
        [
            f"{ct}{suffix}: {charging_types[f'{ct}{suffix}']}"
            if charging_types[f'{ct}{suffix}'] > 0 else ""
            for ct in ['oppb', 'depb']
        ] for suffix in ['', '_neg']
    ]
    ax.bar_label(bars1, labels=labels[0], label_type='center')  # oppb, depb
    ax.bar_label(bars2, labels=labels[1], label_type='center')  # oppb_neg, depb_neg

    ax.set_xlabel("Charging type")
    ax.set_ylabel("Number of rotations")
    ax.yaxis.grid(True)
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.legend(["successful rotations", "negative rotations"])
    ax.set_title("Feasibility of rotations per charging type")
    plt.savefig(extended_plots_path / "charge_types", dpi=300)
    plt.close()


def plot_gc_power_timeseries(extended_plots_path, scenario):
    """Plots the different loads (total, feedin, external) of all grid connectors.

    :param extended_plots_path: directory to save plot to
    :type extended_plots_path: Path
    :param scenario: Provides the data for the grid connectors over time.
    :type scenario: spice_ev.Scenario
    """
    gc_list = list(scenario.components.grid_connectors.keys())

    for gc in gc_list:
        fig, ax = plt.subplots()

        agg_ts = aggregate_timeseries(scenario, gc)
        headers = [
            "grid supply [kW]",
            "fixed load [kW]",
            "local generation [kW]",
            "sum CS power [kW]",
            "battery power [kW]",
            "bat. stored energy [kWh]",
        ]

        has_battery_column = False

        # find time column
        time_index = agg_ts["header"].index("time")
        time_values = [row[time_index] for row in agg_ts["timeseries"]]

        for header_index, header in enumerate(headers):
            try:
                # try to find column with current header
                idx = agg_ts["header"].index(header)
                header_values = [row[idx] for row in agg_ts["timeseries"]]
            except ValueError:
                # column does not exist
                continue

            if header == "bat. stored energy [kWh]":
                has_battery_column = True
                # special plot for battery: same subplot, different y-axis
                ax2 = ax.twinx()
                ax2.set_ylabel("stored battery energy [kWh]")
                # get next color from color cycle (just plotting would start with first color)
                next_color = plt.rcParams['axes.prop_cycle'].by_key()["color"][header_index]
                ax2.plot(
                    time_values, header_values,
                    label=header, c=next_color, linestyle="dashdot")
                ax2.legend()
                fig.set_size_inches(8, 4.8)
            else:
                # normal (non-battery) plot
                ax.plot(time_values, header_values, label=header)

        if has_battery_column:
            # align y axis so that 0 is shared
            # (limits not necessary, as power and energy can't be compared directly)
            ax1_ylims = ax.axes.get_ylim()
            ax1_yratio = ax1_ylims[0] / ax1_ylims[1]
            ax2_ylims = ax2.axes.get_ylim()
            ax2_yratio = ax2_ylims[0] / ax2_ylims[1]
            if ax1_yratio < ax2_yratio:
                ax2.set_ylim(bottom=ax2_ylims[1]*ax1_yratio)
            else:
                ax.set_ylim(bottom=ax1_ylims[1]*ax2_yratio)
            plt.tight_layout()

        ax.legend()
        plt.xticks(rotation=30)
        ax.set_ylabel("Power [kW]")
        ax.set_title(f"Power: {gc}")
        ax.grid(color='gray', linestyle='-')

        # xaxis are datetime strings
        ax.set_xlim(time_values[0], time_values[-1])
        # ax.tick_params(axis='x', rotation=30)
        plt.tight_layout()
        plt.savefig(extended_plots_path / f"{sanitize(gc)}_power_overview.png", dpi=300)
        plt.close(fig)


def ColorGenerator(vehicle_types):
    # generates color according to vehicle_type and charging type of rotation
    colors = util.cycling_generator(plt.rcParams['axes.prop_cycle'].by_key()["color"])
    color_per_vt = {vt: next(colors) for vt in vehicle_types}
    color = None
    while True:
        rotation = yield color
        color = color_per_vt[rotation.vehicle_type]
        rgb = matplotlib.colors.to_rgb(color)
        hsv = matplotlib.colors.rgb_to_hsv(rgb)
        value = hsv[-1]
        if rotation.charging_type == "depb":
            if value >= 0.5:
                value -= 0.3
        elif rotation.charging_type == "oppb":
            if value < 0.5:
                value += 0.3
        else:
            raise NotImplementedError
        hsv[-1] = value
        color = matplotlib.colors.hsv_to_rgb(hsv)


def plot_vehicle_services(schedule, output_path):
    """Plots the rotations serviced by the same vehicle

    :param schedule: Provides the schedule data.
    :type schedule; simba.schedule.Schedule
    :param output_path: Path to the output folder
    :type output_path: pathlib.Path
    """

    # find depots for every rotation
    all_rotations = schedule.rotations
    rotations_per_depot = {r.departure_name: [] for r in all_rotations.values()}
    for rotation in all_rotations.values():
        rotations_per_depot[rotation.departure_name].append(rotation)

    def VehicleIdRowGenerator():
        # generate row ids by using the vehicle_id of the rotation
        vehicle_id = None
        while True:
            rotation = yield vehicle_id
            vehicle_id = rotation.vehicle_id

    rotations_per_depot["All_Depots"] = all_rotations.values()

    output_path_folder = output_path / "vehicle_services"
    output_path_folder.mkdir(parents=True, exist_ok=True)
    for depot, rotations in rotations_per_depot.items():
        # create instances of the generators
        row_generator = VehicleIdRowGenerator()
        color_generator = ColorGenerator(schedule.vehicle_types)

        sorted_rotations = list(
            sorted(rotations, key=lambda x: (x.vehicle_type, x.charging_type, x.departure_time)))
        fig, ax = create_plot_blocks(sorted_rotations, color_generator, row_generator)
        ax.set_ylabel("Vehicle ID")
        # Vehicle ids need small fontsize to fit
        ax.tick_params(axis='y', labelsize=6)
        # add legend
        handles = []
        vts = {f"{rotation.vehicle_type}_{rotation.charging_type}": rotation
               for rotation in sorted_rotations}
        for key, rot in vts.items():
            handles.append(Patch(color=color_generator.send(rot), label=key))
        # Position legend at the top outside of the plot
        ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1),
                  ncol=len(handles)//2+1, prop={"size": 7})
        fig.tight_layout()
        # PDF so Block names stay readable
        fig.savefig(output_path_folder / f"{sanitize(depot)}_vehicle_services.pdf")
        fig.savefig(output_path_folder / f"{sanitize(depot)}_vehicle_services.png", dpi=300)
        plt.close(fig)


def plot_blocks_dense(schedule, output_path):
    """Plots the different loads (total, feedin, external) of all grid connectors.

    :param schedule: Provides the schedule data.
    :type schedule; simba.schedule.Schedule
    :param output_path: Path to the output folder
    :type output_path: pathlib.Path
    """
    # find depots for every rotation
    all_rotations = schedule.rotations
    rotations_per_depot = {r.departure_name: [] for r in all_rotations.values()}
    for rotation in all_rotations.values():
        rotations_per_depot[rotation.departure_name].append(rotation)

    def DenseRowGenerator():
        # Generates row id by trying to find lowes row with a rotation which arrived
        arrival_times = []
        row_nr = None
        while True:
            rotation = yield row_nr
            for i in range(len(arrival_times)):
                at = arrival_times[i]
                if at <= rotation.departure_time:
                    arrival_times[i] = rotation.arrival_time
                    row_nr = i
                    break
            else:
                arrival_times.append(rotation.arrival_time)
                row_nr = len(arrival_times) - 1

    rotations_per_depot["All_Depots"] = all_rotations.values()
    output_path_folder = output_path / "block_distribution"
    output_path_folder.mkdir(parents=True, exist_ok=True)
    for depot, rotations in rotations_per_depot.items():
        # create instances of the generators
        row_generator = DenseRowGenerator()
        color_generator = ColorGenerator(schedule.vehicle_types)
        # sort by departure_time and longest duration
        sorted_rotations = list(
            sorted(rotations,
                   key=lambda r: (r.departure_time, -(r.departure_time - r.arrival_time))))
        fig, ax = create_plot_blocks(sorted_rotations, color_generator, row_generator)
        ax.set_ylabel("Block")
        ax.yaxis.get_major_locator().set_params(integer=True)
        # add legend
        handles = []
        vts = {f"{rotation.vehicle_type}_{rotation.charging_type}": rotation
               for rotation in sorted_rotations}
        for key, rot in vts.items():
            handles.append(Patch(color=color_generator.send(rot), label=key))
        # Position legend at the top outside of the plot
        ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1),
                  ncol=len(handles)//2+1, prop={"size": 7})

        fig.tight_layout()
        # PDF so Block names stay readable
        fig.savefig(output_path_folder / f"{sanitize(depot)}_block_distribution.pdf")
        fig.savefig(output_path_folder / f"{sanitize(depot)}_block_distribution.png", dpi=300)
        plt.close(fig)


def create_plot_blocks(sorted_rotations, color_generator, row_generator):
    # initialize row_generator by yielding None
    next(row_generator)
    next(color_generator)
    y_size = len(sorted_rotations) * 0.02 + 1.5
    fig, ax = plt.subplots(figsize=(7, y_size))
    for rotation in sorted_rotations:
        row_nr = row_generator.send(rotation)
        width = rotation.arrival_time - rotation.departure_time
        artist = ax.barh([row_nr], width=[width], height=0.8, left=rotation.departure_time,
                         label=rotation.id, color=color_generator.send(rotation))
        ax.bar_label(artist, labels=[rotation.id], label_type='center', fontsize=2)

    # Large heights create too much margin with default value of margins: Reduce value to 0.01
    ax.axes.margins(y=0.01)
    ax.grid(axis='x')
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%D %H:%M'))
    ax.set_xlim(min(r.departure_time for r in sorted_rotations) - datetime.timedelta(minutes=30),
                max(r.arrival_time for r in sorted_rotations) + datetime.timedelta(minutes=30))
    ax.tick_params(axis='x', rotation=30)
    return fig, ax


def count_active_rotations(scenario, schedule):
    num_active_rotations = [0] * scenario.n_intervals
    for rotation in schedule.rotations.values():
        ts_start = int((rotation.departure_time - scenario.start_time) / scenario.interval)
        ts_end = ceil((rotation.arrival_time - scenario.start_time) / scenario.interval)
        # ignore rotations after scenario end
        ts_end = min(ts_end, scenario.n_intervals)
        for ts_idx in range(ts_start, ts_end):
            num_active_rotations[ts_idx] += 1
    return num_active_rotations


def plot_active_rotations(extended_plots_path, scenario, schedule):
    """Generate a plot with number of active rotations over time.

    :param extended_plots_path: directory to save plot to
    :type extended_plots_path: Path
    :param scenario: Provides the data for the grid connectors over time.
    :type scenario: spice_ev.Scenario
    :param schedule: Driving schedule for the simulation. schedule.rotations are used
    :type schedule: simba.schedule.Schedule
    """
    ts = [scenario.start_time + scenario.interval * i for i in range(scenario.n_intervals)]
    num_active_rotations = count_active_rotations(scenario, schedule)
    plt.plot(ts, num_active_rotations)
    fig = plt.gcf()
    fig.set_size_inches(8, 4.8)
    ax = plt.gca()
    ax.xaxis_date()
    ax.set_xlim(ts[0], ts[-1])
    plt.xticks(rotation=30)
    plt.ylabel("Number of active rotations")
    ax.yaxis.get_major_locator().set_params(integer=True)
    plt.grid(axis="y")
    plt.title("Active Rotations")
    plt.tight_layout()
    plt.savefig(extended_plots_path / "active_rotations", dpi=300)
    plt.close()
