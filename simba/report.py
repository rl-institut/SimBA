""" Module to generate meaningful output files and/or figures to describe simulation process. """
import csv
import datetime
import logging
from typing import Iterable

import matplotlib.pyplot as plt
from spice_ev.report import aggregate_global_results, plot, generate_reports, aggregate_timeseries


def open_for_csv(filepath):
    """ Create a file handle to write to.

    :param filepath: Path to new file, overwritten if existing
    :type filepath: string or pathlib.Path
    :return: Function to open file
    :rtype: function
    """
    return open(filepath, "w", newline='', encoding='utf-8')


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

    with open_for_csv(args.results_directory / "gc_power_overview_timeseries.csv") as f:
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

    with open_for_csv(args.results_directory / "gc_overview.csv") as f:
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
    """ Save plots as png and pdf.

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
    """ Generates all output files/ plots and saves them in the args.results_directory.

    :param schedule: Driving schedule for the simulation.
    :type schedule: simba.Schedule
    :param scenario: Scenario for with to generate timeseries.
    :type scenario: spice_ev.Scenario
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    :raises Exception: if cost calculation cannot be written to file and
        args.propagate_mode_errors=True
    """

    # generate if needed extended output plots
    if args.extended_output_plots:
        bus_type_distribution_consumption_rotation(args, schedule)
        charge_type_proportion(args, scenario, schedule)
        gc_power_time_overview(args, scenario)
        active_rotations(args, scenario, schedule)

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

        # add active rotations column to rotation_socs.csv
        active_vehicles = [len(scenario.components.vehicles)] * scenario.n_intervals
        depot_stations = [station for station in scenario.components.grid_connectors
                          if schedule.stations[station]["type"] == "deps"]
        for station_name in depot_stations:
            station_ts = scenario.connChargeByTS[station_name]
            for index, step in enumerate(station_ts):
                active_vehicles[index] -= len(step)
        data[0].append("# active_rotations")
        for i in range(len(data)):
            data[i].append(active_vehicles[i-1])

        write_csv(data, args.results_directory / "rotation_socs.csv",
                  propagate_errors=args.propagate_mode_errors)

        with open_for_csv(args.results_directory / "rotation_summary.csv") as f:
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


def bus_type_distribution_consumption_rotation(args, schedule):
    """Plots the distribution of bus types in consumption brackets as a stacked bar chart.

    :param args: Configuration arguments from cfg file. args.results_directory is used
    :type args: argparse.Namespace
    :param schedule: Driving schedule for the simulation. schedule.rotations are used
    :type schedule: eBus-Toolbox.Schedule
    """

    step = 50
    # get bin_number
    max_con = int(max([schedule.rotations[rot].consumption for rot in schedule.rotations]))
    if max_con % step < step/2:
        bin_number = ((max_con // step) * step)
    else:
        bin_number = (max_con // step) * step + step
    labels = [f"{i - step} - {i}" for i in range(step, bin_number+step, step) if i > 0]
    bins = {v_types: [0 for _ in range(bin_number // step)] for v_types in schedule.vehicle_types}

    # fill bins with rotations
    for rot in schedule.rotations:
        for v_type in schedule.vehicle_types:
            if schedule.rotations[rot].vehicle_type == v_type:
                position = int(schedule.rotations[rot].consumption // step)
                if position >= bin_number/step:
                    position -= bin_number // step - 1
                bins[v_type][position] += 1     # index out of range
                break
    # plot
    fig, ax = plt.subplots()
    bar_bottom = [0 for _ in range(bin_number//step)]
    for v_type in schedule.vehicle_types:
        ax.bar(labels, bins[v_type], width=0.9, label=v_type, bottom=bar_bottom)
        for i in range(bin_number//step):
            bar_bottom[i] += bins[v_type][i]
    ax.set_xlabel('Energieverbrauch in kWh')
    ax.set_ylabel('Anzahl der Umläufe')
    ax.set_title('Verteilung der Bustypen über den Energieverbrauch und den Umläufen')
    ax.legend()
    fig.autofmt_xdate()
    ax.yaxis.grid(True)
    plt.savefig(args.results_directory / "distribution_bustypes_consumption_rotations")
    plt.close()


def charge_type_proportion(args, scenario, schedule):
    """Plots the absolute number of rotations distributed by charging types on a bar chart.

    :param args: Configuration arguments from cfg file. args.results_directory is used
    :type args: argparse.Namespace
    :param schedule: Driving schedule for the simulation. schedule.rotations are used
    :type schedule: eBus-Toolbox.Schedule
    """
    # get plotting data
    charging_types = {'oppb': 0, 'oppb_neg': 0, 'depb': 0, 'depb_neg': 0}
    negative_rotations = schedule.get_negative_rotations(scenario)
    for rot in schedule.rotations:
        if schedule.rotations[rot].charging_type == 'oppb':
            if rot in negative_rotations:
                charging_types['oppb_neg'] += 1
            else:
                charging_types['oppb'] += 1
        elif schedule.rotations[rot].charging_type == 'depb':
            if rot in negative_rotations:
                charging_types['depb_neg'] += 1
            else:
                charging_types['depb'] += 1
        else:
            print("Unknown charging type: ", schedule.rotations[rot].charging_type)
    # plot
    fig, ax = plt.subplots()
    bars1 = ax.bar(
        ["Gelegenheitslader", "Depotlader"],
        [charging_types["oppb"], charging_types["depb"]],
        color=["#6495ED", "#6495ED"],
    )
    bars2 = ax.bar(
        ["Gelegenheitslader", "Depotlader"],
        [charging_types["oppb_neg"], charging_types["depb_neg"]],
        color=["#66CDAA", "#66CDAA"],
        bottom=[charging_types["oppb"], charging_types["depb"]],
    )
    ax.bar_label(bars1, labels=[
        f"{charging_types['oppb']} oppb",
        f"{charging_types['depb']} depb",
    ])
    ax.bar_label(bars2, labels=[
        f"{charging_types['oppb_neg']} oppb_neg",
        f"{charging_types['depb_neg']} depb_neg",
    ])

    ax.set_xlabel("Ladetyp")
    ax.set_ylabel("Umläufe")
    ax.set_title("Verteilung von Gelegenheitslader, Depotlader")

    ax.yaxis.grid(True)
    plt.savefig(args.results_directory / "charge_type_proportion")
    plt.close()


def gc_power_time_overview_example(args, scenario):

    gc_list = list(scenario.constants.grid_connectors.keys())
    for gc in gc_list:
        # data
        ts = getattr(scenario, f"{gc}_timeseries")
        time = ts["time"]
        total = ts["grid power [kW]"]
        feed_in = ts["feed-in [kW]"]
        ext_load = ts["ext.load [kW]"]
        cs = ts["sum CS power"]

        # plot
        plt.plot(time, total, label="total")
        plt.plot(time, feed_in, label="feed_in")
        plt.plot(time, ext_load, label="ext_load")
        plt.plot(time, cs, label="CS")
        plt.legend()
        plt.xticks(rotation=45)

        plt.savefig(args.results_directory / f"{gc}_power_time_overview")
        plt.clf()


def gc_power_time_overview(args, scenario):
    """Plots the different loads (total, feedin, external) of all grid connectors.

    :param args: Configuration arguments from cfg file. args.results_directory is used
    :type args: argparse.Namespace
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
        time_index = agg_ts["header"].index("time")
        time_values = [row[time_index] for row in agg_ts["timeseries"]]

        for header in headers:
            try:
                header_index = agg_ts["header"].index(header)
                header_values = [row[header_index] for row in agg_ts["timeseries"]]

                if header == "bat. stored energy [kWh]":
                    ax2 = ax.twinx()  # Instantiate a second Axes that shares the same x-axis
                    ax2.set_ylabel("Power in kWh")  # We already handled the x-label with ax1
                    ax2.plot(time_values, header_values, label=header)
                    ax2.legend()
                    ax2.tick_params(axis='x', rotation=30)
                else:
                    ax.plot(time_values, header_values, label=header)  # Use ax instead of plt
            except ValueError:
                continue

        ax.legend()
        ax.tick_params(axis='x', rotation=30)
        ax.set_ylabel("Power in kW")
        ax.set_title(f"Power: {gc}")

        gc_cleaned = gc.replace("/", "").replace(".", "")
        plt.savefig(args.results_directory / f"{gc_cleaned}_power_time_overview.png")
        plt.close(fig)


def active_rotations(args, scenario, schedule):
    """Generate a plot where the number of active rotations is shown.

    :param args: Configuration arguments from cfg file. args.results_directory is used
    :type args: argparse.Namespace
    :param schedule: Driving schedule for the simulation. schedule.rotations are used
    :type schedule: eBus-Toolbox.Schedule
    :param scenario: Provides the data for the grid connectors over time.
    :type scenario: spice_ev.Scenario
    """
    ts = [
        scenario.start_time if i == 0 else scenario.start_time + scenario.interval * i
        for i in range(scenario.n_intervals)
    ]
    active_vehicles = [len(scenario.components.vehicles)] * scenario.n_intervals
    depot_stations = [
        station for station in scenario.components.grid_connectors
        if schedule.stations[station]["type"] == "deps"
    ]
    for station_name in depot_stations:
        station_ts = scenario.connChargeByTS[station_name]
        for index, step in enumerate(station_ts):
            active_vehicles[index] -= len(step)

    plt.plot(ts, active_vehicles)
    plt.ylabel("Number of active Vehicles")
    plt.title("Active Rotations")
    plt.xticks(rotation=30)
    plt.grid(axis="y")
    plt.savefig(args.results_directory / "active_rotations")
    plt.close()


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
