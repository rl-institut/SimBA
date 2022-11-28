""" Module to generate meaningful output files and/or figures to describe simulation process
"""
import csv
import datetime
import warnings
import json
import matplotlib.pyplot as plt

from src.report import aggregate_timeseries, aggregate_local_results, aggregate_global_results


def generate_vehicle_socs(scenario, args):

    sim_start_time = scenario.start_time
    v_list = []
    for v in scenario.vehicle_socs:
        v_list.append(v)

    with open(args.output_directory / "vehicle_socs.csv", "w+", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["timestep", "time", ] + v_list)
        for i, row in enumerate(zip(*scenario.vehicle_socs.values())):
            t = sim_start_time + i * scenario.interval
            csv_writer.writerow((i, t,) + row)


def generate_station_name_csv(scenario, args):
    for gc in scenario.constants.grid_connectors.keys():
        gc_info = aggregate_timeseries(scenario, gc)
        with open(args.output_directory / f"simulation_spiceEV_{gc}.csv", "w+", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(gc_info["header"])
            for elem in gc_info["timeseries"]:
                csv_writer.writerow(elem)


def generate_station_name_json(scenario, args):
    for gc in scenario.constants.grid_connectors.keys():
        gc_info = aggregate_local_results(scenario, gc)
        with open(args.output_directory / f"simulation_spiceEV_{gc}.json", 'w') as f:
            json.dump(gc_info, f, indent=2)


def generate_cs_power_overview(scenario, args):
    sim_start_time = scenario.start_time
    gc_list = []
    for gc in scenario.constants.grid_connectors.keys():
        gc_list.append(gc)

    with open(args.output_directory / "cs_power_overview.csv", "w+", newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["timestep", "time"] + gc_list)
        for i in range(scenario.n_intervals):
            t = sim_start_time + i * scenario.interval
            row = []
            for gc in gc_list:
                gc_info = aggregate_timeseries(scenario, gc)
                cs_power = gc_info['timeseries']
                sum_cs_power = cs_power[i][9]
                row.append(sum_cs_power)

            csv_writer.writerow((i, t, ) + tuple(row))


def plot(scenario, args):
    print('Done. Create plots...')

    xlabels = []
    for r in scenario.results:
        xlabels.append(r['current_time'])

    # batteries
    if scenario.batteryLevels:
        plots_top_row = 3
        ax = plt.subplot(2, plots_top_row, 3)
        ax.set_title('Batteries')
        ax.set(ylabel='Stored power in kWh')
        for name, values in scenario.batteryLevels.items():
            ax.plot(xlabels, values, label=name)
        ax.legend()
    else:
        plots_top_row = 2

    # vehicles
    ax = plt.subplot(2, plots_top_row, 1)
    ax.set_title('Vehicles')
    ax.set(ylabel='SoC')
    lines = ax.plot(xlabels, scenario.socs)
    # reset color cycle, so lines have same color
    ax.set_prop_cycle(None)

    ax.plot(xlabels, scenario.disconnect, '--')
    if len(scenario.constants.vehicles) <= 10:
        ax.legend(lines, sorted(scenario.constants.vehicles.keys()))

    # charging stations
    ax = plt.subplot(2, plots_top_row, 2)
    ax.set_title('Charging Stations')
    ax.set(ylabel='Power in kW')
    lines = ax.step(xlabels, scenario.sum_cs, where='post')
    if len(scenario.constants.charging_stations) <= 10:
        ax.legend(lines, sorted(scenario.constants.charging_stations.keys()))

    # total power
    ax = plt.subplot(2, 2, 3)
    ax.step(xlabels, list([sum(cs) for cs in scenario.sum_cs]), label="CS", where='post')
    gc_ids = scenario.constants.grid_connectors.keys()
    for gcID in gc_ids:
        for name, values in scenario.loads[gcID].items():
            ax.step(xlabels, values, label=name, where='post')
    # draw schedule
    if scenario.strat.uses_window:
        for gcID, schedule in scenario.gcWindowSchedule.items():
            if all(s is not None for s in schedule):
                # schedule exists
                window_values = [v * int(max(scenario.totalLoad[gcID])) for v in schedule]
                ax.step(xlabels, window_values, label="window {}".format(gcID),
                        linestyle='--', where='post')
    if scenario.strat.uses_schedule:
        for gcID, schedule in scenario.gcPowerSchedule.items():
            if any(s is not None for s in schedule):
                ax.step(xlabels, schedule, label="Schedule {}".format(gcID), where='post')

    ax.step(xlabels, scenario.all_totalLoad, label="Total", where='post')
    ax.set_title('Power')
    ax.set(ylabel='Power in kW')
    ax.legend()
    ax.xaxis_date()  # xaxis are datetime objects

    # price
    ax = plt.subplot(2, 2, 4)
    prices = list(zip(*scenario.prices.values()))
    lines = ax.step(xlabels, prices, where='post')
    ax.set_title('Price for 1 kWh')
    ax.set(ylabel='€')
    if len(gc_ids) <= 10:
        ax.legend(lines, sorted(gc_ids))

    # figure title
    fig = plt.gcf()
    fig.suptitle('Strategy: {}'.format(scenario.strat.description), fontweight='bold')

    # fig.autofmt_xdate()  # rotate xaxis labels (dates) to fit
    # autofmt removes some axis labels, so rotate by hand:
    for ax in fig.get_axes():
        ax.set_xlim(scenario.start_time, scenario.stop_time)
        plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # set size of figure and save it
    plt.gcf().set_size_inches(10, 10)
    plt.savefig(args.output_directory / "run_overview.png")
    plt.savefig(args.output_directory / "run_overview.pdf")
    # plt.show()


def generate(schedule, scenario, args):

    # generate csv out of vehicle's socs
    generate_vehicle_socs(scenario, args)

    # generate csv and json for all stations
    generate_station_name_csv(scenario, args)
    generate_station_name_json(scenario, args)

    # generate cs power overview
    generate_cs_power_overview(scenario, args)

    # save plots as png and pdf
    aggregate_global_results(scenario)
    plot(scenario, args)

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
