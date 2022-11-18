import os
import warnings

from calculate_costs import calculate_costs as calc_costs_spice_ev


def calculate_costs(c_params, scenario, schedule, args):
    """ Calculates annual costs of all necessary vehicles and infrastructure.

    :param c_params: Cost of various infrastructure components.
    :type c_params: dict
    :param scenario: Information about the simulated scenario and all used parameters.
    :type scenario: object
    :param schedule: Information about the whole bus schedule and all used parameters.
    :type schedule: object
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    # initialize dictionary with all costs
    costs = {cost: 0 for cost in ["c_vehicles", "c_vehicles_annual", "c_gcs", "c_gcs_annual",
                                  "c_cs", "c_cs_annual", "c_garage_cs", "c_garage_workstations",
                                  "c_garage", "c_garage_annual", "c_invest", "c_invest_annual",
                                  "c_maint_infrastructure_annual", "c_maint_vehicles_annual",
                                  "c_maint_annual", "c_el_procurement_annual",
                                  "c_el_power_price_annual", "c_el_energy_price_annual",
                                  "c_el_taxes_annual", "c_el_feed_in_remuneration_annual",
                                  "c_el_annual"]}

    # INVESTMENT COSTS #

    # VEHICLES
    v_types = schedule.scenario["constants"]["vehicle_types"]
    for v_type, v_keys in v_types.items():
        if schedule.vehicle_type_counts[v_type] > 0:
            try:
                costs_vehicle = c_params["vehicles"][v_type]["capex"]
            except KeyError:
                warnings.warn("No capex defined for vehicle type " + v_type +
                              ". Unable to calculate investment costs for this vehicle type.")
                continue
            # sum up cost of vehicles and their batteries, depending on how often the battery
            # has to be replaced in the lifetime of the vehicles
            c_vehicles_vt = (schedule.vehicle_type_counts[v_type] *
                             (costs_vehicle + (c_params["vehicles"][v_type]["lifetime"] //
                                               c_params["batteries"]["lifetime_battery"]) *
                             v_keys["capacity"] * c_params["batteries"]["cost_per_kWh"]))
            costs["c_vehicles"] += c_vehicles_vt
            # calculate annual cost of vehicles of this type, depending on their lifetime
            costs["c_vehicles_annual"] += c_vehicles_vt / c_params["vehicles"][v_type]["lifetime"]
    # GRID CONNECTION POINTS
    gcs = schedule.scenario["constants"]["grid_connectors"]
    for gcID, gc_keys in gcs.items():
        try:
            distance_transformer = schedule.stations[gcID]["distance_transformer"]
        except KeyError:
            distance_transformer = c_params["gc"]["default_distance"]
        c_gc = (c_params["gc"]["building_cost_subsidy_per_kW"] * gc_keys["max_power"] +
                c_params["gc"]["capex_gc_fix"] +
                c_params["gc"]["capex_gc_per_meter"] * distance_transformer)
        c_transformer = (c_params["gc"]["capex_transformer_fix"] +
                         c_params["gc"]["capex_transformer_per_kW"] * gc_keys["max_power"])
        costs["c_gcs"] += c_gc + c_transformer
        # calculate annual costs of grid connectors, depending the lifetime of gc and transformer
        costs["c_gcs_annual"] += (c_gc / c_params["gc"]["lifetime_gc"] +
                                  c_transformer / c_params["gc"]["lifetime_transformer"])

    # CHARGING INFRASTRUCTURE
    cs = schedule.scenario["constants"]["charging_stations"]
    for csID in cs.values():
        if csID["type"] == "deps":
            costs["c_cs"] += c_params["cs"]["capex_deps_per_kW"] * csID["max_power"]
        elif csID["type"] == "opps":
            costs["c_cs"] += c_params["cs"]["capex_opps_per_kW"] * csID["max_power"]
    # calculate annual cost of charging stations, depending on their lifetime
    costs["c_cs_annual"] = costs["c_cs"] / c_params["cs"]["lifetime_cs"]

    # GARAGE
    costs["c_garage_cs"] = (c_params["garage"]["n_charging_stations"] *
                            c_params["garage"]["power_cs"] * c_params["cs"]["capex_deps_per_kW"])
    costs["c_garage_workstations"] = -(-sum(schedule.vehicle_type_counts.values()) //
                                       c_params["garage"]["vehicles_per_workstation"] *
                                       c_params["garage"]["cost_per_workstation"])
    costs["c_garage"] = costs["c_garage_cs"] + costs["c_garage_workstations"]
    costs["c_garage_annual"] = (costs["c_garage_cs"] / c_params["cs"]["lifetime_cs"] +
                                costs["c_garage_workstations"] /
                                c_params["garage"]["lifetime_workstations"])
    # MAINTENANCE
    costs["c_maint_infrastructure_annual"] = (costs["c_cs"] * c_params["cs"]["c_maint_cs_per_year"]
                                              + c_transformer *
                                              c_params["gc"]["c_maint_transformer_per_year"])
    # calculate (ceil) number of days in scenario
    drive_days = -(-(schedule.scenario["scenario"]["n_intervals"] *
                     schedule.scenario["scenario"]["interval"]) // (24 * 60))
    for rot in schedule.rotations.values():
        v_type_rot = f"{rot.vehicle_type}_{rot.charging_type}"
        try:
            costs["c_maint_vehicles_annual"] += (rot.distance / 1000 / drive_days * 365) * \
                                                c_params["vehicles"][v_type_rot]["c_maint_per_km"]
        except KeyError:
            warnings.warn("No maintenance costs defined for vehicle type " +
                          c_params["vehicles"][v_type_rot] +
                          ". Unable to calculate maintenance costs for this vehicle type.")
    costs["c_maint_annual"] = (costs["c_maint_infrastructure_annual"] +
                               costs["c_maint_vehicles_annual"])
    costs["c_invest"] = costs["c_vehicles"] + costs["c_cs"] + costs["c_gcs"] + costs["c_garage"]
    costs["c_invest_annual"] = (costs["c_vehicles_annual"] + costs["c_cs_annual"] +
                                costs["c_gcs_annual"] + costs["c_garage_annual"])

    # ELECTRICITY COSTS #

    for gcID, gc in scenario.constants.grid_connectors.items():
        pv = sum([pv.nominal_power for pv in scenario.constants.photovoltaics.values()
                  if pv.parent == gcID])
        timeseries = vars(scenario).get(f"{gcID}_timeseries")

        # Todo: Decide, if costs are saved to json here or later in report.py.
        #  If so, the following three lines can be removed
        # add GC name to results file for cost calculation
        file_name, ext = os.path.splitext(args.save_results)
        save_results = f"{file_name}_{gcID}{ext}"

        # Todo: Decide, if currently unnecessary params should be given to calc_costs_spice_ev
        # calculate costs for electricity
        costs_electricity = calc_costs_spice_ev(
            strategy=args.strategy,
            voltage_level=gc.voltage_level,
            interval=scenario.interval,
            timestamps_list=timeseries.get("time"),
            power_grid_supply_list=timeseries.get("grid power [kW]"),
            price_list=timeseries.get("price [EUR/kWh]"),
            power_fix_load_list=timeseries.get("ext.load [kW]"),
            charging_signal_list=timeseries.get("window"),
            core_standing_time_dict=scenario.core_standing_time,
            price_sheet_json=args.cost_parameters_file,
            results_json=save_results,
            power_pv_nominal=pv,
        )
        # ToDo: Decide if gc-specific costs should be added to scenario object to use in report.py
        # setattr(scenario.constants.grid_connectors[gcID], "costs_electricity", costs_electricity)
        costs["c_el_procurement_annual"] += costs_electricity['power_procurement_per_year']
        costs["c_el_power_price_annual"] += costs_electricity['capacity_costs_eur']
        costs["c_el_energy_price_annual"] += costs_electricity['commodity_costs_eur_per_year']
        costs["c_el_taxes_annual"] += costs_electricity['levies_fees_and_taxes_per_year']
        costs["c_el_feed_in_remuneration_annual"] += costs_electricity[
            'feed_in_remuneration_per_year']
        costs["c_el_annual"] += costs_electricity['total_costs_per_year']

    # round all costs to ct
    for key, v in costs.items():
        costs[key] = round(v, 2)

    print(f"\n"
          f"Total costs: \n"
          f"Investment cost: {costs['c_invest']} €. "
          f"Annual investment costs: {costs['c_invest_annual']} €/a. "
          f"Annual maintenance costs: {costs['c_maint_annual']} €/a. "
          f"Annual costs for electricity: {costs['c_el_annual']} €/a.")

    setattr(scenario, "costs", costs)
