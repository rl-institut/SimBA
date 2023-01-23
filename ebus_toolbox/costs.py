import warnings

from calculate_costs import calculate_costs as calc_costs_spice_ev


def calculate_costs(c_params, scenario, schedule, args):
    """ Calculates annual costs of all necessary vehicles and infrastructure.

    :param c_params: Cost of various infrastructure components.
    :type c_params: dict
    :param scenario: Information about the simulated scenario and all used parameters.
    :type scenario: Scenario
    :param schedule: Information about the whole bus schedule and all used parameters.
    :type schedule: Schedule
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    """

    # initialize dictionary with all costs
    costs = {cost: 0 for cost in [
        # investment costs
        "c_vehicles", "c_gcs", "c_cs", "c_garage_cs", "c_garage", "c_garage_workstations",
        "c_stat_storage", "c_invest",
        # annual investment costs
        "c_vehicles_annual",  "c_gcs_annual",  "c_cs_annual", "c_garage_annual",
        "c_stat_storage_annual", "c_invest_annual",
        # annual maintainance costs
        "c_maint_gc_annual", "c_maint_infrastructure_annual", "c_maint_vehicles_annual",
        "c_maint_stat_storage_annual", "c_maint_annual",
        # annual electricity costs
        "c_el_procurement_annual", "c_el_power_price_annual", "c_el_energy_price_annual",
        "c_el_taxes_annual", "c_el_feed_in_remuneration_annual", "c_el_annual"]}

    # INVESTMENT COSTS #

    # VEHICLES
    v_types = schedule.scenario["components"]["vehicle_types"]
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
            c_vehicles_vt = (
                schedule.vehicle_type_counts[v_type] *
                (costs_vehicle + (c_params["vehicles"][v_type]["lifetime"] //
                                  c_params["batteries"]["lifetime_battery"]) *
                 v_keys["capacity"] * c_params["batteries"]["cost_per_kWh"]))
            costs["c_vehicles"] += c_vehicles_vt
            # calculate annual cost of vehicles of this type, depending on their lifetime
            costs["c_vehicles_annual"] += c_vehicles_vt / c_params["vehicles"][v_type]["lifetime"]

    # GRID CONNECTION POINTS
    gcs = schedule.scenario["components"]["grid_connectors"]
    for gcID in gcs.keys():
        # get max. power of grid connector
        gc_timeseries = getattr(scenario, f"{gcID}_timeseries")
        gc_max_power = -min(gc_timeseries["grid power [kW]"])
        # get voltage_level
        voltage_level = schedule.stations[gcID]["voltage_level"]
        # get distance between grid and grid connector
        distance_to_grid = schedule.stations[gcID].get(
            "distance_to_grid", c_params["gc"][voltage_level]["default_distance"])
        # calculate grid connection costs, typically buidling costs and buiding cost subsidy
        c_gc = (c_params["gc"][voltage_level]["capex_gc_fix"] +
                c_params["gc"][voltage_level]["capex_gc_per_kW"] * gc_max_power +
                c_params["gc"][voltage_level]["capex_gc_per_meter"] * distance_to_grid)
        # calculate transformer costs
        c_transformer = (
            c_params["gc"][voltage_level]["capex_transformer_fix"] +
            c_params["gc"][voltage_level]["capex_transformer_per_kW"] * gc_max_power)
        # calculate total cost of grid connection
        costs["c_gcs"] += c_gc + c_transformer
        # calculate annual costs of grid connection, depending on lifetime of gc and transformer
        costs["c_gcs_annual"] += (c_gc / c_params["gc"]["lifetime_gc"] +
                                  c_transformer / c_params["gc"]["lifetime_transformer"])
        # calculate maintenance cost of grid connection
        costs["c_maint_gc_annual"] += (
            c_transformer * c_params["gc"]["c_maint_transformer_per_year"])
        # STATIONARY STORAGE
        # assume there is a stationary storage
        try:
            # calculate costs of stationary storage
            costs["c_stat_storage"] += (
                c_params["stationary_storage"]["capex_fix"] +
                c_params["stationary_storage"]["capex_per_kWh"] *
                schedule.stations[gcID]["battery"]["capacity"])
        except KeyError:
            # if no stationary storage at grid connector: cost is 0
            pass
    costs["c_stat_storage_annual"] = (
        costs["c_stat_storage"] / c_params["stationary_storage"]["lifetime_stat_storage"])
    costs["c_maint_stat_storage_annual"] = (
        costs["c_stat_storage"] * c_params["stationary_storage"]["c_maint_stat_storage_per_year"])

    # CHARGING INFRASTRUCTURE
    cs = schedule.scenario["components"]["charging_stations"]
    # depot charging stations - each charging bus generates one CS
    for csID in cs.values():
        if csID["type"] == "deps":
            costs["c_cs"] += c_params["cs"]["capex_deps_per_kW"] * csID["max_power"]
    # opportunity charging stations - nr of CS depend on max nr of simultaneously occupied CS
    for gcID in gcs:
        if schedule.stations[gcID]["type"] == "opps":
            gc_timeseries = getattr(scenario, f"{gcID}_timeseries")
            costs["c_cs"] += (
                c_params["cs"]["capex_opps_per_kW"] *
                # get cs power of this opps
                schedule.stations[gcID].get("cs_power_opps", vars(args)["cs_power_opps"]) *
                # get max. nr of occupied CS per grid connector
                max(gc_timeseries["CS in use"]))

    # calculate annual cost of charging stations, depending on their lifetime
    costs["c_cs_annual"] = costs["c_cs"] / c_params["cs"]["lifetime_cs"]

    # GARAGE
    costs["c_garage_cs"] = (
        c_params["garage"]["n_charging_stations"]
        * c_params["garage"]["power_cs"] * c_params["cs"]["capex_deps_per_kW"])
    costs["c_garage_workstations"] = -(
        -sum(schedule.vehicle_type_counts.values())
        // c_params["garage"]["vehicles_per_workstation"]
        * c_params["garage"]["cost_per_workstation"])
    costs["c_garage"] = costs["c_garage_cs"] + costs["c_garage_workstations"]
    costs["c_garage_annual"] = (
        costs["c_garage_cs"] / c_params["cs"]["lifetime_cs"]
        + costs["c_garage_workstations"] / c_params["garage"]["lifetime_workstations"])

    # MAINTENANCE
    costs["c_maint_infrastructure_annual"] = (
        costs["c_cs"] * c_params["cs"]["c_maint_cs_per_year"]
        + costs["c_maint_gc_annual"] + costs["c_maint_stat_storage_annual"])
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

    for gcID, gc in scenario.components.grid_connectors.items():
        pv = sum([pv.nominal_power for pv in scenario.components.photovoltaics.values()
                  if pv.parent == gcID])
        timeseries = vars(scenario).get(f"{gcID}_timeseries")

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
            power_pv_nominal=pv,
        )
        # ToDo: Decide if gc-specific costs should be added to scenario object to use in report.py
        # setattr(scenario.components.grid_connectors[gcID], "costs_electricity", costs_electricity)
        costs["c_el_procurement_annual"] += costs_electricity['power_procurement_costs_per_year']
        costs["c_el_power_price_annual"] += costs_electricity['capacity_costs_eur']
        costs["c_el_energy_price_annual"] += costs_electricity['commodity_costs_eur_per_year']
        costs["c_el_taxes_annual"] += costs_electricity['levies_fees_and_taxes_per_year']
        costs["c_el_feed_in_remuneration_annual"] += costs_electricity[
            'feed_in_remuneration_per_year']
        costs["c_el_annual"] += costs_electricity['total_costs_per_year']

    # round all costs to ct
    for key, v in costs.items():
        costs[key] = round(v, 2)

    print("\n"
          "Total costs: \n"
          f"Investment cost: {costs['c_invest']} €. \n"
          f"Annual investment costs: {costs['c_invest_annual']} €/a. \n"
          f"Annual maintenance costs: {costs['c_maint_annual']} €/a. \n"
          f"Annual costs for electricity: {costs['c_el_annual']} €/a.\n")

    setattr(scenario, "costs", costs)
