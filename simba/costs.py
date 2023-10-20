import logging
import warnings

from spice_ev.costs import calculate_costs as calc_costs_spice_ev


def calculate_costs(c_params, scenario, schedule, args):
    """ Calculates annual costs of all necessary vehicles and infrastructure.

    :param c_params: Infrastructure component costs and specific grid operator costs.
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
        "c_vehicles_annual", "c_gcs_annual", "c_cs_annual", "c_garage_annual",
        "c_stat_storage_annual", "c_invest_annual",
        # annual maintenance costs
        "c_maint_vehicles_annual", "c_maint_gc_annual", "c_maint_cs_annual",
        "c_maint_stat_storage_annual", "c_maint_infrastructure_annual", "c_maint_annual",
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
    gc_in_use = scenario.components.grid_connectors
    for gcID, gc in gc_in_use.items():
        # get dict with costs for specific grid operator
        try:
            c_params_go = c_params[gc.grid_operator]
        except KeyError:
            raise Exception(f"{gcID}: Unknown grid operator {gc.grid_operator}, "
                            "might have to set in electrified stations or cost params")
        # get max. power of grid connector
        gc_timeseries = getattr(scenario, f"{gcID}_timeseries")
        gc_max_power = -min(gc_timeseries["grid supply [kW]"])
        # get voltage_level
        voltage_level = schedule.stations[gcID].get("voltage_level", args.default_voltage_level)
        # get distance between grid and grid connector
        distance_to_grid = schedule.stations[gcID].get(
            "distance_to_grid", c_params_go["gc"][voltage_level]["default_distance"])
        # calculate grid connection costs, typically building costs and building cost subsidy
        c_gc = (c_params_go["gc"][voltage_level]["capex_gc_fix"] +
                c_params_go["gc"][voltage_level]["capex_gc_per_kW"] * gc_max_power +
                c_params_go["gc"][voltage_level]["capex_gc_per_meter"] * distance_to_grid)
        # calculate transformer costs
        c_transformer = (
            c_params_go["gc"][voltage_level]["capex_transformer_fix"] +
            c_params_go["gc"][voltage_level]["capex_transformer_per_kW"] * gc_max_power)
        # calculate total cost of grid connection
        costs["c_gcs"] += c_gc + c_transformer
        # calculate annual costs of grid connection, depending on lifetime of gc and transformer
        costs["c_gcs_annual"] += (c_gc / c_params_go["gc"]["lifetime_gc"] +
                                  c_transformer / c_params_go["gc"]["lifetime_transformer"])
        # calculate maintenance cost of grid connection
        costs["c_maint_gc_annual"] += (
            c_transformer * c_params_go["gc"]["c_maint_transformer_per_year"])

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
    stations = schedule.scenario["components"]["charging_stations"]
    for gcID, station in schedule.stations.items():
        # check if electrified station from list is actually used in simulation
        if gcID in gc_in_use:
            # depot charging station - each charging bus generates one CS
            if station["type"] == "deps":
                cs_at_station = [cs for cs in stations.values() if cs["parent"] == gcID]
                if station["n_charging_stations"] is not None:
                    # get nr of CS in use at station
                    n_cs = min(len(cs_at_station), station["n_charging_stations"])
                    # get max. defined power for deps or else max. power at deps in general
                    defined_max_power = max(
                        station.get("cs_power_deps_depb", vars(args)["cs_power_deps_depb"]),
                        station.get("cs_power_deps_oppb", vars(args)["cs_power_deps_oppb"])
                    )
                    # calculate costs with nr of CS at electrified station and its defined max power
                    costs["c_cs"] += c_params["cs"]["capex_deps_per_kW"] * n_cs * defined_max_power
                else:
                    # get sum of installed power at electrified stations without predefined nr of CS
                    sum_power_at_station = sum(cs["max_power"] for cs in cs_at_station)
                    # calculate costs with sum of installed power at electrified station
                    costs["c_cs"] += c_params["cs"]["capex_deps_per_kW"] * sum_power_at_station
            # opportunity charging station - nr of CS depend on max nr of simultaneously occupied CS
            elif station["type"] == "opps":
                gc_timeseries = getattr(scenario, f"{gcID}_timeseries")
                # get nr of CS in use at station
                n_cs = max(gc_timeseries["# CS in use [-]"])
                # get max. defined power for opps or else max. power at opps in general
                defined_max_power = station.get("cs_power_opps", vars(args)["cs_power_opps"])
                # calculate costs with nr of CS at electrified station and its defined max power
                costs["c_cs"] += c_params["cs"]["capex_opps_per_kW"] * n_cs * defined_max_power
            else:
                warnings.warn(f"Electrified station {station} must be deps or opps. "
                              "Unable to calculate investment costs for this station.")

    # calculate annual cost of charging stations, depending on their lifetime
    costs["c_cs_annual"] = costs["c_cs"] / c_params["cs"]["lifetime_cs"]

    costs["c_maint_cs_annual"] = costs["c_cs"] * c_params["cs"]["c_maint_cs_per_year"]

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

    # TOTAL COSTS
    costs["c_invest"] = (costs["c_vehicles"] + costs["c_cs"] + costs["c_gcs"] + costs["c_garage"] +
                         costs["c_stat_storage"])
    costs["c_invest_annual"] = (costs["c_vehicles_annual"] + costs["c_cs_annual"] +
                                costs["c_gcs_annual"] + costs["c_garage_annual"] +
                                costs["c_stat_storage_annual"])

    # MAINTENANCE COSTS #

    costs["c_maint_infrastructure_annual"] = (costs["c_maint_cs_annual"] +
                                              costs["c_maint_gc_annual"] +
                                              costs["c_maint_stat_storage_annual"])
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
            power_grid_supply_list=timeseries.get("grid supply [kW]"),
            price_list=timeseries.get("price [EUR/kWh]"),
            power_fix_load_list=timeseries.get("fixed load [kW]"),
            power_generation_feed_in_list=timeseries.get("generation feed-in [kW]"),
            power_v2g_feed_in_list=timeseries.get("V2G feed-in [kW]"),
            power_battery_feed_in_list=timeseries.get("battery feed-in [kW]"),
            charging_signal_list=timeseries.get("window"),
            price_sheet_path=args.cost_parameters_file,
            grid_operator=gc.grid_operator,
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

    logging.info(
        "\nTotal costs:\n"
        f"Investment cost: {costs['c_invest']} €. \n"
        f"Annual investment costs: {costs['c_invest_annual']} €/a. \n"
        f"Annual maintenance costs: {costs['c_maint_annual']} €/a. \n"
        f"Annual costs for electricity: {costs['c_el_annual']} €/a.\n")

    setattr(scenario, "costs", costs)
