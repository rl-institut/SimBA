import json


def calculate_costs(args, schedule):
    """ Calculates annual costs of all necessary vehicles and infrastructure.

    :param args: Command line arguments and/or arguments from config file.
    :type args: argparse.Namespace
    :param schedule: Information about the whole bus schedule and all used parameters.
    :type schedule: object
    :return: Investment cost [€], annual investment cost [€/a], annual maintenance cost [€/a]
    :rtype: (float, float, float)
    """
    # general settings and imports for calculation
    ROUND_TO_PLACES = 2
    with open(args.cost_params) as json_file:
        c_params = json.load(json_file)
    with open(args.electrified_stations) as json_file:
        el_stations = json.load(json_file)

    # VEHICLES
    c_vehicles = 0
    c_vehicles_annual = 0
    v_types = schedule.scenario["constants"]["vehicle_types"]
    for v_type, v_keys in v_types.items():
        if schedule.vehicle_type_counts[v_type] > 0:
            try:
                costs_vehicle = c_params["vehicles"][v_type.split("_")[0]]["capex"]
                # future solution -> costs_vehicle= c_params["vehicles"][v_type]["capex"]
            except KeyError:
                print("Warning: No capex defined for vehicle type " + v_type.split("_")[0] +
                      ". Unable to calculate investment costs for this vehicle type.")
            else:
                # sum up cost of vehicles and their batteries, depending on how often the battery
                # has to be replaced in the lifetime of the vehicles
                c_vehicles_vt = schedule.vehicle_type_counts[v_type] * \
                                (costs_vehicle - (-c_params["vehicles"][v_type.split("_")[0]]["lifetime"]
                                                 // c_params["batteries"]["lifetime_battery"]) *
                                v_keys["capacity"] * c_params["batteries"]["cost_per_kWh"])
                c_vehicles += c_vehicles_vt
                # calculate annual cost of vehicles of this type, depending on their lifetime
                c_vehicles_annual +=round(c_vehicles_vt/
                                          c_params["vehicles"][v_type.split("_")[0]]["lifetime"],
                                           ROUND_TO_PLACES)


    # GRID CONNECTION POINTS
    c_gcs = 0
    c_gcs_annual = 0
    gcs = schedule.scenario["constants"]["grid_connectors"]
    for gcID, gc_keys in gcs.items():
        try:
            distance_transformer = el_stations[gcID]["distance_transformer"]
        except KeyError:
            distance_transformer = c_params["gc"]["default_distance"]
        c_gc = c_params["gc"]["building_cost_subsidy_per_kW"] * gc_keys["max_power"] + \
               c_params["gc"]["capex_gc_fix"] + \
               c_params["gc"]["capex_gc_per_meter"] * distance_transformer
        c_transformer = c_params["gc"]["capex_transformer_fix"] + \
                        c_params["gc"]["capex_transformer_per_kW"] * gc_keys["max_power"]
        c_gcs += c_gc + c_transformer

        # calculate annual costs of grid connectors, depending the lifetime of gc and transformator
        c_gcs_annual += round(c_gc / c_params["gc"]["lifetime_gc"] +
                             c_transformer / c_params["gc"]["lifetime_transformer"],
                             ROUND_TO_PLACES)

    # CHARGING INFRASTRUCTURE
    c_cs = 0
    cs = schedule.scenario["constants"]["charging_stations"]
    for csID in cs.values():
        if csID["type"] == "deps":
            c_cs += c_params["cs"]["capex_deps_per_kW"] * csID["max_power"]
        elif csID["type"] == "opps":
            c_cs += c_params["cs"]["capex_opps_per_kW"] * csID["max_power"]
    # calculate annual cost of charging stations, depending on their lifetime
    c_cs_annual = round(c_cs / c_params["cs"]["lifetime_cs"], ROUND_TO_PLACES)

    # GARAGE
    c_garage_cs = (c_params["garage"]["n_charging_stations"] * c_params["garage"]["power_cs"] *
                   c_params["cs"]["capex_deps_per_kW"])
    c_garage_workstations = -(-sum(schedule.vehicle_type_counts.values()) //
                              c_params["garage"]["vehicles_per_workstation"] *
                              c_params["garage"]["cost_per_workstation"])
    c_garage = c_garage_cs + c_garage_workstations
    c_garage_annual = round(c_garage_cs / c_params["cs"]["lifetime_cs"]
                            + c_garage_workstations / c_params["garage"]["lifetime_workstations"]
                            , ROUND_TO_PLACES)
    # MAINTENANCE
    m_infra = c_cs * c_params["cs"]["c_maint_cs_per_year"] + \
              c_transformer * c_params["gc"]["c_maint_transformer_per_year"]
    m_vehicles = 0
    # calculate (ceil) number of days in scenario
    drive_days = -(-(schedule.scenario["scenario"]["n_intervals"] *
                     schedule.scenario["scenario"]["interval"]) // (24 * 60))
    for rot in schedule.rotations.values():
        try:
            m_vehicles += (rot.distance / 1000 / drive_days * 365) * \
                          c_params["vehicles"][rot.vehicle_type]["c_maint_per_km"]
        except KeyError:
            print("Warning: No maintanance costs defined for vehicle type " +
                  c_params["vehicles"][rot.vehicle_type] +
                  ". Unable to calculate maintanance costs for this vehicle type.")

    c_maintenance_annual = round(m_infra + m_vehicles, ROUND_TO_PLACES)
    c_invest = c_vehicles + c_cs + c_gcs + c_garage
    c_invest_annual = round(c_vehicles_annual + c_cs_annual + c_gcs_annual + c_garage_annual,
                            ROUND_TO_PLACES)

    return {"c_invest": c_invest, "c_invest_annual": c_invest_annual,
            "c_maintenance_annual": c_maintenance_annual}
