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
    c_busses = 0
    v_types = schedule.scenario["constants"]["vehicle_types"]
    for v_type, v_keys in v_types.items():
        if schedule.vehicle_type_counts[v_type] > 0:
            if "solo" in v_keys["name"]:
                costs_bus = c_params["bus"]["solo_bus"]
            elif "articulated" in v_keys["name"]:
                costs_bus = c_params["bus"]["articulated_bus"]
            else:
                costs_bus = c_params["bus"]["battery_kwh"] = 0
                print("Warning: Name of vehicle_type doesn't include 'solo' or 'articulated. "
                      "Unable to calculate investment and maintenance costs for busses.'")
            # sum up cost of busses and its batteries, depending on how often the battery has to be
            # replaced in the lifespan of the busses
            c_busses += (schedule.vehicle_type_counts[v_type] *
                         (costs_bus - (-c_params["bus"]["lifespan_bus"]
                                       // c_params["bus"]["lifespan_battery"]) *
                          v_keys["capacity"] * c_params["bus"]["battery_kwh"]))
    # calculate annual cost of busses, depending on its lifespan
    c_busses_annual = round(c_busses / c_params["bus"]["lifespan_bus"], ROUND_TO_PLACES)

    # GRID CONNECTION POINTS
    c_gcs = 0
    gcs = schedule.scenario["constants"]["grid_connectors"]
    for gcID, gc_keys in gcs.items():
        c_gcs += c_params["gc"]["building_cost_subsidy_per_kW"] * gc_keys["max_power"] + \
                 c_params["gc"]["transformer"] + c_params["gc"]["fix"] + \
                 c_params["gc"]["variable"] * el_stations[gcID]["distance_transformer"]
    # calculate annual cost of grid connectors, depending on its lifespan
    c_gcs_annual = round(c_gcs / c_params["gc"]["lifespan_gc"], ROUND_TO_PLACES)

    # CHARGING INFRASTRUCTURE
    c_cs = 0
    cs = schedule.scenario["constants"]["charging_stations"]
    for csID in cs.values():
        if csID["type"] == "deps":
            c_cs += c_params["cs"]["depot_cs_per_kW"] * csID["max_power"]
        elif csID["type"] == "opps":
            c_cs += c_params["cs"]["opportunity_cs_per_kW"] * csID["max_power"]
    # calculate annual cost of charging stations, depending on its lifespan
    c_cs_annual = round(c_cs / c_params["cs"]["lifespan_cs"], ROUND_TO_PLACES)

    # GARAGE
    c_garage_cs = (c_params["garage"]["n_charging_stations"] * c_params["garage"]["power_cs"] *
                   c_params["cs"]["depot_cs_per_kW"])
    c_garage_positions = -(-sum(schedule.vehicle_type_counts.values()) //
                           c_params["garage"]["busses_per_positions"] *
                           c_params["garage"]["cost_per_position"])
    c_garage = c_garage_cs + c_garage_positions
    # ToDo: annual costs of garage?

    # MAINTENANCE
    m_infra = c_cs * c_params["maintenance"]["percentage_of_cs_infra"]
    m_bus = 0
    # calculate (ceil) number of days in scenario
    drive_days = -(-(schedule.scenario["scenario"]["n_intervals"] *
                     schedule.scenario["scenario"]["interval"]) // (24 * 60))
    for rot in schedule.rotations.values():
        vt_ct = f"{rot.vehicle_type}_{rot.charging_type}"
        if "solo" in v_types[vt_ct]["name"]:
            m_bus += (rot.distance / 1000 / drive_days * 365) * \
                     c_params["maintenance"]["solo_bus_per_km"]
        elif "articulated" in v_types[vt_ct]["name"]:
            m_bus += (rot.distance / 1000 / drive_days * 365) * \
                     c_params["maintenance"]["articulated_bus_per_km"]
    c_maintenance_annual = round(m_infra + m_bus, ROUND_TO_PLACES)

    c_invest = c_busses + c_cs + c_gcs + c_garage
    c_invest_annual = c_busses_annual + c_cs_annual + c_gcs_annual

    return {"c_invest": c_invest, "c_invest_annual": c_invest_annual,
            "c_maintenance_annual": c_maintenance_annual}
