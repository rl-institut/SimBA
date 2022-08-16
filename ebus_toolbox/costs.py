import json


class Costs:

    def __init__(self, vehicle_types) -> None:
        self.vehicle_types = vehicle_types

    @staticmethod
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
        round_to_numbers = 2
        with open(args.cost_params) as json_file:
            c_params = json.load(json_file)
        with open(args.electrified_stations) as json_file:
            el_stations = json.load(json_file)

        # VEHICLES
        c_busses = 0
        v_types = schedule.scenario["constants"]["vehicle_types"]
        for i in v_types:
            if schedule.vehicle_type_counts[i] > 0:
                if "solo" in v_types[i]["name"]:
                    costs_bus = c_params["bus"]["solo_bus"]
                elif "articulated" in v_types[i]["name"]:
                    costs_bus = c_params["bus"]["articulated_bus"]
                else:
                    costs_bus = c_params["bus"]["battery_kwh"] = 0
                    print("Warning: Name of vehicle_type doesn't include 'solo' or 'articulated. "
                          "Unable to calculate investment and maintenance costs for busses.'")
                c_busses += (schedule.vehicle_type_counts[i] *
                             (costs_bus - (-c_params["bus"]["lifespan_bus"]
                                           // c_params["bus"]["lifespan_battery"]) *
                              v_types[i]["capacity"] * c_params["bus"]["battery_kwh"]))
        c_busses_annual = round(c_busses / c_params["bus"]["lifespan_bus"], round_to_numbers)

        # GRID CONNECTION POINTS
        c_gcs = 0
        gcs = schedule.scenario["constants"]["grid_connectors"]
        for i in gcs:
            c_gcs += c_params["gc"]["building_cost_subsidy_per_kW"] * gcs[i]["max_power"] + \
                         c_params["gc"]["transformer"] + c_params["gc"]["fix"] + \
                         c_params["gc"]["variable_per_meter"] * \
                         el_stations[i]["distance_transformer_in_meter"]
        c_gcs_annual = round(c_gcs / c_params["gc"]["lifespan_gc"], round_to_numbers)

        # CHARGING INFRASTRUCTURE
        c_cs = 0
        cs = schedule.scenario["constants"]["charging_stations"]
        for i in cs:
            if "deps" in i:
                c_cs += c_params["cs"]["depot_cs_per_kW"] * cs[i]["max_power"]
            elif "opps" in i:
                c_cs += c_params["cs"]["opportunity_cs_per_kW"] * cs[i]["max_power"]
        c_cs_annual = round(c_cs / c_params["cs"]["lifespan_cs"], round_to_numbers)

        # GARAGE

        # MAINTENANCE
        m_infra = c_cs * c_params["maintenance"]["percentage_of_cs_infra"]
        m_bus = 0
        drive_days = -(-(schedule.scenario["scenario"]["n_intervals"] / 60 /
                         schedule.scenario["scenario"]["interval"]) // 24)
        for i in schedule.rotations:
            vc_type = f"{schedule.rotations[i].vehicle_type}_{schedule.rotations[i].charging_type}"
            v_name = v_types[vc_type]["name"]
            if "solo" in v_name:
                m_bus += (schedule.rotations[i].distance / 1000 / drive_days * 365) * \
                         c_params["maintenance"]["solo_bus_per_km"]
            elif "articulated" in v_name:
                m_bus += (schedule.rotations[i].distance / 1000 / drive_days * 365) * \
                         c_params["maintenance"]["articulated_bus_per_km"]
        c_maintenance_annual = round(m_infra + m_bus, round_to_numbers)

        c_invest = c_busses + c_cs + c_gcs
        c_invest_annual = c_busses_annual + c_cs_annual + c_gcs_annual

        return {"c_invest": c_invest, "c_invest_annual": c_invest_annual,
                "c_maintenance_annual": c_maintenance_annual}
