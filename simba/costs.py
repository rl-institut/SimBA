import logging
import warnings

import spice_ev.scenario
from spice_ev.costs import calculate_costs as calc_costs_spice_ev

import simba.schedule


class CostsObject:
    """
    Class representing costs associated with vehicles and infrastructure.

    :param schedule: The simulation schedule.
    :type schedule: simba.schedule.Schedule
    :param scenario: The scenario for the simulation.
    :type scenario: spice_ev.scenario.Scenario
    :param args: Arguments for the simulation.
    :param c_params: Parameters for the simulation costs.
    :type c_params: dict
    """
    CUMULATED = "cumulated"
    GARAGE = "garage"
    # Output is sorted by this parameter
    SORT_COLUMN = "c_invest"

    def __init__(self, schedule: simba.schedule.Schedule, scenario: spice_ev.scenario.Scenario,
                 args,
                 c_params: dict):
        """
        Initialize the CostsObject instance.

        :param schedule: The simulation schedule.
        :type schedule: simba.schedule.Schedule
        :param scenario: The scenario for the simulation.
        :type scenario: spice_ev.scenario.Scenario
        :param args: Arguments for the simulation.
        :param c_params: Parameters for the simulation costs.
        :type c_params: dict
        """

        gc_in_use = scenario.components.grid_connectors
        self.costs_per_gc = {gc: {key: 0 for key in self.get_gc_cost_variables()} for gc in
                             [self.GARAGE] + list(gc_in_use)}
        self.units: dict = None
        self.vehicles_per_gc: dict = None
        self.gc_rotations: dict = None
        self.schedule = schedule
        self.scenario = scenario
        self.args = args
        self.params = c_params
        self.rounding_precision = 2

    def info(self):
        """
        Provide information about the total costs.

        :return: A formatted string containing information about total costs.
        :rtype: str
        """
        cumulated = self.costs_per_gc[self.CUMULATED]
        return ("\nTotal costs:\n"
                f"Investment cost: {cumulated['c_invest']} €. \n"
                f"Annual investment costs: {cumulated['c_invest_annual']} €/a. \n"
                f"Annual maintenance costs: {cumulated['c_maint_annual']} €/a. \n"
                f"Annual costs for electricity: {cumulated['c_el_annual']} €/a.\n")

    def set_gc_rotations(self):
        """
        Set grid connector rotations.

        :return: The updated instance.
        :rtype: CostsObject
        """
        rotations_per_gc = {
            gc: [rot for rot in self.schedule.rotations.values() if rot.departure_name == gc]
            for gc in self.get_gcs()}
        self.gc_rotations = rotations_per_gc
        return self

    def get_unit(self, key):
        """
        Get the unit for a given key.

        :param key: The key to get the unit for.
        :type key: str
        :return: The unit associated with the key.
        :rtype: str
        """

        in_costs_per_gc = any(key in d.keys() for d in self.costs_per_gc.values())
        in_vehicles_per_gc = any(key in d.keys() for d in self.vehicles_per_gc.values())
        assert in_costs_per_gc + in_vehicles_per_gc == 1
        if in_costs_per_gc:
            return self.get_costs_per_gc_unit(key)
        if in_vehicles_per_gc:
            return self.get_vehicles_per_gc_unit()

    def get_costs_per_gc_unit(self, key):
        """
        Get the unit for a cost associated with a grid connector.

        :param key: The key to get the unit for.
        :type key: str
        :return: The unit associated with the key.
        :rtype: str
        """
        return self.get_annual_or_not(key)

    def get_vehicles_per_gc_unit(self):
        """
        Get the unit for vehicles per grid connector.

        :return: The unit for vehicles per grid connector.
        :rtype: str
        """
        return "vehicles"

    def get_annual_or_not(self, key):
        """
        Get the unit for annual or non-annual costs.

        :param key: The key to get the unit for.
        :type key: str
        :return: The unit associated with the key.
        :rtype: str
        """
        if "annual" in key:
            return "€/year"
        else:
            return "€"

    def get_gc_cost_variables(self):
        """
         Get the list of variables representing different cost types.

         :return: The list of cost variables.
         :rtype: list
         """
        return [  # investment costs
            "c_vehicles", "c_gcs", "c_cs", "c_garage_cs", "c_garage",
            "c_garage_workstations",
            "c_stat_storage", "c_feed_in", "c_invest",
            # annual investment costs
            "c_vehicles_annual", "c_gcs_annual", "c_cs_annual", "c_garage_annual",
            "c_stat_storage_annual", "c_feed_in_annual", "c_invest_annual",
            # annual maintenance costs
            "c_maint_vehicles_annual", "c_maint_gc_annual", "c_maint_cs_annual",
            "c_maint_stat_storage_annual", "c_maint_feed_in_annual",
            "c_maint_infrastructure_annual",
            "c_maint_annual",
            # annual electricity costs
            "c_el_procurement_annual", "c_el_power_price_annual",
            "c_el_energy_price_annual",
            "c_el_taxes_annual", "c_el_feed_in_remuneration_annual", "c_el_annual"]

    def get_columns(self):
        """
        Get the sorted columns for cost calculations.

        :return: The list of columns.
        :rtype: list
        """
        return list(dict(sorted(self.costs_per_gc.items(), key=lambda x: x[1][self.SORT_COLUMN],
                                reverse=True)).keys())

    def get_gcs(self):
        """
        Get the grid connectors used in the scenario.

        :return: The grid connectors.
        :rtype: dict
        """
        return self.scenario.components.grid_connectors

    def to_csv_lists(self):
        """
        Convert costs to a list of lists easily converteable to a csv.

        :return: The list of lists of parameters, units and costs per gc.
        :rtype: list
        """
        output = [["parameter", "unit"] + self.get_columns()]
        for key in self.get_vehicle_types():
            row = [key, self.get_unit(key)]
            for col in self.get_columns():
                row.append(self.vehicles_per_gc[col][key])
            output.append(row)

        for key in next(iter(self.costs_per_gc.values())):
            row = [key, self.get_unit(key)]
            for col in self.get_columns():
                row.append(round(self.costs_per_gc[col][key], self.rounding_precision))
            output.append(row)
        return output

    def get_vehicle_types(self):
        """
         Get the types of vehicles used.

         :return: The vehicle types.
         :rtype: list
         """
        return next(iter(self.vehicles_per_gc.values())).keys()

    def cumulate(self):
        """
        Cumulate the costs of vehicles and infrastructure.

        :return: The updated instance.
        :rtype: CostsObject
        """
        v_types = self.get_vehicle_types()
        self.vehicles_per_gc[self.CUMULATED] = dict()
        for v_type in v_types:
            self.vehicles_per_gc[self.CUMULATED][v_type] = sum(
                [self.vehicles_per_gc[gc][v_type] for gc in self.get_gcs()])

        # Cumulate gcs variables
        self.costs_per_gc[self.CUMULATED] = dict()
        for key in next(iter(self.costs_per_gc.values())):
            self.costs_per_gc[self.CUMULATED][key] = 0
            for gc in self.costs_per_gc.keys()-[self.CUMULATED]:
                self.costs_per_gc[self.CUMULATED][key] += self.costs_per_gc[gc][key]

        self.costs_per_gc[self.CUMULATED]["c_invest"] += self.costs_per_gc[self.GARAGE]["c_garage"]
        self.costs_per_gc[self.CUMULATED]["c_invest_annual"] += self.costs_per_gc[self.GARAGE][
            "c_garage_annual"]
        return self

    def get_vehicle_count_per_gc(self):
        """
        Calculate the number of vehicles at each grid connector.

        :raises Exception: If a vehicle from the scenario is not found in any rotation
        :return: The number of vehicles at each grid connector.
        :rtype: dict
        """
        v_types = self.schedule.scenario["components"]["vehicle_types"]
        vehicle_count_per_gc = {gc: {v_type_name: 0 for v_type_name in v_types if
                                     self.schedule.vehicle_type_counts[v_type_name] > 0} for gc in
                                list(self.get_gcs()) + [self.GARAGE]}
        for name, vehicle in self.scenario.components.vehicles.items():
            for rot in self.schedule.rotations.values():
                if rot.vehicle_id == name:
                    vehicle_count_per_gc[rot.departure_name][
                        "_".join(name.split("_")[:-1])] += 1
                    break
            else:
                raise Exception(f"Vehicle {name} not found in rotations")
        return vehicle_count_per_gc

    def set_vehicles_per_gc(self):
        """Calculate vehicle numbers and set vehicles_per_gc before vehicle costs can be calculated.

        :return: self
        :rtype: CostsObject
        """
        self.vehicles_per_gc = self.get_vehicle_count_per_gc()
        return self

    def set_vehicle_costs_per_gc(self):
        """
        Calculate and set the costs associated with vehicles at each grid connector.

        :return: self
        :rtype: CostsObject
        """
        # Iterate over each gc and use the vehicles which start and end at this gc for cost
        # calculation
        for gc, vehicle_types in self.vehicles_per_gc.items():
            for v_type, vehicle_number in vehicle_types.items():
                try:
                    costs_vehicle = self.params["vehicles"][v_type]["capex"]
                except KeyError:
                    warnings.warn("No capex defined for vehicle type " + v_type +
                                  ". Unable to calculate investment costs for this vehicle type.")
                    continue
                vehicle_lifetime = self.params["vehicles"][v_type]["lifetime"]
                battery_lifetime = self.params["batteries"]["lifetime_battery"]
                capacity = self.schedule.scenario["components"]["vehicle_types"][v_type]["capacity"]
                cost_per_kWh = self.params["batteries"]["cost_per_kWh"]
                c_vehicles_vt = (vehicle_number * (costs_vehicle + (
                        vehicle_lifetime // battery_lifetime) * capacity * cost_per_kWh))
                c_vehicles_annual = c_vehicles_vt / vehicle_lifetime
                self.costs_per_gc[gc]["c_vehicles"] += c_vehicles_vt
                self.costs_per_gc[gc]["c_vehicles_annual"] += c_vehicles_annual
        return self

    def set_grid_connection_costs(self):
        """Calculate the costs of each grid connection

        :raises Exception: if grid operator of grid connector can not be found in cost params
        :return: self
        :rtype: CostsObject
        """
        # get dict with costs for specific grid operator
        for gcID, gc in self.get_gcs().items():
            try:
                c_params_go = self.params[gc.grid_operator]
            except KeyError:
                raise Exception(f"{gcID}: Unknown grid operator {gc.grid_operator}, "
                                "might have to set in electrified stations or cost params")
            # get max. power of grid connector

            gc_timeseries = getattr(self.scenario, f"{gcID}_timeseries")
            gc_max_power = -min(gc_timeseries["grid supply [kW]"])
            # skip grid connector, if not actually used
            if gc_max_power == 0:
                continue
            # get voltage_level
            voltage_level = self.schedule.stations[gcID].get("voltage_level",
                                                             self.args.default_voltage_level)
            # get distance between grid and grid connector
            distance_to_grid = self.schedule.stations[gcID].get(
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
            # calculate annual costs of grid connection, depending on lifetime of gc and transformer
            c_gc_annual = (c_gc / c_params_go["gc"]["lifetime_gc"] +
                           c_transformer / c_params_go["gc"]["lifetime_transformer"])
            # calculate maintenance cost of grid connection
            c_maint_gc_annual = (c_transformer * c_params_go["gc"]["c_maint_transformer_per_year"])
            # Store the individual costs of the GC as well
            self.costs_per_gc[gcID]["c_gcs"] = c_gc + c_transformer
            self.costs_per_gc[gcID]["c_gcs_annual"] = c_gc_annual
            self.costs_per_gc[gcID]["c_maint_gc_annual"] = c_maint_gc_annual

            # STATIONARY STORAGE
            # assume there is a stationary storage
            try:
                # calculate costs of stationary storage
                c_stat_storage = (
                        self.params["stationary_storage"]["capex_fix"] +
                        self.params["stationary_storage"]["capex_per_kWh"] *
                        self.schedule.stations[gcID]["battery"]["capacity"])
                self.costs_per_gc[gcID]["c_stat_storage"] = c_stat_storage

                self.costs_per_gc[gcID]["c_stat_storage_annual"] = c_stat_storage / self.params[
                    "stationary_storage"]["lifetime_stat_storage"]
                self.costs_per_gc[gcID]["c_maint_stat_storage_annual"] = (
                        c_stat_storage * self.params["stationary_storage"][
                            "c_maint_stat_storage_per_year"])
            except KeyError:
                # if no stationary storage at grid connector: cost is 0
                pass

            # FEED IN (local renewables)
            # assume there is a feed in
            try:
                # calculate costs of feed in
                c_feed_in = (self.params["feed_in"]["capex_fix"] +
                             self.params["feed_in"]["capex_per_kW"] *
                             self.schedule.stations[gcID]["energy_feed_in"]["nominal_power"])
                self.costs_per_gc[gcID]["c_feed_in"] = c_feed_in
                self.costs_per_gc[gcID]["c_feed_in_annual"] = (
                        c_feed_in / self.params["feed_in"]["lifetime_feed_in"])
                self.costs_per_gc[gcID]["c_maint_feed_in_annual"] = (
                        c_feed_in * self.params["feed_in"]["c_maint_feed_in_per_year"])
            except KeyError:
                # if no feed in at grid connector: cost is 0
                pass
        return self

    def set_charging_infrastructure_costs(self):
        """
        Calculate the costs associated with charging infrastructure at each grid connector.

        :return: self
        :rtype: CostsObject
        """

        all_stations = self.schedule.scenario["components"]["charging_stations"]
        # find all stations from the schedule which are grid connectors
        used_stations = {gcId: station for gcId, station in self.schedule.stations.items() if
                         gcId in self.get_gcs()}

        for gcID, station in used_stations.items():
            # depot charging station - each charging bus generates one CS
            if station["type"] == "deps":
                cs_at_station = [cs for cs in all_stations.values() if cs["parent"] == gcID]
                if station["n_charging_stations"] is not None:
                    # get nr of CS in use at station
                    n_cs = min(len(cs_at_station), station["n_charging_stations"])
                    # get max. defined power for deps or else max. power at deps in general
                    defined_max_power = max(
                        station.get("cs_power_deps_depb", vars(self.args)["cs_power_deps_depb"]),
                        station.get("cs_power_deps_oppb", vars(self.args)["cs_power_deps_oppb"])
                    )
                    # calculate costs with nr of CS at electrified station and its defined max power
                    c_cs = self.params["cs"]["capex_deps_per_kW"] * n_cs * defined_max_power
                else:
                    # get sum of installed power at electrified stations without predefined nr of CS
                    sum_power_at_station = sum(cs["max_power"] for cs in cs_at_station)
                    # calculate costs with sum of installed power at electrified station
                    c_cs = self.params["cs"]["capex_deps_per_kW"] * sum_power_at_station

            # opportunity charging station - nr of CS depend on max nr of simultaneously occupied CS
            elif station["type"] == "opps":
                gc_timeseries = getattr(self.scenario, f"{gcID}_timeseries")
                # get nr of CS in use at station
                n_cs = max(gc_timeseries["# CS in use [-]"])
                # get max. defined power for opps or else max. power at opps in general
                defined_max_power = station.get("cs_power_opps", vars(self.args)["cs_power_opps"])
                # calculate costs with nr of CS at electrified station and its defined max power
                c_cs = self.params["cs"]["capex_opps_per_kW"] * n_cs * defined_max_power
            else:
                warnings.warn(f"Electrified station {station} must be deps or opps. "
                              "Unable to calculate investment costs for this station.")
                c_cs = 0
            self.costs_per_gc[gcID]["c_cs"] = c_cs

            # calculate annual cost of charging stations, depending on their lifetime
            self.costs_per_gc[gcID]["c_cs_annual"] = c_cs / self.params["cs"]["lifetime_cs"]

            self.costs_per_gc[gcID]["c_maint_cs_annual"] = c_cs * self.params["cs"][
                "c_maint_cs_per_year"]

            # TOTAL COSTS
            # Note: Grid Connectors do share the costs of the garage.
            self.costs_per_gc[gcID]["c_invest"] = (
                    self.costs_per_gc[gcID]["c_vehicles"]
                    + self.costs_per_gc[gcID]["c_cs"]
                    + self.costs_per_gc[gcID]["c_gcs"]
                    + self.costs_per_gc[gcID]["c_stat_storage"]
                    + self.costs_per_gc[gcID]["c_feed_in"])
            self.costs_per_gc[gcID]["c_invest_annual"] = (
                    self.costs_per_gc[gcID]["c_vehicles_annual"]
                    + self.costs_per_gc[gcID]["c_cs_annual"]
                    + self.costs_per_gc[gcID]["c_gcs_annual"]
                    + self.costs_per_gc[gcID]["c_stat_storage_annual"]
                    + self.costs_per_gc[gcID]["c_feed_in_annual"])

            # MAINTENANCE COSTS #
            self.costs_per_gc[gcID]["c_maint_infrastructure_annual"] = (
                    self.costs_per_gc[gcID]["c_maint_cs_annual"] +
                    self.costs_per_gc[gcID]["c_maint_gc_annual"] +
                    self.costs_per_gc[gcID]["c_maint_stat_storage_annual"] +
                    self.costs_per_gc[gcID]["c_maint_feed_in_annual"])

            # calculate (ceil) number of days in scenario
            drive_days = -(-(self.schedule.scenario["scenario"]["n_intervals"] *
                             self.schedule.scenario["scenario"]["interval"]) // (24 * 60))
            for rot in self.gc_rotations[gcID]:
                v_type_rot = f"{rot.vehicle_type}_{rot.charging_type}"
                try:
                    self.costs_per_gc[gcID]["c_maint_vehicles_annual"] += (rot.distance / 1000 /
                                                                           drive_days * 365) * \
                                                                          self.params["vehicles"][
                                                                              v_type_rot][
                                                                              "c_maint_per_km"]
                except KeyError:
                    warnings.warn("No maintenance costs defined for vehicle type " +
                                  self.params["vehicles"][v_type_rot] +
                                  ". Unable to calculate maintenance costs for this vehicle type.")
            self.costs_per_gc[gcID]["c_maint_annual"] = (
                        self.costs_per_gc[gcID]["c_maint_infrastructure_annual"] +
                        self.costs_per_gc[gcID]["c_maint_vehicles_annual"])
        return self

    def set_electricity_costs(self):
        """
        Calculate the electricity costs at each grid connector.

        :return: self
        :rtype: CostsObject
        """
        for gcID, gc in self.get_gcs().items():
            station = self.schedule.stations[gcID]
            pv = sum([pv.nominal_power for pv in self.scenario.components.photovoltaics.values()
                      if pv.parent == gcID])
            timeseries = vars(self.scenario).get(f"{gcID}_timeseries")

            # calculate costs for electricity
            costs_electricity = calc_costs_spice_ev(
                strategy=vars(self.args)["strategy_" + station.get("type")],
                voltage_level=gc.voltage_level,
                interval=self.scenario.interval,
                timestamps_list=timeseries.get("time"),
                power_grid_supply_list=timeseries.get("grid supply [kW]"),
                price_list=timeseries.get("price [EUR/kWh]"),
                power_fix_load_list=timeseries.get("fixed load [kW]"),
                power_generation_feed_in_list=timeseries.get("generation feed-in [kW]"),
                power_v2g_feed_in_list=timeseries.get("V2G feed-in [kW]"),
                power_battery_feed_in_list=timeseries.get("battery feed-in [kW]"),
                charging_signal_list=timeseries.get("window"),
                price_sheet_path=self.args.cost_parameters_file,
                grid_operator=gc.grid_operator,
                power_pv_nominal=pv,
            )
            self.costs_per_gc[gcID]["c_el_procurement_annual"] = costs_electricity[
                'power_procurement_costs_per_year']
            self.costs_per_gc[gcID]["c_el_power_price_annual"] = costs_electricity[
                'capacity_costs_eur']
            self.costs_per_gc[gcID]["c_el_energy_price_annual"] = costs_electricity[
                'commodity_costs_eur_per_year']
            self.costs_per_gc[gcID]["c_el_taxes_annual"] = costs_electricity[
                'levies_fees_and_taxes_per_year']
            self.costs_per_gc[gcID]["c_el_feed_in_remuneration_annual"] = costs_electricity[
                'feed_in_remuneration_per_year']
            self.costs_per_gc[gcID]["c_el_annual"] = costs_electricity['total_costs_per_year']
        return self

    def set_garage_costs(self):
        """
        Calculate the costs associated with the garage.

        Note that although a garage can have charging stations, electricity costs are not used.

        :return: self
        :rtype: CostsObject
        """
        c_garage_cs = (
                self.params["garage"]["n_charging_stations"]
                * self.params["garage"]["power_cs"] * self.params["cs"]["capex_deps_per_kW"])
        self.costs_per_gc[self.GARAGE]["c_garage_cs"] = c_garage_cs

        c_garage_workstations = -(
                -sum(self.schedule.vehicle_type_counts.values())
                // self.params["garage"]["vehicles_per_workstation"]
                * self.params["garage"]["cost_per_workstation"])
        self.costs_per_gc[self.GARAGE]["c_garage_workstations"] = c_garage_workstations

        self.costs_per_gc[self.GARAGE]["c_garage"] = c_garage_cs + c_garage_workstations
        self.costs_per_gc[self.GARAGE]["c_garage_annual"] = (
                    c_garage_cs / self.params["cs"]["lifetime_cs"]
                    + c_garage_workstations / self.params["garage"][
                        "lifetime_workstations"])
        return self


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
    cost_object = CostsObject(schedule, scenario, args, c_params)

    # Split the rotations to their home depots as a dict
    cost_object.set_gc_rotations()

    # Count the number of vehicles per gc
    cost_object.set_vehicles_per_gc()

    # Calculate the cost these vehicles
    cost_object.set_vehicle_costs_per_gc()

    # Calculate the cost of the grid connections including stationary and feed in costs
    cost_object.set_grid_connection_costs()

    # costs for charging stations
    cost_object.set_charging_infrastructure_costs()

    cost_object.set_garage_costs()

    # Get electricity costs from spiceEV
    cost_object.set_electricity_costs()

    # Cumulate the costs of all gcs
    cost_object.cumulate()

    logging.info(cost_object.info())

    setattr(scenario, "costs", cost_object)
