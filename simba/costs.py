import logging
import warnings

import spice_ev.scenario
from spice_ev.costs import calculate_costs as calc_costs_spice_ev

import simba.schedule


class CostsObject:
    CUMULATED = "cumulated"
    GARAGE = "garage"
    SORT_COLUMN = "c_invest"

    def __init__(self, schedule: simba.schedule.Schedule, scenario: spice_ev.scenario.Scenario,
                 args,
                 c_params: dict):
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
        cumulated = self.costs_per_gc[self.CUMULATED]
        return ("\nTotal costs:\n"
                f"Investment cost: {cumulated['c_invest']} €. \n"
                f"Annual investment costs: {cumulated['c_invest_annual']} €/a. \n"
                f"Annual maintenance costs: {cumulated['c_maint_annual']} €/a. \n"
                f"Annual costs for electricity: {cumulated['c_el_annual']} €/a.\n")

    def set_gc_rotations(self):
        rotations_per_gc = {
            gc: [rot for rot in self.schedule.rotations.values() if rot.departure_name == gc]
            for gc in self.get_gcs()}
        self.gc_rotations = rotations_per_gc
        return self

    def get_unit(self, key):
        in_costs_per_gc = any(key in d.keys() for d in self.costs_per_gc.values())
        in_vehicles_per_gc = any(key in d.keys() for d in self.vehicles_per_gc.values())
        assert in_costs_per_gc + in_vehicles_per_gc == 1
        if in_costs_per_gc:
            return self.get_costs_per_gc_unit(key)
        if in_vehicles_per_gc:
            return self.get_vehicles_per_gc_unit()

    def get_costs_per_gc_unit(self, key):
        return self.get_annual_or_not(key)

    def get_vehicles_per_gc_unit(self):
        return "vehicles"

    def get_annual_or_not(self, key):
        if "annual" in key:
            return "€/year"
        else:
            return "€"

    def get_gc_cost_variables(self) -> list[str]:
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
        return list(dict(sorted(self.costs_per_gc.items(), key=lambda x: x[1][self.SORT_COLUMN],
                                reverse=True)).keys())

    def get_gcs(self) -> dict:
        return self.scenario.components.grid_connectors

    def to_csv_lists(self):
        output = [["parameter", "unit"] + self.get_columns()]
        for key in self.get_vehicle_types():
            row = [key, self.get_unit(key)]
            for col in self.get_columns():
                row.append(self.vehicles_per_gc[col][key])
            output.append(row)

        for key in self.costs_per_gc[self.CUMULATED]:
            row = [key, self.get_unit(key)]
            for col in self.get_columns():
                row.append(round(self.costs_per_gc[col][key], self.rounding_precision))
            output.append(row)
        return output

    def get_vehicle_types(self):
        return next(iter(self.vehicles_per_gc.values())).keys()

    def cumulate(self):
        # Cumulate vehicle numbers
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

    def get_vehicle_count_per_gc(self) -> dict:
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
        for gcID, gc in self.get_gcs().items():
            # get dict with costs for specific grid operator
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
        all_stations = self.schedule.scenario["components"]["charging_stations"]
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
                    ["cumulated"]

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
            self.costs_per_gc[gcID]["c_invest"] = (
                    self.costs_per_gc[gcID]["c_vehicles"]
                    + self.costs_per_gc[gcID]["c_cs"]
                    + self.costs_per_gc[gcID]["c_gcs"]
                    # costs[key]["c_garage"] +
                    + self.costs_per_gc[gcID]["c_stat_storage"]
                    + self.costs_per_gc[gcID]["c_feed_in"])
            self.costs_per_gc[gcID]["c_invest_annual"] = (
                    self.costs_per_gc[gcID]["c_vehicles_annual"]
                    + self.costs_per_gc[gcID]["c_cs_annual"]
                    + self.costs_per_gc[gcID]["c_gcs_annual"]
                    # + self.costs_per_gc[gcID]["c_garage_annual"]
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

    cost_object.set_grid_connection_costs()

    cost_object.set_grid_connection_costs()

    cost_object.set_charging_infrastructure_costs()

    cost_object.set_garage_costs()

    cost_object.set_electricity_costs()

    cost_object.cumulate()

    logging.info(cost_object.info())

    setattr(scenario, "costs", cost_object)


def calculate_costs_old(c_params, scenario, schedule, args):
    """ Calculates annual costs of all necessary vehicles and infrastructure.

    :param c_params: Infrastructure component costs and specific grid operator costs.
    :type c_params: dict
    :param scenario: Information about the simulated scenario and all used parameters.
    :type scenario: Scenario
    :param schedule: Information about the whole bus schedule and all used parameters.
    :type schedule: Schedule
    :param args: Configuration arguments specified in config files contained in configs directory.
    :type args: argparse.Namespace
    :raises Exception: if grid operator of grid connector can not be found in cost params
    """

    # initialize dictionary with all costs
    gc_in_use = scenario.components.grid_connectors

    costs = {gc: {cost: 0 for cost in [
        # investment costs
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
        "c_el_taxes_annual", "c_el_feed_in_remuneration_annual", "c_el_annual"]} for gc in
             ["cumulated"] + list(gc_in_use)}

    units = dict()
    # INVESTMENT COSTS #

    # VEHICLES
    v_types = schedule.scenario["components"]["vehicle_types"]

    rotations_per_gc = {gc: [rot for rot in schedule.rotations.values() if rot.departure_name == gc]
                        for gc in gc_in_use}
    rotations_per_gc["cumulated"] = [rot for rot in schedule.rotations.values()]

    vehicle_count_per_gc = CostsObject().get_vehicle_count_per_gc()

    for v_type, v_keys in v_types.items():
        if schedule.vehicle_type_counts[v_type] > 0:
            try:
                costs_vehicle = c_params["vehicles"][v_type]["capex"]
            except KeyError:
                warnings.warn("No capex defined for vehicle type " + v_type +
                              ". Unable to calculate investment costs for this vehicle type.")
                continue

            units[v_type] = "vehicles"

            for gc in gc_in_use:
                vehicle_number = vehicle_count_per_gc[gc][v_type]
                # sum up cost of vehicles and their batteries, depending on how often the battery
                # has to be replaced in the lifetime of the vehicles
                c_vehicles_vt = (
                        vehicle_number *
                        (costs_vehicle + (c_params["vehicles"][v_type]["lifetime"] //
                                          c_params["batteries"]["lifetime_battery"]) *
                         v_keys["capacity"] * c_params["batteries"]["cost_per_kWh"]))
                c_vehicles_annual = c_vehicles_vt / c_params["vehicles"][v_type]["lifetime"]
                costs[gc][v_type] = vehicle_number
                costs[gc]["c_vehicles"] += c_vehicles_vt
                costs[gc]["c_vehicles_annual"] += c_vehicles_annual
                try:
                    costs["cumulated"][v_type] += vehicle_number
                except KeyError:
                    costs["cumulated"][v_type] = vehicle_number

                costs["cumulated"]["c_vehicles"] += c_vehicles_vt
                # calculate annual cost of vehicles of this type, depending on their lifetime
                costs["cumulated"]["c_vehicles_annual"] += c_vehicles_annual

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
        # skip grid connector, if not actually used
        if gc_max_power == 0:
            continue
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
        costs["cumulated"]["c_gcs"] += c_gc + c_transformer
        # calculate annual costs of grid connection, depending on lifetime of gc and transformer
        c_gc_annual = (c_gc / c_params_go["gc"]["lifetime_gc"] +
                       c_transformer / c_params_go["gc"]["lifetime_transformer"])

        costs["cumulated"]["c_gcs_annual"] += c_gc_annual
        # calculate maintenance cost of grid connection
        c_maint_gc_annual = (c_transformer * c_params_go["gc"]["c_maint_transformer_per_year"])
        costs["cumulated"]["c_maint_gc_annual"] += c_maint_gc_annual

        # Store the individual costs of the GC as well
        costs[gcID]["c_gcs"] = c_gc + c_transformer
        costs[gcID]["c_gc_annual"] = c_gc_annual
        costs[gcID]["c_maint_gc_annual"] = c_maint_gc_annual

        # STATIONARY STORAGE
        # assume there is a stationary storage
        try:
            # calculate costs of stationary storage
            c_stat_storage = (
                    c_params["stationary_storage"]["capex_fix"] +
                    c_params["stationary_storage"]["capex_per_kWh"] *
                    schedule.stations[gcID]["battery"]["capacity"])
            costs["cumulated"]["c_stat_storage"] += c_stat_storage
            costs[gcID]["c_stat_storage"] = c_stat_storage

            costs[gcID]["c_stat_storage_annual"] = (
                    costs[gcID]["c_stat_storage"] / c_params["stationary_storage"][
                        "lifetime_stat_storage"])
            costs[gcID]["c_maint_stat_storage_annual"] = (
                    costs[gcID]["c_stat_storage"] * c_params["stationary_storage"][
                        "c_maint_stat_storage_per_year"])
        except KeyError:
            # if no stationary storage at grid connector: cost is 0
            pass

        # FEED IN (local renewables)
        # assume there is a feed in
        try:
            # calculate costs of feed in
            c_feed_in = (c_params["feed_in"]["capex_fix"] +
                         c_params["feed_in"]["capex_per_kW"] *
                         schedule.stations[gcID]["energy_feed_in"]["nominal_power"])
            costs["cumulated"]["c_feed_in"] += c_feed_in
            costs[gcID]["c_feed_in"] = c_feed_in

            costs["cumulated"]["c_feed_in_annual"] = (
                    costs["c_feed_in"] / c_params["feed_in"]["lifetime_feed_in"])
            costs["cumulated"]["c_maint_feed_in_annual"] = (
                    costs["c_feed_in"] * c_params["feed_in"]["c_maint_feed_in_per_year"])
        except KeyError:
            # if no feed in at grid connector: cost is 0
            pass

    costs["cumulated"]["c_stat_storage_annual"] = (
            costs["cumulated"]["c_stat_storage"] / c_params["stationary_storage"][
                "lifetime_stat_storage"])
    costs["cumulated"]["c_maint_stat_storage_annual"] = (
            costs["cumulated"]["c_stat_storage"] * c_params["stationary_storage"][
                "c_maint_stat_storage_per_year"])

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
                    c_cs = c_params["cs"]["capex_deps_per_kW"] * n_cs * defined_max_power
                else:
                    # get sum of installed power at electrified stations without predefined nr of CS
                    sum_power_at_station = sum(cs["max_power"] for cs in cs_at_station)
                    # calculate costs with sum of installed power at electrified station
                    c_cs = c_params["cs"]["capex_deps_per_kW"] * sum_power_at_station
                    ["cumulated"]

            # opportunity charging station - nr of CS depend on max nr of simultaneously occupied CS
            elif station["type"] == "opps":
                gc_timeseries = getattr(scenario, f"{gcID}_timeseries")
                # get nr of CS in use at station
                n_cs = max(gc_timeseries["# CS in use [-]"])
                # get max. defined power for opps or else max. power at opps in general
                defined_max_power = station.get("cs_power_opps", vars(args)["cs_power_opps"])
                # calculate costs with nr of CS at electrified station and its defined max power
                c_cs = c_params["cs"]["capex_opps_per_kW"] * n_cs * defined_max_power
            else:
                warnings.warn(f"Electrified station {station} must be deps or opps. "
                              "Unable to calculate investment costs for this station.")
                c_cs = 0
            costs["cumulated"]["c_cs"] += c_cs
            costs[gcID]["c_cs"] = c_cs

    # calculate annual cost of charging stations, depending on their lifetime
    for key in costs:
        costs[key]["c_cs_annual"] = costs[key]["c_cs"] / c_params["cs"]["lifetime_cs"]

        costs[key]["c_maint_cs_annual"] = costs[key]["c_cs"] * c_params["cs"]["c_maint_cs_per_year"]

        # GARAGE
        costs[key]["c_garage_cs"] = (
                c_params["garage"]["n_charging_stations"]
                * c_params["garage"]["power_cs"] * c_params["cs"]["capex_deps_per_kW"])

        costs[key]["c_garage_workstations"] = -(
                -sum(schedule.vehicle_type_counts.values())
                // c_params["garage"]["vehicles_per_workstation"]
                * c_params["garage"]["cost_per_workstation"])
        costs[key]["c_garage"] = costs[key]["c_garage_cs"] + costs[key]["c_garage_workstations"]
        costs[key]["c_garage_annual"] = (
                costs[key]["c_garage_cs"] / c_params["cs"]["lifetime_cs"]
                + costs[key]["c_garage_workstations"] / c_params["garage"]["lifetime_workstations"])

        # TOTAL COSTS
        costs[key]["c_invest"] = (
                costs[key]["c_vehicles"] + costs[key]["c_cs"] + costs[key]["c_gcs"] +
                costs[key]["c_garage"] +
                costs[key]["c_stat_storage"] + costs[key]["c_feed_in"])
        costs[key]["c_invest_annual"] = (
                costs[key]["c_vehicles_annual"] + costs[key]["c_cs_annual"] +
                costs[key]["c_gcs_annual"] + costs[key]["c_garage_annual"] +
                costs[key]["c_stat_storage_annual"] + costs[key]["c_feed_in_annual"])

        # MAINTENANCE COSTS #

        costs[key]["c_maint_infrastructure_annual"] = (costs[key]["c_maint_cs_annual"] +
                                                       costs[key]["c_maint_gc_annual"] +
                                                       costs[key]["c_maint_stat_storage_annual"] +
                                                       costs[key]["c_maint_feed_in_annual"])

        # calculate (ceil) number of days in scenario
        drive_days = -(-(schedule.scenario["scenario"]["n_intervals"] *
                         schedule.scenario["scenario"]["interval"]) // (24 * 60))
        for rot in rotations_per_gc[key]:
            v_type_rot = f"{rot.vehicle_type}_{rot.charging_type}"
            try:
                costs[key]["c_maint_vehicles_annual"] += (rot.distance / 1000 /
                                                          drive_days * 365) * \
                                                         c_params["vehicles"][v_type_rot][
                                                             "c_maint_per_km"]
            except KeyError:
                warnings.warn("No maintenance costs defined for vehicle type " +
                              c_params["vehicles"][v_type_rot] +
                              ". Unable to calculate maintenance costs for this vehicle type.")
        costs[key]["c_maint_annual"] = (costs[key]["c_maint_infrastructure_annual"] +
                                        costs[key]["c_maint_vehicles_annual"])

    # ELECTRICITY COSTS #

    for gcID, gc in scenario.components.grid_connectors.items():
        pv = sum([pv.nominal_power for pv in scenario.components.photovoltaics.values()
                  if pv.parent == gcID])
        timeseries = vars(scenario).get(f"{gcID}_timeseries")

        # calculate costs for electricity
        costs_electricity = calc_costs_spice_ev(
            strategy=vars(args)["strategy_" + station.get("type")],
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
        costs[gcID]["c_el_procurement_annual"] = costs_electricity[
            'power_procurement_costs_per_year']
        costs[gcID]["c_el_power_price_annual"] = costs_electricity['capacity_costs_eur']
        costs[gcID]["c_el_energy_price_annual"] = costs_electricity['commodity_costs_eur_per_year']
        costs[gcID]["c_el_taxes_annual"] = costs_electricity['levies_fees_and_taxes_per_year']
        costs[gcID]["c_el_feed_in_remuneration_annual"] = costs_electricity[
            'feed_in_remuneration_per_year']
        costs[gcID]["c_el_annual"] = costs_electricity['total_costs_per_year']

        # setattr(scenario.components.grid_connectors[gcID], "costs_electricity", costs_electricity)
        costs["cumulated"]["c_el_procurement_annual"] += costs_electricity[
            'power_procurement_costs_per_year']
        costs["cumulated"]["c_el_power_price_annual"] += costs_electricity['capacity_costs_eur']
        costs["cumulated"]["c_el_energy_price_annual"] += costs_electricity[
            'commodity_costs_eur_per_year']
        costs["cumulated"]["c_el_taxes_annual"] += costs_electricity[
            'levies_fees_and_taxes_per_year']
        costs["cumulated"]["c_el_feed_in_remuneration_annual"] += costs_electricity[
            'feed_in_remuneration_per_year']
        costs["cumulated"]["c_el_annual"] += costs_electricity['total_costs_per_year']

    # round all costs to ct
    for key1 in costs:
        for key2, value in costs[key1].items():
            costs[key1][key2] = round(value, 2)

    for key in costs["cumulated"]:
        try:
            units[key]
        except KeyError:
            # if no unit is set use this algorithm
            if "annual" in key:
                units[key] = "€/year"
            else:
                units[key] = "€"

    cost_object = CostsObject(costs, units)

    logging.info(
        "\nTotal costs:\n"
        f"Investment cost: {costs['cumulated']['c_invest']} €. \n"
        f"Annual investment costs: {costs['cumulated']['c_invest_annual']} €/a. \n"
        f"Annual maintenance costs: {costs['cumulated']['c_maint_annual']} €/a. \n"
        f"Annual costs for electricity: {costs['cumulated']['c_el_annual']} €/a.\n")

    setattr(scenario, "costs", cost_object)
