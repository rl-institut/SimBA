import logging
import traceback
import warnings

import spice_ev.scenario
from spice_ev.costs import calculate_costs as calc_costs_spice_ev

import simba.schedule


def calculate_costs(c_params, scenario, schedule, args):
    """ Calculates annual costs of all necessary vehicles and infrastructure.

    :param c_params: Infrastructure component costs and specific grid operator costs
    :type c_params: dict
    :param scenario: Information about the simulated scenario and all used parameters
    :type scenario: Scenario
    :param schedule: Information about the whole bus schedule and all used parameters
    :type schedule: Schedule
    :param args: Configuration arguments specified in config files contained in configs directory
    :type args: argparse.Namespace
    """
    cost_object = Costs(schedule, scenario, args, c_params)

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

    # Get electricity costs from SpiceEV
    cost_object.set_electricity_costs()

    # Calculate the annual costs per kilometer
    cost_object.set_annual_invest_per_km()

    # Cumulate the costs of all gcs
    cost_object.cumulate()

    logging.info(cost_object.info())

    setattr(scenario, "costs", cost_object)


class Costs:
    """ Class representing costs associated with vehicles and infrastructure.

    :param schedule: Simulation schedule
    :type schedule: simba.schedule.Schedule
    :param scenario: Scenario for the simulation
    :type scenario: spice_ev.scenario.Scenario
    :param args: Arguments for the simulation
    :param c_params: Parameters for the simulation costs.
    :type c_params: Dict
    """
    CUMULATED = "cumulated"
    GARAGE = "garage"
    NOT_ELECTRIFIED = "Non_electrified_station"

    DAYS_PER_YEAR = 365.2422

    def __init__(self, schedule: simba.schedule.Schedule, scenario: spice_ev.scenario.Scenario,
                 args, c_params: dict):
        """ Initialize the Costs instance.

        :param schedule: Simulation schedule
        :type schedule: simba.schedule.Schedule
        :param scenario: Scenario for the simulation
        :type scenario: spice_ev.scenario.Scenario
        :param args: Arguments for the simulation
        :type args: Namespace
        :param c_params: Parameters for the simulation costs
        :type c_params: dict
        """

        self.gcs = scenario.components.grid_connectors

        # Make sure station names are not reserved keywords
        for reserved in [self.CUMULATED, self.GARAGE, self.NOT_ELECTRIFIED]:
            assert reserved not in self.gcs, f"{reserved} must not be part of station names"

        self.gcs_and_garage = [self.GARAGE, self.NOT_ELECTRIFIED] + list(self.gcs)
        self.costs_per_gc = {gc: {key: 0 for key in self.get_gc_cost_variables()} for gc in
                             self.gcs_and_garage + [self.CUMULATED]}
        self.units: dict = None
        self.vehicles_per_gc: dict = None
        self.gc_rotations: dict = None
        self.schedule = schedule
        self.scenario = scenario
        self.args = args
        self.params = c_params
        self.rounding_precision = 2

    def info(self):
        """ Provide information about the total costs.

        :return: Formatted string containing information about total costs
        :rtype: str
        """
        cumulated = self.costs_per_gc[self.CUMULATED]
        return ("\nTotal costs:\n"
                f"Investment cost: {cumulated['c_invest']} €. \n"
                f"Annual investment costs: {cumulated['c_invest_annual']} €/a. \n"
                f"Annual investment costs per km: {cumulated['c_invest_annual_per_km']} €/a. \n"
                f"Annual maintenance costs: {cumulated['c_maint_annual']} €/a. \n"
                f"Annual costs for electricity: {cumulated['c_el_annual']} €/a.\n")

    def get_annual_or_not(self, key):
        """ Get the unit for annual or non-annual costs.

        :param key: Key to get the unit for
        :type key: str
        :return: Unit associated with the key
        :rtype: str
        """
        if "annual" in key:
            return "€/year"
        else:
            return "€"

    def get_columns(self):
        """ Get the sorted columns for cost calculations.

        :return: List of columns
        :rtype: list
        """
        first_columns = [self.CUMULATED, self.GARAGE, self.NOT_ELECTRIFIED]

        # Use sorting of electrified_stations.json. Return all stations, even if they have no gc
        # in the scenario to allow for consistent output generation
        middle_columns = list(self.schedule.stations.keys())

        # if there are gcs that are not a station, add them at the end
        trailing_columns = [s for s in self.gcs if s not in self.schedule.stations]
        return first_columns + middle_columns + trailing_columns

    def set_gc_rotations(self):
        """ Set grid connector rotations.

        :return: Updated instance
        :rtype: Costs
        """
        self.gc_rotations = {
            gc: [rot for rot in self.schedule.rotations.values() if rot.departure_name == gc]
            for gc in self.gcs}
        return self

    def get_gc_cost_variables(self):
        """ Get the list of variables representing different cost types.

         :return: List of cost variables
         :rtype: list
         """
        return [  # investment costs
            "c_vehicles", "c_gcs", "c_cs", "c_garage_cs",
            "c_garage_workstations",
            "c_stat_storage", "c_feed_in", "c_invest",
            # annual investment costs
            "c_vehicles_annual", "c_gcs_annual", "c_cs_annual",
            "c_stat_storage_annual", "c_feed_in_annual", "c_invest_annual",
            # annual investment costs per km
            "c_vehicles_annual_per_km", "c_gcs_annual_per_km", "c_cs_annual_per_km",
            "c_stat_storage_annual_per_km", "c_feed_in_annual_per_km", "c_invest_annual_per_km",
            # annual maintenance costs
            "c_maint_vehicles_annual", "c_maint_gc_annual", "c_maint_cs_annual",
            "c_maint_stat_storage_annual", "c_maint_feed_in_annual",
            "c_maint_infrastructure_annual",
            "c_maint_annual",
            # annual electricity costs
            "c_el_procurement_annual", "c_el_power_price_annual",
            "c_el_energy_price_annual",
            "c_el_taxes_annual", "c_el_feed_in_remuneration_annual", "c_el_annual"]

    def get_vehicle_types(self):
        """ Get the types of vehicles used.

         :return: Vehicle types
         :rtype: list
         """
        return [k for k, v in self.schedule.vehicle_type_counts.items() if v > 0]

    def set_charging_infrastructure_costs(self):
        """ Calculate the costs associated with charging infrastructure at each grid connector.

        :return: self
        :rtype: Costs
        """

        all_stations = self.schedule.scenario["components"]["charging_stations"]
        # find all stations from the schedule which are grid connectors
        used_stations = {gcId: station for gcId, station in self.schedule.stations.items()
                         if gcId in self.gcs}

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
                    c_maint = ((rot.distance / 1000 / drive_days * 365) *
                               self.params["vehicles"][v_type_rot]["c_maint_per_km"])
                    self.costs_per_gc[gcID]["c_maint_vehicles_annual"] += c_maint
                except KeyError:
                    warnings.warn("No maintenance costs defined for vehicle type " +
                                  self.params["vehicles"][v_type_rot] +
                                  ". Unable to calculate maintenance costs for this vehicle type.")
            self.costs_per_gc[gcID]["c_maint_annual"] = (
                        self.costs_per_gc[gcID]["c_maint_infrastructure_annual"] +
                        self.costs_per_gc[gcID]["c_maint_vehicles_annual"])
        return self

    def set_electricity_costs(self):
        """ Calculate the electricity costs at each grid connector.

        :return: self
        :rtype: Costs
        """
        for gcID, gc in self.gcs.items():
            station = self.schedule.stations[gcID]
            pv = sum([pv.nominal_power for pv in self.scenario.components.photovoltaics.values()
                      if pv.parent == gcID])
            timeseries = vars(self.scenario).get(f"{gcID}_timeseries")

            # calculate costs for electricity
            try:
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
                    charging_signal_list=timeseries.get("window signal [-]"),
                    price_sheet_path=self.args.cost_parameters_file,
                    grid_operator=gc.grid_operator,
                    power_pv_nominal=pv,
                )
            except Exception:
                costs_electricity = dict()
                logging.warning(f"SpiceEV calculation of costs for {gcID} failed due to "
                                f"{traceback.format_exc()}.")
            error_value = 0
            self.costs_per_gc[gcID]["c_el_procurement_annual"] = costs_electricity.get(
                'power_procurement_costs_per_year', error_value)
            self.costs_per_gc[gcID]["c_el_power_price_annual"] = costs_electricity.get(
                'capacity_costs_eur', error_value)
            self.costs_per_gc[gcID]["c_el_energy_price_annual"] = costs_electricity.get(
                'commodity_costs_eur_per_year', error_value)
            self.costs_per_gc[gcID]["c_el_taxes_annual"] = costs_electricity.get(
                'levies_fees_and_taxes_per_year', error_value)
            self.costs_per_gc[gcID]["c_el_feed_in_remuneration_annual"] = costs_electricity.get(
                'feed_in_remuneration_per_year', error_value)
            self.costs_per_gc[gcID]["c_el_annual"] = costs_electricity.get('total_costs_per_year',
                                                                           error_value)
        return self

    def set_garage_costs(self):
        """ Calculate the costs associated with the garage.

        Note that although a garage can have charging stations, electricity costs are not used.

        :return: self
        :rtype: Costs
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

        self.costs_per_gc[self.GARAGE]["c_invest"] = c_garage_cs + c_garage_workstations
        self.costs_per_gc[self.GARAGE]["c_invest_annual"] = (
                c_garage_cs / self.params["cs"]["lifetime_cs"]
                + c_garage_workstations / self.params["garage"][
                    "lifetime_workstations"])

        return self

    def set_grid_connection_costs(self):
        """ Calculate the costs of each grid connection.

        :raises Exception: If grid operator of grid connector cannot be found in cost params
        :return: self
        :rtype: Costs
        """
        # get dict with costs for specific grid operator
        for gcID, gc in self.gcs.items():
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

    def set_vehicles_per_gc(self):
        """ Calculate the number of vehicles at each grid connector.

        Calculate vehicle numbers and set vehicles_per_gc before vehicle costs can be calculated.

        :return: self
        :rtype: Costs
        """
        v_types = self.schedule.scenario["components"]["vehicle_types"]
        vehicles_per_gc = {gc: {v_type_name: set() for v_type_name in v_types
                                if self.schedule.vehicle_type_counts[v_type_name] > 0}
                           for gc in self.gcs_and_garage}
        for rot in self.schedule.rotations.values():
            # Get the vehicle_type_name including charging type by removing the number of the id
            vehicle_type_name = "_".join(rot.vehicle_id.split("_")[:-1])
            try:
                vehicles_per_gc[rot.departure_name][vehicle_type_name].add(rot.vehicle_id)
            except KeyError:
                # rotation might start at non-electrified station, therefore not found in gc keys.
                vehicles_per_gc[self.NOT_ELECTRIFIED][vehicle_type_name].add(rot.vehicle_id)

        self.vehicles_per_gc = {
            gc: {v_type_name: len(s) for v_type_name, s in vehicles_dict.items()} for
            gc, vehicles_dict in vehicles_per_gc.items()}

        self.vehicles_per_gc[self.CUMULATED] = dict()
        for v_type in v_types:
            self.vehicles_per_gc[self.CUMULATED][v_type] = self.schedule.vehicle_type_counts[v_type]

        return self

    def set_vehicle_costs_per_gc(self):
        """ Calculate and set the costs associated with vehicles at each grid connector.

        :return: self
        :rtype: Costs
        """
        # Iterate over GCs and use the vehicles which start and end at this gc for cost calculation
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

    def set_annual_invest_per_km(self):
        """ Calculate the annual investment costs per total driven distance.

        :return: self
        :rtype: Costs
        """
        total_annual_km = self.get_total_annual_km()
        if total_annual_km:
            attributes = [
                "c_vehicles_annual", "c_gcs_annual", "c_cs_annual",
                "c_stat_storage_annual", "c_feed_in_annual", "c_invest_annual"
            ]
            for gcID in self.gcs:
                for key in attributes:
                    self.costs_per_gc[gcID][f"{key}_per_km"] = (
                            self.costs_per_gc[gcID][key] / total_annual_km)
            self.costs_per_gc[self.GARAGE]["c_invest_annual_per_km"] = (
                    self.costs_per_gc[self.GARAGE]["c_invest_annual"] /
                    total_annual_km)
        return self

    def cumulate(self):
        """ Cumulate the costs of vehicles and infrastructure.

        :return: Updated instance
        :rtype: Costs
        """
        # Cumulate gcs variables
        for key in self.get_gc_cost_variables():
            # Vehicle costs cannot be cumulated since they might be double counted for vehicles
            # with multiple depots. Instead, the vehicle costs were previously calculated in
            # set_vehicle_costs_per_gc()
            if key in ["c_vehicles", "c_vehicles_annual"]:
                continue
            self.costs_per_gc[self.CUMULATED][key] = 0
            for gc in self.costs_per_gc.keys() - [self.CUMULATED]:
                self.costs_per_gc[self.CUMULATED][key] += self.costs_per_gc[gc][key]

        # Since c_vehicles and c_vehicles_annual might have double counting over the other gcs,
        # total invest for the cumulated gc is calculated separately
        # Garage costs are added separately, since their costs come from the garage_workstation
        # and garage_cs which are different to "normal" stations
        self.costs_per_gc[self.CUMULATED]["c_invest"] = (
                self.costs_per_gc[self.CUMULATED]["c_vehicles"]
                + self.costs_per_gc[self.CUMULATED]["c_cs"]
                + self.costs_per_gc[self.CUMULATED]["c_gcs"]
                + self.costs_per_gc[self.CUMULATED]["c_stat_storage"]
                + self.costs_per_gc[self.CUMULATED]["c_feed_in"]
                + self.costs_per_gc[self.GARAGE]["c_invest"]
        )

        self.costs_per_gc[self.CUMULATED]["c_invest_annual"] = (
                self.costs_per_gc[self.CUMULATED]["c_vehicles_annual"]
                + self.costs_per_gc[self.CUMULATED]["c_cs_annual"]
                + self.costs_per_gc[self.CUMULATED]["c_gcs_annual"]
                + self.costs_per_gc[self.CUMULATED]["c_stat_storage_annual"]
                + self.costs_per_gc[self.CUMULATED]["c_feed_in_annual"]
                + self.costs_per_gc[self.GARAGE]["c_invest_annual"]
        )

        total_annual_km = self.get_total_annual_km()
        if total_annual_km:
            self.costs_per_gc[self.CUMULATED]["c_invest_annual_per_km"] = (
                    self.costs_per_gc[self.CUMULATED]["c_invest_annual"] /
                    total_annual_km
            )

        return self

    def get_total_annual_km(self):
        """ Calculate total annual kilometers driven, linearly extrapolated from schedule distance.

        :return: Total annual kilometers driven
        :rtype: float
        """
        # calculate yearly driven kilometers
        total_km = self.schedule.get_total_distance() / 1000
        # use full days to caclulate yearly km, since schedules can overlap with ones from next day
        simulated_days = max(
            round(self.scenario.n_intervals / self.scenario.stepsPerHour / 24, 0),
            1)
        return self.DAYS_PER_YEAR / simulated_days * total_km

    def to_csv_lists(self):
        """ Convert costs to a list of lists easily convertible to a CSV.

        :return: List of lists of parameters, units and costs per gc
        :rtype: list
        """
        output = [["parameter", "unit"] + self.get_columns()]
        for key in self.get_vehicle_types():
            # The first two columns contain the parameter and unit 'vehicles'
            row = [key, "vehicles"]

            for col in self.get_columns():
                # Get the number of vehicles at this gc. Stations which have no gc get 0
                row.append(self.vehicles_per_gc.get(col, {}).get(key, 0))
            output.append(row)

        # Take a single station and take the cost parameters
        cost_parameters = self.get_gc_cost_variables()
        for key in cost_parameters:
            # The first two columns contain the parameter and unit
            row = [key, self.get_annual_or_not(key)]
            for col in self.get_columns():
                # Get the cost at this gc. Stations which have no gc get 0
                num = self.costs_per_gc.get(col, {}).get(key, 0)
                row.append(round(num, self.rounding_precision))
            output.append(row)

        transposed_output = []
        for column, _ in enumerate(output[0]):
            transposed_output.append([row[column] for row in output])
        return transposed_output
