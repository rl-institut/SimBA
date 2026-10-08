import logging
import traceback
import warnings

import simba.schedule

import spice_ev.costs
import spice_ev.report
import spice_ev.scenario
import spice_ev.util


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
    :return: cost object
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

    # Correct electricity costs for storage level drift between start and end of simulation
    cost_object.set_storage_level_cost_correction()

    # Calculate total annual costs
    cost_object.set_total_annual_costs()

    # Calculate the costs per kilometer
    cost_object.set_invest_per_total_km()

    # Set KPIs like total_energy per station, km of rotations per station, station type and more
    cost_object.set_various_kpis()

    # Cumulate the costs of all gcs
    cost_object.cumulate()

    logging.info(cost_object.info())

    return cost_object


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

    DAYS_PER_YEAR = 365

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
                f"Annual maintenance costs: {cumulated['c_maint_annual']} €/a. \n"
                f"Annual costs for electricity: {cumulated['c_el_annual']} €/a.\n"
                f"Costs per km: {cumulated['c_total_per_km']} €/km. \n"
                f"Total kilometers per year: {self.get_total_annual_km()} km. \n")

    def get_unit(self, key):
        """ Get the unit for annual or non-annual costs.

        :param key: Key to get the unit for
        :type key: str
        :return: Unit associated with the key
        :rtype: str
        """
        if "annual" in key:
            return "€/year"
        elif "station_type" == key:
            return "[-]"
        elif "maximum Nr charging stations" == key:
            return "[-]"
        elif "total_km_per_year" == key:
            return "km/year"
        elif "power" in key:
            return "kW"
        elif "per_km" in key:
            return "€/km"
        elif "kWh" in key:
            return "kWh/year"
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
        return [
            # various kpis and attributes
            "station_type", "maximum_gc_power", "max_power_in_plw", "maximum Nr charging stations",
            "total_km_per_year", "annual_kWh_from_grid", "annual_kWh_from_feed_in",

            # investment costs
            "c_vehicles", "c_gcs", "c_cs", "c_garage_cs",
            "c_garage_workstations",
            "c_stat_storage", "c_feed_in", "c_invest",
            # annual investment costs
            "c_vehicles_annual", "c_gcs_annual", "c_cs_annual",
            "c_stat_storage_annual", "c_feed_in_annual", "c_invest_annual",
            # annual maintenance costs
            "c_maint_vehicles_annual", "c_maint_gc_annual", "c_maint_cs_annual",
            "c_maint_stat_storage_annual", "c_maint_feed_in_annual",
            "c_maint_infrastructure_annual",
            "c_maint_annual",
            # annual electricity costs
            "c_el_procurement_annual", "c_el_power_price_annual",
            "c_el_energy_price_annual",
            "c_el_taxes_annual", "c_el_feed_in_remuneration_annual",
            "c_el_storage_level_correction_annual", "c_el_annual",
            # total annual costs
            "c_total_annual",
            # costs per kilometer
            "c_vehicles_per_km", "c_invest_per_km", "c_maint_vehicles_per_km",
            "c_maint_infrastructure_per_km",
            "c_maint_per_km", "c_el_per_km", "c_total_per_km"]

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
            scenario_duration = self.scenario.stop_time - self.scenario.start_time
            drive_days = -(-scenario_duration.total_seconds() // (60*60*24))
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

            # use procurement and commodity costs read from CSV instead of SpiceEV prices, if exist
            prices = station.get("prices", timeseries.get("price [EUR/kWh]"))

            # Get the calculation method from args.
            cost_calculation_name = "cost_calculation_method_" + station.get("type")
            cost_calculation_method = vars(self.args).get(cost_calculation_name)
            if cost_calculation_method is None:
                # not given in config: use the same method as the charging strategy
                default_strategy = vars(self.args)["strategy_" + station.get("type")]
                cost_calculation_method = spice_ev.costs.DEFAULT_COST_CALCULATION[default_strategy]

            # calculate costs for electricity
            try:
                is_peak_load_window = cost_calculation_method.endswith("w_plw")
                if is_peak_load_window and timeseries.get("window signal [-]") is None:
                    logging.info("Generating peak load window signal for cost calculation")
                    timeseries["window signal [-]"] = spice_ev.util.get_time_windows_from_json(
                        self.args.time_windows, gc.grid_operator, gc.voltage_level, self.scenario)
                costs_electricity = spice_ev.costs.calculate_costs(
                    cc_type=cost_calculation_method,
                    voltage_level=gc.voltage_level,
                    interval=self.scenario.interval,
                    timestamps_list=timeseries.get("time"),
                    power_grid_supply_list=timeseries.get("grid supply [kW]"),
                    price_list=prices,
                    power_fix_load_list=timeseries.get("fixed load [kW]"),
                    power_generation_feed_in_list=timeseries.get("generation feed-in [kW]"),
                    power_v2g_feed_in_list=timeseries.get("V2G feed-in [kW]"),
                    power_battery_feed_in_list=timeseries.get("battery feed-in [kW]"),
                    window_signal_list=timeseries.get("window signal [-]"),
                    price_sheet_path=self.args.cost_parameters_path,
                    grid_operator=gc.grid_operator,
                    fee_type=None,  # "RLM" or "SLP" or None (based on energy consumption)
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
            self.costs_per_gc[gcID]["c_el_annual"] = costs_electricity.get(
                'total_costs_per_year', error_value)
            self.costs_per_gc[gcID]["max_power_in_plw"] = costs_electricity.get(
                'peak_power_in_windows', error_value)

        return self

    def set_storage_level_cost_correction(self):
        """ Correct energy-proportional electricity costs for storage level drift.

        If a stationary battery or a vehicle battery ends the simulation with more or less
        stored energy than it started with, part of the energy bought (or still needed) is not
        attributable to bus operation within the simulated period. Since electricity costs are
        annualized based on the simulated period (see :func:`set_various_kpis`), such a one-off
        difference would otherwise be scaled up (or down) as if it recurred every such period
        over the whole year, distorting comparisons between scenarios with different storage
        behavior (e.g. with and without a stationary battery).

        Surplus energy (more stored at the end than at the start) is valued at a moving-average
        price of the device's own charging history, analogous to weighted-average inventory
        costing: already-consumed, historically cheap charging does not count towards the value
        of what is still in the device at the end. A deficit (less stored at the end) is valued
        at the plain volume-weighted average price of the device's own charging events, as an
        estimate of what it would typically cost this device to recharge.

        Only the energy-proportional part of the electricity costs is affected
        (``c_el_storage_level_correction_annual``, folded into ``c_el_annual``); peak-power based
        costs and investment/maintenance costs are untouched. Only has an effect if
        ``args.clean_costs_storage_level`` is set.

        :return: self
        :rtype: Costs
        """
        if not vars(self.args).get("clean_costs_storage_level"):
            return self

        scenario = self.scenario
        scenario_duration_s = scenario.n_intervals * scenario.interval.total_seconds()
        if scenario_duration_s <= 0:
            return self
        annual_factor = (self.DAYS_PER_YEAR * 24 * 60 * 60) / scenario_duration_s

        if not hasattr(scenario, "vehicle_socs"):
            # not yet generated (usually happens later, when report.generate() writes
            # vehicle_socs.csv) - reuse the same logic here instead of duplicating it
            spice_ev.report.generate_soc_timeseries(scenario)

        price_series_cache = {}

        def get_price_series(gc_id):
            """ EUR/kWh at each timestep for a gc. For a "fixed" cost_calculation_method the
            price is time-invariant by definition (a static price sheet, see
            :func:`set_electricity_costs`) - the on-the-fly price SpiceEV sees during the
            simulation is only a placeholder used to give a charging strategy something to
            react to and is not the authoritative billed price. For that case, the flat price
            actually implied by the electricity costs already calculated for this gc is used
            instead. Otherwise (a genuine "variable"/market-based method), the CSV-provided
            price if given, else SpiceEV's own (ct/kWh) price signal, is used. """
            if gc_id in price_series_cache:
                return price_series_cache[gc_id]
            n = scenario.n_intervals
            station = self.schedule.stations.get(gc_id, {})
            station_type = station.get("type", "")
            cost_calculation_method = vars(self.args).get(
                f"cost_calculation_method_{station_type}")
            if cost_calculation_method is None:
                default_strategy = vars(self.args).get(f"strategy_{station_type}")
                cost_calculation_method = spice_ev.costs.DEFAULT_COST_CALCULATION.get(
                    default_strategy)

            if cost_calculation_method and cost_calculation_method.startswith("fixed"):
                annual_kwh = (sum(scenario.totalLoad.get(gc_id, []))
                              / scenario.stepsPerHour * annual_factor)
                energy_costs = self.costs_per_gc.get(gc_id, {})
                annual_cost = (energy_costs.get("c_el_procurement_annual", 0)
                               + energy_costs.get("c_el_energy_price_annual", 0)
                               + energy_costs.get("c_el_taxes_annual", 0))
                flat_price = annual_cost / annual_kwh if annual_kwh else 0
                prices = [flat_price] * n
            else:
                prices = station.get("prices")
                if not prices:
                    prices = [p / 100 for p in scenario.prices.get(gc_id, [0] * n)]
            price_series_cache[gc_id] = prices
            return prices

        correction_per_gc = {gc_id: 0 for gc_id in self.gcs}

        # stationary batteries: fixed to a single gc
        for bat_id, battery in scenario.components.batteries.items():
            gc_id = battery.parent
            if gc_id not in self.gcs:
                continue
            prices = get_price_series(gc_id)
            soc_kwh = scenario.batteryLevels.get(bat_id)
            if not soc_kwh:
                continue
            correction = self._storage_correction_eur(soc_kwh, prices, battery.efficiency)
            if correction:
                correction_per_gc[gc_id] += correction

        # vehicle batteries: attributed to the gc of the vehicle's last started rotation
        last_rotation_gc = self._last_rotation_gc_per_vehicle()
        for v_id, vehicle in scenario.components.vehicles.items():
            soc_fractions = scenario.vehicle_socs.get(v_id)
            if not soc_fractions:
                continue
            capacity = vehicle.battery.capacity
            soc_kwh = [None if s is None else s * capacity for s in soc_fractions]
            prices = self._vehicle_price_series(v_id, get_price_series)
            correction = self._storage_correction_eur(soc_kwh, prices, vehicle.battery.efficiency)
            if correction:
                gc_id = last_rotation_gc.get(v_id, self.NOT_ELECTRIFIED)
                correction_per_gc.setdefault(gc_id, 0)
                correction_per_gc[gc_id] += correction

        # apply, annualize and log the change for every affected gc
        # (self.costs_per_gc[CUMULATED] is only populated later, in cumulate(), so the total
        # before/after is accumulated here instead of read from there)
        cumulated_correction_annual = 0
        total_before = 0
        for gc_id, correction in correction_per_gc.items():
            if not correction or gc_id not in self.costs_per_gc:
                continue
            correction_annual = correction * annual_factor
            before = self.costs_per_gc[gc_id]["c_el_annual"]
            after = before + correction_annual
            self.costs_per_gc[gc_id]["c_el_storage_level_correction_annual"] = correction_annual
            self.costs_per_gc[gc_id]["c_el_annual"] = after
            cumulated_correction_annual += correction_annual
            total_before += before
            logging.info(
                f"{gc_id}: electricity costs corrected for storage level drift: "
                f"{before:.2f} €/a -> {after:.2f} €/a (correction: {correction_annual:+.2f} €/a)")

        if cumulated_correction_annual:
            total_after = total_before + cumulated_correction_annual
            logging.info(
                "Cumulated electricity costs corrected for storage level drift: "
                f"{total_before:.2f} €/a -> {total_after:.2f} €/a "
                f"(correction: {cumulated_correction_annual:+.2f} €/a)")

        return self

    def _last_rotation_gc_per_vehicle(self):
        """ Find the grid connector each vehicle's last started rotation departs from.

        Used to attribute a vehicle's storage level correction to a single grid connector,
        ignoring any intermediate opportunity charging elsewhere.

        :return: mapping of vehicle_id to grid connector id (or NOT_ELECTRIFIED)
        :rtype: dict
        """
        last_departure = {}
        last_gc = {}
        for rot in self.schedule.rotations.values():
            if rot.vehicle_id is None:
                continue
            if (rot.vehicle_id not in last_departure
                    or rot.departure_time > last_departure[rot.vehicle_id]):
                last_departure[rot.vehicle_id] = rot.departure_time
                gc_id = rot.departure_name if rot.departure_name in self.gcs \
                    else self.NOT_ELECTRIFIED
                last_gc[rot.vehicle_id] = gc_id
        return last_gc

    def _vehicle_price_series(self, v_id, get_price_series):
        """ Build a per-timestep price series (EUR/kWh) for whichever gc a vehicle is
        connected to at each timestep.

        :param v_id: vehicle id
        :type v_id: str
        :param get_price_series: callable returning the price series (EUR/kWh) for a gc id
        :type get_price_series: callable
        :return: price at each timestep, None where the vehicle is not connected or the price
            is unknown
        :rtype: list
        """
        scenario = self.scenario
        charging_stations = scenario.components.charging_stations
        prices = []
        for conn in scenario.connected:
            cs_id = conn.get(v_id)
            price = None
            if cs_id is not None:
                gc_id = charging_stations[cs_id].parent
                gc_prices = get_price_series(gc_id)
                idx = len(prices)
                if gc_prices is not None and idx < len(gc_prices):
                    price = gc_prices[idx]
            prices.append(price)
        return prices

    @staticmethod
    def _storage_correction_eur(soc_kwh, prices, efficiency):
        """ Electricity-cost correction (EUR, simulated period, not annualized) for one storage
        device's change in stored energy between the start and end of the simulation.

        :param soc_kwh: energy content of the device at each timestep (kWh), None where unknown
            (e.g. leading/trailing timesteps of a vehicle not yet or no longer tracked)
        :type soc_kwh: list
        :param prices: price (EUR/kWh) at each timestep, aligned with soc_kwh. Only prices at
            timesteps where the device actually charges are used.
        :type prices: list
        :param efficiency: charging efficiency of the device
        :type efficiency: float
        :return: correction in EUR for the simulated period (positive: cost increase, negative:
            cost decrease), or None if it could not be determined (e.g. device never charged)
        :rtype: float or None
        """
        known_idx = [i for i, v in enumerate(soc_kwh) if v is not None]
        if len(known_idx) < 2:
            return None
        start_idx, end_idx = known_idx[0], known_idx[-1]

        # device's own volume-weighted average price of its charging events (used directly for
        # a deficit, and as the seed value of the moving average for a surplus)
        charged_energy = 0
        charged_cost = 0
        for t in range(start_idx, end_idx):
            delta = soc_kwh[t + 1] - soc_kwh[t]
            if delta > 0 and prices[t] is not None:
                grid_energy = delta / efficiency
                charged_energy += grid_energy
                charged_cost += grid_energy * prices[t]
        if charged_energy <= 0:
            # device never charged during the simulation: no price information available
            return None
        own_avg_price = charged_cost / charged_energy

        delta_e = soc_kwh[end_idx] - soc_kwh[start_idx]
        if delta_e == 0:
            return 0

        if delta_e < 0:
            # deficit: estimate cost of recharging at the device's own average price
            return (-delta_e / efficiency) * own_avg_price

        # surplus: moving average ("gleitender Durchschnittspreis") of what remains, tracked as
        # EUR per kWh of energy actually stored in the device (not per grid kWh, unlike
        # own_avg_price above) - charging efficiency is priced in once, when a charging event
        # is blended in, and must not be applied again when valuing the final surplus.
        # Seeded with the device's own average price, converted to the same per-stored-kWh
        # basis, as the assumed cost of the energy present at the start.
        stock = soc_kwh[start_idx]
        avg_price = own_avg_price / efficiency
        for t in range(start_idx, end_idx):
            delta = soc_kwh[t + 1] - soc_kwh[t]
            if delta > 0:
                price = prices[t] if prices[t] is not None else own_avg_price
                price_per_stored_kwh = price / efficiency
                new_stock = stock + delta
                avg_price = (stock * avg_price + delta * price_per_stored_kwh) / new_stock
                stock = new_stock
            elif delta < 0:
                stock += delta  # withdrawal: avg_price unchanged
        return -delta_e * avg_price

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

    def set_various_kpis(self):
        """ Set several kpis.

        :return: self
        :rtype: Costs
        """
        stepsPerHour = self.scenario.stepsPerHour
        # factor for scaling to a year
        scenario_duration_s = self.scenario.n_intervals * self.scenario.interval.total_seconds()
        if scenario_duration_s > 0:
            annual_factor = (self.DAYS_PER_YEAR * 24 * 60 * 60) / scenario_duration_s
        else:
            annual_factor = 0
        # get dict with costs for specific grid operator
        for gcID, gc in self.gcs.items():
            station = self.schedule.stations[gcID]
            self.costs_per_gc[gcID]["station_type"] = station["type"]
            timeseries = getattr(self.scenario, f"{gcID}_timeseries")
            self.costs_per_gc[gcID]["maximum_gc_power"] = abs(min(timeseries["grid supply [kW]"]))
            self.costs_per_gc[gcID]["annual_kWh_from_grid"] = \
                (sum(timeseries["grid supply [kW]"]) / stepsPerHour) * annual_factor
            self.costs_per_gc[gcID]["annual_kWh_from_feed_in"] = \
                -sum(timeseries.get("local generation [kW]", [0])) / stepsPerHour * annual_factor
            cs_at_station = max(getattr(self.scenario, f"{gcID}_timeseries")["# CS in use [-]"])
            self.costs_per_gc[gcID]["maximum Nr charging stations"] = cs_at_station

        # total_km_per_year
        rotations_per_station = {gcid: [] for gcid in self.gcs_and_garage}
        for rotation in self.schedule.rotations.values():
            try:
                rotations_per_station[rotation.departure_name].append(rotation)
            except KeyError:
                # rotation starts at a non electrified station
                rotations_per_station[self.NOT_ELECTRIFIED].append(rotation)
        for gcID, rotations in rotations_per_station.items():
            self.costs_per_gc[gcID]["total_km_per_year"] = \
                (sum([r.distance for r in rotations]) * annual_factor) / 1000

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

    def set_total_annual_costs(self):
        """ Calculate the total annual investment.

        :return: self
        :rtype: Costs
        """
        attributes = [
            "c_invest_annual", "c_maint_annual", "c_el_annual"
        ]
        for gcID in self.gcs:
            for key in attributes:
                self.costs_per_gc[gcID]["c_total_annual"] += self.costs_per_gc[gcID][key]
        for key in attributes:
            self.costs_per_gc[self.GARAGE]["c_total_annual"] += self.costs_per_gc[self.GARAGE][key]
        return self

    def set_invest_per_total_km(self):
        """ Calculate the annual investment costs per total driven distance.

        :return: self
        :rtype: Costs
        """
        total_annual_km = self.get_total_annual_km()
        if total_annual_km:
            attributes = [
                "c_vehicles", "c_invest", "c_maint_vehicles", "c_maint_infrastructure",
                "c_maint", "c_el", "c_total"
            ]
            for gcID in self.gcs:
                for key in attributes:
                    self.costs_per_gc[gcID][f"{key}_per_km"] = (
                            self.costs_per_gc[gcID][f"{key}_annual"] / total_annual_km)
            self.costs_per_gc[self.GARAGE]["c_invest_per_km"] = (
                    self.costs_per_gc[self.GARAGE]["c_invest_annual"] /
                    total_annual_km)
            self.costs_per_gc[self.GARAGE]["c_total_per_km"] = (
                    self.costs_per_gc[self.GARAGE]["c_total_annual"] /
                    total_annual_km)
        return self

    def cumulate(self):
        """ Cumulate the costs of vehicles and infrastructure.

        :return: Updated instance
        :rtype: Costs
        """
        # Cumulate gcs variables
        for key in self.get_gc_cost_variables():
            # skip special attributes
            if key in ["station_type", "maximum_gc_power", "max_power_in_plw"]:
                continue
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
            self.costs_per_gc[self.CUMULATED]["c_invest_per_km"] = \
                self.costs_per_gc[self.CUMULATED]["c_invest_annual"] / total_annual_km
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
            row = [key, self.get_unit(key)]
            for col in self.get_columns():
                # Get the cost at this gc. Stations which have no gc get 0
                num = self.costs_per_gc.get(col, {}).get(key, 0)
                try:
                    row.append(round(num, self.rounding_precision))
                except TypeError:
                    row.append(num)
            output.append(row)
        return output
