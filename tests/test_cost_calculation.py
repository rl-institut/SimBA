import logging
from datetime import datetime
from types import SimpleNamespace

import pytest

from tests.test_schedule import BasicSchedule
from simba.util import uncomment_json_file
from simba.costs import Costs, calculate_costs


class TestCostCalculation:
    def test_cost_calculation(self):
        schedule, scenario, args = BasicSchedule().basic_run()
        file = args.cost_parameters_path
        with open(file, "r") as file:
            cost_params = uncomment_json_file(file)

        assert args.strategy_deps == "balanced"
        assert args.strategy_opps == "greedy"
        assert args.cost_calculation

        args.cost_calculation_method_deps = None
        args.cost_calculation_method_opps = None

        costs_vanilla = calculate_costs(cost_params, scenario, schedule, args)

        assert args.strategy_deps == "balanced"
        assert args.strategy_opps == "greedy"

        args.cost_calculation_method_deps = "fixed_wo_plw"
        args.cost_calculation_method_opps = "fixed_wo_plw"
        costs_with_same_strat = calculate_costs(cost_params, scenario, schedule, args)

        # assert all costs are the same
        for station in costs_vanilla.costs_per_gc:
            for key in costs_vanilla.costs_per_gc[station]:
                assert (costs_vanilla.costs_per_gc[station][key] ==
                        costs_with_same_strat.costs_per_gc[station][key]), station

        # test with different method (balanced_market costs for balanced strategy)
        args.cost_calculation_method_deps = "balanced_market"
        args.cost_calculation_method_opps = "balanced_market"
        costs_with_other_strat = calculate_costs(cost_params, scenario, schedule, args)
        station = "cumulated"
        cost_vanilla = costs_vanilla.costs_per_gc[station]
        cost_other = costs_with_other_strat.costs_per_gc[station]
        assert cost_vanilla["c_el_power_price_annual"] != cost_other["c_el_power_price_annual"], (
            "electric power costs should differ")

        # PLW: will create window time series before cost calculation
        args.cost_calculation_method_deps = "fixed_w_plw"
        args.cost_calculation_method_opps = "fixed_w_plw"
        costs_with_other_strat = calculate_costs(cost_params, scenario, schedule, args)
        """
        station = "cumulated"
        for key in costs_vanilla.costs_per_gc[station]:
            if "el_energy" not in key:
                continue
            assert (costs_vanilla.costs_per_gc[station][key] !=
                costs_with_other_strat.costs_per_gc[station][key]), key
        """


class TestStorageCorrectionUnit:
    """ Pure unit tests for Costs._storage_correction_eur, independent of any simulation. """

    def test_never_charged_returns_none(self):
        # SoC only ever decreases: no charging event, no price information available
        assert Costs._storage_correction_eur([100, 80, 60], [0.2, 0.2], 1.0) is None

    def test_zero_net_delta_returns_zero(self):
        # charges up, then fully discharges back to the starting level: no net drift
        assert Costs._storage_correction_eur([0, 500, 0], [0.2, None], 0.9) == 0

    def test_surplus_single_charge_efficiency_not_double_counted(self):
        # regression test for a bug where the surplus branch divided by efficiency twice:
        # charging 5000 kWh from empty at a constant price must cost exactly
        # (5000 / efficiency) * price, not (5000 / efficiency**2) * price
        efficiency = 0.95
        price = 0.20
        correction = Costs._storage_correction_eur([0, 5000], [price], efficiency)
        expected = -(5000 / efficiency) * price
        assert correction == pytest.approx(expected)
        # sanity check that the old (buggy) formula would have given a different result
        buggy = -(5000 / efficiency) * (price / efficiency)
        assert correction != pytest.approx(buggy)

    def test_deficit_uses_own_average_charging_price(self):
        # two charging events at different prices, ending below the starting level
        soc_kwh = [1000, 1300, 700, 900]
        prices = [0.10, None, 0.40]
        efficiency = 0.8
        correction = Costs._storage_correction_eur(soc_kwh, prices, efficiency)
        # own_avg_price = (300/eff*0.10 + 200/eff*0.40) / (300/eff + 200/eff) = 0.22
        own_avg_price = (300 / efficiency * 0.10 + 200 / efficiency * 0.40) / (
            300 / efficiency + 200 / efficiency)
        expected = (100 / efficiency) * own_avg_price
        assert correction == pytest.approx(expected)
        assert correction > 0  # deficit increases cost

    def test_surplus_moving_average_discounts_already_consumed_cheap_charge(self):
        # charge cheap, discharge most of it, then charge expensive: the remaining surplus
        # should be valued closer to the expensive price, not a blend of both
        soc_kwh = [0, 1000, 400, 1500]
        prices = [0.10, None, 0.30]
        efficiency = 0.9
        correction = Costs._storage_correction_eur(soc_kwh, prices, efficiency)
        assert correction == pytest.approx(-411.11111111111114)

        # a naive average over the whole charging history would have given a different
        # (too optimistic, less negative) result - the moving average must not match it
        own_avg_price = (1000 / efficiency * 0.10 + 1100 / efficiency * 0.30) / (
            1000 / efficiency + 1100 / efficiency)
        naive = -(1500 / efficiency) * own_avg_price
        assert correction != pytest.approx(naive)
        assert correction < naive  # moving average must value the surplus higher

    def test_leading_and_trailing_unknown_socs_are_ignored(self):
        # None at the edges (e.g. a vehicle not yet or no longer tracked) must be trimmed
        # before computing delta_e, not treated as part of the known range
        with_none = Costs._storage_correction_eur([None, 0, 5000, None], [None, 0.2], 0.95)
        without_none = Costs._storage_correction_eur([0, 5000], [0.2], 0.95)
        assert with_none == pytest.approx(without_none)

    def test_missing_price_at_a_charging_step_falls_back_to_own_average(self):
        # three charging steps, the middle one has no price information
        soc_kwh = [0, 100, 300, 400]
        prices = [0.20, None, 0.40]
        efficiency = 1.0
        correction = Costs._storage_correction_eur(soc_kwh, prices, efficiency)
        assert correction == pytest.approx(-120.0)

    def test_interior_gap_in_soc_kwh_raises(self):
        # a None strictly between the first and last known soc (not merely at the edges) would
        # silently be skipped by soc_kwh[t + 1] - soc_kwh[t] otherwise - must fail loudly instead
        soc_kwh = [0, 100, None, 300]
        prices = [0.20, 0.20, 0.20]
        with pytest.raises(AssertionError):
            Costs._storage_correction_eur(soc_kwh, prices, 1.0)


class TestStorageCorrectionHelpers:
    """ Unit tests for the attribution/price-lookup helpers, using lightweight stand-ins
    instead of a full simulation. """

    def test_last_rotation_gc_per_vehicle(self):
        early = datetime(2024, 1, 1, 8, 0)
        late = datetime(2024, 1, 1, 14, 0)
        rotations = {
            "r1": SimpleNamespace(vehicle_id="v1", departure_time=early, departure_name="A"),
            "r2": SimpleNamespace(vehicle_id="v1", departure_time=late, departure_name="B"),
            "r3": SimpleNamespace(vehicle_id="v2", departure_time=early, departure_name="Z"),
            "r4": SimpleNamespace(vehicle_id=None, departure_time=early, departure_name="A"),
        }
        costs = Costs.__new__(Costs)
        costs.schedule = SimpleNamespace(rotations=rotations)
        costs.gcs = {"A": None, "B": None}

        result = costs._last_rotation_gc_per_vehicle()

        # v1's later rotation (departing at B) wins over the earlier one (departing at A)
        assert result["v1"] == "B"
        # v2's only rotation departs at a non-electrified station
        assert result["v2"] == Costs.NOT_ELECTRIFIED
        # rotations without a vehicle_id are ignored
        assert "v3" not in result

    def test_vehicle_price_series_follows_connected_gc(self):
        costs = Costs.__new__(Costs)
        connected = [{"v1": "cs_A"}, {}, {"v1": "cs_B"}]
        costs.scenario = SimpleNamespace(
            connected=connected,
            n_intervals=len(connected),
            components=SimpleNamespace(charging_stations={
                "cs_A": SimpleNamespace(parent="A"),
                "cs_B": SimpleNamespace(parent="B"),
            })
        )
        price_series = {"A": [10, 20, 30], "B": [1, 2, 3]}

        def get_price_series(gc_id):
            return price_series[gc_id]

        result = costs._vehicle_price_series("v1", get_price_series)

        # connected to cs_A (gc A) at t=0, disconnected at t=1, connected to cs_B (gc B) at t=2
        assert result == [10, None, 3]

    def test_vehicle_price_series_batches_multiple_runs(self):
        # exercises several contiguous runs in one series: two timesteps at cs_A, a disconnect,
        # three timesteps at cs_B, then a run at cs_C that extends to the very last timestep
        # (regression test for an off-by-one that would index scenario.connected out of range
        # when a run reaches the end of the simulation)
        costs = Costs.__new__(Costs)
        connected = [
            {"v1": "cs_A"}, {"v1": "cs_A"}, {}, {"v1": "cs_B"}, {"v1": "cs_B"},
            {"v1": "cs_B"}, {"v1": "cs_C"},
        ]
        costs.scenario = SimpleNamespace(
            connected=connected,
            n_intervals=len(connected),
            components=SimpleNamespace(charging_stations={
                "cs_A": SimpleNamespace(parent="A"),
                "cs_B": SimpleNamespace(parent="B"),
                "cs_C": SimpleNamespace(parent="C"),
            })
        )
        price_series = {
            "A": [10, 11, 12, 13, 14, 15, 16],
            "B": [20, 21, 22, 23, 24, 25, 26],
            # "C" deliberately has no price series, to check the "gc_prices is None" guard
        }

        def get_price_series(gc_id):
            return price_series.get(gc_id)

        result = costs._vehicle_price_series("v1", get_price_series)

        assert result == [10, 11, None, 23, 24, 25, None]


class TestStorageCorrectionIntegration:
    """ Integration tests wiring set_storage_level_cost_correction into calculate_costs, using
    a real (schedule, scenario) pair. Exact correction amounts are already covered by
    TestStorageCorrectionUnit; these tests check that the flag/logging/CSV plumbing works. """

    @staticmethod
    def _basic_setup():
        schedule, scenario, args = BasicSchedule().basic_run()
        with open(args.cost_parameters_path, "r") as f:
            cost_params = uncomment_json_file(f)
        return schedule, scenario, args, cost_params

    def test_flag_off_leaves_correction_at_zero(self):
        schedule, scenario, args, cost_params = self._basic_setup()
        args.clean_costs_storage_level = False

        costs = calculate_costs(cost_params, scenario, schedule, args)

        for gc, values in costs.costs_per_gc.items():
            assert values["c_el_storage_level_correction_annual"] == 0, gc

    def test_flag_on_wires_correction_into_totals_and_csv(self):
        schedule, scenario, args, cost_params = self._basic_setup()

        args.clean_costs_storage_level = False
        costs_off = calculate_costs(cost_params, scenario, schedule, args)

        # The "balanced" strategy used by basic.cfg never actually cycles the stationary
        # battery at Station-0 (it stays at its initial SoC throughout the whole run), so a
        # synthetic surplus is injected directly into the recorded battery levels here to
        # deterministically exercise the correction, independent of charging-strategy behavior.
        bat_id = "Station-0 storage"
        n = scenario.n_intervals
        scenario.batteryLevels[bat_id] = [0.0] * (n - 1) + [123.0]

        args.clean_costs_storage_level = True
        costs_on = calculate_costs(cost_params, scenario, schedule, args)

        gc_id = "Station-0"
        correction = costs_on.costs_per_gc[gc_id]["c_el_storage_level_correction_annual"]
        assert correction != 0

        # wiring: corrected c_el_annual must equal the uncorrected value plus the correction
        assert costs_on.costs_per_gc[gc_id]["c_el_annual"] == pytest.approx(
            costs_off.costs_per_gc[gc_id]["c_el_annual"] + correction)

        # a gc without any storage drift stays untouched
        assert costs_on.costs_per_gc["Station-3"]["c_el_storage_level_correction_annual"] == 0

        # the correction is listed as its own row in the CSV output
        csv_rows = costs_on.to_csv_lists()
        keys = [row[0] for row in csv_rows]
        assert "c_el_storage_level_correction_annual" in keys

    def test_negative_vehicle_soc_is_clamped_to_zero(self, monkeypatch):
        schedule, scenario, args, cost_params = self._basic_setup()
        args.clean_costs_storage_level = True
        n = scenario.n_intervals
        v_id, vehicle = next(iter(scenario.components.vehicles.items()))
        capacity = vehicle.battery.capacity
        efficiency = vehicle.battery.efficiency

        # only this vehicle drifts: charges -100 kWh -> 0 at 0.20 EUR/kWh, then 0 -> 50 kWh at
        # 0.40 EUR/kWh (without clamping, the moving average would divide by zero)
        scenario.vehicle_socs = {
            v_id: [-100 / capacity] * (n - 2) + [0, 50 / capacity]}
        prices = [0.20] * (n - 2) + [0.40, 0.40]
        monkeypatch.setattr(Costs, "_vehicle_price_series", lambda self, v, f: prices)

        costs = calculate_costs(cost_params, scenario, schedule, args)

        # only the 50 kWh stored above 0 count as surplus, valued at the price they were bought
        annual_factor = (Costs.DAYS_PER_YEAR * 24 * 3600) / (
            n * scenario.interval.total_seconds())
        expected = -(50 / efficiency) * 0.40 * annual_factor
        total = sum(values["c_el_storage_level_correction_annual"]
                    for gc, values in costs.costs_per_gc.items() if gc != Costs.CUMULATED)
        assert total == pytest.approx(expected)

    def test_logging_reports_before_after_and_delta(self, caplog):
        schedule, scenario, args, cost_params = self._basic_setup()

        bat_id = "Station-0 storage"
        n = scenario.n_intervals
        scenario.batteryLevels[bat_id] = [0.0] * (n - 1) + [123.0]
        args.clean_costs_storage_level = True

        with caplog.at_level(logging.INFO):
            costs = calculate_costs(cost_params, scenario, schedule, args)

        correction = costs.costs_per_gc["Station-0"]["c_el_storage_level_correction_annual"]
        assert correction != 0
        assert any(
            "Station-0: electricity costs corrected for storage level drift" in r.message
            for r in caplog.records)
        assert any(
            "Cumulated electricity costs corrected for storage level drift" in r.message
            for r in caplog.records)


class TestAnnualKwhFromGrid:
    """ Regression tests for the fixed-tariff billed-energy calculation: a gc's flat price must
    not depend on how much it exports, since exported energy is never billed (only remunerated
    separately as feed-in). See :func:`Costs._annual_kwh_from_grid` for the bug this guards
    against (adding generation back onto totalLoad instead of clipping it). """

    def test_export_periods_are_excluded_not_subtracted(self):
        total_load = [100, 100, -500, 50, -1000]
        result = Costs._annual_kwh_from_grid(total_load, steps_per_hour=60, annual_factor=100)
        expected = (100 + 100 + 0 + 50 + 0) / 60 * 100
        assert result == pytest.approx(expected)

    def test_result_is_unaffected_by_the_size_of_export_periods(self):
        # a station exporting a little vs. a lot must be billed identically, as long as the
        # import (positive) periods are the same
        base = [100, 100, 50]
        with_small_export = base + [-1]
        with_large_export = base + [-100000]
        assert Costs._annual_kwh_from_grid(
            with_small_export, 60, 100) == pytest.approx(
            Costs._annual_kwh_from_grid(with_large_export, 60, 100))

    def test_matches_plain_sum_when_never_exporting(self):
        # the old (buggy) and new formula coincide exactly when totalLoad never goes negative,
        # which is why the bug only showed up once a gc started exporting
        total_load = [100, 50, 200, 0, 75]
        assert Costs._annual_kwh_from_grid(total_load, 60, 100) == pytest.approx(
            sum(total_load) / 60 * 100)


class TestPvAttributionUnit:
    """ Pure unit tests for the PV/grid merit-order split used for stationary batteries (see
    docs/source/simba_features.rst, "PV attribution for stationary batteries"), independent of
    any simulation. """

    def test_pv_surplus_series_converts_power_to_energy_and_subtracts_other_load(self):
        costs = Costs.__new__(Costs)
        costs.scenario = SimpleNamespace(
            n_intervals=3,
            stepsPerHour=60,
            localGenerationPower={"A": [120, 0, 300]},
            fixedLoads={"A": [{"x": 20}, {}, {"x": 50}]},
            connChargeByTS={"A": [{"cs1": 10}, {}, {"cs1": 100}]},
        )
        # other_load = [30, 0, 150]; raw surplus (kW) = [90, 0, 150]; /stepsPerHour (60) -> kWh
        assert costs._pv_surplus_series("A") == pytest.approx([1.5, 0.0, 2.5])

    def test_pv_surplus_series_is_zero_without_any_pv(self):
        costs = Costs.__new__(Costs)
        costs.scenario = SimpleNamespace(
            n_intervals=2, stepsPerHour=60, localGenerationPower={"A": [0, 0]},
            fixedLoads={}, connChargeByTS={})
        assert costs._pv_surplus_series("A") == [0, 0]

    def test_feed_in_remuneration_price_picks_matching_kwp_tier(self):
        costs = Costs.__new__(Costs)
        costs.scenario = SimpleNamespace(components=SimpleNamespace(
            photovoltaics={"pv1": SimpleNamespace(parent="A", nominal_power=50)}))
        costs.gcs = {"A": SimpleNamespace(grid_operator="sh-netz")}
        costs.params = {"sh-netz": {"feed-in_remuneration": {
            "PV": {"kWp": [10, 100, 1000], "remuneration": [8, 5, 3]}}}}

        # 50 kWp falls into the second tier (<= 100 kWp) -> 5 ct/kWh -> 0.05 EUR/kWh
        assert costs._feed_in_remuneration_price("A") == pytest.approx(0.05)

    def test_feed_in_remuneration_price_is_zero_without_pv_plant(self):
        costs = Costs.__new__(Costs)
        costs.scenario = SimpleNamespace(components=SimpleNamespace(photovoltaics={}))
        costs.gcs = {"A": SimpleNamespace(grid_operator="sh-netz")}
        costs.params = {}
        assert costs._feed_in_remuneration_price("A") == 0.0

    def test_feed_in_remuneration_price_warns_and_is_zero_above_highest_tier(self, caplog):
        costs = Costs.__new__(Costs)
        costs.scenario = SimpleNamespace(components=SimpleNamespace(
            photovoltaics={"pv1": SimpleNamespace(parent="A", nominal_power=5000)}))
        costs.gcs = {"A": SimpleNamespace(grid_operator="sh-netz")}
        costs.params = {"sh-netz": {"feed-in_remuneration": {
            "PV": {"kWp": [10, 100, 1000], "remuneration": [8, 5, 3]}}}}
        with caplog.at_level(logging.WARNING):
            assert costs._feed_in_remuneration_price("A") == 0.0
        assert any("no feed-in remuneration defined" in r.message for r in caplog.records)

    def test_blend_battery_prices_with_pv_splits_charging_by_merit_order(self):
        # charge 0 -> 50 (step 0), idle (step 1), charge 50 -> 150 (step 2)
        soc_kwh = [0, 50, 50, 150]
        grid_prices = [0.30, 0.30, 0.30, 0.30]
        pv_surplus = [30, 0, 40]  # kWh available at each step
        price_pv = 0.05
        efficiency = 1.0

        blended = Costs._blend_battery_prices_with_pv(
            soc_kwh, grid_prices, pv_surplus, price_pv, efficiency)

        # step 0: 30 kWh of the 50 kWh charge from PV, 20 from grid
        #   -> (20*0.30 + 30*0.05) / 50 = 0.15
        # step 1: no charging, price is irrelevant and left at the grid price
        # step 2: 40 kWh of the 100 kWh charge from PV, 60 from grid
        #   -> (60*0.30 + 40*0.05) / 100 = 0.20
        assert blended == pytest.approx([0.15, 0.30, 0.20, 0.30])

    def test_blend_battery_prices_with_pv_is_a_no_op_without_pv_surplus(self):
        soc_kwh = [0, 50, 150]
        grid_prices = [0.30, 0.30, 0.30]
        pv_surplus = [0, 0]
        blended = Costs._blend_battery_prices_with_pv(soc_kwh, grid_prices, pv_surplus, 0.05, 1.0)
        assert blended == grid_prices


class TestPvAttributionIntegration:
    """ Integration test for the PV/grid merit-order split wired into
    set_storage_level_cost_correction, using a real (schedule, scenario) pair. """

    @staticmethod
    def _basic_setup():
        return TestStorageCorrectionIntegration._basic_setup()

    def test_pv_surplus_is_wired_into_the_battery_correction(self, monkeypatch):
        # Station-0's gc also receives vehicle corrections unrelated to its stationary battery
        # (whichever vehicles' last rotation departs from there), so the per-gc correction
        # cannot be compared directly against a hand-computed battery-only expectation. Instead,
        # this runs the real pipeline twice, changing only the mocked PV surplus available to
        # the battery between the two runs - the (identical in both runs) vehicle contribution
        # then cancels out, isolating the effect of the merit-order split on the battery.
        schedule, scenario, args, cost_params = self._basic_setup()
        bat_id = "Station-0 storage"
        gc_id = "Station-0"
        n = scenario.n_intervals
        # battery charges from 0 to 123 kWh in the last step, i.e. during the transition at
        # index n - 2 (soc_kwh[n-2] -> soc_kwh[n-1])
        soc_kwh = [0.0] * (n - 1) + [123.0]
        args.clean_costs_storage_level = True
        monkeypatch.setattr(Costs, "_feed_in_remuneration_price", lambda self, gc: 0.05)

        scenario.batteryLevels[bat_id] = list(soc_kwh)
        monkeypatch.setattr(Costs, "_pv_surplus_series", lambda self, gc: [0.0] * n)
        costs_no_pv = calculate_costs(cost_params, scenario, schedule, args)
        correction_no_pv = costs_no_pv.costs_per_gc[gc_id]["c_el_storage_level_correction_annual"]

        surplus_with_pv = [0.0] * n
        surplus_with_pv[n - 2] = 1000.0  # abundant PV exactly when the battery charges
        scenario.batteryLevels[bat_id] = list(soc_kwh)
        monkeypatch.setattr(Costs, "_pv_surplus_series", lambda self, gc: surplus_with_pv)
        costs_with_pv = calculate_costs(cost_params, scenario, schedule, args)
        correction_with_pv = (
            costs_with_pv.costs_per_gc[gc_id]["c_el_storage_level_correction_annual"])

        # covering the battery's charging with (cheaper) PV instead of grid energy must shrink
        # the cost reduction the surplus otherwise provides
        assert correction_with_pv > correction_no_pv
