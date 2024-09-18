from tests.test_schedule import BasicSchedule
from simba.util import uncomment_json_file
from simba.costs import calculate_costs


class TestCostCalculation:
    def test_cost_calculation(self):
        schedule, scenario, args = BasicSchedule().basic_run()
        file = args.cost_parameters_path
        with open(file, "r") as file:
            cost_params = uncomment_json_file(file)

        assert args.strategy_deps == "balanced"
        assert args.strategy_opps == "greedy"

        args.cost_calculation_strategy_deps = None
        args.cost_calculation_strategy_opps = None

        costs_vanilla = calculate_costs(cost_params, scenario, schedule, args)

        assert args.strategy_deps == "balanced"
        assert args.strategy_opps == "greedy"

        args.cost_calculation_strategy_deps = "balanced"
        args.cost_calculation_strategy_opps = "greedy"
        costs_with_same_strat = calculate_costs(cost_params, scenario, schedule, args)

        # assert all costs are the same
        for station in costs_vanilla.costs_per_gc:
            for key in costs_vanilla.costs_per_gc[station]:
                assert (costs_vanilla.costs_per_gc[station][key] ==
                        costs_with_same_strat.costs_per_gc[station][key]), station

        args.cost_calculation_strategy_opps = "balanced_market"
        args.cost_calculation_strategy_deps = "balanced_market"
        costs_with_other_strat = calculate_costs(cost_params, scenario, schedule, args)
        print(costs_vanilla.costs_per_gc["cumulated"]["c_total_annual"])
        print(costs_with_other_strat.costs_per_gc["cumulated"]["c_total_annual"])
        station = "cumulated"
        for key in costs_vanilla.costs_per_gc[station]:
            if "el_energy" not in key:
                continue
            assert (costs_vanilla.costs_per_gc[station][key] !=
                    costs_with_other_strat.costs_per_gc[station][key]), key

        args.cost_calculation_strategy_opps = "peak_load_window"
        args.cost_calculation_strategy_deps = "peak_load_window"
        costs_with_other_strat = calculate_costs(cost_params, scenario, schedule, args)
        station = "cumulated"
        for key in costs_vanilla.costs_per_gc[station]:
            if "el_energy" not in key:
                continue
            assert (costs_vanilla.costs_per_gc[station][key] !=
                    costs_with_other_strat.costs_per_gc[station][key]), key

        args.cost_calculation_strategy_opps = "peak_shaving"
        args.cost_calculation_strategy_deps = "peak_shaving"
        costs_with_other_strat = calculate_costs(cost_params, scenario, schedule, args)
        # assert all costs are the same
        for station in costs_vanilla.costs_per_gc:
            for key in costs_vanilla.costs_per_gc[station]:
                assert (costs_vanilla.costs_per_gc[station][key] ==
                        costs_with_other_strat.costs_per_gc[station][key]), station
