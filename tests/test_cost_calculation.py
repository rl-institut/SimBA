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
        print(costs_vanilla.costs_per_gc["cumulated"]["c_total_annual"])
        print(costs_with_other_strat.costs_per_gc["cumulated"]["c_total_annual"])
        station = "cumulated"
        for key in costs_vanilla.costs_per_gc[station]:
            if "el_energy" not in key:
                continue
            assert (costs_vanilla.costs_per_gc[station][key] !=
                    costs_with_other_strat.costs_per_gc[station][key]), key

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
