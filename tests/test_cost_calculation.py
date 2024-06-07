from tests.test_schedule import BasicSchedule
from simba.util import uncomment_json_file
from simba.costs import calculate_costs

class TestCostCalculation:
    def test_cost_calculation(self):
        schedule, scenario, args = BasicSchedule().basic_run()
        file = args.cost_parameters_file
        with open(file, "r") as file:
            cost_params = uncomment_json_file(file)

        costs = calculate_costs(cost_params, scenario, schedule, args)

        assert args.strategy_deps == "balanced"
        assert args.strategy_opps == "greedy"

        args.cost_calculation_strategy_opps == "balanced"
