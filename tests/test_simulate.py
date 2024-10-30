from argparse import Namespace
from csv import DictReader
import json
import logging
from pathlib import Path
import pytest
import warnings

from simba import util
from simba.simulate import simulate, modes_simulation

root_path = Path(__file__).parent.parent
example_path = root_path / "data/examples"


class TestSimulate:
    # Add propagate_mode_errors as developer setting to raise Exceptions.
    NON_DEFAULT_VALUES = {
        "vehicle_types_path": example_path / "vehicle_types.json",
        "electrified_stations_path": example_path / "electrified_stations.json",
        "station_data_path": example_path / "all_stations.csv",
        "cost_parameters_path": example_path / "cost_params.json",
        "outside_temperature_over_day_path": example_path / "default_temp_summer.csv",
        "level_of_loading_over_day_path": example_path / "default_level_of_loading_over_day.csv",
        "schedule_path": example_path / "trips_example.csv",
        "mode": [],
        "interval": 15,
        "propagate_mode_errors": True,
        "preferred_charging_type": "oppb"
    }

    MANDATORY_ARGS = {
        "vehicle_types_path": example_path / "vehicle_types.json",
        "electrified_stations_path": example_path / "electrified_stations.json",
        "cost_parameters_path": example_path / "cost_params.json",
        "outside_temperature_over_day_path": example_path / "default_temp_summer.csv",
        "level_of_loading_over_day_path": example_path / "default_level_of_loading_over_day.csv",
        "schedule_path": example_path / "trips_example.csv",
        "mode": [],
        "min_recharge_deps_oppb": 1,
        "min_recharge_deps_depb": 1,
        "gc_power_opps": 100000,
        "gc_power_deps": 100000,
        "cs_power_opps": 300,
        "cs_power_deps_depb": 150,
        "cs_power_deps_oppb": 150,
        "interval": 15,
        "days": None,
        "signal_time_dif": 10,
        "include_price_csv": None,
        "rotation_filter_variable": None,
        "seed": 1,
        "default_buffer_time_opps": 0,
        "default_buffer_time_deps": 0,
        "desired_soc_opps": 1,
        "desired_soc_deps": 1,
        "min_charging_time": 0,
        "default_voltage_level": "MV",
    }

    def get_args(self):
        # try to run a mode that does not exist
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        parser = util.get_parser()
        # Set the parser defaults to the specified non default values
        parser.set_defaults(**self.NON_DEFAULT_VALUES)
        # get all args with default values
        args, _ = parser.parse_known_args()
        args = util.replace_deprecated_arguments(args)

        return args

    def test_basic(self):
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        simulate(args)

    def test_mandatory_missing(self):
        values = self.MANDATORY_ARGS.copy()

        for k, v in self.MANDATORY_ARGS.items():
            del values[k]
            with pytest.raises(Exception):
                simulate(Namespace(**values))
            # reset
            values[k] = v

        # required file missing
        for fpath in ["vehicle_types_path", "electrified_stations_path", "cost_parameters_path"]:
            values[fpath] = ""
            with pytest.raises(Exception):
                simulate(Namespace(**values))
            # reset
            values[fpath] = self.MANDATORY_ARGS[fpath]

    def test_unknown_mode(self, caplog):
        # try to run a mode that does not exist
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        args.mode = "foo"
        with caplog.at_level(logging.ERROR):
            simulate(args)
        assert caplog.record_tuples == [('root', logging.ERROR, 'Unknown mode foo ignored')]

    def test_late_sim(self, caplog):
        # sim mode has no function, just produces a log info later
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        args.mode = ["sim", "sim"]
        with caplog.at_level(logging.INFO):
            simulate(args)
        # also captures INFO about running SpiceEV, so only compare second element
        assert caplog.record_tuples[1] == ('root', logging.INFO, 'Intermediate sim ignored')

    def test_error_handling_in_modes(self, caplog):
        args = self.get_args()
        args.mode = "station_optimization"
        args.propagate_mode_errors = False
        del args.optimizer_config_path
        with caplog.at_level(logging.ERROR):
            simulate(args)
        assert len(caplog.record_tuples) > 0

    def test_modes(self, caplog, tmp_path):
        args = self.get_args()
        args.output_path = tmp_path
        schedule, scenario = simulate(args)
        with caplog.at_level(logging.ERROR):
            args.mode = "sim"
            modes_simulation(schedule, scenario, args)

            args.mode = "sim_greedy"
            modes_simulation(schedule, scenario, args)

            args.mode = "neg_depb_to_oppb"
            modes_simulation(schedule, scenario, args)

            args.mode = "split_negative_depb"
            modes_simulation(schedule, scenario, args)

            args.mode = "station_optimization_single_step"
            modes_simulation(schedule, scenario, args)

            args.mode = "station_optimization"
            modes_simulation(schedule, scenario, args)

        assert len(caplog.record_tuples) == 0

    def test_mode_service_opt(self):
        # basic run
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        args.mode = "service_optimization"
        simulate(args)
        # all rotations remain negative
        args.desired_soc_deps = 0
        args.desired_soc_opps = 0
        args.ALLOW_NEGATIVE_SOC = True
        simulate(args)

    def test_mode_change_charge_type(self):
        # all rotations remain negative
        args = self.get_args()
        args.mode = "neg_oppb_to_depb"
        args.desired_soc_deps = 0
        args.desired_soc_opps = 0
        args.ALLOW_NEGATIVE_SOC = True
        simulate(args)

    def test_mode_remove_negative(self):
        args = self.get_args()
        args.mode = "remove_negative"
        args.desired_soc_deps = 0
        args.ALLOW_NEGATIVE_SOC = True
        simulate(args)

    def test_mode_report(self, tmp_path):
        # report with cost calculation, write to tmp
        args = self.get_args()
        args.mode = "report"
        args.cost_calculation = True
        args.output_path = tmp_path
        args.strategy_deps = "balanced"
        args.strategy_opps = "greedy"
        args.show_plots = False

        # tuned so that some rotations don't complete
        args.days = .33
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)

        # get power procurement price from price sheet
        with open(args.cost_parameters_path) as file:
            cost_params = json.load(file)
            procurement_price = cost_params["default_grid_operator"]["power_procurement"]["charge"]
            # given in ct/kWh, convert to €/kWh with sign change (direction of grid power)
            procurement_price = -procurement_price / 100
        # read out vehicle costs, procurement price must match
        with (tmp_path / "report_1/summary_vehicles_costs.csv").open() as csvfile:
            reader = DictReader(csvfile)
            for row in reader:
                if row["parameter"] == "unit":
                    continue
                energy = float(row["annual_kWh_from_grid"])
                price = float(row["c_el_procurement_annual"])
                if energy == 0:
                    assert price == 0
                else:
                    assert pytest.approx(price / energy) == procurement_price

    def test_empty_report(self, tmp_path):
        # report with no rotations
        args = self.get_args()
        args.mode = ["remove_negative", "report"]
        args.desired_soc_deps = 0
        args.ALLOW_NEGATIVE_SOC = True
        args.cost_calculation = True
        args.output_path = tmp_path
        args.show_plots = False
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)

    def test_extended_plot(self, tmp_path):
        args = self.get_args()
        args.mode = "report"
        args.output_path = tmp_path
        args.show_plots = False
        args.extended_output_plots = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)
            for expected_file in [
                    "active_rotations.png", "charge_types.png", "Station-0_power_overview.png",
                    "distribution_consumption.png", "distribution_distance.png"
            ]:
                assert (tmp_path / f"report_1/extended_plots/{expected_file}").exists(), (
                    f"{expected_file} not found in extended plots")

    def test_create_trips_in_report(self, tmp_path):
        # create_trips_in_report option: must generate valid input trips.csv
        args_dict = vars(self.get_args())
        update_dict = {
            "mode": ["report"],
            "desired_soc_deps": 0,
            "ALLOW_NEGATIVE_SOC": True,
            "cost_calculation": False,
            "output_path": tmp_path,
            "show_plots": False,
            "create_trips_in_report": True,
        }
        args_dict.update(update_dict)

        # simulate base scenario, report generates new trips.csv in (tmp) output
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(Namespace(**args_dict))
        # new simulation with generated trips.csv
        args_dict = vars(self.get_args())
        args_dict["schedule_path"] = tmp_path / "report_1/trips.csv"
        simulate(Namespace(**(args_dict)))

    def test_mode_recombination(self):
        args = self.get_args()
        args.mode = "recombination"
        simulate(args)
