from argparse import Namespace
import logging
from pathlib import Path
import pytest
import warnings

from simba import util
from simba.simulate import simulate


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
        "input_schedule": example_path / "trips_example.csv",
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
        "input_schedule": example_path / "trips_example.csv",
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
        args.output_directory = tmp_path
        args.strategy_deps = "balanced"
        args.strategy_opps = "greedy"
        args.show_plots = False

        # tuned so that some rotations don't complete
        args.days = .33
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)

    def test_empty_report(self, tmp_path):
        # report with no rotations
        args = self.get_args()
        args.mode = ["remove_negative", "report"]
        args.desired_soc_deps = 0
        args.ALLOW_NEGATIVE_SOC = True
        args.cost_calculation = True
        args.output_directory = tmp_path
        args.show_plots = False
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)

    def test_extended_plot(self, tmp_path):
        args = self.get_args()
        args.mode = "report"
        args.output_directory = tmp_path
        args.show_plots = False
        args.extended_output_plots = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)

    def test_create_trips_in_report(self, tmp_path):
        # create_trips_in_report option: must generate valid input trips.csv
        args_dict = vars(self.get_args())
        update_dict = {
            "mode": ["report"],
            "desired_soc_deps": 0,
            "ALLOW_NEGATIVE_SOC": True,
            "cost_calculation": False,
            "output_directory": tmp_path,
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
        args_dict["input_schedule"] = tmp_path / "report_1/trips.csv"
        simulate(Namespace(**(args_dict)))

    def test_mode_recombination(self):
        args = self.get_args()
        args.mode = "recombination"
        simulate(args)
