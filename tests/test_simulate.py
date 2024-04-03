from argparse import Namespace
import logging
from pathlib import Path
import pytest
import warnings

from simba.simulate import simulate


root_path = Path(__file__).parent.parent
example_path = root_path / "data/examples"


class TestSimulate:
    # Add propagate_mode_errors as developer setting to raise Exceptions.
    DEFAULT_VALUES = {
        "vehicle_types_path": example_path / "vehicle_types.json",
        "electrified_stations": example_path / "electrified_stations.json",
        "cost_parameters_file": example_path / "cost_params.json",
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
        "desired_soc_opps": 1,
        "desired_soc_deps": 1,
        "min_charging_time": 0,
        "default_voltage_level": "MV",
        "propagate_mode_errors": True,
    }

    def test_basic(self):
        args = Namespace(**(self.DEFAULT_VALUES))
        simulate(args)

    def test_missing(self):
        # every value in DEFAULT_VALUES is expected to be set, so omitting one should raise an error
        values = self.DEFAULT_VALUES.copy()
        # except propagate_modes_error, since this is only checked when error needs propagation
        del self.DEFAULT_VALUES["propagate_mode_errors"]
        for k, v in self.DEFAULT_VALUES.items():
            del values[k]
            with pytest.raises(Exception):
                simulate(Namespace(**values))
            # reset
            values[k] = v
        # restore the setting for further testing
        self.DEFAULT_VALUES["propagate_mode_errors"] = values["propagate_mode_errors"]

        # required file missing
        for file_type in ["vehicle_types_path", "electrified_stations", "cost_parameters_file"]:
            values[file_type] = ""
            with pytest.raises(Exception):
                simulate(Namespace(**values))
            # reset
            values[file_type] = self.DEFAULT_VALUES[file_type]

    def test_unknown_mode(self, caplog):
        # try to run a mode that does not exist
        args = Namespace(**(self.DEFAULT_VALUES))
        args.mode = "foo"
        with caplog.at_level(logging.ERROR):
            simulate(args)
        assert caplog.record_tuples == [('root', logging.ERROR, 'Unknown mode foo ignored')]

    def test_late_sim(self, caplog):
        # sim mode has no function, just produces a log info later
        args = Namespace(**(self.DEFAULT_VALUES))
        args.mode = ["sim", "sim"]
        with caplog.at_level(logging.INFO):
            simulate(args)
        # also captures INFO about running SpiceEV, so only compare second element
        assert caplog.record_tuples[1] == ('root', logging.INFO, 'Intermediate sim ignored')

    def test_mode_service_opt(self):
        # basic run
        values = self.DEFAULT_VALUES.copy()
        values["mode"] = "service_optimization"
        simulate(Namespace(**values))
        # all rotations remain negative
        values["desired_soc_deps"] = 0
        values["desired_soc_opps"] = 0
        values["ALLOW_NEGATIVE_SOC"] = True
        simulate(Namespace(**values))

    def test_mode_change_charge_type(self):
        # all rotations remain negative
        values = self.DEFAULT_VALUES.copy()
        values["mode"] = "neg_oppb_to_depb"
        values["desired_soc_deps"] = 0
        values["desired_soc_opps"] = 0
        values["ALLOW_NEGATIVE_SOC"] = True
        simulate(Namespace(**values))

    def test_mode_remove_negative(self):
        values = self.DEFAULT_VALUES.copy()
        values["mode"] = "remove_negative"
        values["desired_soc_deps"] = 0
        # values["desired_soc_opps"] = 0
        values["ALLOW_NEGATIVE_SOC"] = True
        simulate(Namespace(**values))

    def test_mode_report(self, tmp_path):
        # report with cost calculation, write to tmp
        values = self.DEFAULT_VALUES.copy()
        values["mode"] = "report"
        values["cost_calculation"] = True
        values["output_directory"] = tmp_path
        values["strategy"] = "distributed"
        values["strategy_deps"] = "balanced"
        values["strategy_opps"] = "greedy"

        values["show_plots"] = False
        # tuned so that some rotations don't complete
        values["days"] = .33
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(Namespace(**values))

    def test_empty_report(self, tmp_path):
        # report with no rotations
        values = self.DEFAULT_VALUES.copy()
        values.update({
            "mode": ["remove_negative", "report"],
            "desired_soc_deps": 0,
            "ALLOW_NEGATIVE_SOC": True,
            "cost_calculation": True,
            "output_directory": tmp_path,
            "strategy": "distributed",
            "show_plots": False,
        })
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(Namespace(**values))
