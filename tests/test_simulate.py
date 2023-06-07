from argparse import Namespace
from pathlib import Path
import pytest
import warnings

from ebus_toolbox.simulate import simulate


root_path = Path(__file__).parent.parent
example_path = root_path / "data/examples"


class TestSimulate:

    DEFAULT_VALUES = {
        "vehicle_types": example_path / "vehicle_types.json",
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
        "seed": None,
        "default_buffer_time_opps": 0,
        "desired_soc_opps": 1,
        "desired_soc_deps": 1,
        "min_charging_time": 0,
        "default_voltage_level": "MV",
    }

    def test_basic(self):
        args = Namespace(**(self.DEFAULT_VALUES))
        simulate(args)

    def test_missing(self):
        # every value in DEFAULT_VALUES is expected to be set, so omitting one should raise an error
        values = self.DEFAULT_VALUES.copy()
        for k, v in self.DEFAULT_VALUES.items():
            del values[k]
            with pytest.raises(Exception):
                simulate(Namespace(**values))
            # reset
            values[k] = v

        # required file missing
        for file_type in ["vehicle_types", "electrified_stations", "cost_parameters_file"]:
            values[file_type] = ""
            with pytest.raises(Exception):
                simulate(Namespace(**values))
            # reset
            values[file_type] = self.DEFAULT_VALUES[file_type]

    def test_unknown_mode(self):
        # try to run a mode that does not exist
        args = Namespace(**(self.DEFAULT_VALUES))
        args.mode = "foo"
        with pytest.warns(UserWarning, match="Unknown mode"):
            simulate(args)

    def test_late_sim(self):
        # sim mode has no function, just produces a warning later
        args = Namespace(**(self.DEFAULT_VALUES))
        args.mode = ["sim", "sim"]
        with pytest.warns(UserWarning, match="Intermediate sim"):
            simulate(args)

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
