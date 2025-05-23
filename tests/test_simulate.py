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
        "vehicle_types_path": example_path / "vehicle_types/vehicle_types.json",
        "electrified_stations_path":
            example_path / "electrified_stations/electrified_stations.json",
        "station_data_path": example_path / "consumption/all_stations.csv",
        "cost_parameters_path": example_path / "costs/cost_params.json",
        "outside_temperature_over_day_path": example_path / "consumption/default_temp_summer.csv",
        "level_of_loading_over_day_path":
            example_path / "consumption/default_level_of_loading_over_day.csv",
        "schedule_path": example_path / "trips/trips_example.csv",
        "mode": [],
        "interval": 15,
        "propagate_mode_errors": True,
        "preferred_charging_type": "oppb"
    }

    def get_args(self):
        # try to run a mode that does not exist
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        parser = util.get_parser()
        # Set the parser defaults to the specified non default values
        parser.set_defaults(**self.NON_DEFAULT_VALUES)
        # get all args with default values
        args, _ = parser.parse_known_args()
        util.mutate_args_for_spiceev(args)
        args = util.replace_deprecated_arguments(args)

        return args

    def test_basic(self):
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        simulate(args)

    def test_mandatory_missing(self):
        mandatory_args = {
            "vehicle_types_path": example_path / "vehicle_types/vehicle_types.json",
            "electrified_stations_path":
                example_path / "electrified_stations/electrified_stations.json",
            "cost_parameters_path": example_path / "costs/cost_params.json",
            "outside_temperature_over_day_path":
                example_path / "consumption/default_temp_summer.csv",
            "level_of_loading_over_day_path":
                example_path / "consumption/default_level_of_loading_over_day.csv",
            "schedule_path": example_path / "trips/trips_example.csv",
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

        values = mandatory_args.copy()
        for k, v in mandatory_args.items():
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
            values[fpath] = mandatory_args[fpath]

    def test_unknown_mode(self, caplog):
        # try to run a mode that does not exist
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        args.mode = ["foo"]
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
        # also captures INFO about running SpiceEV, so only compare last element
        assert caplog.record_tuples[-1] == ('root', logging.INFO, 'Intermediate sim ignored')

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
        args.propagate_mode_errors = True
        schedule, scenario = simulate(args)
        with caplog.at_level(logging.ERROR):
            args.mode = ["sim"]
            modes_simulation(schedule, scenario, args)

            args.mode = ["sim_greedy"]
            modes_simulation(schedule, scenario, args)

            args.mode = ["neg_depb_to_oppb"]
            modes_simulation(schedule, scenario, args)

            args.mode = ["split_negative_depb"]
            modes_simulation(schedule, scenario, args)

            args.mode = ["station_optimization_single_step"]
            modes_simulation(schedule, scenario, args)

            args.mode = ["station_optimization"]
            modes_simulation(schedule, scenario, args)

        # warning: no column for Station-3 price CSV
        assert len(caplog.record_tuples) == 1

    def test_mode_service_opt(self):
        # basic run
        # Get the parser from util. This way the test is directly coupled to the parser arguments
        args = self.get_args()
        args.mode = ["service_optimization"]
        simulate(args)
        # all rotations remain negative
        args.desired_soc_deps = 0
        args.desired_soc_opps = 0
        simulate(args)

    def test_mode_change_charge_type(self):
        # all rotations remain negative
        args = self.get_args()
        args.mode = ["neg_oppb_to_depb"]
        args.desired_soc_deps = 0
        args.desired_soc_opps = 0
        simulate(args)

    def test_mode_remove_negative(self):
        args = self.get_args()
        args.mode = ["remove_negative"]
        args.desired_soc_deps = 0
        simulate(args)

    def test_mode_report(self, tmp_path):
        # report with cost calculation, write to tmp
        args = self.get_args()
        args.mode = ["report"]
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
            # save procurement costs and annual energy for each Station
            station_data = {s: {"energy_costs": None, "energy_amount": None}
                            for s in reader.fieldnames if s.startswith("Station")}
            for row in reader:
                if row["parameter"] == "c_el_procurement_annual":
                    for name in station_data:
                        station_data[name]["energy_costs"] = float(row[name])
                if row["parameter"] == "annual_kWh_from_grid":
                    for name in station_data:
                        station_data[name]["energy_amount"] = float(row[name])
            # check quotient: no energy must mean no cost, otherwise procurement price
            for name in station_data:
                energy = station_data[name]["energy_amount"]
                cost = station_data[name]["energy_costs"]
                if energy == 0:
                    assert cost == 0, f"{name} has costs without energy"
                else:
                    assert pytest.approx(cost/energy) == procurement_price, f"{name}: procurement"

    def test_empty_report(self, tmp_path):
        # report with no rotations
        args = self.get_args()
        args.mode = ["remove_negative", "report"]
        args.desired_soc_deps = 0
        args.cost_calculation = True
        args.output_path = tmp_path
        args.show_plots = False
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)

    def test_extended_plot(self, tmp_path):
        args = self.get_args()
        args.mode = ["report"]
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
        args = self.get_args()
        args.mode = ["report"]
        args.desired_soc_deps = 0
        args.cost_calculation = False
        args.output_path = tmp_path
        args.show_plots = False
        args.create_trips_in_report = True

        # simulate base scenario, report generates new trips.csv in (tmp) output
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            simulate(args)
        # new simulation with generated trips.csv
        args = self.get_args()
        args.input_schedule = tmp_path / "report_1/trips.csv"
        simulate(args)

    def test_pickle(self, tmp_path):
        # create pickle in report
        args = self.get_args()
        args.mode = ["report"]
        args.show_plots = False
        args.cost_calculation = False
        args.output_path = tmp_path
        args.create_pickle_in_report = True
        simulate(args)
        pickle_path = tmp_path / "report_1/scenario.pkl"
        assert pickle_path.exists()

        # read in pickle for new simulation. load_pickle can't be only mode, add dummy sim
        args.mode = ["load_pickle", "sim"]
        args.output_path = None
        args.load_pickle_path = pickle_path
        args.cost_parameters_path = None  # keep original cost parameters
        schedule, _ = simulate(args)
        assert "foo" not in schedule.data_container.cost_parameters_data

        # replace cost parameters after loading pickle
        args.mode = ["load_pickle", "sim"]
        args.cost_parameters_path = tmp_path / "cost_params.json"
        with open(args.cost_parameters_path, "w") as f:
            f.write('{"foo": 1}')
        schedule, _ = simulate(args)
        assert schedule.data_container.cost_parameters_data["foo"] == 1

    def test_mode_recombination(self):
        args = self.get_args()
        args.mode = ["recombination"]
        simulate(args)
