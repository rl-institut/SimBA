import pathlib
import pytest
import sys

import ebus_toolbox.optimizer_util as opt_util
import ebus_toolbox.util as util
from ebus_toolbox.station_optimization import run_optimization
from ebus_toolbox.trip import Trip
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from spice_ev.report import generate_soc_timeseries

file_root = pathlib.Path(__file__).parent / "test_input_files/optimization"


class TestStationOptimization:

    # Using the pytest fixture tmp_path to have access to a temporary directory as class attribute
    @pytest.fixture(autouse=True)
    def initialize_tmp_path(self, tmp_path):
        self.tmp_path = tmp_path

    mandatory_args = {
        "min_recharge_deps_oppb": 0,
        "min_recharge_deps_depb": 0,
        "gc_power_opps": 1000,
        "gc_power_deps": 1000,
        "cs_power_opps": 100,
        "cs_power_deps_depb": 50,
        "cs_power_deps_oppb": 150
    }
    with open(file_root / "vehicle_types.json", "r", encoding='utf-8') as file:
        vehicle_types = util.uncomment_json_file(file)

    with open(file_root / "electrified_stations.json", "r", encoding='utf-8') as file:
        electrified_stations = util.uncomment_json_file(file)

    def test_basic_run(self, trips_file_name="trips.csv"):
        """ Check if running a basic example works and if a scenario object is returned.

        :param trips_file_name: file name of the trips file. Has to be inside the test_input_file
            folder
        :type trips_file_name: str
        :return: schedule, scenario"""
        path_to_trips = file_root / trips_file_name
        sys.argv = ["foo", "--config", str(file_root / "ebus_toolbox.cfg")]
        args = util.get_args()

        Trip.consumption = Consumption(self.vehicle_types,
                                       outside_temperatures=None,
                                       level_of_loading_over_day=None)

        path_to_all_station_data = file_root / "all_stations.csv"
        generated_schedule = Schedule.from_csv(path_to_trips, self.vehicle_types,
                                               self.electrified_stations,
                                               **self.mandatory_args,
                                               station_data_path=path_to_all_station_data)

        scen = generated_schedule.run(args)
        # optimization depends on vehicle_socs, therefore they need to be generated
        generate_soc_timeseries(scen)
        args.output_directory = self.tmp_path
        args.results_directory = self.tmp_path
        assert self.tmp_path
        return generated_schedule, scen, args

    def test_basic_optimization(self):
        """ Test if the base optimization finishes without raising errors"""
        trips_file_name = "trips_for_optimizer.csv"
        sched, scen, args = self.test_basic_run(trips_file_name)
        config_path = file_root / "optimizer.cfg"
        conf = opt_util.read_config(config_path)

        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
        conf.solver = "spiceev"
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)

    def test_deep_optimization(self):
        """ Check if deep analysis finds the prepared optimal solution for the test case.

        The Test case is a a 3 star like network with two rotations which are negative without
        electrification.
        Rotation 1:
        Depot --> Station-1 ---> Station-2 -->Station-1-->Depot
        Rotation 2:
        Depot --> Station-1 ---> Station-3 -->Station-1-->Depot

        Greedy optimization will electrify Station-1 first since it is helpful for both rotations.
        Since electrifying Station-1 is not enough to electrify either rotations, greedy
        optimization will electrify Station-2 and Station-3 as well. Deep Analysis will expand
        the checked nodes --> It should find that electrifying Station-2 and Station-3 is enough
        without electrifying Station-1

        """
        trips_file_name = "trips_for_optimizer_deep.csv"
        sched, scen, args = self.test_basic_run(trips_file_name=trips_file_name)
        args.input_schedule = file_root / trips_file_name
        config_path = file_root / "optimizer.cfg"
        conf = opt_util.read_config(config_path)

        solvers = ["quick", "spiceev"]
        node_choices = ["step-by-step", "brute"]
        conf.opt_type = "deep"
        for solver in solvers:
            for node_choice in node_choices:
                conf.solver = solver
                conf.node_choice = node_choice
                opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
                assert "Station-1" not in opt_sched.stations
                assert "Station-2" in opt_sched.stations
                assert "Station-3" in opt_sched.stations

        trips_file_name = "trips_extended.csv"
        sched, scen, args = self.test_basic_run(trips_file_name=trips_file_name)
        args.input_schedule = file_root / trips_file_name
        config_path = file_root / "optimizer.cfg"
        conf = opt_util.read_config(config_path)

        solvers = ["quick", "spiceev"]
        node_choices = ["step-by-step", "brute"]
        conf.opt_type = "deep"
        opt_stat = None
        for solver in solvers:
            for node_choice in node_choices:
                conf.solver = solver
                conf.node_choice = node_choice
                opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
                if opt_stat is None:
                    opt_stat = {stat for stat in opt_sched.stations}
                    assert len(opt_stat) > 0
                else:
                    assert opt_stat == set(opt_sched.stations)

    def test_critical_stations_optimization(self, caplog):
        """Test if station 2 and 3 are correctly recognized as critical stations.

        :param caplog: pytest fixture, which is automatically created,
            to have access to logging data
        """
        trips_file_name = "trips_for_optimizer_deep.csv"
        sched, scen, args = self.test_basic_run(trips_file_name=trips_file_name)
        args.input_schedule = file_root / trips_file_name
        config_path = file_root / "optimizer.cfg"
        conf = opt_util.read_config(config_path)

        conf.check_for_must_stations = True
        conf.solver = "quick"
        conf.node_choice = "step-by-step"
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
        assert ("must stations {'Station-3', 'Station-2'}" in caplog.text or
                "must stations {'Station-2', 'Station-3'}" in caplog.text)
