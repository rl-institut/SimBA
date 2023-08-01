import pathlib
import re
import shutil
from copy import copy

import pytest
import sys
import json
from tests.conftest import example_root
import ebus_toolbox.optimizer_util as opt_util
import ebus_toolbox.util as util
from ebus_toolbox.station_optimization import run_optimization
from ebus_toolbox.trip import Trip
from ebus_toolbox.consumption import Consumption
from ebus_toolbox.schedule import Schedule
from spice_ev.report import generate_soc_timeseries

file_root = pathlib.Path(__file__).parent / "test_input_files/optimization"


class TestStationOptimization:
    tmp_path = None

    # this format allows calling setup path directly, when not using pytest
    # for debugging purposes
    @pytest.fixture(autouse=True)
    def setup_test_fixture(self, tmp_path):
        self.setup_test(tmp_path)

    # setup for test including setting up tmp_path and file creation
    def setup_test(self, tmp_path):
        # This will guarantee proper file setup before every test. Generates files and inputs
        # by adjusting example files and storing them in a temporary directory.

        self.tmp_path = tmp_path

        # Create a temporary config file as copy from the example configuration.
        src = example_root / "ebus_toolbox.cfg"
        src_text = src.read_text()

        # don't show plots. spaces are optional, so use regex
        src_text = re.sub(r"show_plots\s*=\s*true", "show_plots = false", src_text)

        vehicles_dest = tmp_path / "vehicle_types.json"
        # store vehicles temporarily and use them in the config file
        shutil.copy(example_root / "vehicle_types.json", vehicles_dest)

        # opens the vehicle type file, and adjust capacity and mileage of all vehicles
        self.vehicle_types = adjust_vehicle_file(vehicles_dest, capacity=50, mileage=10)

        # remove escape characters from string
        vehicles_dest_str = str(vehicles_dest).replace('\\', '/')
        # replace line which defines vehicle_types up to line break. line break is concatenated in
        # the replacement, to keep format
        src_text = re.sub(
            r"(vehicle_types\s=.*)(:=\r\n|\r|\n)",
            "vehicle_types = " + vehicles_dest_str + r"\g<2>", src_text)

        # Use the default electrified stations from example folder but change some values
        with open(example_root / "electrified_stations.json", "r", encoding='utf-8') as file:
            self.electrified_stations = util.uncomment_json_file(file)
        # only keep Station-0 electrified and remove the other staitons
        self.electrified_stations = {"Station-0": self.electrified_stations["Station-0"]}
        del self.electrified_stations["Station-0"]["external_load"]
        del self.electrified_stations["Station-0"]["battery"]
        del self.electrified_stations["Station-0"]["energy_feed_in"]

        # store the adjusted electrified_stations temporarily and use them in the config file
        electrified_stations_dest = tmp_path / "electrified_stations.json"
        with open(electrified_stations_dest, "w", encoding='utf-8') as file:
            json.dump(self.electrified_stations, file)

        # remove escape characters from string. \1 refers to the replacement of the first group
        # in the regex expression, i.e. not replacing the newline characters
        electrified_stations_dest_str = str(electrified_stations_dest).replace('\\', '/')
        src_text = re.sub(
            r"(electrified_stations\s=.*)(:=\r\n|\r|\n)",
            "electrified_stations = " + electrified_stations_dest_str + r"\g<2>", src_text)

        src_text = re.sub(
            r"(preferred_charging_type\s=.*)(:=\r\n|\r|\n)",
            "preferred_charging_type = oppb"r"\g<2>", src_text)

        # change config file with adjusted temporary paths to vehicles and electrified stations
        dst = tmp_path / "ebus_toolbox.cfg"
        dst.write_text(src_text)

    def basic_run(self, trips_file_name="trips.csv"):
        """ Check if running a basic example works and if a scenario object is returned.

        :param trips_file_name: file name of the trips file. Has to be inside the test_input_file
            folder
        :type trips_file_name: str
        :return: schedule, scenario"""
        path_to_trips = file_root / trips_file_name
        sys.argv = ["foo", "--config", str(self.tmp_path / "ebus_toolbox.cfg")]
        args = util.get_args()
        args.input_schedule = path_to_trips
        Trip.consumption = Consumption(self.vehicle_types,
                                       outside_temperatures=None,
                                       level_of_loading_over_day=None)
        args2 = copy(args)
        del args2.vehicle_types
        generated_schedule = Schedule.from_csv(path_to_trips, self.vehicle_types,
                                               self.electrified_stations,
                                               **vars(args2))

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
        sched, scen, args = self.basic_run(trips_file_name)
        config_path = example_root / "default_optimizer.cfg"
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
        sched, scen, args = self.basic_run(trips_file_name=trips_file_name)
        args.input_schedule = file_root / trips_file_name
        config_path = example_root / "default_optimizer.cfg"
        conf = opt_util.read_config(config_path)

        solvers = ["quick", "spiceev"]
        node_choices = ["step-by-step", "brute"]
        conf.opt_type = "deep"
        for solver in solvers:
            for node_choice in node_choices:
                conf.solver = solver
                conf.node_choice = node_choice
                opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
                assert len(opt_sched.get_negative_rotations(opt_scen)) == 0
                assert "Station-1" not in opt_sched.stations
                assert "Station-2" in opt_sched.stations
                assert "Station-3" in opt_sched.stations

        trips_file_name = "trips_extended.csv"
        # adjust mileage so scenario is not possible without adding electrification
        self.vehicle_types = adjust_vehicle_file(args.vehicle_types, mileage=2, capacity=150)
        sched, scen, args = self.basic_run(trips_file_name=trips_file_name)
        # optimization can only be properly tested if negative rotations exist
        assert len(sched.get_negative_rotations(scen)) > 0
        args.input_schedule = file_root / trips_file_name

        opt_stat = None
        for solver in solvers:
            for node_choice in node_choices:
                conf.solver = solver
                conf.node_choice = node_choice
                opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
                neg_rots = opt_sched.get_negative_rotations(opt_scen)
                assert len(neg_rots) == 0
                if opt_stat is None:
                    opt_stat = {stat for stat in opt_sched.stations}
                    assert len(opt_stat) > 0
                else:
                    assert len(opt_stat) == len(set(opt_sched.stations))

    def test_critical_stations_optimization(self, caplog):
        """Test if station 2 and 3 are correctly recognized as critical stations.

        :param caplog: pytest fixture, which is automatically created,
            to have access to logging data
        """
        trips_file_name = "trips_for_optimizer_deep.csv"
        sched, scen, args = self.basic_run(trips_file_name=trips_file_name)
        args.input_schedule = file_root / trips_file_name
        config_path = example_root / "default_optimizer.cfg"
        conf = opt_util.read_config(config_path)

        conf.check_for_must_stations = True
        conf.solver = "quick"
        conf.node_choice = "step-by-step"
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
        assert ("must stations {'Station-3', 'Station-2'}" in caplog.text or
                "must stations {'Station-2', 'Station-3'}" in caplog.text)


def adjust_vehicle_file(source, capacity=None, mileage=0):
    # use the default vehicles from example folder
    with open(source, "r", encoding='utf-8') as file:
        vehicle_types = util.uncomment_json_file(file)

    for vehicle in vehicle_types:
        for vtype in vehicle_types[vehicle]:
            if capacity is not None:
                vehicle_types[vehicle][vtype]["capacity"] = capacity
            if mileage is not None:
                vehicle_types[vehicle][vtype]["mileage"] = mileage
    with open(source, "w", encoding='utf-8') as file:
        json.dump(vehicle_types, file)
    return vehicle_types
