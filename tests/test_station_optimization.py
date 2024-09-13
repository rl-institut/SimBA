import logging
from copy import copy, deepcopy
import json
from pathlib import Path
import pytest
import random
import re
import sys
import shutil

from simba import station_optimizer
import simba.optimizer_util as opt_util
from simba.data_container import DataContainer
from simba.simulate import pre_simulation
from simba.station_optimization import run_optimization
import simba.util as util
from spice_ev.report import generate_soc_timeseries
from tests.conftest import example_root

file_root = Path(__file__).parent / "test_input_files/optimization"


def slow_join_all_subsets(subsets):
    def join_subsets(subsets):
        subsets = [s.copy() for s in subsets]
        for i in range(len(subsets)):
            for ii in range(i + 1, len(subsets)):
                intersec = subsets[i].intersection(subsets[ii])
                if len(intersec) > 0:
                    subsets[i] = subsets[i].union(subsets[ii])
                    subsets.remove(subsets[ii])
                    return True, subsets
        return False, subsets

    joined_subset = True
    while joined_subset:
        joined_subset, subsets = join_subsets(subsets)
    return subsets


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
        src = example_root / "simba.cfg"
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
            r"(vehicle_types_path\s=.*)(:=\r\n|\r|\n)",
            "vehicle_types_path = " + vehicles_dest_str + r"\g<2>", src_text)

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
        dst = tmp_path / "simba.cfg"
        dst.write_text(src_text)

    def generate_datacontainer_args(self, trips_file_name="trips.csv"):
        """ Check if running a basic example works and if a scenario object is returned.

        :param trips_file_name: file name of the trips file. Has to be inside the test_input_file
            folder
        :type trips_file_name: str
        :return: schedule, scenario, args"""

        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.output_directory = self.tmp_path
        args.results_directory = self.tmp_path
        assert self.tmp_path

        args.input_schedule = file_root / trips_file_name
        data_container = DataContainer().fill_with_args(args)
        return data_container, args

    def generate_schedule_scenario(self, args, data_container):
        generated_schedule, args = pre_simulation(args, data_container)
        scen = generated_schedule.run(args)

        # optimization depends on vehicle_socs, therefore they need to be generated
        generate_soc_timeseries(scen)
        return generated_schedule, scen

    def test_join_all_subsets(self):
        subsets = [{1, 2, 3}, {3, 4, 6, 7}, {7, 8}, {20, 21}, {21, 22}, {6}]
        assert len(slow_join_all_subsets(deepcopy(subsets))) == 2
        assert len(opt_util.join_all_subsets(deepcopy(subsets))) == 2
        random.seed(42)
        subsets = []
        for _ in range(15):
            subset = set()
            for rnd in range(random.randint(0, 10)):
                subset.add(random.randint(0, 30))
            subsets.append(subset)

        joined_subsets1 = slow_join_all_subsets(deepcopy(subsets))
        joined_subsets2 = opt_util.join_all_subsets(deepcopy(subsets))
        assert len(joined_subsets1) == len(joined_subsets2)
        for subset in joined_subsets1:
            assert subset in joined_subsets2

    def test_fast_calculations_and_events(self):
        """ Test if the base optimization finishes without raising errors"""
        trips_file_name = "trips_for_optimizer.csv"
        data_container, args = self.generate_datacontainer_args(trips_file_name)
        data_container.stations_data = {}
        args.preferred_charging_type = "oppb"
        sched, scen = self.generate_schedule_scenario(args, data_container)
        t = sched.rotations["2"].trips[0]
        t.distance = 100000
        sched.calculate_consumption()
        sched.stations["Station-1"] = {"type": "opps", "n_charging_stations": None}
        sched.stations["Station-2"] = {"type": "opps", "n_charging_stations": None}
        scen = sched.run(args)
        generate_soc_timeseries(scen)

        config = opt_util.OptimizerConfig()
        sopt = station_optimizer.StationOptimizer(sched, scen, args, config=config,
                                                  logger=logging.getLogger())

        # create charging dicts which contain soc over time, which is numerically calculated
        sopt.create_charging_curves()
        # remove none values from socs in the vehicle_socs
        sopt.replace_socs_from_none_to_value()
        vehicle_socs_fast = sopt.timeseries_calc(list(sched.stations.keys()))
        for vehicle, socs in scen.vehicle_socs.items():
            assert vehicle_socs_fast[vehicle][-1] == pytest.approx(socs[-1], 0.01, abs=0.01)

        events = sopt.get_low_soc_events(soc_data=vehicle_socs_fast, rel_soc=True)
        assert len(events) == 1
        e = events[0]
        e1 = copy(e)
        vehicle_socs_reduced = {vehicle: [soc - 1 for soc in socs] for vehicle, socs in
                                scen.vehicle_socs.items()}
        events = sopt.get_low_soc_events(soc_data=vehicle_socs_reduced, rel_soc=True)
        assert len(events) == 1
        e2 = events[0]
        assert e1.start_idx == e2.start_idx
        assert e1.end_idx == e2.end_idx
        assert e1.min_soc == e2.min_soc

        vehicle_socs_increased = {vehicle: [min(soc + abs(e1.min_soc) - 0.1, 1) for soc in socs] for
                                  vehicle, socs in scen.vehicle_socs.items()}
        events = sopt.get_low_soc_events(soc_data=vehicle_socs_increased, rel_soc=True)
        e3 = events[0]
        assert e1.start_idx != e3.start_idx
        assert e1.end_idx == e3.end_idx
        assert e1.min_soc != e3.min_soc

        vehicle_socs_more_increased = {
            vehicle: [min(soc + abs(e1.min_soc) + 0.1, 1) for soc in socs] for vehicle, socs in
            scen.vehicle_socs.items()}
        events = sopt.get_low_soc_events(soc_data=vehicle_socs_more_increased, rel_soc=True)
        assert len(events) == 0

    def test_basic_optimization(self):
        """ Test if the base optimization finishes without raising errors"""
        trips_file_name = "trips_for_optimizer.csv"
        data_container, args = self.generate_datacontainer_args(trips_file_name)
        data_container.stations_data = {}
        args.preferred_charging_type = "oppb"
        sched, scen = self.generate_schedule_scenario(args, data_container)
        config_path = example_root / "default_optimizer.cfg"
        conf = opt_util.read_config(config_path)

        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)
        conf.solver = "spiceev"
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, args=args)

    def test_schedule_consistency(self):
        """ Test if the optimization returns all rotations even when some filters are active"""
        trips_file_name = "trips_for_optimizer.csv"
        data_container, args = self.generate_datacontainer_args(trips_file_name)
        args.preferred_charging_type = "oppb"
        sched, scen = self.generate_schedule_scenario(args, data_container)
        config_path = example_root / "default_optimizer.cfg"
        conf = opt_util.read_config(config_path)

        sched_impossible = deepcopy(sched)
        # Make a single trip impossible
        next(iter(sched_impossible.rotations.values())).trips[0].distance = 99999999
        sched_impossible.calculate_consumption()
        scen = sched_impossible.run(args)

        amount_rotations = len(sched_impossible.rotations)
        conf.remove_impossible_rotations = True
        conf.run_only_oppb = False
        opt_sched, opt_scen = run_optimization(conf, sched=sched_impossible, scen=scen, args=args)
        assert len(opt_sched.rotations) == amount_rotations

        # set a single rotation to depot
        next(iter(sched.rotations.values())).set_charging_type("depb")
        # and run again to make everything consistent
        scen = sched.run(args)

        amount_rotations = len(sched.rotations)
        conf.run_only_oppb = True
        opt_sched, opt_scen = run_optimization(conf, sched=deepcopy(sched), scen=scen, args=args)
        assert len(opt_sched.rotations) == amount_rotations

        conf.run_only_oppb = False
        opt_sched, opt_scen = run_optimization(conf, sched=deepcopy(sched), scen=scen, args=args)
        assert len(opt_sched.rotations) == amount_rotations

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
        data_container, args = self.generate_datacontainer_args(trips_file_name)
        data_container.stations_data = {}
        args.preferred_charging_type = "oppb"
        for trip_d in data_container.trip_data:
            trip_d["distance"] *= 15
        sched, scen = self.generate_schedule_scenario(args, data_container)
        config_path = example_root / "default_optimizer.cfg"
        conf = opt_util.read_config(config_path)

        assert len(sched.get_negative_rotations(scen)) == 2

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

    def test_deep_optimization_extended(self):
        trips_file_name = "trips_extended.csv"
        data_container, args = self.generate_datacontainer_args(trips_file_name)
        data_container.stations_data = {}
        args.preferred_charging_type = "oppb"
        sched, scen = self.generate_schedule_scenario(args, data_container)
        # optimization can only be properly tested if negative rotations exist
        assert len(sched.get_negative_rotations(scen)) > 0
        args.input_schedule = file_root / trips_file_name
        config_path = example_root / "default_optimizer.cfg"
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
        data_container, args = self.generate_datacontainer_args(trips_file_name)
        for trip_d in data_container.trip_data:
            trip_d["distance"] *= 15
        data_container.stations_data = {}
        args.preferred_charging_type = "oppb"
        sched, scen = self.generate_schedule_scenario(args, data_container)
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
