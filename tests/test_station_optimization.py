from tests.helpers import generate_basic_schedule
import ebus_toolbox.optimizer_util as util
from ebus_toolbox.station_optimization import run_optimization
import tests.test_schedule
import pathlib

test_root = pathlib.Path(__file__).parent
file_root = test_root / "test_input_files"

class TestStationOptimization:

    def test_basic_optimization(self):
        trips_file_name = "trips_for_optimizer.csv"
        schedule_test = tests.test_schedule.TestSchedule()

        sched, scen, args = schedule_test.test_basic_run(trips_file_name)

        config_path = "./tests/test_input_files/optimizer.cfg"
        conf = util.read_config(config_path)
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, this_args=args)




    def test_deep_optimization(self):
        """ Check if deep analysis finds the prepared optimal solution for the test case.
        The Test case is a a 3 star like network with two rotations which are negative without
        electrification.
        Rotation 1:
        Depot --> Station-1 ---> Station-2 -->Station-1-->Depot
        Rotation 2:
        Depot --> Station-1 ---> Station-3 -->Station-1-->Depot

        Greedy optimization will electrify Station-1 first since its helpful for both rotations.
        Since electrifiying Station-1 is not enough to electrify either rotations, greedy
        optimization will electrify Station-2 and Station-3 as well. Deep Analysis will expand
        the checked nodes --> It should find that electrifiying Station-2 and Station-3 is enough

        """
        trips_file_name = "trips_for_optimizer_deep.csv"
        schedule_test = tests.test_schedule.TestSchedule()

        sched, scen, args = schedule_test.test_basic_run(trips_file_name)
        args.input_schedule = file_root / trips_file_name
        config_path = "./tests/test_input_files/optimizer.cfg"
        conf = util.read_config(config_path)
        conf.solver = "quick"
        conf.opt_type = "deep"
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, this_args=args)
        assert "Station-1" not in opt_sched.stations
        assert "Station-2" in opt_sched.stations
        assert "Station-3" in opt_sched.stations

        conf.solver = "spiceev"
        conf.opt_type = "deep"
        # sched = deepcopy(schedb)
        # scen = deepcopy(scenb)
        # args = deepcopy(argsb)
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, this_args=args)
        assert "Station-1" not in opt_sched.stations
        assert "Station-2" in opt_sched.stations
        assert "Station-3" in opt_sched.stations



t = TestStationOptimization()
t.test_deep_optimization()