from tests.helpers import generate_basic_schedule
import ebus_toolbox.optimizer_util as util
from ebus_toolbox.station_optimization import run_optimization
import tests.test_schedule

class TestStationOptimization:
    trips_file_name = "trips_for_optimizer.csv"

    def test_spice_ev_optimization(self):
        schedule_test = tests.test_schedule.TestSchedule()

        sched, scen, args = schedule_test.test_basic_run(self.trips_file_name)

        config_path = "./tests/test_input_files/optimizer.cfg"
        conf = util.read_config(config_path)
        opt_sched, opt_scen = run_optimization(conf, sched=sched, scen=scen, this_args=args)


t = TestStationOptimization()
t.test_spice_ev_optimization()