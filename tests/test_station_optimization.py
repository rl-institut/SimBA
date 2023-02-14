from tests.helpers import generate_basic_schedule
import ebus_toolbox.optimizer_util as util
from ebus_toolbox.station_optimization import run_optimization
class TestStationOptimization:

    def test_spice_ev_optimization(self):
        util.print_time()
        config_path = "./tests/test_input_files/optimizer.cfg"
        conf = util.read_config(config_path)
        sched, scen = 123
        opt_sched, opt_scen = run_optimization(conf, sched=None, scen=None, this_args=None)
        util.print_time()
        import pickle
        with open("schedule_opt.pickle", "wb") as f:
            pickle.dump(opt_sched, f)
        with open("scenario_opt.pickle", "wb") as f:
            pickle.dump(opt_scen, f)
