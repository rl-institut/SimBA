import schedule
from src import scenario
from optimizer_config import read_config, OptimizerConfig

class ScheduleOptimizer:
    def __init__(self, sched:schedule.Schedule(), scen:scenario.Scenario(), args,
                 config:OptimizerConfig(), **kwargs):
        self.schedule = sched
        self.scenario = scen
        self.args = args
        self.config = config
        self.logger = kwargs.get("logger", None)

        if config.reduce_rots:
            self.sched.rotations = {rot: sched.rotations[rot] for rot in config.rots}
