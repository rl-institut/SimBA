""" Reusable functions that support tests
"""
import sys
from simba import util
from simba.data_container import DataContainer
from simba.simulate import pre_simulation
from tests.conftest import example_root


def generate_basic_schedule():
    sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
    args = util.get_args()
    args.preferred_charging_type = "oppb"
    mandatory_args = {
        "min_recharge_deps_oppb": 0,
        "min_recharge_deps_depb": 0,
        "gc_power_opps": 1000,
        "gc_power_deps": 1000,
        "cs_power_opps": 100,
        "cs_power_deps_depb": 50,
        "cs_power_deps_oppb": 150,
        "desired_soc_deps": 1,
    }
    vars(args).update(mandatory_args)
    data_container = DataContainer().fill_with_args(args)
    generated_schedule, args = pre_simulation(args, data_container)
    return generated_schedule, args
