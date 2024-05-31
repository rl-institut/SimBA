""" Reusable functions that support tests
"""
import sys
from copy import deepcopy

from simba import trip, consumption, util
from simba.data_container import DataContainer
from simba.simulate import pre_simulation
from tests.conftest import example_root


def generate_basic_schedule(cache=[None]):
    if cache[0] is not None:
        return deepcopy(cache[0])
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
    cache[0] = deepcopy((generated_schedule, args))
    return generated_schedule, args


def initialize_consumption(vehicle_types):
    data_container = DataContainer()
    data_container.add_vehicle_types(vehicle_types)
    data_container.add_consumption_data_from_vehicle_type_linked_files()
    trip.Trip.consumption = consumption.Consumption(vehicle_types)
    for name, df in data_container.consumption_data.items():
        trip.Trip.consumption.set_consumption_interpolation(name, df)
    return trip.Trip.consumption
