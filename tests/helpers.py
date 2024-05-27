""" Reusable functions that support tests
"""
from argparse import Namespace

from simba import schedule, trip, consumption, util
from simba.data_container import DataContainer


def generate_basic_schedule():
    schedule_path = 'data/examples/trips_example.csv'
    station_path = 'data/examples/electrified_stations.json'
    temperature_path = 'data/examples/default_temp_winter.csv'
    lol_path = 'data/examples/default_level_of_loading_over_day.csv'
    with open("data/examples/vehicle_types.json", 'r', encoding='utf-8') as f:
        vehicle_types = util.uncomment_json_file(f)

    initialize_consumption(vehicle_types)

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
    generated_schedule = schedule.Schedule.from_csv(
        schedule_path, vehicle_types, station_path, **mandatory_args,
        outside_temperature_over_day_path=temperature_path, level_of_loading_over_day_path=lol_path)
    generated_schedule.assign_vehicles(Namespace(**mandatory_args))
    return generated_schedule


def initialize_consumption(vehicle_types):
    data_container = DataContainer()
    data_container.add_vehicle_types(vehicle_types)
    data_container.add_consumption_data_from_vehicle_type_linked_files()
    trip.Trip.consumption = consumption.Consumption(vehicle_types)
    for name, df in data_container.consumption_data.items():
        trip.Trip.consumption.set_consumption_interpolation(name, df)
    return trip.Trip.consumption
