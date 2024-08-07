""" Reusable functions that support tests
"""
from argparse import Namespace
from simba import schedule, trip, consumption, util


def generate_basic_schedule():
    schedule_path = 'data/examples/trips_example.csv'
    station_path = 'data/examples/electrified_stations.json'
    temperature_path = 'data/examples/default_temp_winter.csv'
    lol_path = 'data/examples/default_level_of_loading_over_day.csv'
    with open("data/examples/vehicle_types.json", 'r', encoding='utf-8') as f:
        vehicle_types = util.uncomment_json_file(f)
    with open(station_path, 'r', encoding='utf-8') as f:
        stations = util.uncomment_json_file(f)
    trip.Trip.consumption = consumption.Consumption(
        vehicle_types, Namespace(
            outside_temperature_over_day_path=temperature_path,
            level_of_loading_over_day_path=lol_path
        )
    )

    mandatory_args = {
        "min_recharge_deps_oppb": 0,
        "min_recharge_deps_depb": 0,
        "gc_power_opps": 1000,
        "gc_power_deps": 1000,
        "cs_power_opps": 100,
        "cs_power_deps_depb": 50,
        "cs_power_deps_oppb": 150
    }
    generated_schedule =\
        schedule.Schedule.from_csv(schedule_path, vehicle_types, stations, **mandatory_args)
    generated_schedule.calculate_consumption()

    return generated_schedule
