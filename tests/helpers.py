""" Reusable functions that support tests
"""
import json
from ebus_toolbox import schedule, trip, consumption


def generate_basic_schedule():
    schedule_path = 'data/examples/trips_example.csv'
    station_path = 'data/examples/electrified_stations.json'
    temperature_path = 'data/examples/default_temp_winter.csv'
    lol_path = 'data/examples/default_level_of_loading_over_day.csv'
    with open("data/examples/vehicle_types.json", 'r', encoding='utf-8') as f:
        vehicle_types = json.load(f)
    trip.Trip.consumption = consumption.Consumption(vehicle_types,
                                                    outside_temperatures=temperature_path,
                                                    level_of_loading_over_day=lol_path)

    mandatory_args = {
        "min_recharge_deps_oppb": 0,
        "min_recharge_deps_depb": 0,
        "gc_power_opps": 1000,
        "gc_power_deps": 1000,
        "cs_power_opps": 100,
        "cs_power_deps_depb": 50,
        "cs_power_deps_oppb": 150
    }

    return schedule.Schedule.from_csv(schedule_path, vehicle_types, station_path, **mandatory_args)
