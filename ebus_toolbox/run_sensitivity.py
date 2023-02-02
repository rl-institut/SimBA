from ebus_toolbox.sensitivity_analysis import get_buffer_times
from ebus_toolbox.sensitivity_analysis import get_temperature
from ebus_toolbox.sensitivity_analysis import get_default_hpc
from ebus_toolbox.sensitivity_analysis import get_battery_aging
from ebus_toolbox.sensitivity_analysis import get_depot_delay
from ebus_toolbox.sensitivity_analysis import get_reduced_power
from pathlib import Path


def prepare_sensitivity(args):

    temperature = get_temperature()
    temperature.to_csv('data/bvg/default_temp_sensitivity.csv')

    args.outside_temperature_over_day_path = 'data/bvg/default_temp_sensitivity.csv'
    temperature.to_csv(args.output_directory / 'default_temp_sensitivity.csv')

    return args


def run_sensitivity (args, scenario_id):
    print('Run sensitivity analysis')

    buffertimes = get_buffer_times()
    args.default_buffer_time_opps = buffertimes

    default_hpc = get_default_hpc()
    args.electrified_stations = default_hpc

    battery_aging = get_battery_aging()
    args.vehicle_types = battery_aging

    reduced_power_opps, reduced_power_depb, reduced_power_oppb = get_reduced_power()
    args.cs_power_opps = reduced_power_opps
    args.cs_power_deps_depb = reduced_power_depb
    args.cs_power_deps_oppb = reduced_power_oppb

    # depot_delay = get_depot_delay()
    # args.depot_delay = delay_dep

    # metadaten(args) speichern
    export_args = {'scenario': args}
    export_str = 'meta_scenario_' + str(scenario_id) + '.txt'
    export_file = open(args.output_directory / export_str, 'wt')
    export_file.write(str(export_args))
    export_file.close()
    return args