import pandas as pd
from ebus_toolbox.sensitivity_analysis import get_buffer_times
from ebus_toolbox.sensitivity_analysis import get_default_hpc
from ebus_toolbox.sensitivity_analysis import get_temperature
from ebus_toolbox.sensitivity_analysis import get_battery_aging
from ebus_toolbox.sensitivity_analysis import get_reduced_power

def prepare_sensitivity(args):

    buffertimes = get_buffer_times()
    args.default_buffer_time_opps = buffertimes

    default_hpc = get_default_hpc()
    args.electrified_stations = default_hpc

    battery_aging = get_battery_aging()
    args.vehicle_types = battery_aging

    temperature = get_temperature()
    temperature.to_csv('data/bvg_test/default_temp_sensitivity.csv')

    reduced_power = get_reduced_power()
    args.cs_power_opps = reduced_power[0]
    args.cs_power_deps_depb = reduced_power[1]
    args.cs_power_deps_oppb = reduced_power[2]

def run_sensitivity (args, scenario_id):

    print('Run sensitivity analysis')

    # true/false in config für arten der störgrößen, z.b. delay = true -> run sensitivity_analysis.delay -> ändere args
    # verspätung immer wieder neu aufrufen, dann kommen jedes mal neue werte
    # metadaten(args) speichern
    export_args = {'scenario': args}
    export_str = 'meta_scenario_' + str(scenario_id) + '.txt'
    export_file = open(export_str, 'wt')
    export_file.write(str(export_args))
    export_file.close()
    return args

# run_sensitivity in 2 funktionen aufteilen -> prepare und run
