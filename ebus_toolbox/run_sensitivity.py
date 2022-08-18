import pandas as pd
from ebus_toolbox.sensitivity_analysis import get_buffer_times
from ebus_toolbox.sensitivity_analysis import get_default_hpc
from ebus_toolbox.sensitivity_analysis import get_temperature

def run_sensitivity (args, scenario_id):
    print('Run sensitivity analysis')
    buffertimes = get_buffer_times()
    args.default_buffer_time_opps = buffertimes
    default_hpc = get_default_hpc()
    print(default_hpc)
    args.electrified_stations = default_hpc
    temperature = get_temperature()
    temperature.to_csv('data/bvg_test/default_temp_sensitivity.csv') #todo: consumption wird vorher aufgerufen --> wie bekommen wir die daten vorher da hin?
    # true/false in config für arten der störgrößen, z.b. delay = true -> run sensitivity_analysis.delay -> ändere args
    # verspätung immer wieder neu aufrufen, dann kommen jedes mal neue werte
    # metadaten(args) speichern
    export_args = {'scenario': args}
    export_str = 'meta_scenario_' + str(scenario_id) + '.txt'
    export_file = open(export_str, 'wt')
    export_file.write(str(export_args))
    export_file.close()
    return args


