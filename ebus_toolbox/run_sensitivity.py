from ebus_toolbox.sensitivity_analysis import get_buffer_times

def run_sensitivity (args, scenario_id):
    print('Run sensitivity analysis')
    buffertimes = get_buffer_times()
    args.default_buffer_time_opps = buffertimes
    # true/false in config für arten der störgrößen, z.b. delay = true -> run sensitivity_analysis.delay -> ändere args
    # verspätung immer wieder neu aufrufen, dann kommen jedes mal neue werte
    # metadaten(args) speichern
    export_args = {'scenario': args}
    export_str = 'meta_scenario_' + str(scenario_id) + '.txt'
    export_file = open(args.output_directory / export_str, 'wt')
    export_file.write(str(export_args))
    export_file.close()
    return args