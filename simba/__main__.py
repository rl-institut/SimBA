from datetime import datetime
import logging
from pathlib import Path
import shutil

from simba import simulate, util

if __name__ == '__main__':
    args = util.get_args()

    time_str = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    if vars(args).get("scenario_name"):
        dir_name = time_str + '_' + args.scenario_name
    else:
        dir_name = time_str
    if args.output_directory is not None:
        args.output_directory = Path(args.output_directory) / dir_name
        # create subfolder for specific sim results with timestamp.
        # if folder doesn't exist, create folder.
        # needs to happen after set_options_from_config since
        # args.output_directory can be overwritten by config
        args.output_directory_input = args.output_directory / "input_data"
        try:
            args.output_directory_input.mkdir(parents=True, exist_ok=True)
        except NotADirectoryError:
            # can't create new directory (may be write protected): no output
            args.output_directory = None
    if args.output_directory is not None:
        # copy input files to output to ensure reproducibility
        copy_list = [args.config, args.electrified_stations_path, args.vehicle_types_path]
        if "station_optimization" in args.mode:
            copy_list.append(args.optimizer_config)

        # only copy cost params if they exist
        if args.cost_parameters_path is not None:
            copy_list.append(args.cost_parameters_path)
        for c_file in map(Path, copy_list):
            shutil.copy(c_file, args.output_directory_input / c_file.name)

        util.save_version(args.output_directory_input / "program_version.txt")

    util.setup_logging(args, time_str)

    try:
        simulate.simulate(args)
    except Exception as e:
        logging.error(e)
        raise
    finally:
        logging.shutdown()
