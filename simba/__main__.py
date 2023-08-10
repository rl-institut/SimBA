from datetime import datetime
import logging
from pathlib import Path
import shutil

from simba import simulate, util

if __name__ == '__main__':
    args = util.get_args()

    time_str = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    args.output_directory = Path(args.output_directory) / time_str
    # create subfolder for specific sim results with timestamp.
    # if folder doesn't exist, create folder.
    # needs to happen after set_options_from_config since
    # args.output_directory can be overwritten by config
    args.output_directory_input = args.output_directory / "input_data"
    args.output_directory_input.mkdir(parents=True, exist_ok=True)

    # copy input files to output to ensure reproducibility
    copy_list = [args.config, args.electrified_stations, args.vehicle_types]
    if "station_optimization" in args.mode:
        copy_list.append(args.optimizer_config)

    # only copy cost params if they exist
    if args.cost_parameters_file is not None:
        copy_list.append(args.cost_parameters_file)
    for c_file in map(Path, copy_list):
        shutil.copy(c_file, args.output_directory_input / c_file.name)

    util.save_version(args.output_directory_input / "program_version.txt")
    util.setup_logging(args, time_str)

    try:
        simulate.simulate(args)
    except Exception as e:
        logging.error(e)
        raise
