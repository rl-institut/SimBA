from datetime import datetime
import logging
from pathlib import Path

from simba import simulate, util

if __name__ == '__main__':
    args = util.get_args()

    time_str = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    if args.output_directory is not None:
        args.output_directory = Path(args.output_directory) / time_str
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

    # copy basic input files to output to ensure reproducibility
    if args.output_directory is not None:
        copy_list = [
            args.config, args.input_schedule,
            args.electrified_stations, args.vehicle_types,
            args.cost_parameters_file]
        for input_file in copy_list:
            util.save_input_file(input_file, args)
        util.save_version(args.output_directory_input / "program_version.txt")

    util.setup_logging(args, time_str)

    try:
        simulate.simulate(args)
    except Exception as e:
        logging.error(e)
        raise
    finally:
        logging.shutdown()
