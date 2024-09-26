import argparse
import csv
import json
import logging
from pathlib import Path
import shutil
import subprocess
from datetime import datetime, timedelta

from spice_ev.strategy import STRATEGIES
from spice_ev.util import set_options_from_config


def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()


def save_version(file_path):
    with open(file_path, "w", encoding='utf-8') as f:
        f.write("Git Hash SimBA:" + get_git_revision_hash())


def save_input_file(file_path, args):
    """ Copy given file to output folder, to ensure reproducibility.

    *file_path* must exist and *output_directory_input* must be set in *args*.
    If either condition is not met or the file has already been copied, nothing is done.

    :param file_path: source file
    :type file_path: string or Path
    :param args: general info, output_directory_input is required
    :type args: Namespace
    """
    if file_path is None:
        return
    output_directory_input = vars(args).get("output_directory_input", None)
    if output_directory_input is None:
        # input directory was not created
        return
    source_path = Path(file_path)
    target_path = output_directory_input / source_path.name
    if not source_path.exists():
        # original file missing
        return
    if target_path.exists():
        # already saved
        return
    shutil.copy(source_path, target_path)


def uncomment_json_file(f, char='//'):
    """ Remove comments from JSON file.

    Wouldn't it be nice to have comments in JSON data? Now you can!
    Both full-line comments and trailing lines are supported.

    :param f: file to read
    :type f: JSON file handle
    :param char: char sequence used for commenting, defaults to '//'
    :type char: string
    :return: JSON file content
    :rtype: dict
    """
    uncommented_data = ""
    for line in f:
        comment_idx = line.find(char)
        if comment_idx == -1:
            # no comment in line
            uncommented_data += line
        else:
            # remove comment from line
            uncommented_data += line[:comment_idx]
    return json.loads(uncommented_data)


def get_mean_from_hourly_dict(hourly_dict: dict, start: datetime, end: datetime) -> float:
    """ Get the mean value from hourly data.

    Use daterange from start to end for calculating a mean value by looking up hourly data.
    :param hourly_dict: dictionary with hourly keys and data
    :type hourly_dict: dict
    :param start: start of the range for interpolation
    :type start: datetime
    :param end: end of the range for interpolation
    :type end: datetime
    :return: mean value
    :rtype: float
    """
    total_duration = end - start
    # special case for shared hour of the same day
    if total_duration < timedelta(hours=1) and start.hour == end.hour:
        return hourly_dict.get(start.hour)

    timestep = timedelta(hours=1)
    # proportionally add the start value until the next hour
    next_full_hour = (start+timestep).replace(hour=0, minute=0, second=0, microsecond=0)
    value = hourly_dict.get(start.hour) * (next_full_hour - start).total_seconds()
    start = next_full_hour
    for dt in daterange(start, end, timestep):
        # proportionally apply value according to seconds inside the current hour.
        step_duration = min((end - dt).total_seconds(), 3600)
        value += (hourly_dict.get(dt.hour) * step_duration)
    # divide by total seconds to get mean value
    value /= (total_duration.total_seconds())
    return value


def get_csv_delim(path, other_delims=set()):
    """ Get the delimiter of a character separated file.
     Checks the file for ",", "tabulator" and ";" as well as optional other characters
     or strings.
     In case no clear delimiter is found "," as delimiter is returned

    :param path: Path to file to be checked
    :type path: str
    :param other_delims: Other delimiters besides the default delimiters to be checked
    :type other_delims: set() of str
    :return: delimiter
    :rtype: str
    """

    # create union of default and optional other delimiters
    possible_delims = {",", ";", "\t"}.union(other_delims)
    # create a dict which counts the occurrences of the delimiter per row

    with open(path, "r", encoding='utf-8') as f:
        # count delimiters in first line
        line = f.readline()
        counters = {d: line.count(d) for d in possible_delims if line.count(d) > 0}

        # count delimiters in all following lines, rejecting those with different count than before
        for line_nr, line in enumerate(f):
            # for every delimiter in the dictionary. Casting to set creates new instance
            # needed in case of counter changing during the iteration.
            possible_delims = set(counters.keys())
            for delim in possible_delims:
                # compare the counted amount with the first row values
                amount = line.count(delim)
                # delete the counter if it is different to the first row
                if counters[delim] != amount:
                    del counters[delim]
            # if only one delimiter is remaining
            if len(counters) == 1:
                # take the last item and return the key
                return counters.popitem()[0]
            # if not even a single delimiter is remaining
            elif not counters:
                logging.warning("Warning: Delimiter could not be found.\n"
                                "Returning standard Delimiter ','")
                return ","
    #  multiple delimiters are possible. Every row was checked but more than 1 delimiter
    # has the same amount of occurrences (>0) in every row.
    logging.warning("Warning: Delimiter could not be found.\n"
                    "Returning standard delimiter ','")
    return ","


def get_dict_from_csv(column, file_path, index):
    """ Get a dictonary with the key of a numeric index and the value of a numeric column

    :param column: column name for dictionary values. Content needs to be castable to float
    :type column: str
    :param file_path: file path
    :type file_path: str or Path
    :param index: column name of the index / keys of the dictionary.
        Content needs to be castable to float
    :return: dictionary with numeric keys of index and numeric values of column
    """
    output = dict()
    with open(file_path, "r") as f:
        delim = get_csv_delim(file_path)
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            output[float(row[index])] = float(row[column])
    return output


def nd_interp(input_values, lookup_table):
    """ Interpolates a value from a table.

    Get the interpolated output value for a given input value of n dimensions
    and a given lookup table with n+1 columns.
    Input values outside of the lookup table are clipped to the bounds.

    :param input_values: tuple with n input values
    :type input_values: tuple(floats)
    :param lookup_table: Table with n+1 columns, with the output values in the (n+1)th column
    :type lookup_table: list() of lists()
    :return: interpolated value
    :rtype: float
    """
    # find all unique values in table per column
    dim_sets = [set() for _ in input_values]
    for row in lookup_table:
        for i, v in enumerate(row[:-1]):
            dim_sets[i].add(v)
    dim_values = [sorted(s) for s in dim_sets]
    # find nearest value(s) per column
    # go through sorted column values until last less / first greater
    lower = [None] * len(input_values)
    upper = [None] * len(input_values)
    for i, v in enumerate(input_values):
        # initialize for out of bound values -> Constant value since lower and upper will both
        # be the same boundary value. Still allows for interpolation in other dimensions
        # forcing lower<upper could be implemented for extrapolation beyond the bounds.
        lower[i] = dim_values[i][0]
        upper[i] = dim_values[i][-1]
        for c in dim_values[i]:
            if v >= c:
                lower[i] = c
            if v <= c:
                upper[i] = c
                break
    # find rows in table made up of only lower or upper values
    points = []
    for row in lookup_table:
        for i, v in enumerate(row[:-1]):
            if lower[i] != v and upper[i] != v:
                break
        else:
            points.append(row)

    # Make points unique
    points = [tuple(p) for p in points]
    points = list(set(points))
    # interpolate between points that differ only in current dimension
    for i, x in enumerate(input_values):
        new_points = []
        # find points that differ in just that dimension
        for j, p1 in enumerate(points):
            for p2 in points[j + 1:]:
                for k in range(len(input_values)):
                    if p1[k] != p2[k] and i != k:
                        break
                else:
                    # differing row found
                    x1 = p1[i]
                    y1 = p1[-1]
                    x2 = p2[i]
                    y2 = p2[-1]
                    dx = x2 - x1
                    dy = y2 - y1
                    m = dy / dx
                    n = y1 - m * x1
                    y = m * x + n
                    # generate new point at interpolation
                    p = [v for v in p1]
                    p[i] = x
                    p[-1] = y
                    new_points.append(p)
                    # only couple
                    break
            else:
                # no matching row (singleton dimension?)
                new_points.append(p1)
        points = new_points

    return points[0][-1]


def daterange(start_date, end_date, time_delta):
    """ Iterate over a datetime range using a time_delta step.

    Like range(), the end_value is excluded.
    :param start_date: first value of iteration
    :type start_date: datetime.datetime
    :param end_date: excluded end value of iteration
    :type end_date: datetime.datetime
    :param time_delta: step size of iteration
    :type time_delta: datetime.timedelta
    :yields: iterated value
    :rtype: Iterator[datetime.datetime]
    """
    while start_date < end_date:
        yield start_date
        start_date += time_delta


def cast_float_or_none(val: any) -> any:
    """ Cast a value to float. If a ValueError or TypeError is raised, None is returned

    :param val: value to cast
    :type val: any
    :return: casted value
    """
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def cycling_generator(cycle: []):
    """Generator to loop over lists
    :param cycle: list to cycle through
    :type cycle: list()
    :yields: iterated value
    :rtype:  Iterator[any]
    """
    i = 0
    while True:
        yield cycle[i % len(cycle)]
        i += 1


def setup_logging(args, time_str):
    """ Setup logging.

    :param args: command line arguments. Used: logfile, loglevel, output_path
    :type args: argparse.Namespace
    :param time_str: log file name if args.logfile is not given
    :type time_str: str
    """
    # always to console
    log_level = vars(logging)[args.loglevel.upper()]
    console = logging.StreamHandler()
    console.setLevel(log_level)
    log_handlers = [console]

    if args.logfile is not None and args.output_path is not None:
        # optionally to file in output dir
        if args.logfile:
            log_name = args.logfile
        else:
            log_name = f"{time_str}.log"
        log_path = args.output_path / log_name
        print(f"Writing log to {log_path}")
        file_logger = logging.FileHandler(log_path, encoding='utf-8')
        log_level_file = vars(logging).get((args.loglevel_file or args.loglevel).upper())
        file_logger.setLevel(log_level_file)
        log_handlers.append(file_logger)
        log_level = min(log_level, log_level_file)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=log_handlers
    )
    logging.captureWarnings(True)


def get_buffer_time(trip, default=0):
    """ Get buffer time at arrival station of a trip.

    Buffer time is an abstraction of delays like
    docking procedures and is added to the planned arrival time.

    :param trip: trip to calculate buffer time for
    :type trip: simba.Trip
    :param default: Default buffer time if no station specific buffer time is given. [minutes]
    :type default: dict, numeric
    :return: buffer time in minutes
    :rtype: dict or int

    NOTE: Buffer time dictionaries map hours of the day to a buffer time.
    Keys are ranges of hours and corresponding values provide buffer time in
    minutes for that time range.
    An entry with key "else" is a must if not all hours of the day are covered.
    Example: ``buffer_time = {"10-22": 2, "22-6": 3, "else": 1}``
    """

    schedule = trip.rotation.schedule
    buffer_time = schedule.stations.get(trip.arrival_name, {}).get('buffer_time', default)

    # distinct buffer times depending on time of day can be provided
    # in that case buffer time is of type dict instead of int
    if isinstance(buffer_time, dict):
        # sort dict to make sure 'else' key is last key
        buffer_time = {key: buffer_time[key] for key in sorted(buffer_time)}
        current_hour = trip.arrival_time.hour
        for time_range, buffer in buffer_time.items():
            if time_range == 'else':
                buffer_time = buffer
                break
            else:
                start_hour, end_hour = [int(t) for t in time_range.split('-')]
                if end_hour < start_hour:
                    if current_hour >= start_hour or current_hour < end_hour:
                        buffer_time = buffer
                        break
                else:
                    if start_hour <= current_hour < end_hour:
                        buffer_time = buffer
                        break
    return buffer_time


def replace_deprecated_arguments(args):
    # handling of args with default values
    # Pairs of deprecated names and new names
    deprecated_names = [
        ("input_schedule", "schedule_path"),
        ("vehicle_types", "vehicle_types_path"),
        ("output_directory", "output_path"),
    ]
    for old, new in deprecated_names:
        if vars(args)[old] is not None:
            logging.warning(
                f"Parameter '{old}' is deprecated. Use '{new}' instead. The value of args.{new}: "
                f"{args.__getattribute__(new)} is replaced with {args.__getattribute__(old)} .")
            # Replace value of current name with value of deprecated name
            args.__setattr__(new, args.__getattribute__(old))
        # delete deprecated name
        args.__delattr__(old)

    # Pairs of deprecated names and new names and verbose name
    deprecated_names = [
        ("electrified_stations", "electrified_stations_path", "electrified stations paths"),
        ("cost_parameters_file", "cost_parameters_path", "costs parameter paths"),
        ("rotation_filter", "rotation_filter_path", "rotation filter paths"),
        ("optimizer_config", "optimizer_config_path", "station optimizer config paths"),
    ]

    for old, new, description in deprecated_names:
        if vars(args)[old] is not None:
            assert vars(args)[new] is None, \
                f"Multiple {description} are not supported. Found values for {old} and {new}."
            logging.warning(f"The parameter '{old}' is deprecated. "
                            f"Use '{new}' instead.")
            # Replace value of current name with value of deprecated name
            args.__setattr__(new, args.__getattribute__(old))
        # delete deprecated name
        args.__delattr__(old)
    return args


def mutate_args_for_spiceev(args):
    # arguments relevant to SpiceEV, setting automatically to reduce clutter in config
    args.margin = 1
    args.ALLOW_NEGATIVE_SOC = True
    args.PRICE_THRESHOLD = -100  # ignore price for charging decisions


def get_args():
    parser = get_parser()

    args = parser.parse_args()

    # arguments relevant to SpiceEV, setting automatically to reduce clutter in config
    mutate_args_for_spiceev(args)

    # If a config is provided, the config will overwrite previously parsed arguments
    set_options_from_config(args, check=parser, verbose=False)

    # Check if deprecated arguments were given and change them accordingly
    args = replace_deprecated_arguments(args)

    # rename special options
    args.timing = args.eta

    mandatory_arguments = ["schedule_path", "electrified_stations_path"]
    missing = [a for a in mandatory_arguments if vars(args).get(a) is None]
    if missing:
        raise Exception("The following arguments are required: {}".format(", ".join(missing)))

    return args


def get_parser():
    parser = argparse.ArgumentParser(
        description='SimBA - Simulation toolbox for Bus Applications.')
    parser.add_argument('--scenario-name', help='Identifier of scenario, appended to results')

    # #### Paths #####
    parser.add_argument('--schedule-path',
                        help='Path to CSV file containing all trips of schedule to be analyzed')
    parser.add_argument('--output-path', default="data/sim_outputs",
                        help='Location where all simulation outputs are stored')
    parser.add_argument('--electrified-stations-path', help='include electrified_stations json')
    parser.add_argument('--vehicle-types-path', default="data/examples/vehicle_types.json",
                        help='location of vehicle type definitions')
    parser.add_argument('--station-data-path', default=None,
                        help='Use station data to back calculation of consumption with height\
                         information of stations')
    parser.add_argument('--outside_temperature_over_day_path', default=None,
                        help="Use csv. data with 'hour' and temperature' columns to set \
                        temperatures in case they are not in trips.csv")
    parser.add_argument('--level_of_loading_over_day_path', default=None,
                        help="Use csv. data with 'hour' and level_of_loading' columns to set \
                        level of loading in case they are not in trips.csv")
    parser.add_argument('--cost-parameters-path', default=None,
                        help='include cost parameters json, needed if cost_calculation==True')
    parser.add_argument('--rotation-filter-path', default=None,
                        help='Use json data with rotation ids')

    # #### Modes #####
    mode_choices = [
        'sim', 'neg_depb_to_oppb', 'neg_oppb_to_depb', 'service_optimization',
        'station_optimization', 'remove_negative', 'split_negative_depb', 'report']
    parser.add_argument('--mode', default=['sim', 'report'], nargs='*', choices=mode_choices,
                        help=f"Specify what you want to do. Choose one or more from \
                        {', '.join(mode_choices)}. \
                        sim runs a single simulation with the given inputs. \
                        neg_depb_to_oppb changes charging type of negative depb rotations. \
                        neg_oppb_to_depb changes charging type of negative oppb rotations. \
                        service optimization finds the largest set of electrified rotations. \
                        station_optimization finds the smallest set of electrified stations.\
                        remove_negative removes all negative rotations.\
                        split_negative_depb splits and merges negative depb rotations. \
                        report generates simulation output files.")

    # #### Flags #####
    parser.add_argument('--cost-calculation', '-cc', action='store_true',
                        help='Calculate costs')
    parser.add_argument('--check-rotation-consistency', action='store_true',
                        help='Check rotation assumptions when building schedule.')
    parser.add_argument('--skip_inconsistent_rotations', action='store_true',
                        help='Remove rotations from schedule that violate assumptions. ')
    parser.add_argument('--show-plots', action='store_true',
                        help='show plots for users to view in "report" mode')
    parser.add_argument('--extended-output-plots', action='store_true',
                        help='show extended plots')
    parser.add_argument('--propagate-mode-errors', action='store_true',
                        help='Re-raise errors instead of continuing during simulation modes')
    parser.add_argument('--create-scenario-file', help='Write scenario.json to file')
    parser.add_argument('--create-trips-in-report', action='store_true',
                        help='Write a trips.csv during report mode')
    parser.add_argument('--optimizer-config-path', default=None,
                        help="For station_optimization an optimizer_config is needed. \
                        Input a path to an .cfg file or use the default_optimizer.cfg")
    parser.add_argument('--rotation-filter-variable', default=None,
                        choices=[None, 'include', 'exclude'],
                        help='set mode for filtering schedule rotations')
    parser.add_argument('--zip-output', '-z', action='store_true',
                        help='compress output folder after simulation')

    # #### Charging strategy #####
    parser.add_argument('--preferred-charging-type', '-pct', default='depb',
                        choices=['depb', 'oppb'], help="Preferred charging type. Choose one\
                        from {depb, oppb}. opp stands for opportunity.")
    parser.add_argument('--strategy-deps', default='balanced', choices=STRATEGIES,
                        help='strategy to use in depot')
    parser.add_argument('--strategy-opps', default='greedy', choices=STRATEGIES,
                        help='strategy to use at station')

    # #### Cost calculation strategy #####
    parser.add_argument('--cost-calculation-strategy-deps', choices=STRATEGIES,
                        help='Strategy for cost calculation to use in depot')
    parser.add_argument('--cost-calculation-strategy-opps', choices=STRATEGIES,
                        help='Strategy for cost calculation to use at station')

    parser.add_argument('--strategy-options-deps', default={},
                        type=lambda s: s if type(s) is dict else json.loads(s),
                        help='special strategy options to use in depot')
    parser.add_argument('--strategy-options-opps', default={},
                        type=lambda s: s if type(s) is dict else json.loads(s),
                        help='special strategy options to use at electrified station')

    # #### Physical setup of environment #####
    parser.add_argument('--gc-power-opps', metavar='POPP', type=float, default=100000,
                        help='max power of grid connector at opp stations')
    parser.add_argument('--gc-power-deps', metavar='PDEP', type=float, default=100000,
                        help='max power of grid connector at depot stations')
    parser.add_argument('--cs-power-opps', metavar='CSPOPP', type=float, default=300,
                        help='max power of charging station at opp stations')
    parser.add_argument('--cs-power-deps-depb', metavar='CSPDEPDEP', type=float, default=150,
                        help='max power of charging station at depot stations for depot busses')
    parser.add_argument('--cs-power-deps-oppb', metavar='CSPDEPOPP', type=float, default=150,
                        help='max power of charging station at depot stations for opp busses')
    parser.add_argument('--desired-soc-deps', metavar='SOC', type=float, default=1.0,
                        help='set minimum desired SOC (0 - 1) for depot charging')
    parser.add_argument('--desired-soc-opps', metavar='SOC', type=float, default=1.0,
                        help='set minimum desired SOC (0 - 1) for opportunity charging')
    # Disposition in depots
    parser.add_argument('--min-recharge-deps-oppb', default=1,
                        help='Minimum fraction of capacity for recharge when leaving the depot.')
    parser.add_argument('--min-recharge-deps-depb', default=1,
                        help='Minimum fraction of capacity for recharge when leaving the depot.')
    parser.add_argument('--min-charging-time', help='define minimum time of charging',
                        default=0)
    parser.add_argument('--assign-strategy', default='adaptive',
                        choices=['adaptive', 'fixed_recharge'],
                        help='Strategy for vehicle disposition.')
    parser.add_argument('--default-buffer-time-opps', help='time to subtract from of standing time '
                        'at opp station to simulate docking procedure.', default=0)
    parser.add_argument('--default-buffer-time-deps', default=0,
                        help='time to subtract from of standing time '
                        'at depot station to simulate docking procedure.')
    parser.add_argument('--default-voltage-level', help='Default voltage level for '
                        'charging stations if not set in electrified_stations file',
                        default='MV', choices=['HV', 'HV/MV', 'MV', 'MV/LV', 'LV'])
    parser.add_argument('--peak-load-window-power-opps', type=float, default=1000,
                        help='reduced power of opp stations during peak load windows')
    parser.add_argument('--peak-load-window-power-deps', type=float, default=1000,
                        help='reduced power of depot stations during peak load windows')
    parser.add_argument('--default-mean-speed', type=float, default=30,
                        help='Default assumed mean speed for busses in km/h')
    parser.add_argument('--default-depot-distance', type=float, default=5,
                        help='Default assumed average distance from any station to a depot in km')
    # #### Simulation Parameters #####
    parser.add_argument('--days', metavar='N', type=int, default=None,
                        help='set duration of scenario as number of days')
    parser.add_argument('--interval', metavar='MIN', type=int, default=1,
                        help='set number of minutes for each timestep (Δt)')
    parser.add_argument('--signal-time-dif', default=10,
                        help='time difference between signal time and actual '
                             'start time of a vehicle event im min.')
    parser.add_argument('--eta', action='store_true',
                        help='Show estimated time to finish simulation after each step, '
                             'instead of progress bar. Not recommended for fast computations.')
    parser.add_argument('--skip-flex-report', action='store_true',
                        help='Skip flex band creation when generating reports.')

    # #### LOGGING PARAMETERS #### #
    parser.add_argument('--loglevel', default='INFO', type=str.upper,
                        choices=logging._nameToLevel.keys(), help='Log level.')
    parser.add_argument('--logfile', default='', help='Log file suffix. null: no log file.')
    parser.add_argument('--loglevel_file', default='', type=str.upper,
                        choices=list(logging._nameToLevel.keys()) + [''],
                        help='Log level for file logger.')

    # #### SpiceEV PARAMETERS ONLY DEFAULT VALUES NOT IN SimBA CONFIG #####
    parser.add_argument('--seed', default=1, type=int, help='set random seed')
    parser.add_argument('--include-price-csv',
                        help='include CSV for energy price. \
                            You may define custom options with --include-price-csv-option')
    parser.add_argument('--include-price-csv-option', '-po', metavar=('KEY', 'VALUE'),
                        nargs=2, default=[], action='append',
                        help='append additional argument to price signals')
    parser.add_argument('--time-windows', metavar='FILE',
                        help='use peak load windows to force lower power '
                        'during times of high grid load')

    # Deprecated options for downwards compatibility
    parser.add_argument('--input-schedule',  default=None,
                        help='Deprecated use "schedule-path" instead')
    parser.add_argument('--electrified-stations', default=None,
                        help='Deprecated use "electrified-stations-path" instead')
    parser.add_argument('--vehicle-types', default=None,
                        help='Deprecated use "vehicle-types-path" instead')
    parser.add_argument('--cost-parameters-file', default=None,
                        help='Deprecated use "cost-parameters-path" instead')
    parser.add_argument('--rotation-filter', default=None,
                        help='Deprecated use "rotation-filter-path" instead')
    parser.add_argument('--output-directory', default=None,
                        help='Deprecated use "output-path" instead')
    parser.add_argument('--optimizer_config', default=None,
                        help='Deprecated use "optimizer-config-path" instead')

    parser.add_argument('--config', help='Use config file to set arguments')
    return parser
