import argparse
import json
import logging
import subprocess

from spice_ev.strategy import STRATEGIES
from spice_ev.util import set_options_from_config


def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()


def save_version(file_path):
    with open(file_path, "w", encoding='utf-8') as f:
        f.write("Git Hash SimBA:" + get_git_revision_hash())


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


def setup_logging(args, time_str):
    """ Setup logging.

    :param args: command line arguments. Used: logfile, loglevel, output_directory
    :type args: argparse.Namespace
    :param time_str: log file name if args.logfile is not given
    :type time_str: str
    """
    # always to console
    log_level = vars(logging)[args.loglevel.upper()]
    console = logging.StreamHandler()
    console.setLevel(log_level)
    log_handlers = [console]
    if args.logfile is not None and args.output_directory is not None:
        # optionally to file in output dir
        if args.logfile:
            log_name = args.logfile
        else:
            log_name = f"{time_str}.log"
        log_path = args.output_directory / log_name
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

def get_args():
    parser = get_parser()

    args = parser.parse_args()

    # arguments relevant to SpiceEV, setting automatically to reduce clutter in config
    mutate_args_for_spiceev(args)

    # If a config is provided, the config will overwrite previously parsed arguments
    set_options_from_config(args, check=parser, verbose=False)

    # rename special options
    args.timing = args.eta

    missing = [a for a in ["input_schedule", "electrified_stations"] if vars(args).get(a) is None]
    if missing:
        raise Exception("The following arguments are required: {}".format(", ".join(missing)))

    return args


def mutate_args_for_spiceev(args):
    # arguments relevant to SpiceEV, setting automatically to reduce clutter in config
    args.strategy = 'distributed'
    args.margin = 1
    args.ALLOW_NEGATIVE_SOC = True
    args.PRICE_THRESHOLD = -100  # ignore price for charging decisions



def get_parser():
    parser = argparse.ArgumentParser(
        description='SimBA - Simulation toolbox for Bus Applications.')

    # #### Paths #####
    parser.add_argument('--input-schedule',
                        help='Path to CSV file containing all trips of schedule to be analyzed.')
    parser.add_argument('--output-directory', default="data/sim_outputs",
                        help='Location where all simulation outputs are stored')
    parser.add_argument('--electrified-stations', help='include electrified_stations json')
    parser.add_argument('--vehicle-types-path', default="data/examples/vehicle_types.json",
                        help='location of vehicle type definitions')
    parser.add_argument('--station_data_path', default=None,
                        help='Use station data to back calculation of consumption with height\
                         information of stations')
    parser.add_argument('--outside_temperature_over_day_path', default=None,
                        help="Use csv. data with 'hour' and temperature' columns to set \
                        temperatures in case they are not in trips.csv")
    parser.add_argument('--level_of_loading_over_day_path', default=None,
                        help="Use csv. data with 'hour' and level_of_loading' columns to set \
                        level of loading in case they are not in trips.csv")
    parser.add_argument('--cost-parameters-file', default=None,
                        help='include cost parameters json, needed if cost_calculation==True')
    parser.add_argument('--rotation-filter', default=None,
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
    parser.add_argument('--propagate-mode-errors', default=False,
                        help='Re-raise errors instead of continuing during simulation modes')
    parser.add_argument('--create-scenario-file', help='Write scenario.json to file')
    parser.add_argument('--rotation-filter-variable', default=None,
                        choices=[None, 'include', 'exclude'],
                        help='set mode for filtering schedule rotations')

    # #### Charging strategy #####
    parser.add_argument('--preferred-charging-type', '-pct', default='depb',
                        choices=['depb', 'oppb'], help="Preferred charging type. Choose one\
                        from {depb, oppb}. opp stands for opportunity.")
    parser.add_argument('--strategy-deps', default='balanced', choices=STRATEGIES,
                        help='strategy to use in depot')
    parser.add_argument('--strategy-opps', default='greedy', choices=STRATEGIES,
                        help='strategy to use at station')
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
    parser.add_argument('--optimizer_config', default=None,
                        help="For station_optimization an optimizer_config is needed. \
                        Input a path to an .cfg file or use the default_optimizer.cfg")
    parser.add_argument('--time-windows', metavar='FILE',
                        help='use peak load windows to force lower power '
                        'during times of high grid load')

    parser.add_argument('--config', help='Use config file to set arguments')
    return parser

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