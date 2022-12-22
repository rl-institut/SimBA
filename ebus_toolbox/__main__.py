import argparse
import shutil
from ebus_toolbox import simulate, util
from pathlib import Path
from datetime import datetime

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='eBus-Toolbox - \
        simulation program for electric bus fleets.')
    parser.add_argument('--input-schedule', nargs='?',
                        help='Path to CSV file containing all trips of schedule to be analyzed.')
    parser.add_argument('--mode', default='sim', choices=['sim', 'service_optimization'],
                        help='Specify what you want to do. Choose one from {sim, \
                        service_optimization}. sim runs a single simulation with the given inputs. \
                        service optimization finds the largest set of electrified rotations.')
    parser.add_argument('--output-directory', default="./data/sim_outputs", nargs='?',
                        help='Location where all simulation outputs are stored')
    parser.add_argument('--preferred-charging-type', '-pct', default='depb',
                        choices=['depb', 'oppb'], help="Preferred charging type. Choose one\
                        from {depb, oppb}. opp stands for opportunity.")
    parser.add_argument('--vehicle-types', default="./data/examples/vehicle_types.json",
                        help='location of vehicle type definitions')
    parser.add_argument('--min-recharge-deps-oppb', default=1,
                        help='Minimum fraction of capacity for recharge when leaving the depot.')
    parser.add_argument('--min-recharge-deps-depb', default=1,
                        help='Minimum fraction of capacity for recharge when leaving the depot.')
    parser.add_argument('--check-rotation-consistency', action='store_true',
                        help='Check rotation assumptions when building schedule.')
    parser.add_argument('--ignore-inconsistent-rotations', action='store_true',
                        help='Remove rotations from schedule that violate assumptions. ')
    parser.add_argument('--days', metavar='N', type=int, default=None,
                        help='set duration of scenario as number of days')
    parser.add_argument('--interval', metavar='MIN', type=int, default=15,
                        help='set number of minutes for each timestep (Δt)')
    parser.add_argument('--desired-soc-deps', metavar='SOC', type=float, default=1.0,
                        help='set minimum desired SOC (0 - 1) for depot charging')
    parser.add_argument('--desired-soc-opps', metavar='SOC', type=float, default=0.8,
                        help='set minimum desired SOC (0 - 1) for opportunity charging')
    parser.add_argument('--gc-power-opps', metavar='POPP', type=float, default=None,
                        help='max power of grid connector at opp stations')
    parser.add_argument('--gc-power-deps', metavar='PDEP', type=float, default=None,
                        help='max power of grid connector at depot stations')
    parser.add_argument('--cs-power-opps', metavar='CSPOPP', type=float, default=450,
                        help='max power of charging station at opp stations')
    parser.add_argument('--cs-power-deps-depb', metavar='CSPDEPDEP', type=float, default=100,
                        help='max power of charging station at depot stations for depot busses')
    parser.add_argument('--cs-power-deps-oppb', metavar='CSPDEPOPP', type=float, default=150,
                        help='max power of charging station at depot stations for opp busses')
    parser.add_argument('--battery', '-b', default=[], nargs=2, type=float, action='append',
                        help='add battery with specified capacity in kWh and C-rate \
                            (-1 for variable capacity, second argument is fixed power))')
    parser.add_argument('--seed', default=None, type=int, help='set random seed')
    parser.add_argument('--iterations', default=1, type=int, help='iterations for optimization')
    parser.add_argument('--include-ext-load-csv',
                        help='include CSV for external load. \
                            You may define custom options with --include-ext-csv-option')
    parser.add_argument('--include-ext-csv-option', '-eo', metavar=('KEY', 'VALUE'),
                        nargs=2, action='append',
                        help='append additional argument to external load')
    parser.add_argument('--include-feed-in-csv',
                        help='include CSV for energy feed-in, e.g., local PV. \
                            You may define custom options with --include-feed-in-csv-option')
    parser.add_argument('--include-feed-in-csv-option', '-fo', metavar=('KEY', 'VALUE'),
                        nargs=2, action='append', help='append additional argument to feed-in load')
    parser.add_argument('--include-price-csv',
                        help='include CSV for energy price. \
                            You may define custom options with --include-price-csv-option')
    parser.add_argument('--include-price-csv-option', '-po', metavar=('KEY', 'VALUE'),
                        nargs=2, default=[], action='append',
                        help='append additional argument to price signals')
    parser.add_argument('--start-date', default='2018-01-02',
                        help='Provide start date of simulation in format YYYY-MM-DD.E.g. '
                             '2018-01-31')
    parser.add_argument('--electrified-stations', help='include electrified_stations json',
                        default='data/examples/electrified_stations.json')
    parser.add_argument('--cost-calculation', '-cc', action='store_true',
                        help='Calculate costs')
    parser.add_argument('--cost-parameters-file', help='include cost parameters json', default=None)
    parser.add_argument('--pv-power', type=int, default=0, help='set nominal power for local '
                                                                'photovoltaic power plant in kWp')
    parser.add_argument('--min-charging-time', help='define minimum time of charging',
                        default=0)
    parser.add_argument('--default-buffer-time-opps', help='time to subtract off of standing time '
                        'at opp station to simulate docking procedure.', default=0)
    parser.add_argument('--signal-time-dif', help='time difference between signal time and actual '
                                                  'start time of a vehicle event im min.',
                        default=10)
    parser.add_argument('--visual', '-v', action='store_true', help='Show plots of the results')
    parser.add_argument('--eta', action='store_true',
                        help='Show estimated time to finish simulation after each step, \
                        instead of progress bar. Not recommended for fast computations.')
    parser.add_argument('--save-timeseries', help='Write timesteps to file')
    parser.add_argument('--save-results', help='Write general info to file')
    parser.add_argument('--save-soc', help='Write SOC info to file')
    parser.add_argument('--strategy', '-s', default='greedy',
                        help='Specify the charging strategy. One of {}. You may define \
                        custom options with --strategy-option.'.format('greedy, balanced'))
    parser.add_argument('--margin', '-m', metavar='X', type=float, default=0.05,
                        help=('Add margin for desired SOC [0.0 - 1.0].\
                        margin=0.05 means the simulation will not abort if vehicles \
                        reach at least 95%% of the desired SOC before leaving. \
                        margin=1 -> the simulation continues with every positive SOC value.'))
    parser.add_argument('--strategy-option', '-so', metavar=('KEY', 'VALUE'),
                        nargs=2, action='append',
                        help='Append additional options to the charging strategy.')
    parser.add_argument('--config', help='Use config file to set arguments')
    parser.add_argument('--station_data_path', help='Use station data to back calculation       \
                                                    of consumption with height information of   \
                                                    stations')
    parser.add_argument('--outside_temperature_over_day_path', default=None,
                        help="Use csv. data with 'hour' and temperature' columns to set \
                        temperatures in case they are not in trips.csv")

    parser.add_argument('--level_of_loading_over_day_path', default=None,
                        help="Use csv. data with 'hour' and level_of_loading' columns to set \
                        level of loading in case they are not in trips.csv")

    args = parser.parse_args()
    # arguments relevant to SpiceEV, setting automatically to reduce clutter in config
    args.ALLOW_NEGATIVE_SOC = True
    args.attach_vehicle_soc = True

    util.set_options_from_config(args, check=True, verbose=False)

    args.output_directory = Path(args.output_directory) / \
        str(datetime.now().strftime("%Y-%m-%d-%H-%M-%S") + "_eBus_results")

    # create subfolder for specific sim results with timestamp.
    # if folder doesnt exists, create folder.
    # needs to happen after set_options_from_config since
    # args.output_directory can be overwritten by config
    args.output_directory.mkdir(parents=True, exist_ok=True)

    if not args.save_timeseries:
        args.save_timeseries = args.output_directory / "simulation_spiceEV.csv"
    if not args.save_results:
        args.save_results = args.output_directory / "simulation_spiceEV.json"
    if not args.save_soc:
        args.save_soc = args.output_directory / "simulation_soc_spiceEV.csv"

    # copy input files to output to ensure reproducibility
    copy_list = [args.config, args.electrified_stations, args.vehicle_types]

    # only copy cost params if they exist
    if args.cost_parameters_file is not None:
        copy_list.append(args.cost_parameters_file)
    for c_file in copy_list:
        shutil.copy(str(c_file), str(args.output_directory / Path(c_file).name))

    util.save_version(Path(args.output_directory / "program_version.txt"))

    # rename special options
    args.timing = args.eta

    if args.input_schedule is None:
        raise SystemExit("The following argument is required: input_schedule")

    simulate.simulate(args)
