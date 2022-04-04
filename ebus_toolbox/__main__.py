import argparse
from ebus_toolbox import simulate, util

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ebus Toolbox - \
        Simulation Program for Ebus Fleets.')
    parser.add_argument('input_schedule', nargs='?', help='Set the scenario JSON file')
    parser.add_argument('input', nargs='?', help='output file name (example.json)')
    parser.add_argument('--preferred_charging_type', '-pct', default='depot',
                        choices=['depot', 'opp'], help="Preferred charging type. Choose one\
                        from <depb> and <oppb>. opp stands for opportunity.")
    parser.add_argument('--vehicle-types', default=None,
                        help='location of vehicle type definitions')
    parser.add_argument('--min-standing-time_depot', default=6,
                        help='Minimum standing time in hours for all buses in the depot before \
                        they become available for dispatch again.')
    parser.add_argument('--days', metavar='N', type=int, default=None,
                        help='set duration of scenario as number of days')
    parser.add_argument('--interval', metavar='MIN', type=int, default=15,
                        help='set number of minutes for each timestep (Δt)')
    parser.add_argument('--desired-soc', metavar='SOC', type=float, default=0.8,
                        help='set minimum desired SOC (0 - 1) for each charging process')
    parser.add_argument('--gc_power_opps', metavar='POPP', type=float, default=1500,
                        help='max power of grid connector at opp stations')
    parser.add_argument('--gc_power_deps', metavar='PDEP', type=float, default=1500,
                        help='max power of grid connector at depot stations')
    parser.add_argument('--cs_power_opps', metavar='CSPOPP', type=float, default=150,
                        help='max power of charging station at opp stations')
    parser.add_argument('--cs_power_deps_depb', metavar='CSPDEPDEP', type=float, default=150,
                        help='max power of charging station at depot stations for depot busses')
    parser.add_argument('--cs_power_deps_oppb', metavar='CSPDEPOPP', type=float, default=150,
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
    parser.add_argument('--start_date', default='2018-01-02',
                        help='Provide start date of simulation in format YYYY-MM-DD.E.g. '
                             '2018-01-31')
    parser.add_argument('--electrified_stations', help='include electrified_stations json',
                        default='examples/electrified_stations.json')
    parser.add_argument('--vehicle_types', help='include vehicle_types json',
                        default='examples/vehicle_types.json')
    parser.add_argument('--min_charging_time_opp', help='define minimum time of charging',
                        default=2)
    parser.add_argument('--visual', '-v', action='store_true', help='Show plots of the results')
    parser.add_argument('--desired_soc', default=1, help='desired_soc of vehicles')
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
    args = parser.parse_args()

    util.set_options_from_config(args, check=True, verbose=False)

    simulate.simulate(args)
