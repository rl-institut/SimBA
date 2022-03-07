import argparse
from ebus_toolbox import simulate, util

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ebus Toolbox - \
        Simulation Program for Ebus Fleets.')
    parser.add_argument('input', nargs='?', help='Set the scenario JSON file')
    parser.add_argument('--preferred_charging_type', '-pct', default='depot',
                        choices=['depot', 'opp'], help="Preferred charging type. Choose one\
                        from <depot> and <opp>. opp stands for opportunity.")
    parser.add_argument('--vehicle-types', default=None,
                        help='location of vehicle type definitions')
    parser.add_argument('--min-standing-time', default=6,
                        help='Minimum standing time in hours for all buses in the depot before \
                        they become available for dispatch again.')
    parser.add_argument('--visual', '-v', action='store_true', help='Show plots of the results')
    parser.add_argument('--desired_soc', default = 1, help='desired_soc of vehicles')
    parser.add_argument('--eta', action='store_true',
                        help='Show estimated time to finish simulation after each step, \
                        instead of progress bar. Not recommended for fast computations.')
    parser.add_argument('--save-timeseries', help='Write timesteps to file')
    parser.add_argument('--save-results', help='Write general info to file')
    parser.add_argument('--config', help='Use config file to set arguments')
    args = parser.parse_args()

    util.set_options_from_config(args, check=True, verbose=False)

    simulate.simulate(args)
