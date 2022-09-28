import json


def set_options_from_config(args, check=False, verbose=True):
    """Read options from config file, update given args, try to parse options
    , ignore comment lines (begin with #)

    :param args: input arguments
    :type args: argparse.Namespace
    :param check: raise ValueError on unknown options
    :type check: bool
    :param verbose: gives final overview of arguments
    :type bool

    :raises ValueError: Raised if unknown options are given.
    """

    if "config" in args and args.config is not None:
        # read options from config file
        with open(args.config, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    # comment
                    continue
                if len(line) == 0:
                    # empty line
                    continue
                k, v = line.split('=')
                k = k.strip()
                v = v.strip()
                try:
                    # option may be special: number, array, etc.
                    v = json.loads(v)
                except ValueError:
                    # or not
                    pass
                # known option?
                if (k not in args) and check:
                    raise ValueError("Unknown option {}".format(k))
                # set option
                vars(args)[k] = v
        # Give overview of options
        if verbose:
            print("Options: {}".format(vars(args)))


def get_buffer_time(trip, default=0):
    """ Get buffer time at arrival station of a trip.
        Buffer_time is an abstraction of delays like docking procedures and
        is added to the planned arrival time

    :param trip: The of buffer time of this trips arrival is returned.
    :type trip: ebus_toolbox.Trip
    :param default: Default buffer time if no station specific buffer time is given. [minutes]
    :type default: dict, numeric
    :return: Buffer time
    :rtype: numeric

    NOTE: Buffertime dictionaries map hours of the day to a buffer time.
    Keys are ranges of hours and corresponding values provide buffer time in
    minutes for that time range.
    An entry with key "else" is a must if not all
    hours of the day are covered.
    E.g.
        buffer_time = {
            "10-22": 2,
            "22-6": 3,
            "else": 1
        }
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
