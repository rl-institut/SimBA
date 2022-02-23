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
