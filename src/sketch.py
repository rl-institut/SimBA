# imports
from spice_ev import args_from_config
import optimizer

def checkInput():
    pass

def add_charging_types():
    pass

def calc_trip_rotation_consumption():
    pass

def assign_vehicles():
    pass

def filter_rotations():
    pass


def evaluate(args):
    # read CSV input file
    # check for correct data format / standard
    trips = read(csv_file)
    if not checkInput(trips):
        # exit sim
        pass

    gc_info[-1]["target"] = event.target if event.target is not None else gc_info[-1]["targejsdhajkshdksdsdashdksajhdkhaskdhkashkdhaskhdkhashdklast"]
    trips = filter_trips(args.rotations_to_be_ignored)
    opt = Optimizer(args.opt_config) if args.opt_config else None

    # may be updated by optimizer later
    create_configs_spiceev()

    while(True):
        # construct szenario and simulate in spice ev
        # stop if no optimization wanted (simple yes or no response)
        # stop if optimi
        args_from_config()
        vehicle_types = read(vt.json)

        add_charging_types()
        calc_trip_rotation_consumption()
        assign_vehicles()

        # RUN SPICE EV
   
        if opt.optimize(inputs, outputs) == 'satisfied':
            break

    generate_report()


if __name__ == 'main':
    # parse args from cmd or file
    evaluate(args)