""" Optimization that tries minimizing the amount of electrified stations to achieve full
electrification.

"""

from datetime import datetime
import json
from pathlib import Path
import logging
import shutil
import ebus_toolbox.station_optimizer
from ebus_toolbox.station_optimizer import util
config = util.OptimizerConfig()


def setup_logger(conf):
    """ Setup file and stream logging by config and args arguments
    :param conf: configuration object
    :return: logger
    :rtype: Logger
    """
    this_logger = logging.getLogger(__name__)
    this_logger.setLevel(logging.DEBUG)

    # logging to one file which keeps track of optimization over many runs
    file_handler_all_opts = logging.FileHandler('optimizer.log')
    file_handler_all_opts.setLevel(conf.debug_level)

    # and logging to a file which is put in the folder with the other optimizer results
    file_handler_this_opt = logging.FileHandler(Path(conf.optimizer_output_dir) /
                                                Path('optimizer.log'))
    file_handler_this_opt.setLevel(conf.debug_level)

    formatter = logging.Formatter('%(asctime)s:%(message)s',
                                  "%m%d %H%M%S")

    file_handler_all_opts.setFormatter(formatter)
    file_handler_this_opt.setFormatter(formatter)

    formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(conf.console_level)
    this_logger.addHandler(file_handler_this_opt)
    this_logger.addHandler(file_handler_all_opts)
    this_logger.addHandler(stream_handler)
    return this_logger


def main():
    """ main call"""
    util.print_time()
    config_path = "./data/examples/default_optimizer.cfg"
    conf = util.read_config(config_path)
    opt_sched, opt_scen = run_optimization(conf)
    util.print_time()
    import pickle
    with open("schedule_opt.pickle", "wb") as f:
        pickle.dump(opt_sched, f)
    with open("scenario_opt.pickle", "wb") as f:
        pickle.dump(opt_scen, f)


def prepare_filesystem(this_args, conf):
    """ Prepare files and folders in the optimization results folder
    :param conf: configuration object
    :param this_args: Namespace object of arguments for ebus toolbox
    """
    now = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    this_args.output_directory = Path(this_args.output_directory)
    conf.optimizer_output_dir = Path(this_args.output_directory) / (now + "_optimizer")
    # create sub folder for specific sim results with timestamp.
    # if folder doesnt exists, create folder.
    conf.optimizer_output_dir.mkdir(parents=True, exist_ok=True)

    # copy paste the config file
    copy_file = Path(conf.path)
    destination = conf.optimizer_output_dir / Path("optimizer_config.cfg")
    shutil.copy(copy_file, destination)


def run_optimization(conf, sched=None, scen=None, args=None):
    """    Optimizes scenario by adding electrified stations sparingly
    until scenario has no below 0 soc events.

    The optimizer config file is used for setting the optimization up.

    :param conf: Configuration object of optimization
    :type conf: OptimizerConfig

    :param sched: Simulation schedule containing buses, rotations etc.
    :type sched: ebus_toolbox.Schedule

    :param scen: Simulation scenario containing simulation results
                     including the SoC of all vehicles over time
    :type scen: spice_ev.Scenario
    :param args: Simulation arguments for manipulation or generated outputs
    :type args: object

    :return: (Schedule,Scenario) optimized schedule and Scenario
    :rtype: tuple(ebus_toolbox.Schedule, spice_ev.Scenario)
    """

    # load pickle files only in the case they are not given as arguments
    if sched is None or scen is None or args is None:
        # if no schedule was given as argument, make sure no scenario
        # and args input was given as well.
        assert sched == scen == args is None
        sched, scen, args = util.toolbox_from_pickle(conf.schedule,
                                                     conf.scenario, conf.args)

    # prepare Filesystem with folders and paths and copy config
    prepare_filesystem(args, conf)

    if args.save_soc:
        args.save_soc = args.output_directory / "simulation_soc_spiceEV.csv"
    new_ele_stations_path = conf.optimizer_output_dir / Path("optimized_stations" + ".json")

    logger = setup_logger(conf)

    if args.desired_soc_deps != 1 and conf.solver == "quick":
        logger.error("Fast calc is not yet optimized for desired socs unequal to 1")

    optimizer = ebus_toolbox.station_optimizer.StationOptimizer(sched, scen, args, conf, logger)

    # set battery and charging curves through config file if wished for
    optimizer.set_battery_and_charging_curves()

    # filter out depot chargers if option is set
    if conf.run_only_oppb:
        sched.rotations = {r: sched.rotations[r] for r in sched.rotations
                           if "oppb" in sched.rotations[r].vehicle_id}
        assert len(sched.rotations) > 0, "Removing depot chargers led to a schedule without any " \
                                         "rotations."

    # rebasing the scenario meaning simulating it again with the given conditions of
    # included and excluded stations and rotations and maybe changed battery sizes
    if conf.rebase_scenario:
        must_include_set, ele_stations = optimizer.rebase_spice_ev()
    else:
        # no new spice ev calculation will take place but some variables need to be adjusted.
        must_include_set, ele_stations = optimizer.rebase_simple()

    # create charging dicts which contain soc over time, which is numerically calculated
    optimizer.create_charging_curves()

    # remove none values from socs in the vehicle_socs
    optimizer.remove_none_socs()

    # check if rotations cant be operated electric, even with all stations electrified.
    if conf.remove_impossible_rotations:
        neg_rots = optimizer.get_negative_rotations_all_electrified()
        optimizer.config.exclusion_rots.update(neg_rots)
        optimizer.schedule.rotations = {r: optimizer.schedule.rotations[r]
                                        for r in optimizer.schedule.rotations if
                                        r not in optimizer.config.exclusion_rots}
        logger.warning("%s negative rotations %s were removed from schedule",
                       len(neg_rots), neg_rots)
        assert len(optimizer.schedule.rotations) > 0, "Schedule cant be optimized, since" \
                                                      "rotations cant be electrified."

    # if the whole network cant be fully electrified if even just a single station is not
    # electrified, this station must be included in a fully electrified network
    # this can make solving networks much simpler. Some information in the electrification path
    # gets lost though
    if conf.check_for_must_stations:
        must_stations = optimizer.get_must_stations_and_rebase(relative_soc=False)
        logger.warning("%s must stations %s", len(must_stations), must_stations)

    logger.debug("Starting greedy station optimization")
    print("Starting greedy station optimization")
    ele_stations, ele_station_set = optimizer.loop()
    ele_station_set = ele_station_set.union(must_include_set)
    logger.debug("%s electrified stations : %s", len(ele_station_set), ele_station_set)
    logger.debug("%s total stations", len(ele_stations))
    logger.debug("These rotations could not be electrified: %s", optimizer.could_not_be_electrified)

    vehicle_socs = optimizer.timeseries_calc()

    new_events = optimizer.get_low_soc_events(soc_data=vehicle_socs)

    if len(new_events) > 0:
        logger.debug("Still not electrified with abs. soc with fast calc")
        for event in new_events:
            logger.debug(event.rotation.id)
    with open(new_ele_stations_path, "w", encoding="utf-8") as file:
        json.dump(ele_stations, file, ensure_ascii=False, indent=2)

    logger.debug("Spice EV is calculating optimized case as a complete scenario")
    _, __ = optimizer.preprocessing_scenario(
        electrified_stations=ele_stations, run_only_neg=False)

    logger.warning("Still negative rotations: %s", optimizer.schedule.
                   get_negative_rotations(optimizer.scenario))
    print("Station optimization finished after ", end="")
    util.print_time()

    return optimizer.schedule, optimizer.scenario


if __name__ == "__main__":
    main()
