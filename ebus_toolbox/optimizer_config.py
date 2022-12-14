from datetime import datetime

class OptimizerConfig():
    def __init__(self):
        pass

def read_config(config_path):
    import configparser
    import json
    config_path = config_path
    config_parser = configparser.ConfigParser()
    config_parser.sections()
    config_parser.read(config_path)

    conf = OptimizerConfig()
    default = config_parser["DEFAULT"]
    conf.debug_level = int(default.get("debug_level", 0))
    sce = config_parser["SCENARIO"]
    conf.exclusion_rots = set(json.loads(sce.get("exclusion_rots", "[]")))
    conf.exclusion_stations = set(json.loads(sce.get("exclusion_stations", "[]")))
    conf.inclusion_stations = set(json.loads(sce.get("inclusion_stations", "[]")))
    conf.standard_opp_station = dict(json.loads(sce.get("standard_opp_station", "{}")))

    pi = config_parser["PICKLE"]

    conf.schedule = pi.get("schedule", "")
    conf.scenario = pi.get("scenario", "")
    conf.args = pi.get("args", "")

    vehicle = config_parser["VEHICLE"]
    conf.charge_eff = float(vehicle.get("charge_eff", 0.95))
    conf.battery_capacity = float(vehicle.get("battery_capacity", 0))
    if conf.battery_capacity==0:
        conf.battery_capacity=None
    conf.charging_curve = json.loads(vehicle.get("charging_curve", "[]"))
    if conf.charging_curve==[]:
        conf.charging_curve=None
    conf.charging_power = float(vehicle.get("charging_power", 0))
    if conf.charging_power==0:
        conf.charging_power=None
    conf.min_soc = float(vehicle.get("min_soc", 0.0))



    optimizer = config_parser["OPTIMIZER"]
    conf.solver = optimizer.get("solver", "spiceev")
    conf.rebase_scenario = optimizer.getboolean("rebase_scenario", True)
    conf.pickle_rebased = optimizer.getboolean("pickle_rebased", False)
    conf.pickle_rebased_name = optimizer.get("pickle_rebased_name",
            "rebased_" +str(datetime.now().strftime("%Y-%m-%d-%H-%M-%S")))
    conf.opt_type = optimizer.get("opt_type", "greedy")
    conf.remove_impossible_rots = optimizer.getboolean("remove_impossible_rots", False)
    conf.node_choice = optimizer.get("node_choice", "step-by-step")
    conf.max_brute_loop = int(optimizer.get("max_brute_loop", 200))
    conf.run_only_neg = optimizer.getboolean("run_only_neg", False)
    conf.estimation_threshold = float(optimizer.get("estimation_threshold", 0.8))
    conf.output_path= optimizer.get("output_path")
    conf.check_for_must_stations = optimizer.getboolean("check_for_must_stations", True)

    special = config_parser["SPECIAL"]
    conf.decision_tree_path = special.get("decision_tree_path", None)
    if conf.decision_tree_path in ["", '""', "''"] :
        conf.decision_tree_path=None
    conf.save_decision_tree = special.getboolean("save_decision_tree", False)
    conf.reduce_rots = special.getboolean("reduce_rots", False)
    conf.rots = json.loads(special.get("rots", []))

    return conf