# eBus Toolbox

The eBus Toolbox assists the user in analysing and optimising electrified bus fleets and schedules.

### Usage

At the current stage, only a single functionality is implemented, which is processing a bus schedule stored in a specific CSV format (see `data/examples/trips_examples.csv`) and run it through a module called SpiceEV for an indepth SOC analysis.

To try it out, run the following command from the root directory of this repository after installing the dependencies from `requirements.txt`. 

`python -m ebus_toolbox --config data/configs/ebus_toolbox.cfg`

The repo provides an example for each necessary input file, so the example case can be executed without the need for the user to provide any data themselves.