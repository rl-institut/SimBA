Getting Started
===============

The eBus-Toolbox assists the user in analysing and optimising electrified bus fleets and schedules.

.. _installation:

Installation
------------

At the current stage, only a single functionality is implemented, which is processing a bus schedule stored in a specific CSV format (see `data/examples/trips_examples.csv`) and run it through a module called SpiceEV for an in-depth SOC analysis.

To try it out, first clone this repository and then install the required packages to your current environment by running

`pip install -r requirements.txt`

Now you can start the eBus Toolbox module with all configurations stored at `data/configs/ebus_toolbox.cfg` via the command

``python -m ebus_toolbox --config data/configs/ebus_toolbox.cfg``

The repo provides an example for each necessary input file, so the example case can be executed without the need for the user to provide any data themselves.

To run the eBus Toolbox with your own `schedule.csv` (see details [below](#input-data)) file and default configurations run

`python -m ebus_toolbox --input_schedule path/to/schedule.csv`

Default configurations are detailed at `data/configs/ebus_toolbox.cfg`.


General Concept
---------------
