.. image:: https://user-images.githubusercontent.com/104760879/217226792-4297d3c8-8a7c-45ad-894f-5efd03031f49.png
    :alt: ebus_toolbox_logo

Getting Started
===============

The eBus-Toolbox assists the user in analysing and optimising electrified bus fleets and schedules.

.. Without creating links like in the line below, subpages go missing from the sidebar

.. _installation_label:

Installation
------------
To try it out, first clone `this repository <https://github.com/rl-institut/eBus-Toolbox>`_ and then install the required packages to your current environment by running

`pip install -r requirements.txt`

Now you can start the eBus Toolbox module with all configurations stored at `data/configs/ebus_toolbox.cfg` via the command

``python -m ebus_toolbox --config data/configs/ebus_toolbox.cfg``

The repo provides an example for each necessary input file, so the example case can be executed without the need for the user to provide any data themselves.

To run the eBus Toolbox with your own `schedule.csv` (see details [below](#input-data)) file and default configurations run

`python -m ebus_toolbox --input_schedule path/to/schedule.csv`

Default configurations are detailed at `data/configs/ebus_toolbox.cfg`.


General Concept
---------------
At the current stage, only a single functionality is implemented, which is processing a bus schedule stored in a specific CSV format (see `data/examples/trips_examples.csv`) and run it through a module called SpiceEV for an in-depth SOC analysis.

.. figure:: https://user-images.githubusercontent.com/104760879/217225545-5e6858c1-d056-4519-beea-6274d06533c7.png
    :alt:  ebus_toolbox_modules
    :width: 800

    Modules of the eBus-Toolbox

.. figure:: https://user-images.githubusercontent.com/104760879/217226800-647956c5-9d63-4988-a710-f3326a8304d5.png
    :alt:  ebus_toolbox_default_plot
    :width: 800

    Default output plot for a single simulation.

More text
