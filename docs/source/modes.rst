# Without creating links like in the line below, subpages go missing from the sidebar

.. _sim_modes:

Modes of the eBus-Toolbox
=========================

The eBus-Toolbox assists the user in analysing and optimising electrified bus fleets and schedules. Besides a simple simulation run, several
different modes support the user in finding optimal solutions for their eBus-System. Supported Modes are

* simple simulation
* negative depot to opportunity charger
* negative opportunity to depot charger
* service optimization
* station optimization
* report

While the default mode is the simple simulation modes can be chained together to achieve the desired results. The chain of modes is defined in the config file (default: ebus_toolbox.cfg) under the keyword *mode*:
::
    mode = ["sim", "report"]

This results in a simple simulation with a following report. To run a simulation the ebus-toolbox creates a schedule which contains all information about how the bus system is supposed to run. Some modes are allowed to mutate the schedule in some way, which makes chaining of modes especially useful. An extended simple use case would be
    mode = ["sim", "report" ,"neg_depb_to_oppb", "report]




simple simulation
-----------------
The simple simulation case is the default mode. Its usage is explained in :ref:`Getting Started`. The simulation takes the scenario as is. No parameters will be adjusted, optimized or changed in any way. The report will include information about the expected socs, power loads at the charging stations or depots, default plots for the scenario and other useful data


negative depot to opportunity charger
-------------------------------------
| This mode is the first kind of optimization provided by the eBus-Toolbox. It takes a scenario and uses the provided depot charger vehicle data to check if the scenario can run with only depot charger vehicles without socs falling below 0. Vehicles with a soc below 0 are changed to opportunity chargers with the provided opportunity charger vehicle data and the simulation is run again.
| NOTE: Charging types are only switched by the eBus-Toolbox if the corresponding vehicle type exists.
| TBC



negative opportunity to depot charger
-------------------------------------
| This mode is the similar to the previous optimization. It takes a scenario and uses the provided opportunity charger vehicle data to check if the scenario can run with only opportunity charger vehicles without socs falling below 0. Vehicles with a soc below 0 are changed to opportunity chargers with the provided depot charger vehicle data and the simulation is run again.
| NOTE: Charging types are only switched by the eBus-Toolbox if the corresponding vehicle type exists.
| TBC


service optimization
--------------------
This mode optimizes a scenario by creating sub-scenarios, so the given sub-scenario runs without socs falling below 0. The sub-scenario has the same input has the scenario but only consists of a sub-sets of rotations. These sub-sets are are optimized to be as big as possible. Previous negative rotations can become positive since effects like blocked charging points by other rotations can be reduced.

station optimization
--------------------
.. _optimization_loop:
.. figure:: https://user-images.githubusercontent.com/104760879/217225177-66201146-d31a-4127-9ca0-4d6e6e5a3cc4.png
    :width: 600
    :alt: optimization_loop

    Caption

:numref:`optimization_loop` shows xxxxxxxxxxxxxxxxxxx

.. figure:: https://user-images.githubusercontent.com/104760879/217225588-abfad83d-9d2a-463a-8597-584e29f5f885.png
    :width: 600
    :alt: below_0_soc_event

    Caption


chained modes
-------------
Lorem ipsum ...