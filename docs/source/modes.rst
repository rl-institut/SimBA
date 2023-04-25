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

chained modes
-------------
While the default mode of the ebus-toolbox is the simple simulation together with a report, modes can be chained together differently to achieve the desired results. The chain of modes is defined in the config file (default: ebus_toolbox.cfg) under the keyword *mode*:
::
    mode = ["sim", "report"]

This results in a simple simulation with a following report. To run a simulation the ebus-toolbox creates a schedule which contains all information about how the bus system is supposed to run. Some modes are allowed to mutate the schedule in some way, which makes chaining of modes especially useful. Their output describes the simulation outcome of this mutated schedule. An extended simple use case would be
::
    mode = ["sim", "report" ,"neg_depb_to_oppb", "report]

Where the scenario is run as is, a report is generated, the schedule is changed and simulated again and a second report is generated. The description what the modes do
can be found below.

simple simulation
-----------------
The simple simulation case is the default mode. Its usage is explained in :ref:`Getting Started`. Every chain of modes starts with a simple simulation, even if it is not explicitly called. The simulation takes the scenario as is. No parameters will be adjusted, optimized or changed in any way. The charging type for each rotation is read from the trips.csv if this data is included. If not the *preferred_charging_type* from the config file is used, as long as the provided vehicles data provides the preferred_charging_type for the specified vehicle type.


negative depot to opportunity charger
-------------------------------------

This mode is the first kind of optimization provided by the eBus-Toolbox and is called by
::
    mode = ["neg_depb_to_oppb"]
It takes the results of the previous simulation, attains all rotations which had a negative soc, and changes their vehicle type to depot chargers.
| NOTE: Charging types are only switched by the eBus-Toolbox if the corresponding vehicle type as depot charger exists in the provided vehicles_data.json.



negative opportunity to depot charger
-------------------------------------
|This mode is analogous to *neg_depb_to_oppb*
This mode is the second kind of optimization provided by the eBus-Toolbox and is called by
::
    mode = ["neg_oppb_to_depb"]
It takes the results of the previous simulation, attains all rotations which had a negative soc, and changes their vehicle type to opportunity chargers.
| NOTE: Charging types are only switched by the eBus-Toolbox if the corresponding vehicle type as opportunity charger exists in the provided vehicles_data.json.


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


report
-----------------
The  The report will include information about the expected socs, power loads at the charging stations or depots, default plots for the scenario and other useful data.


chained modes
-------------
While the default mode of the ebus-toolbox is the simple simulation, modes can be chained together to achieve the desired results. The chain of modes is defined in the config file (default: ebus_toolbox.cfg) under the keyword *mode*:
::
    mode = ["sim", "report"]

This results in a simple simulation with a following report. To run a simulation the ebus-toolbox creates a schedule which contains all information about how the bus system is supposed to run. Some modes are allowed to mutate the schedule in some way, which makes chaining of modes especially useful. An extended simple use case would be
::
    mode = ["sim", "report" ,"neg_depb_to_oppb", "report]

Where the scenario is run as is, a report is generated, the schedule is changed and simulated again and a second report is generated. The description what the modes do
can be found below.
=======
Lorem ipsum ...


Report
-------------
Cost calculation
~~~~~~~~~~~~~

This mode calculates investment and maintenance costs of the infrastructure as well as energy costs in the scenario. The costs are calculated based on the price sheet, given as input in the ``costs_params.json``.
The following costs are calculated as both total and annual, depending on the lifetime of each component. See `SpiceEV <https://spice-ev.readthedocs.io/en/latest/charging_strategies_incentives.html#incentive-scheme>`_ for the calculation of electricity costs.

* Investment
    * **Busses**: Costs for all busses used in the simulation. Costs include battery swaps, depending on the lifetime of both busses and batteries.
    * **Charging infrastructure**: Costs for all depot and opportunity charging stations, depending on the number of actually used charging stations at each grid connector.
    * **Grid connectors**: Costs for grid connectors and transformers, depending on the voltage level and the distance to the grid.
    * **Garages**: Costs for workstations and charging infrastructure at garages.
    * **Stationary storages**: Costs for stationary batteries at depot and opportunity stations, depending on its capacity.
* Maintenance
    * Depending on the lifetime of each component maintenance costs are calculated for busses, charging infrastructure, grid connectors and stationary storages.
* Electricity
    * **Power procurement**: Costs for the procurement of energy.
    * **Grid fees**: Costs for power and energy price, depending on the voltage level and the utilization time per year.
    * **Taxes**: Taxes like electricity taxes, depending on given taxes by price sheet.
    * **Feed-in remuneration**: Remuneration for electricity fed into the grid.

As result the following table is saved as CSV:

+---------------------------------+----------+-----------------------------------------------------------------------+
|**parameter**                    | **unit** | **description**                                                       |
+=================================+==========+=======================================================================+
|c_vehicles                       | EUR      | Investment costs of all busses                                        |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_gcs                            | EUR      | Investment costs of all grid connectors                               |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_cs                             | EUR      | Investment costs of all charging stations                             |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_garage_cs                      | EUR      | Investment costs of charging stations at garages                      |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_garage                         | EUR      | Investment costs of garages itself                                    |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_garage_workstations            | EUR      | Investment costs of working stations at garages                       |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_stat_storage                   | EUR      | Investment costs of stationary storages                               |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_invest                         | EUR      | Sum of all investment costs                                           |
+---------------------------------+----------+-----------------------------------------------------------------------+
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_vehicles_annual                | EUR/year | Annual investment costs of all busses                                 |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_gcs_annual                     | EUR/year | Annual investment costs of all grid connectors                        |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_cs_annual                      | EUR/year | Annual investment costs of all charging stations                      |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_garage_annual                  | EUR/year | Sum of annual investment costs of garages                             |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_stat_storage_annual            | EUR/year | Annual investment costs of all stationary storages                    |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_invest_annual                  | EUR/year | Sum of all annual investment costs                                    |
+---------------------------------+----------+-----------------------------------------------------------------------+
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_gc_annual                | EUR/year | Annual maintenance costs of grid connectors                           |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_infrastructure_annual    | EUR/year | Annual maintenance costs of charging stations and stationary storages |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_vehicles_annual          | EUR/year | Annual maintenance costs of busses                                    |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_stat_storage_annual      | EUR/year | Annual maintenance costs of stationary storages                       |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_annual                   | EUR/year | Sum of annual maintenance costs                                       |
+---------------------------------+----------+-----------------------------------------------------------------------+
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_el_procurement_annual          | EUR/year | Annual costs of power procurement                                     |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_el_power_price_annual          | EUR/year | Annual grid fee for highest load peak                                 |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_el_energy_price_annual         | EUR/year | Annual grid fee for drawn energy                                      |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_el_taxes_annual                | EUR/year | Annual costs for all electricity related taxes                        |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_el_feed_in_remuneration_annual | EUR/year | Annual feed-in remuneration                                           |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_el_annual                      | EUR/year | Sum of all annual electricity costs                                   |
+---------------------------------+----------+-----------------------------------------------------------------------+












