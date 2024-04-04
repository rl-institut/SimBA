
.. _simba_features:

Features of SimBA
=================

.. _consumption_analysis:

Consumption analysis
--------------------

The consumption can be calculated in two ways: Either with a constant average specific consumption or with a trip specific consumption, where the consumption depends on the temperature, the incline, the level of loading and the average speed.

To use a constant consumption, this value can be defined in the :ref:`vehicle_types` as "mileage" in the unit of [kWh/km]. To use the trip specific consumption, a consumption table needs to be created with data for each vehicle type. The path to the consumption table is assigned to the variable "mileage" in :ref:`vehicle_types`. The consumption table should have the columns "vehicle_type", "level_of_loading", "incline" ,"mean_speed_kmh", "t_amb", "consumption_kwh_per_km" and its creation is at ths point of development not part of SimBA. The consumption for each trip is then calculated by interpolation of the datapoints provided in the table.

The level_of_loading describes the share between an empty vehicle (0) and a fully loaded vehicle (1) and can be handed over in two ways: Either directly as a column in :ref:`schedule` with specific values for each individual trip or as an input file "default_level_of_loading.csv" (path defined in :ref:`config`) containing values for each hour of the day. SimBA will first look into schedule and take this value and if it is not defined take the one from the .csv file. The temperature information is obtained in the same way as the level of loading and should be given in °C. In order to calculate the consumption based on the incline, an extra "all_stations.csv" file (path defined in :ref:`config`) has to be provided containing information on the elevation in m of each station. The incline is then calculated as the average incline between start and stop of the trip, which is an assumption that creates good results for vehicles with recuperation technology, as most electric vehicles have. The mean speed is calculated using the information provided in :ref:`schedule` about duration and length of each trip.


.. _vehicle_dispatch:

Vehicle Dispatch
----------------

To allocate the rotations to vehicles, vehicles of the needed type to fulfil the rotation are used. If no suitable vehicle is available, a new vehicle is created. A vehicle is defined as "available" if it is currently not serving another rotation, has the same depot and if it had enough time after return to the depot to be charged. This "minimum standing time" at the depot is calculated using the variable min_recharge_deps_oppb or min_recharge_deps_depb from the :ref`config` together with the respective battery capacity of the vehicle and assuming the maximum available power of the depot charging stations.

Charging simulation
-------------------

The charging simulation is carried out in the open source software `SpiceEV <https://github.com/rl-institut/spice_ev>`_, that is included in SimBA as a package. SimBA therefore uses SpiceEVs charging strategy "distributed", that allows to separate the charging strategy depending on the station type. The station types can be either depot charging station (deps) or opportunity charging station (opps). The charging strategy that is applied at depot and opportunity charging stations can be chosen in the config as "strategy_deps" respectively "strategy_opps", together with the strategy options. See the SpiceEV documentation for the description of all possible `strategies <https://spice-ev.readthedocs.io/en/latest/charging_strategies_incentives.html>`_ and  `options <https://spice-ev.readthedocs.io/en/latest/simulating_with_spiceev.html#strategy-options>`_.

.. _generate_report:

Generate report
---------------

The generation of the report is implemented as a mode, that can be activated with the keyword "report" in modes (:ref:`report`). The report generates most of the output files. The report can be called any number of times e.g. mode = ["sim", "report", "neg_depb_to_oppb", "report", "service_optimization", "report"]. For each report, a sub-folder is created in the output directory (as defined in the :ref:`config`)  named "report_[nr]" with the respective number.

The generation of the report can be modified using the flag "cost_calculation" in :ref:`config`. If this flag is set to true, each report will also generate the file "summary_vehicles_costs.csv".

Default outputs
###############

| **Grid Connector Overview (gc_overview.csv)**
| Contains information about charging stations, including their names, types, maximum power, maximum number of charging stations, total energy usage, and use factors for the least, second least, and third least utilized charging stations.

| **Grid Connector Time Series (gc_power_overview_timeseries.csv)**
| Time series of power flow in kW for every grid connector

| **Rotation SoC Data (rotation_socs.csv)**
| Time series of SoC for each rotation.

| **Vehicle SoC Data (vehicle_socs.csv)**
| Time series of SoC for each vehicle.

| **Rotation Summary (rotation_summary.csv)**
| Contains data related to the rotation of vehicles, including the start and end times of each rotation, the type and ID of the vehicle, the depot name, the lines the vehicle traveled, total energy consumption in kWh, distance traveled in m, and various charging-related metrics such as charging type and SoC at arrival, minimum SoC and if the rotation had negative SoC values.

| **Overview Plots (run_overview.pdf and run_overview.png)**
| Contains plots for SoCs for every vehicle, power at each charging station, batteries, external loads and feed-ins as well as price time series for each station.

| **Station Data Summary (simulation_station_xy.json)**
| Contains information about the simulation interval, grid connector, photovoltaics, charging strategy, average flexible power range per time window, total drawn energy from the grid, average duration of standing events, maximum drawn power, total energy fed into the grid, maximum stored energy in each battery, number of load cycles for stationary batteries and vehicles, and number of times vehicle SoC was below the desired SoC on departure.

| **Station Data Time Series (simulation_timeseries_station_xy.csv)**
| Contains station specific time series including price of electricity, grid supply, fixed loads, battery power, energy stored in battery, flex band boundaries, battery feed, charging station power use, occupied charging stations and charging stations in use as well as vehicles which are at the station.

| **Overview on costs and vehicles (summary_vehicles_costs.csv)**
| If cost_calculation is activated, this file contains the cost report as described below in :ref:`cost_calculation`.

| **Used modes (used_modes.txt)**
| This text file lists all modes executed until the report is created.

.. _cost_calculation:

Cost calculation
################
| **Cost calculation (summary_vehicles_costs.csv)**
| This is an optional output which calculates investment and maintenance costs of the infrastructure as well as energy costs in the scenario. The costs are calculated based on the price sheet, given as input in the :ref:`cost_params`.
| The energy costs and the grid connector costs are specific for each grid operator, as given by the :ref:`cost_params`.
| The following costs are calculated as both total and annual, depending on the lifetime of each component. See `SpiceEV documentation <https://spice-ev.readthedocs.io/en/latest/charging_strategies_incentives.html#incentive-scheme>`_ for the calculation of electricity costs.

* Investment
    * **Buses**: Costs for all buses used in the simulation. Costs include battery swaps, depending on the lifetime of both buses and batteries.
    * **Charging infrastructure**: Costs for all depot and opportunity charging stations, depending on the number of actually used charging stations at each grid connector.
    * **Grid connectors**: Costs for grid connectors and transformers, depending on the voltage level and the distance to the grid.
    * **Garages**: Costs for workstations and charging infrastructure at garages.
    * **Stationary storages**: Costs for stationary batteries at depot and opportunity stations, depending on its capacity [kWh].
    * **Feed-in**: Costs for feed-in (power generation plant such as PV modules) at depot and opportunity stations, depending on its capacity [kW].
* Maintenance
    * Depending on the lifetime of each component maintenance costs are calculated for buses, charging infrastructure, grid connectors,  stationary storages and feed-in.
* Electricity
    * **Power procurement**: Costs for the procurement of energy.
    * **Grid fees**: Costs for power and energy price, depending on the voltage level and the utilization time per year.
    * **Taxes**: Taxes like electricity taxes, depending on given taxes by price sheet.
    * **Feed-in remuneration**: Remuneration for electricity fed into the grid.

As result the following table is saved as CSV. The first two columns of the csv file are "parameter" and "unit" as descibed below. The third column "cumulated" returns the results for the whole scenario. In the consecutive columns, the results for each electrified station are displayed separately. Be aware that:

* Busses, that start or stop at more then one depot are displayed in each column of the respective electrified station, but only once in the column "cumulated", so the sum of all depots in not necessarily equal to the sum of all electrified stations. This is valid for both number of vehicles and investment costs.
* Costs are only calculated for :ref:`electrified_stations`. Vehicles for rotations starting at non electrified stations get added to the column "Non_elctrified_station".


+---------------------------------+----------+-----------------------------------------------------------------------+
|**parameter**                    | **unit** | **description**                                                       |
+=================================+==========+=======================================================================+
|c_vehicles                       | EUR      | Investment costs of all buses                                         |
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
|c_feed_in                        | EUR      | Investment costs of feed in (power generation plant)                  |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_invest                         | EUR      | Sum of all investment costs                                           |
+---------------------------------+----------+-----------------------------------------------------------------------+
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_vehicles_annual                | EUR/year | Annual investment costs of all buses                                  |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_gcs_annual                     | EUR/year | Annual investment costs of all grid connectors                        |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_cs_annual                      | EUR/year | Annual investment costs of all charging stations                      |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_garage_annual                  | EUR/year | Sum of annual investment costs of garages                             |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_stat_storage_annual            | EUR/year | Annual investment costs of all stationary storages                    |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_feed_in                        | EUR/year | Annual investment costs of all feed-ins                               |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_invest_annual                  | EUR/year | Sum of all annual investment costs                                    |
+---------------------------------+----------+-----------------------------------------------------------------------+
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_gc_annual                | EUR/year | Annual maintenance costs of grid connectors                           |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_infrastructure_annual    | EUR/year | Annual maintenance costs of charging stations and stationary storages |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_vehicles_annual          | EUR/year | Annual maintenance costs of buses                                     |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_stat_storage_annual      | EUR/year | Annual maintenance costs of stationary storages                       |
+---------------------------------+----------+-----------------------------------------------------------------------+
|c_maint_feed_in_annual           | EUR/year | Annual maintenance costs of feed-in                                   |
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

Optimization
------------

There are several options for optimizations that are implemented as :ref:`sim_modes`. These options are currently:

* :ref:`neg_depb_to_oppb`
* :ref:`neg_oppb_to_depb`
* :ref:`Service Optimization`
* :ref:`Station Optimization`

.. _consistency_check:

Consistency check
-----------------

SimBA makes certain assumption, that have to be valid to trust the results. these are:

* The trips inside a rotation are chronologically sorted
* The trip time is not negative, so the arrival of the trip is later or equal to its departure.
* The break time between trips is not negative, so the departure of the consecutive trip is later or equal to the arrival of the preceding trip.
* Each rotation has a defined and fixed depot, so the rotation starts and ends at the same station
* Every trip within a rotation starts where the previous trip ended

In order to test these assumptions, the flag "check_rotation_consistency" can be activated in the :ref:`config`, which will result in the display of cases where assumptions are broken in the console and in the log file. Additionally, the inconsistent rotations can be filtered out of the simulation by setting the "skip_inconsistent_rotations" flag to true.


.. _rotation_filter:

Rotation filter
---------------

Before all rotations specified in the :ref:`schedule` are simulated, there is the option to filter only the ones relevant to for the actual analysis. This is activated by setting the "rotation_filter_variable" in the :ref:`config` to either "include" to consider only certain rotations from the schedule, or to "exclude" to exclude certain rotations from the analysis. The list of rotations for both options is specified as "rotation_filter" in the Path paragraph of the :ref:`config`.

Logging
-------

SimBA uses the "logging" package for logging. All logging messages are both displayed in the Terminal and written to a .log file. The filepath and the loglevel can be defined in the :ref:`config`. Four log levels are available in the following order: DEBUG, INFO, WARN and ERROR. INFO includes INFO, WARN and ERROR but excludes DEBUG.