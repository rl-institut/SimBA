..
    # Without creating links like in the line below, subpages go missing from the sidebar

.. _sim_modes:

Modes of the eBus-Toolbox
=========================

The eBus-Toolbox assists the user in analyzing and optimising electrified bus fleets and schedules. Besides a simple simulation run, several
different modes support the user in finding optimal solutions for their eBus-System. Supported Modes are

* simple simulation
* negative depot to opportunity charger
* negative opportunity to depot charger
* service optimization
* station optimization
* report

Chained Modes
-------------
While the default mode of the ebus-toolbox is the simple simulation together with a report, modes can be chained together differently to achieve the desired results. The chain of modes is defined in the config file (default: ebus_toolbox.cfg) under the keyword *mode*:

::

        mode = ["sim", "report"]

This results in a simple simulation with a following report. To run a simulation the ebus-toolbox creates a schedule which contains all information about how the bus system is supposed to run. Some modes are allowed to mutate the schedule in some way, which makes chaining of modes especially useful. Their output describes the simulation outcome of this mutated schedule. An extended simple use case would be:

::

    mode = ["sim", "report" ,"neg_depb_to_oppb", "report]

Where the scenario is run as is, a report is generated, the schedule is changed and simulated again and a second report is generated. The description what the modes do
can be found below.

Simple Simulation
-----------------
The simple simulation case is the default mode. Its usage is explained in :ref:`Getting Started`. Every chain of modes starts with a simple simulation, even if it is not explicitly listed in the modes. The simulation takes the scenario as is. No parameters will be adjusted, optimized or changed in any way. The charging type for each vehicle is read from the rotation information from the trips.csv if this data is included. If the data is not included *preferred_charging_type* from the config file is used, as long as the provided vehicles data provides the preferred_charging_type for the specified vehicle type.


Negative Depot to Opportunity Charger
-------------------------------------
This mode is the first kind of optimization provided by the eBus-Toolbox and is called by:

::

    mode = ["neg_depb_to_oppb"]

It takes the results of the previous simulation, attains all rotations which had a negative soc, and changes their vehicle type to depot chargers.

.. note:: Charging types are only switched by the eBus-Toolbox if the corresponding vehicle type as depot charger exists in the provided vehicles_data.json.



Negative Opportunity to Depot Charger
-------------------------------------
This mode is analogous to *neg_depb_to_oppb*.
This mode is the second kind of optimization provided by the eBus-Toolbox and is called by

::

    mode = ["neg_oppb_to_depb"]

It takes the results of the previous simulation, attains all rotations which had a negative soc, and changes their vehicle type to opportunity chargers.

.. note:: Charging types are only switched by the eBus-Toolbox if the corresponding vehicle type as opportunity charger exists in the provided vehicles_data.json.

Service Optimization
--------------------
This mode optimizes a scenario by creating sub-scenarios, so the given sub-scenario runs without socs falling below 0. The sub-scenario has the same input has the scenario but only consists of a sub-sets of rotations. These sub-sets are are optimized to be as big as possible. Previous negative rotations can become positive since effects like blocked charging points by other rotations can be reduced.

Station Optimization
--------------------
This mode optimizes a scenario by electrifying as few opportunity stations as possible using a greedy approach. The network with no opportunity charging station is first analyzed to find rotations which fail at the current stage and to estimate the potential of electrifying each station by its own. *Step-by-step* new opportunity stations are electrified until full electrification is reached. The optimization assumes that at every electrified station unlimited charging points exist, i.e. the number of simultaneously charging buses is not limited. In between each electrification a simulation is run and the network is analyzed again. The first run called the **base optimization** leads to a scenario which often times is better than extensively optimizing the scenario by hand. Since a greedy approach can not guarantee a global optimum a second extensive optimization can be chained to this base optimization. This *deep* optimization can make use of a *step-by-step* decision tree expansion which evaluates new combinations of electrified stations starting with the most promising combinations **OR** use a *brute* force approach trying to reduce the amount of electrified stations by one in comparison to the base optimization. The step-by-step process of the optimization follows :numref:`optimization_loop`

.. _optimization_loop:
.. figure:: https://user-images.githubusercontent.com/104760879/217225177-66201146-d31a-4127-9ca0-4d6e6e5a3cc4.png
    :width: 600
    :alt: optimization_loop

    Steps of the optimization loop until full electrification is reached.

After a single simulation is run the rotations are analyzed. Any time a vehicle goes below an soc of zero (or a self defined value) a low soc event is triggered. This event saves information about when the soc reached its minimal value and the history before that up to a point of an upper soc threshold, with the default value being 1. Stations inside of this time span are potentially able to mitigate the low soc and are stored with other information about the event. :numref:`low_soc_event` shows a possible soc history with a low soc event.

.. _low_soc_event:
.. figure:: https://user-images.githubusercontent.com/104760879/217225588-abfad83d-9d2a-463a-8597-584e29f5f885.png
    :width: 600
    :alt: below_0_soc_event

    Low soc event and classification of stations.

The next step groups low soc events based on the stations which were found earlier. Events which share at least one station could possibly interact with each other, e.g. vehicles could share a charging station. Therefore groups are build which do not share any stations in between groups. This speeds up the optimization process since for every electrification and simulation only rotations are calculated which could be impacted by the change.

Since greedy approaches execute the step which seems most promising in the current situation an evaluation function is needed. One possible approach could be to simulate each scenario, meaning simulating every case in which one of all possible stations is electrified and continuing with the best case. The optimizer does not use this approach. Instead an approximation function is used to evaluate the potential of electrifying a station. This approximation function analyzes the duration at each stop, the possible charging time, the soc and resulting possible charging power (battery with high socs are charged at a lower rate) as well as the upper soc threshold and minimal soc of the event. While this methodology is not accurate in all cases, e.g. a station could exist multiple times inside of a low soc event, therefore charging the first time at this station would alter the soc and charging power the vehicle has the second time it reaches the station, it seems well suited as heuristic for choosing the most promising station. The objective function of choosing what the *best* station is, is the mitigation of missing charge, i.e. what is the minimal amount of energy that needs to be inserted into the battery, so that no soc is below 0.

After the evaluation selected a station to be electrified the scenario input data is altered so that vehicles at this station are charged without limitation of charging points. This is followed up by a detailed simulation which can make use of a highly accurate solver for charging events called *SpiceEV* or a less accurate but faster solver. Now the resulting system has less missing charge and the potentials of stations might be decreased. Also a single group might have been split up into several smaller groups which can be analyzed even quicker. Therefore the loop repeats up until the point the missing charge in the system is zero or in other words the system is fully electrified.

At the current stage the scenario to be optimized needs depot charging stations at the start and end of each rotation. The scenario should not contain any opportunity charging stations. If for a given scenario opportunity charging stations are predefined, i.e. the scenario should contain a specific electrification and is set in the *electrified_station.json* the solver type *spice_ev* should be used in the *optimizer.cfg*. If the *quick* solver is supposed to be used the station can be listed in *inclusion_stations* while the *electrified_stations.json* should only contain depot stations. Stations can be also excluded from optimization by adding their name to *exclusion_stations*.

Deep Optimization
####################
The greedy algorithm in the base optimization can not guarantee that the solution is the global optimum. This is why the use of the *deep* mode is recommended for systems with high requirements. After the first run, instead of electrifying the station with the highest potential the second best station is electrified. This is similar to a decision tree, where every node is a set of electrified stations, with the first node being zero stations electrified and the last node being all stations electrified. The nodes in between correlate with every possible state of electrification. Each branch therefore represents an additional electrification of a single station. . The algorithm continues electrifying the best station, as long as this node has not been evaluated yet. This way gradually all possible nodes are checked. The search stops whenever the number of stations surpasses the number of the current optimal solution. If several options with the same optimal number of stations arise, they can be found in the log file of the optimizer, but only one file with optimized stations is produced.

**Pruning** is used to stop evaluation of branches, whenever foresight predicts that no better solution will be reached. This is done through the simple heuristic of checking the sum of potential of the n remaining stations with the highest potentials, with n being the number until the number of stations of the current optimal solution is reached.

| **Example:**
| The base optimization found a set of 5 stations to fully electrify the scenario. These stations are *A*, *B*, *C*, *D* and *E* which were chosen in the same order. The whole scenario consists of the whole alphabet of stations. The deep optimization starts with evaluating a scenario without any electrified opportunity stations. Depot stations are electrified. The first evaluation gives a sorted list of potentials by

========  =====
Station   Potential
========  =====
*A*       85
*X*       75
*B*       30
*E*       25
...       ...
========  =====

In the base optimization Station *A* was chosen since it showed the highest potential. The deep optimization ignores this node since it has been evaluated already and chooses station *X* instead. After a detailed simulation with *X* electrified, the remaining stations are evaluated again.

========  =====
Station   Potential
========  =====
*B*       28
*E*       25
*C*       20
*G*       18
...       ...
========  =====

For every vehicle the amount of missing energy is calculated and summed up. In this example case the missing energy is 85. Since 4 stations are remaining until the current optimum of 5 stations is reached, the 4 stations with the highest potential are evaluated in this case

.. math::

   Pot = Pot_B + Pot_E + Pot_C + Pot_G = 28 + 25 + 20 +18 = 91

In this case the potential is high enough to continue the exploration of this branch. If the potential would have been below 85 the branch would have been pruned, meaning it would not be explored any further and labeled as *not promising*. It is not promising since it will not lead to a better solution than the current one. This is the case since on one hand the evaluation by approximation tends to overestimate the potential while the missing energy is accurately calculated and on the other hand electrification of stations can reduce the potential of other stations, for example if 2 stations charge the same rotation, electrifying one station might fully electrify the rotation meaning the potential of the other station drops to zero.
This concept can reduce the amount of nodes which have to be checked.

**Mandatory stations** can be attained to increase the optimization process. Mandatory stations are defined by being stations which are needed for a fully electrified system. To check if a station *Y* is a mandatory station can be easily attained by simulating the network with every station electrified except *Y*. If the system has vehicle socs which drop below the minimal soc (default value is 0) in this scenario, the station is mandatory. In the later exploration of best combinations of stations this station will be included in any case.

**Impossible rotations** are rotations which given the settings are not possible to be run as opportunity chargers, given the vehicle properties, even when every station is electrified. Before starting an optimization it is recommended to remove these rotations from the optimization, since the optimizer will not reach the goal of full electrification.

**Continuing optimizations** can be useful in cases where simulation of the base case is slow or considerable effort was put into optimization before. The user might want to continue the optimization from the state where they left off. To speed up multiple optimizations or split up a big optimization in multiple smaller calculations two features are in early development. Experienced users can use these features on their own accord with a few minor implementation steps. To skip a potentially long simulation, with the simulation of the scenario being the first step of every ebus-toolbox run, the optimizer.config allows for using pickle files for the three major objects args, schedule and scenario. After pickling the resulting objects, the optimizer can be prompted to use them instead of using whatever other input is fed into the optimizer. This is done by giving the paths to the pickle files in the optimizer.cfg.

::

    args = data/args.pickle
    schedule = data/schedule.pickle
    scenario = data/scenario.pickle

If they are provided they are used automatically. All three pickle files need to be set.

If a deep optimization takes to long to run it in one go, it is possible to save the state of the decision tree as pickle file as well. Reloading of the state is possible and will lead to a continuation of the previous optimization.
To make use of this feature the parameters in the optimizer.cfg have to be set.

::

    decision_tree_path = data/last_optimization.pickle
    save_decision_tree = True


The functionality of the optimizer is controlled through the optimizer.config













Report
-------------
The report will include information about the expected socs, power loads at the charging stations or depots, default plots for the scenario and other useful data.
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












