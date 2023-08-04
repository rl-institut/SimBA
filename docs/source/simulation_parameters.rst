.. _simulation_parameters:

Simulation Parameters
=====================

The simulation of an eBus-System relies on a variety of simulation parameters.
SimBA provides most of them as default values. Depending on specific needs adjusting
these values can increase the accuracy of the simulation outputs. The SimBA input files are described
in detail in the following subsections as well as their default parameters.
When providing the user defined input files, the user should make sure the files are either 'utf-8'
encoded or not contain regional characters.

.. _config:

Configuration
-------------
The configuration file config.cfg is provided as example in ./examples/ and provides the user with most of the functionality surrounding the settings and boundary conditions of a simulation. The config file is divided into several sections:

| Paths: Here all in- and output paths are defined
| Modes: Here the modes of execution can be defined
| Flags: Here, optional functionality can be activated or deactivated
| Physical setup of environment: Here, the physical setup is characterized
| Simulation Parameters: The simulation can be adjusted using these parameters

The example (data/simba.cfg) contains parameter descriptions which are explained here in more detail:

.. list-table:: config.cfg parameters
   :header-rows: 1

   * - Parameter
     - Default value
     - Expected values
     - Description
   * - input_schedule
     - Mandatory: no default given
     - Path as string
     - Input file containing :ref:`schedule` information
   * - Output_directory
     - Data/sim_outputs
     - Path as string
     - Output files are stored here
   * - electrified_stations
     - ./data/examples/vehicle_types.json
     - Path as string
     - Path to Electrified stations data
   * - vehicle_types
     - ./data/examples/vehicle_types.json
     - Path as string
     - Path to :ref:`vehicle_types`
   * - station_data_path
     - Optional: no default given
     - Path as string
     - Path to :ref:`station_geo_data`
   * - outside_temperature_over_day_path
     - Optional: no default given
     - Path as string
     - Path to :ref:`temperature_data`
   * - level_of_loading_over_day_path
     - Optional: no default given
     - Path as string
     - Path to :ref:`level_of_loading`
   * - cost_parameters_file
     - Optional: no default given
     - Path as string
     - Path to :ref:`cost_params`
   * - mode
     - ['sim', 'report']
     - List of modes is any order in range of ['sim', 'neg_depb_to_oppb', 'neg_oppb_to_depb', 'service_optimization', 'report']
     - The order of :ref:`sim_modes` is defined here
   * - cost_calculation
     - false
     - Boolean
     - Activates the :ref:`cost_calculation`
   * - check_rotation_consistency
     - false
     - Boolean
     - Activates the :ref:`consistency_check`
   * - skip_inconsistent_rotations
     - false
     - Boolean
     - If check_rotation_consistency is active, rotations that don't comply with the checked assumptions are removed from the schedule if skip_inconsistent_rotations is set true
   * - show_plots
     - false
     - Boolean
     - If activated, plots are displayed with every run of :ref:`report` mode

   * - preferred_charging_type
     - depb
     - depb, oppb
     - All rotations that have no specification of charging type in :ref:`Schedule` are assigned the charging type defined here
   * - gc_power_opps
     - 100000
     - Numeric
     -  Default max. power [kW] of grid connectors at opportunity charging stations, Individual gc_power per gc can be defined in :ref:`electrified_stations`
   * - gc_power_deps
     - 100000
     - Numeric
     -  Default max. power [kW] of grid connectors at depot charging stations, Individual gc_power per gc can be defined in :ref:`electrified_stations`
   * - cs_power_opps
     - 300
     - Numeric
     - Default max. power [kW] of opportunity charging stations
   * - cs_power_deps_depb
     - 300
     - Numeric
     - Default max. power [kW] of depot charging stations for depot charging buses. Individual cs_power per gc and cs type can be defined in :ref:`electrified_stations`
   * - cs_power_deps_oppb
     - 300
     - Numeric
     - Default max. power [kW] of depot charging stations for opportunity charging buses. Individual cs_power per gc and cs type can be defined in :ref:`electrified_stations`
   * - desired_soc_deps
     - 1
     - 0...1
     - Minimum allowed state of charge when leaving a depot station after charging. Also used to initialize the vehicles SoCs at the beginning of the simulation.
   * - desired_soc_opps
     - 1
     - 0...1
     - Minimum allowed state of charge when leaving an opportunity station after charging
   * - min_recharge_deps_oppb
     - 1
     - 0...1
     - This value is used to calculate the minimum standing time of opportunity charging busses at the depot, which is needed for the :ref:`vehicle_dispatch`
   * - min_recharge_deps_depb
     - 1
     - 0...1
     - This value is used to calculate the minimum standing time of depot charging busses at the depot, which is needed for the :ref:`vehicle_dispatch`
   * - min_charging_time
     - 0
     - Numeric
     - Only stops that are longer than the time defined here are used for charging
   * - default_buffer_time_opps
     - 0
     - Numeric or dict e.g. {"10-22": 5, "else": 2} (else clause is a must if using the dict definition)
     - The buffer time is deducted off of the planned standing time at each opportunity station. It can be used to model things like delays and/or docking procedures. This value is used if no specific buffer is defined per station in :ref:`electrified_stations`. It can either be given as constant or depending on the time of the day using a dict.
   * - default_voltage_level
     - MV
     - HV, HV/MV, MV, MV/LV, LV
     - The default voltage level is used, if no specific voltage level is defined per station in :ref:`electrified_stations`. It is used to calculate the costs. Choices describe high voltage (HV), transformer between high and medium voltage (HV/MV), medium voltage MV, transformer between medium and low voltage (MV/LV) and low voltage (LV)

   * - days
     - Optional: no default given
     - Numeric
     - If this value is defined only the first number of 'days' of the schedule are simulated
   * - interval
     - 1
     - Numeric
     - Timestep in minutes
   * - signal_time_dif
     - 10
     - Numeric
     - Some strategies use limited foresight. E.g. prioritization of vehicles at limited number of charging stations is carried out only for this time ahead of actual time step. Also used in spiceEV as time difference between signal time and actual start time of a vehicle event in min.
   * - eta
     - false
     - Boolean
     - Show estimated time to finish simulation after each step. Not recommended for fast computations


Schedule
--------

To analyze your own electric bus schedule, the data needs to be provided as a .csv file where each row contains the details of a single trip of that schedule. Find the details about the various columns in this file below. The first table lists the **mandatory** columns while the second one (tbd) lists optional parameters. Refer to `data/examples/trips.csv` for an example.

.. list-table:: schedule mandatory input
   :widths: 150 300 150
   :header-rows: 1

   * - Column Name
     - Description
     - Example
   * - rotation_id
     - Unique alphanumeric ID to identify rotations
     - 27312
   * - departure_name
     - Name of the station the trip starts at
     - Warschauer Straße
   * - departure_time
     - Date and time at which bus starts trip (ISO-Format)
     - 2022-03-13T10:25
   * - arrival_name
     - Name of the station the trip ends at
     - Ostbahnhof Berlin
   * - arrival_time
     - Date and Time at which bus completes trip (ISO-Format) (e.g. yyyy-mm-ddThh:mm[:ss])
     - 2022-03-13T10:30
   * - distance
     - Distance traveled in **m**
     - 1340
   * - vehicle_type
     - | ID of vehicle type defined in vehicle types file. Set path of this file in config
       | (see default for reference: `data/examples/vehicle_types.json`)
     - some_bus_type

.. list-table:: schedule optional input
   :widths: 150 300 150
   :header-rows: 1

   * - Column Name
     - Description
     - Example
   * - line
     - The bus line
     - 512, M10, X11 etc.
   * - charging_type
     - | The preferred charging type for this trip.
       | NOTE: All trips of a rotation need to have the same charging type.
       | If omitted, charging type is set according to preferred charging type provided in the config file.
     - Options: **depb**,  **oppb**
   * - temperature
     - Temperature of the trip in **degC**
     - 25
   * - level_of_loading
     - The level of loading of the bus on this trip in between 0 and 1
     - 0.5.

This is how a schedule file might look like.

+--------+----------------+---------------------+---------------------+--------------+----------+-------------+--------------+-------------+------------------+
| line   | departure_name | departure_time      | arrival_time        | arrival_name | distance | rotation_id | vehicle_type | temperature | level_of_loading |
+========+================+=====================+=====================+==============+==========+=============+==============+=============+==================+
| LINE_0 | Station-0      | 2022-03-07 21:28:00 | 2022-03-07 21:31:00 | Station-1    | 1530     | 1           | 12m_bus      | 20          | 0                |
+--------+----------------+---------------------+---------------------+--------------+----------+-------------+--------------+-------------+------------------+
| LINE_0 | Station-1      | 2022-03-07 21:31:00 | 2022-03-07 22:04:00 | Station-3    | 14519    | 1           | 12m_bus      | -5          | 0.9              |
+--------+----------------+---------------------+---------------------+--------------+----------+-------------+--------------+-------------+------------------+
| LINE_0 | Station-3      | 2022-03-07 22:08:00 | 2022-03-07 22:43:00 | Station-1    | 13541    | 1           | 12m_bus      |             |                  |
+--------+----------------+---------------------+---------------------+--------------+----------+-------------+--------------+-------------+------------------+
| LINE_0 | Station-1      | 2022-03-07 22:51:00 | 2022-03-07 23:24:00 | Station-2    | 14519    | 1           | 12m_bus      |             |                  |
+--------+----------------+---------------------+---------------------+--------------+----------+-------------+--------------+-------------+------------------+


.. _vehicle_types:

Vehicle types
-------------
vehicle_type.json
tbc

.. _electrified_stations:

Electrified stations
--------------------
Stations which are electrified. TBC

.. _cost_params:

Cost parameters
---------------
TBC



TBC

.. _station_geo_data:

Station data
------------
Geodata. TBV


.. _level_of_loading:

Level of loading
----------------
TBC

.. _temperature_data:

Temperatures
------------

TBC

.. _consumption_table:

Consumption table
-----------------
The consumption table can be referenced in the :ref:`vehicle_types` file. Instead of constant consumption SimBA uses provided temperatures, level of loadings, mean speeds, average inclines and the vehicle type to interpolate the consumption value from this data table. Level of loading and temperatures are read from the :ref:`schedule` if the trips provide them. If they are missing from the schedule, they are looked up from the files :ref:`level_of_loading` and :ref:`temperature_data`. The average incline is calculated from :ref:`station_geo_data` and the mean speed is calculated by using the departure and arrival time and distance provided by the schedule.