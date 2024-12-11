.. _simulation_parameters:

Simulation Parameters
=====================

The simulation of an eBus-System relies on a variety of simulation parameters.
SimBA provides most of them as default values. Depending on specific needs adjusting
these values can increase the accuracy of the simulation outputs. The SimBA input files are described in detail in the following subsections as well as their default parameters. When providing the user defined input files, the user should make sure the files are either 'utf-8' encoded or not contain regional characters. Except for the ``electrified_stations.json``, all .json files can be enriched with in-line comments starting with :code:`\\\\`.

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
   * - scenario_name
     - Optional: no default given
     - string
     - scenario identifier, appended to output directory name and report file names
   * - schedule_path
     - Mandatory: no default given
     - Path as string
     - Input file containing :ref:`schedule` information
   * - output_path
     - Data/sim_outputs
     - Path as string
     - Output files are stored here; set to null to deactivate
   * - electrified_stations_path
     - ./data/examples/vehicle_types.json
     - Path as string
     - Path to Electrified stations data
   * - vehicle_types_path
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
   * - cost_parameters_path
     - Optional: no default given
     - Path as string
     - Path to :ref:`cost_params`
   * - optimizer_config_path
     - Optional: no default given
     - Path as string
     - Path to station optimizer config :ref:`optimizer_config`
   * - rotation_filter_path
     - Optional: no default given
     - Path as string
     - Path to rotation filter json
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
   * - rotation_filter_variable
     - null
     - string
     - How to filter rotations according to file 'rotation_filter': options are "include" (whitelist), "exclude" (blacklist), null (ignore)
   * - create_trips_in_report
     - false
     - Boolean
     - Write a new trips.csv during report mode to output directory?
   * - create_pickle_in_report
     - false
     - Boolean
     - Pickle current schedule and scenario during report mode
   * - load_pickle
     - Optional, no default given
     - Path to pickle file
     - Load schedule and scenario from this pickle file, expects load_pickle as first mode
   * - extended_output_plots
     - false
     - Boolean
     - If set, create additional plots when running :ref:`report` mode
   * - strategy_deps
     - balanced
     - SpiceEV Strategies (greedy, balanced, peak_shaving, peak_load_windows, balanced_market)
     - Charging strategy used in depots.
   * - strategy_opps
     - greedy
     - SpiceEV Strategies (greedy, balanced, peak_shaving, peak_load_windows, balanced_market)
     - Charging strategy used in opportunity stations.
   * - cost_calculation_method_deps
     - fixed_wo_plw
     - SpiceEV cost calculation type (fixed_wo_plw, fixed_w_plw, variable_wo_plw, variable_w_plw, balanced_market, flex_window)
     - Method for cost calculation at depots.
   * - cost_calculation_method_opps
     - fixed_wo_plw
     - SpiceEV cost calculation type, same choices as in depot
     - Method for cost calculation at opportunity stations.
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
     - The buffer time in minutes is subtracted from of the planned standing time at each opportunity station. It can be used to model things like delays and/or docking procedures. This value is used if no specific buffer is defined per station in :ref:`electrified_stations`. It can either be given as constant or depending on the time of the day using a dict.
   * - default_buffer_time_deps
     - 0
     - Numeric
     - The buffer time in minutes is subtracted from of the planned standing time at each depot station. It can be used to model things like delays and/or docking procedures. This value is used for every depot station
   * - assign_strategy
     - adaptive
     - adaptive, min_recharge
     - The value of assign_strategy sets the algorithm of vehicle disposition. "adaptive" uses vehicles to service rotations with the lowest soc, without the rotation getting negative. "min_recharge" only uses vehicles which are above the charge type specific threshold (see min_recharge_deps_oppb, min_recharge_deps_depb)
   * - default_voltage_level
     - MV
     - HV, HV/MV, MV, MV/LV, LV
     - The default voltage level is used, if no specific voltage level is defined per station in :ref:`electrified_stations`. It is used to calculate the costs. Choices describe high voltage (HV), transformer between high and medium voltage (HV/MV), medium voltage MV, transformer between medium and low voltage (MV/LV) and low voltage (LV)
   * - loglevel
     - INFO
     - DEBUG, INFO, WARN or ERROR
     - Log level. All logging messages are both displayed in the console and written to a log file
   * - logfile
     - <datetime>.log
     - String
     - Log file name. Set to null to disable logging to file
   * - loglevel_file
     - (same as loglevel)
     - String
     - Log level for file logger
   * - default_mean_speed
     - 30
     - numeric
     - Default assumed mean speed for busses in km/h. Used in split_negative_depb for generating depot trips.
   * - default_depot_distance
     - 5
     - numeric
     - Default assumed average distance from any station to a depot in km. Used in split_negative_depb for generating depot trips.
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
   * - skip_flex_report
     - false
     - Boolean
     - Skip generation of flex_report in SpiceEV. Activating can save time as this feature is rarely used


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

The vehicle types that can be used are defined in the "vehicle_type.json". The path to this file has to be defined in the :ref:`config` and an example is given at "data/examples/vehicle_types.json".

The data is structured as a .json where the top level key represents the vehicle_type, that needs to correspond to the "vehicle_type" defined in the :ref:`schedule`. The next level key defines the charging_type ("oppb" or "depb"). For one vehicle type either one or both charging types can be defined and for each given charging type the specifications in the third level of the .json have to be given. In this level, the parameters for the specified vehicle are be defined. The specification of one vehicle with the vehicle_type "AB" and the charging_types "depb" and "oppb" is given as follows:

.. code-block:: json

    {
        "AB": {  // vehicle_type
            "depb": {  // charging_type
                "name": "articulated bus - depot charging",  // long name
                "capacity": 250,  // battery capacity in kWh
                "charging_curve": [[0, 150], [0.8, 150], [1, 15]],  // charging curve [SoC, kW]
                "min_charging_power": 0,  // min charging power in KW
                "v2g": false,  // Is vehicle capable of vehicle to grid?
                "mileage": "data/examples/energy_consumption_example.csv",  // mileage in kWh/km or link to consumption.csv
                "battery_efficiency": 0.95  // optional. default: 0.95
            },
            "oppb": {
                "name": "articulated bus - opportunity charging",
                "capacity": 150,
                "charging_curve": [[0, 250], [0.8, 250], [1, 25]],
                "min_charging_power": 0,
                "v2g": false,
                "mileage": "data/examples/energy_consumption_example.csv"
            }
        }
    }

.. _electrified_stations:

Electrified stations
--------------------

All stations, that are or could be equipped with charging infrastructure have to be parameterized in the "electrified_stations.json" together with their grid connection, charging infrastructure and local energy systems. The path to this file has to be defined in the :ref:`config`.

The data is structured as a .json where the top level key represents the station name, that needs to correspond to the "departure_name", respectively "arrival_name" defined in the :ref:`schedule`. Each station has two mandatory arguments: "type" defines if a station is a depot ("deps") or an opportunity charging station ("opps") and "n_charging_stations" limits the amount of vehicles, that can simultaneously charge at one station.

Furthermore, the energy system at each station can be characterized in terms of local power generation ("energy_feed_in"), local external loads ("external_load") or local stationary electric power storage ("battery"). An example that displays all further parameters and the specification of the local energy systems is given at "data/examples/electrified_stations.json".


.. _cost_params:

Cost parameters
---------------
In order to run the :ref:`cost_calculation`, all cost parameters are to be defined in the ``cost_params.json``. The file is used as input for both, SimBA and SpiceEV, as both tools do part of the cost calculation and therefore no comments are allowed here. If not otherwise specified the investments/costs are gross prices. A commented example is given below, for a working example please refer to "data/examples/cost_params.json".

.. code-block:: json

    {
        "vehicles": {  // all vehicles and charging types have to be defined here
            "SB_debp": {  // all combinations of vehicle types and charging types have a separate cost definition, the name is to be given as [vehicle_type]_[charging_type]
                "capex": 500000,  // investment cost for one vehicle without vehicle battery
                "c_maint_per_km": 0.24,  // maintenance cost per km
                "lifetime": 14  // lifetime of the vehicle in years
            }
        },
        "batteries": {  // vehicle battery
            "lifetime_battery": 7,   // lifetime of the vehicle battery in years
            "cost_per_kWh": 250  // investment cost for vehicle battery per kWh
        },
        "stationary_storage": {    // stationary electric energy storage
            "capex_fix": 1,  // fix investment cost for stationary storage
            "capex_per_kWh": 1,  //  investment cost for stationary storage per kWh
            "c_maint_stat_storage_per_year": 0.02,  // annual maintenance costs in % of capex
            "lifetime_stat_storage": 20  // lifetime in years
        },
        "cs":{  // charging stations
            "capex_opps_per_kW": 877.5,  //  investment cost for opportunity charging stations per kW
            "capex_deps_per_kW": 1000,  //  investment cost for depot charging stations per kW
            "lifetime_cs": 20,  // lifetime of charging stations in years
            "c_maint_cs_per_year": 0.02  // annual maintenance costs in % of capex
        },
        "garage": {
            "n_charging_stations": 1,  // number of charging stations for the garage
            "power_cs": 50,  // power of the charging stations for the garage
            "vehicles_per_workstation": 20,  // how many vehicles share one workstation
            "cost_per_workstation": 245000,  //  investment cost for one workstation
            "lifetime_workstations": 20  // lifetime in years
        },
        "grid operator": {
            "gc": {  // grid connection
                "LV": {  // grid connection in specific voltage level. Options are "HV", "HV/MV", "MV", "MV/LV", "LV" and all relevant voltage levels have to be defined here
                    "default_distance": 50,  // Used if not specified individually in electrified_stations.json
                    "capex_gc_fix": 100,  // fix investment cost for establishing a grid connection
                    "capex_gc_per_meter": 16.85,  // investment cost per meter
                    "capex_gc_per_kW": 24.14,  // investment cost per kW
                    "capex_transformer_fix": 0,  // fix investment cost for a transformer
                    "capex_transformer_per_kW": 0  // fix investment cost for a transformer per kW
                },
                "lifetime_gc": 50,  // lifetime of the grid connection in years
                "c_maint_transformer_per_year": 0.02,  // annual maintenance costs in % of capex
                "lifetime_transformer": 20  // lifetime in years
            }
        }
    }

All remaining parameters such as grid fees or energy taxes are described in the example file.


.. _station_geo_data:

Station data
------------
The file "all_stations.csv" contains information that is relevant for all stations regardless of their status of electrification. At this stage of development this reduces to the information of station height that is relevant only if a trip specific :ref:`consumption_analysis` is employed. See the example at "data/examples/all_stations.csv" for the required structure.


.. _level_of_loading:

Level of loading
----------------

If a trip specific :ref:`consumption_analysis` is employed, the level of loading for each trip is required. This information can be detailed in the :ref:`schedule`. If not specified there, a default value for every hour of the day can be specified in this file. See the example at "data/examples/default_level_of_loading_over_day.csv" for the required structure.


.. _temperature_data:

Temperatures
------------
If a trip specific :ref:`consumption_analysis` is employed, the temperature for each trip is required. This information can be detailed in the :ref:`schedule`. If not specified there, a default value for every hour of the day can be specified in this file. See the example at "data/examples/default_temp_summer.csv" for the required structure.

.. _consumption_table:

Consumption table
-----------------
The consumption table can be referenced in the :ref:`vehicle_types` file. Instead of constant consumption SimBA uses provided temperatures, level of loadings, mean speeds, average inclines and the vehicle type to interpolate the consumption value from this data table. Level of loading and temperatures are read from the :ref:`schedule` if the trips provide them. If they are missing from the schedule, they are looked up from the files :ref:`level_of_loading` and :ref:`temperature_data`. The average incline is calculated from :ref:`station_geo_data` and the mean speed is calculated by using the departure and arrival time and distance provided by the schedule.
