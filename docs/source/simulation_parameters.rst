.. _simulation_parameters:

Simulation Parameters
=====================

The simulation of an eBus-System relies on a variety of simulation parameters.
The eBus-Toolbox provides most of them as default values. Depending on specific needs adjusting
these values can increase the accuracy of the simulation outputs. The eBus-Toolbox input files are described
in detail in the following subsections as well as their default parameters.
When providing the user defined input files, the user should make sure the files are either 'utf-8'
encoded or not contain regional characters.

Configuration
-------------
The configuration file config.cfg is provided as example in ./examples/ and provides the user with most of the functionality surrounding the settings and boundary conditions of a simulation. The example contains parameter descriptions which are explained here in more detail:

.. _schedule:

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
     - Date and Time at which bus starts trip
     - 2022-03-13T10:25
   * - arrival_name
     - Name of the station the trip ends at
     - Ostbahnhof Berlin
   * - arrival_time
     - Date and Time at which bus completes trip (e.g. yyyy-mm-ddThh:mm[:ss])
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


.. _vehicle_types:

Vehicle types
-------------
vehicle_type.json
tbc

Electrified stations
--------------------
Stations which are electrified. TBC

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
The consumption table can be referenced in the :ref:`vehicle_types` file. Instead of constant consumption the eBus-Toolbox uses provided temperatures, level of loadings, mean speeds, average inclines and the vehicle type to interpolate the consumption value from this data table. Level of loading and temperatures are read from the :ref:`schedule` if the trips provide them. If they are missing from the schedule, they are looked up from the files :ref:`level_of_loading` and :ref:`temperature_data`. The average incline is calculated from :ref:`station_geo_data` and the mean speed is calculated by using the departure and arrival time and distance provided by the schedule.