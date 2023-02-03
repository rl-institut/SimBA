Simulation Parameters
=====================

The simulation of an eBus-Sytem relies on a variety of simulation parameters.
The eBus-Toolbox provides most of them as default values. Depending on specific needs adjusting
these values can increase the accuracy of the simulation outputs. The eBus-Toolbox input files are described
in detail in the following subsections as well as their default parameters.
When providing the user defined input files, the user should make sure the files are either 'utf-8'
encoded or not contain regional characters.

Configuration
-------------

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




