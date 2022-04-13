# eBus Toolbox

The eBus Toolbox assists the user in analysing and optimising electrified bus fleets and schedules.

### Usage

At the current stage, only a single functionality is implemented, which is processing a bus schedule stored in a specific CSV format (see `data/examples/trips_examples.csv`) and run it through a module called SpiceEV for an in-depth SOC analysis.

To try it out, first clone this repository and then install the required packages to your current environment by running

`pip install -r requirements.txt` 

Now you can start the eBus Toolbox with all configurations stored at `data/configs/ebus_toolbox.cfg` module via the command.

``python -m ebus_toolbox --config data/configs/ebus_toolbox.cfg``

The repo provides an example for each necessary input file, so the example case can be executed without the need for the user to provide any data themselves.

To run the eBus Toolbox with your own `schedule.csv` (see details [below](#input-data)) file and default configurations run

`python -m ebus_toolbox --input_schedule path/to/schedule.csv`

Default configurations are detailed at `data/configs/ebus_toolbox.cfg`.



### Input Data

To analyze your own electric bus schedule, the data needs to be provided as a CSV file where each row contains the details of a single trip of that schedule. Find the details about the various columns in this file below. The first table lists the **mandatory** columns while the second one (tbd) lists optional parameters. Refer to `data/examples/trips.csv` for an example.

| Column Name    | Description                                                  | Example           |
| -------------- | ------------------------------------------------------------ | ----------------- |
| rotation_id    | Unique alphanumeric ID to identify rotations                 | 27312             |
| departure_name | Name of the station the trip starts at                       | Warschauer Straße |
| departure_time | Date and Time at which bus starts trip                       | 2022-03-13T10:25  |
| arrival_name   | Name of the station the trip ends at                         | Ostbahnhof Berlin |
| arrival_time   | Date and Time at which bus completes trip (e.g. yyyy-mm-ddThh:mm[:ss]) | 2022-03-13T10:30  |
| distance       | Distance traveled in **m**                                   | 1340              |
| vehicle_type   | ID of vehicle type defined in vehicle types file. Set path of this file in config.<br />(see default for reference: `data/examples/vehicle_types.json`) | some_bus_type     |

