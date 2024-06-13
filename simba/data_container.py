"""Module to handle data access, by the varying SimBA modules."""
import argparse
import csv
import logging
import datetime
from pathlib import Path
from typing import Dict

import pandas as pd

from simba import util, ids
from simba.ids import INCLINE, LEVEL_OF_LOADING, SPEED, T_AMB, TEMPERATURE, CONSUMPTION


class DataContainer:
    def __init__(self):
        # Dictionary of dict[VehicleTypeName][ChargingType] containing the vehicle_info dictionary
        self.vehicle_types_data: Dict[str, any] = {}
        # Dictionary of dict[consumption_lookup_name] containing a consumption lookup
        self.consumption_data: Dict[str, pd.DataFrame] = {}
        # Dictionary of dict[hour_of_day] containing temperature in °C
        self.temperature_data: Dict[int, float] = {}
        # Dictionary of dict[hour_of_day] containing level of loading [-]
        self.level_of_loading_data: Dict[int, float] = {}
        # Dictionary of dict[station_name] containing information about electrification
        self.stations_data: Dict[str, dict] = {}
        # Dictionary containing various infos about investment costs and grid operator
        self.cost_parameters_data: Dict[str, dict] = {}
        # Dictionary containing all stations and their geo location (lng,lat,elevation)
        self.station_geo_data: Dict[str, dict] = {}
        # List of trip dictionaries containing trip information like arrival time and station
        # departure time and station, distance and more
        self.trip_data: [dict] = []

    def fill_with_args(self, args: argparse.Namespace) -> 'DataContainer':
        """ Fill self with data from file_paths defined in args

        :param args: Arguments containing paths for input_schedule, vehicle_types_path,
            electrified_stations_path, cost_parameters_path, outside_temperature_over_day_path,
            level_of_loading_over_day_path, station_data_path
        :return: self
        """

        return self.fill_with_paths(
            trips_file_path=args.input_schedule,
            vehicle_types_path=args.vehicle_types_path,
            electrified_stations_path=args.electrified_stations_path,
            cost_parameters_path=args.cost_parameters_path,
            outside_temperature_over_day_path=args.outside_temperature_over_day_path,
            level_of_loading_over_day_path=args.level_of_loading_over_day_path,
            station_data_path=args.station_data_path,
        )

    def add_trip_data_from_csv(self, file_path: Path) -> 'DataContainer':
        """ Add trip data from csv file to DataContainer

        :param file_path: csv file path
        :return: self with trip data
        """

        self.trip_data = []
        with open(file_path, 'r', encoding='utf-8') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                trip_d = dict(trip)
                trip_d["arrival_time"] = datetime.datetime.fromisoformat(trip["arrival_time"])
                trip_d["departure_time"] = datetime.datetime.fromisoformat(trip["departure_time"])
                trip_d[LEVEL_OF_LOADING] = util.cast_float_or_none(trip.get(LEVEL_OF_LOADING))
                trip_d[TEMPERATURE] = util.cast_float_or_none(trip.get(TEMPERATURE))
                trip_d["distance"] = float(trip["distance"])
                self.trip_data.append(trip_d)
        return self

    def add_station_geo_data(self, data: dict) -> None:
        """Add station_geo data to the data container.

        Used when adding station_geo to a data container from any source
        :param data: data containing station_geo
        :type data: dict
        """
        self.station_geo_data = data

    def fill_with_paths(self,
                        trips_file_path,
                        vehicle_types_path,
                        electrified_stations_path,
                        outside_temperature_over_day_path=None,
                        level_of_loading_over_day_path=None,
                        cost_parameters_path=None,
                        station_data_path=None,
                        ):
        """ Fill self with data from file_paths
        :param trips_file_path: csv path to trips
        :param vehicle_types_path: json file path to vehicle_types
        :param electrified_stations_path: json file path to electrified stations
        :param outside_temperature_over_day_path: csv path to temperatures over the hour of day
        :param level_of_loading_over_day_path: csv path to level of loading over the hour of day
        :param cost_parameters_path: json file path to cost_parameters
        :param station_data_path: csv file path to station geo_data
        :return: self
        """

        # Add the vehicle_types from a json file
        self.add_vehicle_types_from_json(vehicle_types_path)

        # Add consumption data, which is found in the vehicle_type data
        self.add_consumption_data_from_vehicle_type_linked_files()

        if outside_temperature_over_day_path is not None:
            self.add_temperature_data_from_csv(outside_temperature_over_day_path)

        if level_of_loading_over_day_path is not None:
            self.add_level_of_loading_data_from_csv(level_of_loading_over_day_path)

        # Add electrified_stations data
        self.add_stations_from_json(electrified_stations_path)

        # Add cost_parameters_data
        if cost_parameters_path:
            self.add_cost_parameters_from_json(cost_parameters_path)

        # Add station geo data
        if station_data_path is not None:
            self.add_station_geo_data_from_csv(station_data_path)

        self.add_trip_data_from_csv(trips_file_path)
        return self

    def add_station_geo_data_from_csv(self, file_path: Path) -> 'DataContainer':
        # find the temperature and elevation of the stations by reading the .csv file.
        # this data is stored in the schedule and passed to the trips, which use the information
        # for consumption calculation. Missing station data is handled with default values.
        self.station_geo_data = dict()
        line_num = None
        try:
            with open(file_path, "r", encoding='utf-8') as f:
                delim = util.get_csv_delim(file_path)
                reader = csv.DictReader(f, delimiter=delim)
                for line_num, row in enumerate(reader):
                    self.station_geo_data[str(row['Endhaltestelle'])] = {
                        ids.ELEVATION: float(row[ids.ELEVATION]),
                        ids.LATITUDE: float(row.get(ids.LATITUDE, 0)),
                        ids.LONGITUDE: float(row.get(ids.LONGITUDE, 0)),
                    }
        except (FileNotFoundError, KeyError):
            logging.log(msg=f"External csv file {file_path} not found or not named properly. "
                            "(Needed column names are 'Endhaltestelle' and 'elevation')",
                        level=100)
            raise
        except ValueError:
            line_num += 2
            logging.log(msg=f"Can't parse numeric data in line {line_num} from file {file_path}.",
                        level=100)
            raise
        return self

    def add_level_of_loading_data(self, data: dict) -> 'DataContainer':
        """Add level_of_loading data to the data container.

        Used when adding level_of_loading to a data container from any source
        :param data: data containing hour and level_of_loading
        :type data: dict
        :return: DataContainer containing level of loading data
        """
        self.level_of_loading_data = data
        return self

    def add_level_of_loading_data_from_csv(self, file_path: Path) -> 'DataContainer':
        index = "hour"
        column = "level_of_loading"
        level_of_loading_data_dict = util.get_dict_from_csv(column, file_path, index)
        self.add_level_of_loading_data(level_of_loading_data_dict)
        return self

    def add_temperature_data(self, data: dict) -> 'DataContainer':
        """Add temperature data to the data container.

        Used when adding temperature to a data container from any source
        :param data: data containing temperature
        :type data: dict
        :return: DataContainer containing temperature data
        """
        self.temperature_data = data
        return self

    def add_temperature_data_from_csv(self, file_path: Path) -> 'DataContainer':
        """ Get temperature data from a path and raise verbose error if file is not found.

        :param file_path: csv path to temperature over hour of day
        :return: DataContainer containing cost parameters
        """
        index = "hour"
        column = "temperature"
        temperature_data_dict = util.get_dict_from_csv(column, file_path, index)
        self.add_temperature_data(temperature_data_dict)
        return self

    def add_cost_parameters(self, data: dict) -> 'DataContainer':
        """Add cost_parameters data to the data container. cost_parameters will be stored in the
        args instance

        Used when adding cost_parameters to a data container from any source
        :param data: data containing cost_parameters
        :type data: dict
        :return: DataContainer containing station data / electrified stations
        """
        self.cost_parameters_data = data
        return self

    def add_cost_parameters_from_json(self, file_path: Path) -> 'DataContainer':
        """ Get cost parameters data from a path and raise verbose error if file is not found.

        :param file_path: file path
        :return: DataContainer containing cost parameters
        """
        cost_parameters = self.get_json_from_file(file_path, "cost parameters")
        self.add_cost_parameters(cost_parameters)
        return self

    def add_stations(self, data: dict) -> 'DataContainer':
        """Add station data to the data container. Stations will be stored in the
        Schedule instance

        Used when adding stations to a data container from any source
        :param data: data containing stations
        :type data: dict
        :return: DataContainer containing station data / electrified stations
        """
        self.stations_data = data
        return self

    def add_stations_from_json(self, file_path: Path) -> 'DataContainer':
        """ Get electrified_stations data from a file_path
        :param file_path: file path
        :return: DataContainer containing station data / electrified stations

        """
        stations = self.get_json_from_file(file_path, "electrified stations")
        self.add_stations(stations)
        return self

    def add_vehicle_types(self, data: dict) -> None:
        """Add vehicle_type data to the data container. Vehicle_types will be stored in the
        Schedule instance

        Used when adding new vehicle types to a data container from any source
        :param data: data containing vehicle_types
        :type data: dict
        :return: DataContainer containing vehicle types
        """
        self.vehicle_types_data = data
        return self

    def add_vehicle_types_from_json(self, file_path: Path):
        """ Get vehicle_types from a json file_path
        :param file_path: file path
        :return: DataContainer containing vehicle types
        """
        vehicle_types = self.get_json_from_file(file_path, "vehicle types")
        self.add_vehicle_types(vehicle_types)
        return self

    @staticmethod
    def get_json_from_file(file_path: Path, data_type: str) -> any:
        """ Get json data from a file_path and raise verbose error if it fails.
        :raises FileNotFoundError: if file does not exist
        :param file_path: file path
        :param data_type: data type used for error description
        :return: json data
        """
        try:
            with open(file_path, encoding='utf-8') as f:
                return util.uncomment_json_file(f)
        except FileNotFoundError:
            raise FileNotFoundError(f"Path to {file_path} for {data_type} "
                                    "does not exist. Exiting...")

    def add_consumption_data_from_vehicle_type_linked_files(self):
        """ Add mileage data from files linked in the vehicle_types to the container.

        :return: DataContainer containing consumption data
        """
        assert self.vehicle_types_data, "No vehicle_type data in the data_container"
        mileages = list(get_values_from_nested_key("mileage", self.vehicle_types_data))
        mileages = list(filter(lambda x: isinstance(x, str) or isinstance(x, Path), mileages, ))
        # Cast mileages to set, since files only have to be read once
        for mileage_path in set(mileages):
            # creating the interpol function from csv file.
            delim = util.get_csv_delim(mileage_path)
            df = pd.read_csv(mileage_path, sep=delim)
            self.add_consumption_data(mileage_path, df)
        return self

    def add_consumption_data(self, data_name, df: pd.DataFrame) -> 'DataContainer':
        """Add consumption data to the data container. Consumption data will be used by the
        Consumption instance

        data_name must be equal to the mileage attribute in vehicle_types. E.g. Vehicle type
        '12m_opp' has a mileage of 'consumption_12m.csv' -> data_name ='consumption_12m.csv'
        :param data_name: name of the data, linked with vehicle_type
        :type data_name: str
        :param df: dataframe with consumption data and various expected columns
        :type df: pd.DataFrame
        :return: DatacContainer instance with added consumption data
        """
        missing_cols = [c for c in [INCLINE, T_AMB, LEVEL_OF_LOADING, SPEED, CONSUMPTION] if
                        c not in df.columns]
        assert not missing_cols, f"Consumption data is missing {', '.join(missing_cols)}"
        assert data_name not in self.consumption_data, f"{data_name} already exists in data"
        self.consumption_data[data_name] = df

        return self


def get_values_from_nested_key(key, data: dict) -> list:
    """Get all the values of the specified key in a nested dict

    :param key: key to find
    :param data: data to search through
    :yields: values of key in data
    """

    try:
        yield data[key]
    except KeyError:
        pass
    for value in data.values():
        if isinstance(value, dict):
            yield from get_values_from_nested_key(key, value)
