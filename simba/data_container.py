"""Module to handle data access, by the varying SimBA modules."""
import argparse
import csv
import logging
import datetime
from pathlib import Path
from typing import Dict

import pandas as pd

from simba import util
from simba.ids import INCLINE, LEVEL_OF_LOADING, SPEED, T_AMB, TEMPERATURE, CONSUMPTION


class DataContainer:
    def __init__(self):
        self.vehicle_types_data: Dict[str, any] = {}
        self.consumption_data: Dict[str, pd.DataFrame] = {}
        self.temperature_data: Dict[int, float] = {}
        self.level_of_loading_data: Dict[int, float] = {}
        self.stations_data: Dict[str, dict] = {}
        self.cost_parameters_data: Dict[str, dict] = {}
        self.station_geo_data: Dict[str, dict] = {}

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
                trip_d[LEVEL_OF_LOADING] = cast_float_or_none(trip.get(LEVEL_OF_LOADING))
                trip_d[TEMPERATURE] = cast_float_or_none(trip.get(TEMPERATURE))
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
        try:
            with open(file_path, "r", encoding='utf-8') as f:
                delim = util.get_csv_delim(file_path)
                reader = csv.DictReader(f, delimiter=delim)
                for row in reader:
                    self.station_geo_data[str(row['Endhaltestelle'])] = {
                        "elevation": float(row['elevation']),
                        "lat": float(row.get('lat', 0)),
                        "long": float(row.get('long', 0)),
                    }
        except (FileNotFoundError, KeyError):
            logging.warning("External csv file '{}' not found or not named properly "
                          "(Needed column names are 'Endhaltestelle' and 'elevation')".
                          format(file_path),
                          stacklevel=100)
        except ValueError:
            logging.warning("External csv file '{}' should only contain numeric data".
                          format(file_path),
                          stacklevel=100)
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
        level_of_loading_data_dict = get_dict_from_csv(column, file_path, index)
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
        index = "hour"
        column = "temperature"
        temperature_data_dict = get_dict_from_csv(column, file_path, index)
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
            raise FileNotFoundError(f"Path to {data_type} ({file_path}) "
                                    "does not exist. Exiting...")

    def add_consumption_data_from_vehicle_type_linked_files(self):
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

        for expected_col in [INCLINE, T_AMB, LEVEL_OF_LOADING, SPEED, CONSUMPTION]:
            assert expected_col in df.columns, f"Consumption data is missing {expected_col}"
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


def get_dict_from_csv(column, file_path, index):
    output = dict()
    with open(file_path, "r") as f:
        delim = util.get_csv_delim(file_path)
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            output[float(row[index])] = float(row[column])
    return output


def cast_float_or_none(val: any) -> any:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None
