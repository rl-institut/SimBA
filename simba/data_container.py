"""Module to handle data access, by the varying SimBA modules."""
import argparse
import csv
import logging
import datetime
from pathlib import Path
from typing import Dict

import pandas as pd

from simba import util
from simba.ids import INCLINE, LEVEL_OF_LOADING, SPEED, T_AMB, CONSUMPTION


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

    def fill_with_args(self, args: argparse.Namespace):
        return self.fill_with_paths(args.vehicle_types_path,
                                    args.electrified_stations_path,
                                    args.cost_parameters_path,
                                    args.input_schedule,
                                    args.outside_temperature_over_day_path,
                                    args.level_of_loading_over_day_path,
                                    args.station_data_path,
                                    )

    def fill_with_paths(self,
                        vehicle_types_path,
                        electrified_stations_path,
                        cost_parameters_file,
                        trips_file_path,
                        outside_temperature_over_day_path,
                        level_of_loading_over_day_path,
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

        # Add station geo data
        if station_data_path is not None:
            self.add_station_geo_data_from_csv(station_data_path)

        # Add cost_parameters_data
        self.add_cost_parameters_from_json(cost_parameters_file)

        self.add_trip_data_from_csv(trips_file_path)
        return self

    def add_trip_data_from_csv(self, file_path: Path) -> 'DataContainer':
        """ Add trip data from csv file to DataContainer"""

        self.trip_data = []
        with open(file_path, 'r', encoding='utf-8') as trips_file:
            trip_reader = csv.DictReader(trips_file)
            for trip in trip_reader:
                trip_d = dict(trip)
                trip_d["arrival_time"] = datetime.datetime.fromisoformat(trip["arrival_time"])
                trip_d["departure_time"] = datetime.datetime.fromisoformat(trip["departure_time"])
                trip_d["level_of_loading"] = cast_float_or_none(trip.get("level_of_loading"))
                trip_d["temperature"] = cast_float_or_none(trip.get("temperature"))
                trip_d["distance"] = float(trip["distance"])
                self.trip_data.append(trip_d)




    def add_station_geo_data(self, data: dict) -> None:
        """Add station_geo data to the data container.

        Used when adding station_geo to a data container from any source
        :param data: data containing station_geo
        :type data: dict
        """
        self.station_geo_data = data

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
                    self.station_geo_data[str(row['Endhaltestelle'])]= {
                        "elevation": float(row['elevation']),
                        "lat": float(row.get('lat', 0)),
                        "long": float(row.get('long', 0)),
                    }
        except FileNotFoundError or KeyError:
            logging.warning("Warning: external csv file '{}' not found or not named properly "
                          "(Needed column names are 'Endhaltestelle' and 'elevation')".
                          format(file_path),
                          stacklevel=100)
        except ValueError:
            logging.warning("Warning: external csv file '{}' does not contain numeric "
                          "values in the column 'elevation'. Station data is discarded.".
                          format(file_path),
                          stacklevel=100)

        return self

    def add_level_of_loading_data(self, data: dict) -> None:
        """Add level_of_loading data to the data container.

        Used when adding level_of_loading to a data container from any source
        :param data: data containing level_of_loading
        :type data: dict
        """
        self.level_of_loading_data = data

    def add_level_of_loading_data_from_csv(self, file_path: Path) -> 'DataContainer':
        index = "hour"
        column = "level_of_loading"
        level_of_loading_data_dict = get_dict_from_csv(column, file_path, index)
        self.add_level_of_loading_data(level_of_loading_data_dict)
        return self

    def add_temperature_data(self, data: dict) -> None:
        """Add temperature data to the data container.

        Used when adding temperature to a data container from any source
        :param data: data containing temperature
        :type data: dict
        """
        self.temperature_data = data

    def add_temperature_data_from_csv(self, file_path: Path) -> 'DataContainer':
        index = "hour"
        column = "temperature"
        temperature_data_dict = get_dict_from_csv(column, file_path, index)
        self.add_temperature_data(temperature_data_dict)
        return self

    def add_cost_parameters(self, data: dict) -> None:
        """Add cost_parameters data to the data container. cost_parameters will be stored in the
        args instance

        Used when adding cost_parameters to a data container from any source
        :param data: data containing cost_parameters
        :type data: dict
        """
        self.cost_parameters_data = data

    def add_cost_parameters_from_json(self, file_path: Path) -> 'DataContainer':
        """ Get json data from a file_path"""
        try:
            with open(file_path, encoding='utf-8') as f:
                cost_parameters = util.uncomment_json_file(f)
        except FileNotFoundError:
            raise Exception(f"Path to cost parameters ({file_path}) "
                            "does not exist. Exiting...")
        self.add_cost_parameters(cost_parameters)
        return self

    def add_stations(self, data: dict) -> None:
        """Add station data to the data container. Stations will be stored in the
        Schedule instance

        Used when adding stations to a data container from any source
        :param data: data containing stations
        :type data: dict
        """
        self.stations_data = data

    def add_stations_from_json(self, file_path: Path) -> 'DataContainer':
        """ Get json data from a file_path"""
        try:
            with open(file_path, encoding='utf-8') as f:
                stations = util.uncomment_json_file(f)
        except FileNotFoundError:
            raise Exception(f"Path to electrified stations ({file_path}) "
                            "does not exist. Exiting...")
        self.add_stations(stations)
        return self

    def add_vehicle_types(self, data: dict) -> None:
        """Add vehicle_type data to the data container. Vehicle_types will be stored in the
        Schedule instance

        Used when adding new vehicle types to a data container from any source
        :param data: data containing vehicle_types
        :type data: dict
        """
        self.vehicle_types_data = data

    def add_vehicle_types_from_json(self, file_path: Path):
        """ Get json data from a file_path"""
        try:
            with open(file_path, encoding='utf-8') as f:
                vehicle_types = util.uncomment_json_file(f)
        except FileNotFoundError:
            raise Exception(f"Path to vehicle types ({file_path}) "
                            "does not exist. Exiting...")
        self.add_vehicle_types(vehicle_types)
        return self

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