"""Module to handle data access, by the varying SimBA modules."""
from pathlib import Path
from typing import Dict

import pandas as pd

from simba import util
from simba.consumption import Consumption
from simba.ids import INCLINE, LEVEL_OF_LOADING, SPEED, T_AMB, CONSUMPTION


class DataContainer:
    def __init__(self):
        self.vehicle_types_data: Dict[str, any] = {}
        self.consumption_data: Dict[str, pd.DataFrame] = {}

    def add_vehicle_types(self, data: dict) -> None:
        """Add vehicle_type data to the data container. Vehicle_types will be stored in the
        Schedule instance

        :param data: data containing vehicle_types
        :type data: dict
        """
        self.vehicle_types_data = data

    def add_vehicle_types_from_json(self, file_path: Path):
        try:
            with open(file_path, encoding='utf-8') as f:
                vehicle_types = util.uncomment_json_file(f)
        except FileNotFoundError:
            raise Exception(f"Path to vehicle types ({file_path}) "
                            "does not exist. Exiting...")
        self.add_vehicle_types(vehicle_types)

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

    def add_consumption_data(self, data_name, df: pd.DataFrame) -> None:
        """Add consumption data to the data container. Consumption data will be used by the
        Consumption instance

        data_name must be equal to the mileage attribute in vehicle_types. E.g. Vehicle type
        '12m_opp' has a mileage of 'consumption_12m.csv' -> data_name ='consumption_12m.csv'
        :param data_name: name of the data, linked with vehicle_type
        :type data_name: str
        :param df: dataframe with consumption data and various expected columns
        :type df: pd.DataFrame"""
        for expected_col in [INCLINE, T_AMB, LEVEL_OF_LOADING, SPEED, CONSUMPTION]:
            assert expected_col in df.columns, f"Consumption data is missing {expected_col}"
        assert data_name not in self.consumption_data, f"{data_name} already exists in data"
        self.consumption_data[data_name] = df

    def to_consumption(self) -> Consumption:
        """Build a consumption instance from the stored data
        :returns: Consumption instance
        """
        # setup consumption calculator that can be accessed by all trips
        consumption = Consumption(self.vehicle_types_data)
        for name, df in self.consumption_data.items():
            consumption.set_consumption_interpolation(name, df)
        return consumption


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
