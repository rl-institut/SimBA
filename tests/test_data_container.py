import sys

import simba.data_container
from simba import util
from tests.conftest import example_root


class TestDataContainer:
    def test_get_values_from_nested_key(self):
        nested_dict = {"foo": {"bar": "baz1"}, "bob": {"bar": "baz2"}}
        gen = simba.data_container.get_values_from_nested_key("bar", nested_dict)
        assert next(gen) == "baz1"
        assert next(gen) == "baz2"

    def test_data_container(self):
        data_container = simba.data_container.DataContainer()
        sys.argv = ["foo", "--config", str(example_root / "configs/minimal.cfg")]
        args = util.get_args()
        data_container.fill_with_args(args)

    def test_resolve_csv_files(self, tmp_path):
        stations = {"Station-1": {
            "fixed_load": {"load_1": {"csv_file": "loads/load.csv"}},
            "price_csv": {"csv_file": str(tmp_path / "prices.csv")},
            "other": [{"csv_file": "pv.csv"}],
        }}
        simba.data_container.resolve_csv_files(stations, tmp_path)
        station = stations["Station-1"]
        assert station["fixed_load"]["load_1"]["csv_file"] == str(tmp_path / "loads/load.csv")
        assert station["price_csv"]["csv_file"] == str(tmp_path / "prices.csv")
        assert station["other"][0]["csv_file"] == str(tmp_path / "pv.csv")
        # no data_dir: unchanged
        stations = {"Station-1": {"price_csv": {"csv_file": "prices.csv"}}}
        simba.data_container.resolve_csv_files(stations, None)
        assert stations["Station-1"]["price_csv"]["csv_file"] == "prices.csv"
