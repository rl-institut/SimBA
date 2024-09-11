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
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        data_container.fill_with_args(args)
