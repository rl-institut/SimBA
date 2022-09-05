"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""
from ebus_toolbox import schedule
import json


def test_schedule_from_csv():
    with open('./data/examples/vehicle_types.json') as f:
        vehicle_types = json.load(f)
    empty_schedule = schedule.Schedule(vehicle_types, './data/examples/electrified_stations.json')
    assert empty_schedule.rotations == {}
