"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""
from ebus_toolbox import schedule


def test_schedule_from_csv():
    s = schedule.Schedule.from_csv('./data/private_examples/trips_example-bvg_datetime.csv',
                                   './data/examples/vehicle_types.json')
    assert len(s.rotations) == 6
