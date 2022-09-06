"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""
from tests.mocks import mockSchedule


def test_set_charging_type():
    s = mockSchedule()
    rot = list(s.rotations.values())[0]

    # set charging type to oppb
    rot.set_charging_type('oppb')
    assert rot.charging_type == 'oppb'
    assert rot.consumption == 420

    # try setting charging type to depb
    # not possible since capacity < consumption
    rot.set_charging_type('depb')
    assert rot.charging_type == 'oppb'
    assert rot.consumption == 420

    # increase capacity of depb and try again
    s.vehicle_types["CKB_depb"]["capacity"] = 700
    rot.set_charging_type("depb")
    assert rot.charging_type == 'depb'
    assert rot.consumption == 630
