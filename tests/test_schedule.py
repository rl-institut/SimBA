"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""
from datetime import timedelta

from tests.helpers import generate_basic_schedule
from ebus_toolbox.schedule import Schedule


def test_schedule_from_csv():
    schedule = generate_basic_schedule()
    assert len(schedule.rotations) == 1


def test_consistency():
    sched = generate_basic_schedule()
    # check if no error is thrown in the basic case
    assert len(Schedule.check_consistency(sched)) == 0

    error = "Trip time is negative"
    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    faulty_trip = faulty_rot.trips[0]
    # create error through moving trip arrival 1 day before departure
    faulty_trip.arrival_time = faulty_trip.departure_time - timedelta(days=1)
    assert Schedule.check_consistency(sched)["1"] == error

    error = "Break time is negative"
    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    faulty_trip = faulty_rot.trips[1]
    # create error through moving trip departure before last arrival
    faulty_trip.departure_time = faulty_rot.trips[0].arrival_time-timedelta(minutes=1)
    assert Schedule.check_consistency(sched)["1"] == error

    error = "Trips are not sequential"
    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    faulty_rot.trips[1].arrival_name = "foo"
    faulty_rot.trips[0].departure_name = "bar"
    assert Schedule.check_consistency(sched)["1"] == error

    error = "Start and end of rotation differ"
    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    departure_trip = list(faulty_rot.trips)[0]
    departure_trip.departure_name = "foo"
    faulty_rot.trips[-1].arrival_name = "bar"
    assert Schedule.check_consistency(sched)["1"] == error

    error = "Rotation data differs from trips data"

    # check arrival data in rotation
    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    faulty_rot.trips[-1].arrival_name = "foo"
    faulty_rot.trips[0].departure_name = "foo"
    faulty_rot.arrival_name = "bar"
    assert Schedule.check_consistency(sched)["1"] == error

    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    arrival_trip = faulty_rot.trips[-1]
    faulty_rot.arrival_time = arrival_trip.arrival_time - timedelta(minutes=1)
    assert Schedule.check_consistency(sched)["1"] == error

    # check departure data in rotation
    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    faulty_rot.trips[-1].arrival_name = "foo"
    faulty_rot.trips[0].departure_name = "foo"
    faulty_rot.departure_name = "bar"
    assert Schedule.check_consistency(sched)["1"] == error

    sched = generate_basic_schedule()
    faulty_rot = list(sched.rotations.values())[0]
    departure_trip = faulty_rot.trips[0]
    faulty_rot.departure_time = departure_trip.departure_time - timedelta(minutes=1)
    assert Schedule.check_consistency(sched)["1"] == error
