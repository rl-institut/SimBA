from simba.rotation import Rotation
from tests.helpers import generate_basic_schedule


def test_set_charging_type():
    s = generate_basic_schedule()
    rot = list(s.rotations.values())[0]

    # set different mileages for different charging types to make sure consumption is properly
    # calculated
    for vehicle_key, vehicle_class in s.vehicle_types.items():
        for charging_key, vehicle in vehicle_class.items():
            if charging_key == "depb":
                vehicle["mileage"] = 10
            else:
                vehicle["mileage"] = 20
    rot.consumption = rot.calculate_consumption()

    # set charging type to oppb
    rot.set_charging_type('oppb')
    assert rot.charging_type == 'oppb'
    # save the consumption of this type of charger
    consumption_oppb = rot.consumption

    # set charging type to depb
    rot.set_charging_type('depb')
    assert rot.charging_type == 'depb'
    # check that the consumption changed due to the change in charging type. The proper calculation
    # of consumption is tested in test_consumption
    assert rot.consumption != consumption_oppb


def test_calculate_consumption():
    s = generate_basic_schedule()
    vt = next(iter(s.vehicle_types))

    r = Rotation("my_rot", vt, s)
    # consumption without trips
    assert r.calculate_consumption() == 0

    some_rot = next(iter(s.rotations.values()))
    first_trip = some_rot.trips[0]
    first_trip.charging_type = "depb"
    del first_trip.rotation
    r.add_trip(vars(first_trip))
    r.trips[0].consumption = None
    r.calculate_consumption()
    assert r.trips[0].consumption is not None
    second_trip = some_rot.trips[1]
    del second_trip.rotation
    r.add_trip(vars(second_trip))
    for trip in r.trips:
        trip.consumption = None
    r.calculate_consumption()
    for trip in r.trips:
        assert trip.consumption is not None
