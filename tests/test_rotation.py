import pytest

from simba.rotation import Rotation
from tests.helpers import generate_basic_schedule


def test_set_charging_type():
    s, _ = generate_basic_schedule()
    rot = list(s.rotations.values())[0]
    # set different mileages for different charging types to make sure consumption is properly
    # calculated
    for vehicle_key, vehicle_class in s.vehicle_types.items():
        for charging_key, vehicle in vehicle_class.items():
            if charging_key == "depb":
                vehicle["mileage"] = 10
            else:
                vehicle["mileage"] = 20

    # set charging type to oppb
    rot.set_charging_type('oppb')
    assert rot.charging_type == 'oppb'
    # save the consumption of this type of charger
    consumption_oppb = s.calculate_rotation_consumption(rot)

    # set charging type to depb
    rot.set_charging_type('depb')
    assert rot.charging_type == 'depb'
    # check that the consumption changed due to the change in charging type. The proper calculation
    # of consumption is tested in test_consumption
    consumption_depb = s.calculate_rotation_consumption(rot)
    assert consumption_depb * 2 == pytest.approx(consumption_oppb)

