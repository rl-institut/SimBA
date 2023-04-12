from tests.helpers import generate_basic_schedule


def test_set_charging_type():
    s = generate_basic_schedule()
    rot = list(s.rotations.values())[0]

    # set charging type to oppb
    rot.set_charging_type('oppb')
    assert rot.charging_type == 'oppb'
    assert rot.consumption == 420

    # set charging type to depb
    rot.set_charging_type('depb')
    assert rot.charging_type == 'depb'
    assert rot.consumption == 630
