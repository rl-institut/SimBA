import pytest

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


def test_min_standing_time_charging_curve():
    s, _ = generate_basic_schedule()
    # rotation ending at a depot
    rot = next(r for r in s.rotations.values()
               if s.stations.get(r.arrival_name, {}).get("type") == "deps")
    ct = rot.charging_type
    v_info = s.vehicle_types[rot.vehicle_type][ct]
    charge_power = s.stations[rot.arrival_name].get(
        f"cs_power_deps_{ct}", vars(s)[f"cs_power_deps_{ct}"])
    setattr(s, f"min_recharge_deps_{ct}", 1)
    # consumption must be below capacity to recharge consumption completely
    v_info["capacity"] = max(v_info["capacity"], 2 * rot.consumption)
    # standing time with constant charging power and without losses
    constant_time = rot.consumption / charge_power

    # constant charging curve without losses: same standing time
    v_info["charging_curve"] = [[0, charge_power], [1, charge_power]]
    v_info["battery_efficiency"] = 1
    s.recharge_curves = s.get_charge_curves(s.get_deps_charge_levels(), 1, time_step_min=1)
    assert rot.min_standing_time == pytest.approx(constant_time, rel=1e-2)

    # charging losses: longer standing time
    v_info["battery_efficiency"] = 0.9
    s.recharge_curves = s.get_charge_curves(s.get_deps_charge_levels(), 1, time_step_min=1)
    assert rot.min_standing_time == pytest.approx(constant_time / 0.9, rel=1e-2)

    # charging power decreasing at high SoC: even longer standing time
    v_info["charging_curve"] = [[0, charge_power], [0.8, charge_power], [1, charge_power / 10]]
    s.recharge_curves = s.get_charge_curves(s.get_deps_charge_levels(), 1, time_step_min=1)
    assert rot.min_standing_time > constant_time / 0.9

    # recharge limited by min_recharge_deps: shorter standing time
    setattr(s, f"min_recharge_deps_{ct}", 0.1)
    assert rot.min_standing_time < constant_time
