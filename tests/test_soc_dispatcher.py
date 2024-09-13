from copy import deepcopy
from datetime import timedelta
import sys
import pandas as pd
import pytest

from simba.simulate import pre_simulation
from tests.conftest import example_root
from simba import util
from simba.data_container import DataContainer


class TestSocDispatcher:
    def basic_run(self):
        """Returns a schedule, scenario and args after running SimBA.
        :return: schedule, scenario, args
        """
        # set the system variables to imitate the console call with the config argument.
        # first element has to be set to something or error is thrown
        sys.argv = ["foo", "--config", str(example_root / "simba.cfg")]
        args = util.get_args()
        args.seed = 5
        args.attach_vehicle_soc = True
        data_container = DataContainer().fill_with_args(args)
        sched, args = pre_simulation(args, data_container=data_container)

        # Copy an opportunity rotation twice, so dispatching can be tested
        assert sched.rotations["1"].charging_type == "oppb"
        sched.rotations["11"] = deepcopy(sched.rotations["1"])
        sched.rotations["12"] = deepcopy(sched.rotations["1"])
        sched.rotations["11"].id = "11"
        sched.rotations["12"].id = "12"

        # Mutate the first copy, so it ends later and has higher consumption
        sched.rotations["11"].trips[-1].arrival_time += timedelta(minutes=5)
        sched.rotations["11"].arrival_time += timedelta(minutes=5)

        sched.rotations["11"].trips[-1].distance += 10000
        # Mutate the second copy, so it starts later than "1" ends.
        # This way a prev. vehicle can be used.
        dt = sched.rotations["1"].arrival_time - sched.rotations["12"].departure_time + timedelta(
            minutes=20)
        for t in sched.rotations["12"].trips:
            t.arrival_time += dt
            t.departure_time += dt
        sched.rotations["12"].departure_time += dt
        sched.rotations["12"].arrival_time += dt

        # Copy a depot rotation, so a vehicle can be used again
        assert sched.rotations["2"].charging_type == "depb"

        sched.rotations["21"] = deepcopy(sched.rotations["2"])
        sched.rotations["21"].id = "21"
        dt = sched.rotations["2"].arrival_time - sched.rotations["21"].departure_time + timedelta(
            minutes=30)
        for t in sched.rotations["21"].trips:
            t.arrival_time += dt
            t.departure_time += dt
        sched.rotations["21"].departure_time += dt
        sched.rotations["21"].arrival_time += dt

        for v_type in sched.vehicle_types.values():
            for charge_type in v_type.values():
                charge_type["mileage"] = 1

        # calculate consumption of all trips
        sched.calculate_consumption()

        # Create soc dispatcher
        sched.init_soc_dispatcher(args)

        sched.assign_vehicles(args)
        scen = sched.run(args)

        for rot in sched.rotations.values():
            print(rot.id, rot.vehicle_id)

        return sched, scen, args

    @pytest.fixture
    def custom_vehicle_assignment(self):
        custom_vehicle_assignment = []

        custom_vehicle_assignment.append(dict(rot="41", v_id="AB_depb_1", soc=1))
        custom_vehicle_assignment.append(dict(rot="4", v_id="AB_depb_1", soc=1))
        custom_vehicle_assignment.append(dict(rot="3", v_id="AB_depb_2", soc=0.8))
        custom_vehicle_assignment.append(dict(rot="31", v_id="AB_depb_2", soc=0.8))
        custom_vehicle_assignment.append(dict(rot="21", v_id="AB_depb_3", soc=0.69))
        custom_vehicle_assignment.append(dict(rot="2", v_id="AB_depb_3", soc=1))
        custom_vehicle_assignment.append(dict(rot="1", v_id="AB_oppb_1", soc=1))
        custom_vehicle_assignment.append(dict(rot="11", v_id="AB_oppb_2", soc=0.6))
        custom_vehicle_assignment.append(dict(rot="12", v_id="AB_oppb_1", soc=0.945))
        return custom_vehicle_assignment

    def test_basic_dispatching(self, custom_vehicle_assignment):
        """Returns a schedule, scenario and args after running SimBA.
        :param custom_vehicle_assignment: list of assignments
        """
        sched, scen, args = self.basic_run()
        pd.DataFrame(scen.vehicle_socs).plot()

        sched.assign_vehicles_custom(custom_vehicle_assignment)
        scen = sched.run(args)

        pd.DataFrame(scen.vehicle_socs).plot()

    def test_basic_missing_rotation(self, custom_vehicle_assignment):
        """Test if missing a rotation throws an error
        :param custom_vehicle_assignment: list of assignments
        """
        sched, scen, args = self.basic_run()
        # delete data for a single rotation but keep the rotation_id
        missing_rot_id = custom_vehicle_assignment[-1]["rot"]
        del custom_vehicle_assignment[-1]

        # if data for a rotation is missing an error containing the rotation id should be raised
        with pytest.raises(Exception, match=missing_rot_id):
            sched.assign_vehicles_custom(custom_vehicle_assignment)
