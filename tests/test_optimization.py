from argparse import Namespace
from copy import deepcopy

from simba import optimization
from tests.helpers import generate_basic_schedule


class TestOptimization:
    def test_prepare_trips(self):
        schedule = generate_basic_schedule()
        for r in schedule.rotations.values():
            r.set_charging_type("depb")
        trips, depot_trips = optimization.prepare_trips(schedule)

        # trips must not start/end in depot
        for trip_list in trips.values():
            assert len(trip_list) > 0
            for trip in trip_list:
                departure_station = schedule.stations.get(trip.departure_name)
                arrival_station = schedule.stations.get(trip.arrival_name)
                assert departure_station is None or departure_station["type"] != "deps"
                assert arrival_station is None or arrival_station["type"] != "deps"

        # depot_trips must only contain depot trips
        for station_trips in depot_trips.values():
            # station in depot_trips must have associated trips
            assert len(station_trips) > 0
            for trip in station_trips.values():
                departure_station = schedule.stations.get(trip.departure_name)
                arrival_station = schedule.stations.get(trip.arrival_name)
                if departure_station is None:
                    assert arrival_station["type"] == "deps"
                elif arrival_station is None:
                    assert departure_station["type"] == "deps"
                else:
                    assert arrival_station["type"] == "deps" or departure_station["type"] == "deps"

        # only depb rotations, no filter: trips must contain all trips except initial/final
        for rot_id, rot in schedule.rotations.items():
            assert len(trips[rot_id]) == len(rot.trips) - 2

        # Station-0 only depot
        for station_name in schedule.stations.keys():
            try:
                assert "Station-0" in depot_trips[station_name]
            except KeyError:
                # no trip from this station to depot
                pass

    def test_generate_depot_trip_data_dict(self):
        schedule = generate_basic_schedule()
        trip1 = schedule.rotations["1"].trips[0]
        trip2 = schedule.rotations["1"].trips[-1]
        trip2.distance = trip1.distance / 2
        depot_trips = {
            "opps-1": {
                "Station-0": trip1,  # is name of depot in example schedule
            },
            "opps-2": {
                "Station-0": trip1,
                "Station-3": trip2,
            }
        }
        DEFAULT_DISTANCE = 5  # km
        DEFAULT_SPEED = 30  # km/h
        # station known, only one trip
        trip_dict = optimization.generate_depot_trip_data_dict(
            "opps-1", depot_trips, DEFAULT_DISTANCE, DEFAULT_SPEED)
        assert trip_dict["name"] == "Station-0" and trip_dict["distance"] == trip1.distance

        # station known, two depot trips: use trip with lower distance
        trip_dict = optimization.generate_depot_trip_data_dict(
            "opps-2", depot_trips, DEFAULT_DISTANCE, DEFAULT_SPEED)
        assert trip_dict["name"] == "Station-3" and trip_dict["distance"] == trip2.distance

        # station unknown: use defaults for any depot
        trip_dict = optimization.generate_depot_trip_data_dict(
            "station-unknown", depot_trips, DEFAULT_DISTANCE, DEFAULT_SPEED)
        depot_station = schedule.stations[trip_dict["name"]]
        assert depot_station["type"] == "deps"
        assert trip_dict["distance"] == DEFAULT_DISTANCE * 1000

    def test_recombination(self):
        schedule = generate_basic_schedule()
        for r in schedule.rotations.values():
            r.set_charging_type("depb")
        args = Namespace(**({
            "desired_soc_deps": 1,
            "default_depot_distance": 5,
            "default_mean_speed": 30
        }))
        trips, depot_trips = optimization.prepare_trips(schedule)
        original_trips = deepcopy(trips)
        original_schedule = deepcopy(schedule)
        recombined_schedule = optimization.recombination(schedule, args, trips, depot_trips)
        # recombined rotations must be possible
        for rotation in recombined_schedule.rotations.values():
            vehicle_type = schedule.vehicle_types[rotation.vehicle_type][rotation.charging_type]
            if rotation.charging_type == "deps":
                assert rotation.consumption <= vehicle_type["capacity"]
        # all trips must be accounted for
        for rot_id, trip_list in original_trips.items():
            new_rot_name = f"{rot_id}_r"
            counter = 0
            trip_counter = 0
            while new_rot_name in recombined_schedule.rotations:
                trip_counter += len(recombined_schedule.rotations[new_rot_name].trips) - 2
                counter += 1
                new_rot_name = f"{rot_id}_r_{counter}"
            assert len(trip_list) == trip_counter

        # make one trip impossible
        schedule = original_schedule
        trips = deepcopy(original_trips)
        trips["1"][0].delta_soc = -2
        recombined_schedule = optimization.recombination(schedule, args, trips, depot_trips)
        # add Ein/Aussetzfahrt, subtract one impossible trip
        assert len(original_trips["1"]) + 2 - 1 == len(recombined_schedule.rotations["1_r"].trips)
