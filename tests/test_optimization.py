from argparse import Namespace

from simba import optimization
from tests.helpers import generate_basic_schedule


class TestOptimization:
    def test_prepare_trips(self):
        schedule = generate_basic_schedule()
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

    def test_generate_depot_trip(self):
        schedule = generate_basic_schedule()
        trip1 = schedule.rotations["1"].trips[0]
        trip2 = schedule.rotations["1"].trips[-1]
        trip2.distance = trip1.distance / 2
        depot_trips = {
            "station-0": {
                "depot-1": trip1,
            },
            "station-1": {
                "depot-1": trip1,
                "depot-2": trip2,
            }
        }
        # station known, only one trip
        trip_dict = optimization.generate_depot_trip("station-0", depot_trips)
        assert trip_dict["name"] == "depot-1" and trip_dict["distance"] == trip1.distance

        # station known, two depot trips: use trip with lower distance
        trip_dict = optimization.generate_depot_trip("station-1", depot_trips)
        assert trip_dict["name"] == "depot-2" and trip_dict["distance"] == trip2.distance

        # station unknown: use defaults for any depot
        trip_dict = optimization.generate_depot_trip("station-unknown", depot_trips)
        assert trip_dict["name"] in depot_trips.keys()
        assert trip_dict["distance"] == optimization.DEFAULT_DEPOT_DISTANCE * 1000

    def test_recombination(self):
        schedule = generate_basic_schedule()
        args = Namespace(**({"desired_soc_deps": 1}))
        trips, depot_trips = optimization.prepare_trips(schedule)
        recombined_schedule = optimization.recombination(schedule, args, trips, depot_trips)
        # recombined rotations must be possible
        for rotation in recombined_schedule.rotations.values():
            vehicle_type = schedule.vehicle_types[rotation.vehicle_type][rotation.charging_type]
            assert rotation.consumption <= vehicle_type["capacity"]
