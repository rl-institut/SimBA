from datetime import datetime, timedelta
import json
from simba import util
from tests.test_schedule import BasicSchedule


class TestUtil:
    def test_uncomment_json_file(self, tmp_path):
        p = tmp_path / "test.json"
        with open(p, 'w+', encoding='utf-8') as f:
            # no comment
            json.dump({"test1": 1}, f)
            f.seek(0)
            assert util.uncomment_json_file(f)["test1"] == 1

        with open(p, 'w+', encoding='utf-8') as f:
            # multiline, no comment
            f.seek(0)
            f.write('{\n"test1.1": 1.1,\n"test1.2": 1.2\n}')
            f.seek(0)
            assert sum(util.uncomment_json_file(f).values()) == 2.3

        with open(p, 'w+', encoding='utf-8') as f:
            # single char comment in single line
            f.write('# this is a comment\n{"test2": 2}\n#and another comment')
            f.seek(0)
            assert util.uncomment_json_file(f, '#')["test2"] == 2

        with open(p, 'w+', encoding='utf-8') as f:
            # single char comment at end of line
            f.write('{"test3": 3} % this is a comment')
            f.seek(0)
            assert util.uncomment_json_file(f, '%')["test3"] == 3

        with open(p, 'w+', encoding='utf-8') as f:
            # multi-char comment in own line
            f.write('// {"test4": "comment"}\n{"test4": 4}')
            f.seek(0)
            assert util.uncomment_json_file(f, '//')["test4"] == 4

        with open(p, 'w+', encoding='utf-8') as f:
            # multi-char comment at end of line
            f.write('{"test5": 5} """ {test5: "comment"} """')
            f.seek(0)
            assert util.uncomment_json_file(f, '"""')["test5"] == 5

    def test_get_git_revision_hash(self):
        git_hash = util.get_git_revision_hash()
        assert type(git_hash) is str

    def test_get_buffer_time(self):
        schedule, scenario, _ = BasicSchedule().basic_run()
        trip = next(iter(schedule.rotations.values())).trips.pop(0)
        util.get_buffer_time(trip)
        buffer_time = {"10-22": 2,
                       "22-6": 3,
                       "else": 1
                       }

        trip.arrival_time = datetime(year=2023, month=1, day=1, hour=10)
        assert util.get_buffer_time(trip, default=buffer_time) == 2

        trip.arrival_time = datetime(year=2023, month=1, day=1, hour=21, minute=59, second=59)
        assert util.get_buffer_time(trip, default=buffer_time) == 2

        trip.arrival_time = datetime(year=2023, month=1, day=1, hour=22)
        assert util.get_buffer_time(trip, default=buffer_time) == 3

        trip.arrival_time = datetime(year=2023, month=1, day=1, hour=22, second=1)
        assert util.get_buffer_time(trip, default=buffer_time) == 3

        trip.arrival_time = datetime(year=2023, month=1, day=2, hour=0, second=1)
        assert util.get_buffer_time(trip, default=buffer_time) == 3

        trip.arrival_time = datetime(year=2023, month=1, day=2, hour=6, second=1)
        assert util.get_buffer_time(trip, default=buffer_time) == 1

    def test_get_mean_from_hourly_dict(self):
        # Dict with values 0-23
        hourly_dict = {x: x for x in range(0, 24)}

        # Range off a whole day should give 23/2 = 11.5
        start = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
        end = (start + timedelta(days=1))
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 11.5

        # First hour of the day has value of 0
        end = (start + timedelta(hours=1))
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 0

        # Second hour of the day has value of 1. Mean of [0,1] is 0.5
        end = (start + timedelta(hours=2))
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 0.5

        # 1 minute of the first hour and second hour with values of 0 and 1 should result in 0.5
        start = start.replace(hour=0, minute=59, second=0, microsecond=0)
        end = start.replace(hour=1, minute=1, second=0, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 0.5

        # 1 minute of the first hour and 2 minutes of the second hour with values of 0 and 1
        # should result in (0+1+1)/3 = 2/3
        start = start.replace(hour=0, minute=59, second=0, microsecond=0)
        end = start.replace(hour=1, minute=2, second=0, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 2 / 3

        # 1 minute of the first hour and 59 minutes of the second hour with values of 0 and 1
        # should result in (1*0+59*1)/60 = 59/60
        start = start.replace(hour=0, minute=59, second=0, microsecond=0)
        end = start.replace(hour=1, minute=59, second=0, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 59 / 60

        # Daterange which covers part of the first hour should result in the value of the first hour
        start = start.replace(hour=0, minute=0, second=0, microsecond=0)
        end = start.replace(hour=0, minute=10, second=0, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 0

        # No duration should work. 1 hour results in the lookup 1.
        start = start.replace(hour=1, minute=0, second=0, microsecond=0)
        end = start.replace(hour=1, minute=0, second=0, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 1

        # Check lookup of second hour. Make sure 1 is returned.
        start = start.replace(hour=1, minute=0, second=0, microsecond=0)
        end = start.replace(hour=1, minute=1, second=0, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 1

        # Check if seconds are properly handled. 1s of 0 and 2s of 1 -> 2/3
        start = start.replace(hour=0, minute=59, second=59, microsecond=0)
        end = start.replace(hour=1, minute=0, second=2, microsecond=0)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 2 / 3

        # Check if seconds are properly handled. 1ms of 0 and 2ms of 1 -> 2/3
        start = start.replace(hour=0, minute=59, second=59, microsecond=int(1e6-1))
        end = start.replace(hour=1, minute=0, second=0, microsecond=2)
        assert util.get_mean_from_hourly_dict(hourly_dict, start, end) == 2 / 3
