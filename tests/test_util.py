from datetime import datetime
import json
from ebus_toolbox import util
from tests.test_schedule import TestSchedule


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
        schedule, scenario = TestSchedule().basic_run()
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
