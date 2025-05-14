from pathlib import Path
import pytest
import re
import subprocess

ROOT_PATH = Path(__file__).parent.parent
EXAMPLE_PATH = ROOT_PATH / "data/examples"


def adjust_paths(cfg_text, tmp_path):
    # adjust paths in config text to use absolute paths
    # provide path to input data
    cfg_text = cfg_text.replace("data/examples", str(EXAMPLE_PATH))
    # write output to tmp
    cfg_text = re.sub(r"output_path.+", f"output_path = {str(tmp_path.as_posix())}", cfg_text)
    # reduce horizon
    cfg_text = re.sub(r"HORIZON\":\d+", "HORIZON\":1", cfg_text)
    # don't show plots. spaces are optional, so use regex
    cfg_text = re.sub(r"show_plots\s*=\s*true", "show_plots = false", cfg_text)
    return cfg_text


class TestExampleConfigs:

    # create own test for every example config (can run tests in parallel)
    @pytest.mark.parametrize("cfg_path",  (EXAMPLE_PATH / "configs").glob("*.cfg"))
    def test_eval(self, cfg_path, tmp_path):
        # test single example config: does it run successfully?
        if cfg_path.stem == "minimal_pickle":
            # no example pickle file: skip this
            return
        # copy cfg to tmp, adjust paths
        src_text = cfg_path.read_text()
        dst = tmp_path / "simba.cfg"
        dst.write_text(adjust_paths(src_text, tmp_path))

        assert subprocess.call(["python", "-m", "simba", "--config", dst]) == 0, (
            cfg_path.stem + " failed")

    def test_required_files(self, tmp_path):
        # test if all necessary / expected files are copied into output folder
        # copy cfg to tmp, adjust paths
        src = EXAMPLE_PATH / "configs/extensive.cfg"
        dst = tmp_path / "simba.cfg"
        src_text = src.read_text()
        dst.write_text(adjust_paths(src_text, tmp_path))

        # call toolbox from shell
        subprocess.call(["python", "-m", "simba", "--config", dst])

        # make sure all required input files have been copied to output folder
        expected = [
            'simba.cfg',
            'program_version.txt',
            'trips_example.csv',
            'electrified_stations.json',
            'vehicle_types.json',

            'cost_params.json',
            'default_level_of_loading_over_day.csv',
            'default_temp_winter.csv',
            'energy_consumption_example.csv',
            'example_pv_feedin.csv',
            'example_external_load.csv',
            'price_timeseries.csv',
            'time_windows.json',
        ]
        input_dir = next(tmp_path.glob('*/input_data'))
        missing = list()
        for file_name in expected:
            if not (input_dir / file_name).exists():
                missing.append(file_name)
        if missing:
            raise Exception("Missing input files in output directory: " + ', '.join(missing))
