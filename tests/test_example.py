from pathlib import Path
import re
import subprocess

ROOT_PATH = Path(__file__).parent.parent
EXAMPLE_PATH = ROOT_PATH / "data/examples"


class TestExampleConfigs:

    def test_example_cfg(self, tmp_path):
        # copy cfg to tmp, adjust paths
        src = EXAMPLE_PATH / "simba.cfg"
        src_text = src.read_text()
        # provide path to input data
        src_text = src_text.replace("data/examples", str(EXAMPLE_PATH))
        # write output to tmp
        src_text = re.sub(r"output_path.+", f"output_path = {str(tmp_path.as_posix())}",
                          src_text)
        dst = tmp_path / "simba.cfg"
        # don't show plots. spaces are optional, so use regex
        src_text = re.sub(r"show_plots\s*=\s*true", "show_plots = false", src_text)
        dst.write_text(src_text)

        # call toolbox from shell
        assert subprocess.call([
            "python", "-m", "simba", "--config", dst
        ]) == 0

        # make sure all required files have been copied to output folder
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
