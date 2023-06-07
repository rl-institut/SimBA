from pathlib import Path
import re
import subprocess

ROOT_PATH = Path(__file__).parent.parent
EXAMPLE_PATH = ROOT_PATH / "data/examples"


class TestExampleConfigs:

    def test_example_cfg(self, tmp_path):
        # copy cfg to tmp, adjust paths
        src = EXAMPLE_PATH / "ebus_toolbox.cfg"
        src_text = src.read_text()
        # provide path to input data
        src_text = src_text.replace("data/examples", str(EXAMPLE_PATH))
        # write output to tmp
        src_text = src_text.replace("data/sim_outputs", str(tmp_path))
        dst = tmp_path / "ebus_toolbox.cfg"
        # don't show plots. spaces are optional, so use regex
        src_text = re.sub(r"show_plots\s*=\s*true", "show_plots = false", src_text)
        dst.write_text(src_text)

        # call toolbox from shell
        assert subprocess.call([
            "python", "-m", "ebus_toolbox", "--config", dst
        ]) == 0
