import matplotlib
from pathlib import Path

matplotlib.use("Agg")
test_root = Path(__file__).parent
file_root = test_root / "test_input_files"
example_root = Path(__file__).parent.parent / "data/examples"
