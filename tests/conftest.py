import matplotlib
# Force the non-interactive Agg backend for the whole test session. Without
# this, generate_plots() (simba/report.py) calls plt.ion(), which makes
# matplotlib pick an interactive GUI backend (TkAgg on Windows) on first use.
# That spins up real Tk/ttk widgets, which can fail with
# "_tkinter.TclError: Can't find a usable tk.tcl" depending on the local
# Tcl/Tk setup, and shouldn't be creating GUI windows during automated tests
# anyway. Must be set before matplotlib.pyplot is imported anywhere, so this
# has to live in conftest.py, which pytest always loads first.
matplotlib.use("Agg")

from pathlib import Path

test_root = Path(__file__).parent
file_root = test_root / "test_input_files"
example_root = Path(__file__).parent.parent / "data/examples"
