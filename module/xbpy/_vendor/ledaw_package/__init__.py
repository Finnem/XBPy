import os
import glob
import importlib
import warnings

# Get the directory where this __init__.py is located
module_dir = os.path.dirname(__file__)

# List all the Python files in this directory (except __init__.py)
module_files = glob.glob(os.path.join(module_dir, "*.py"))

# GUI/Qt-only modules.  They depend on PySide6, which is a heavy optional
# dependency only useful for the LEDAW desktop GUI.  Headless callers
# (batch CovaLED extraction over many ORCA outputs) never need them, so
# skip them up front to avoid a hard ImportError when PySide6 is missing.
_GUI_ONLY_MODULES = {"gui_engine", "job_engine"}

# Import all functions from all modules.  Skip GUI-only modules (require
# PySide6); for everything else, surface optional-dep failures as warnings
# so a single missing extra cannot kill the whole vendored package.
for module_file in module_files:
    module_name = os.path.basename(module_file)[:-3]  # remove ".py" from the name
    if module_name == "__init__":
        continue
    if module_name in _GUI_ONLY_MODULES:
        continue
    try:
        module = importlib.import_module(f'.{module_name}', package=__name__)
    except ImportError as e:
        warnings.warn(
            f"xbpy._vendor.ledaw_package: skipping optional submodule "
            f"'{module_name}' because of missing dependency: {e}",
            RuntimeWarning,
            stacklevel=2,
        )
        continue
    globals().update({name: getattr(module, name) for name in dir(module) if not name.startswith('_')})