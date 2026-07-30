import logging
import logging.config
import os
from pathlib import Path


# The returned handle must remain alive for the whole application lifetime.
_UMFPACK_DLL_HANDLE = None


def configure_umfpack_runtime():
    """Register the SuiteSparse DLL directory on Windows."""
    global _UMFPACK_DLL_HANDLE

    if os.name != "nt":
        return

    dll_dir = os.environ.get("OPENSC2_UMFPACK_DLL_DIR")

    if not dll_dir:
        return

    dll_path = Path(dll_dir)

    if not dll_path.is_dir():
        raise RuntimeError(
            "OPENSC2_UMFPACK_DLL_DIR does not point to a valid directory: "
            f"{dll_path}"
        )

    _UMFPACK_DLL_HANDLE = os.add_dll_directory(str(dll_path))


# This must run before importing OPENSC2 modules that may load UMFPACK.
configure_umfpack_runtime()

from scipy.sparse.linalg import use_solver

solver = os.environ.get("OPENSC2_SPARSE_SOLVER", "superlu").lower()

if solver == "umfpack":
    use_solver(useUmfpack=True)
elif solver == "superlu":
    use_solver(useUmfpack=False)
else:
    raise ValueError(
        "OPENSC2_SPARSE_SOLVER must be either 'umfpack' or 'superlu', "
        f"not {solver!r}."
    )

print(f"OPENSC2 sparse solver requested: {solver}")

from opensc2_gui import OPENSC2_GUI

# from opensc2_gui_simpl import OPENSC2_GUI


logging.config.fileConfig(
    fname="logging_opensc2.conf",
    disable_existing_loggers=True,
)

# Get the logger specified in the file, this will be the parent logger.
logger = logging.getLogger("opensc2Logger")

# Make an instance of class OPENSC2_GUI (cdp, 12/2020).
gui = OPENSC2_GUI()

# Infinite loop of the main_window to start the program (cdp, 11/2020).
gui.main_window.mainloop()