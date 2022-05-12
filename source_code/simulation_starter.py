import logging
import logging.config

from opensc2_gui import OPENSC2_GUI

# from opensc2_gui_simpl import OPENSC2_GUI

logging.config.fileConfig(
    fname="logging_opensc2.conf", disable_existing_loggers=True
)

# Get the logger specified in the file, this will be the parent logger.
logger = logging.getLogger("opensc2Logger")

# make an instance of class OPENSC2_GUI (cdp, 12/2020)
gui = OPENSC2_GUI()
# Infinite loop of the main_window to start the program (cdp, 11/2020)
gui.main_window.mainloop()
