import re
import logging
# ---------------- logging setup ----------------
logger = logging.getLogger("hyperedges")
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler()
_handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
logger.addHandler(_handler)

def set_debug(enabled: bool = False) -> None:
    """Turn verbose debugging on/off globally."""
    logger.setLevel(logging.DEBUG if enabled else logging.INFO)

# ---------------- parsing / counting ----------------
DIVIDER_RE = re.compile(r'^\s*%+\s*$')                   # %%%... separators
LAYER_LINE_RE = re.compile(r'^\s*(\d+)\s*:\s*(.+?)\s*$') # "12: 82, 71, ..."

N_LAYERS = 60  # fixed number of layers
K = 8  # size of the expert groups
MIN_FREQ = 1 # minimum  frequency required to store the hyperedge