import logging

from rich.logging import RichHandler

logger = logging.getLogger("nanoflux")
logger.setLevel(logging.DEBUG)

# configure handler + formatter once
handler = RichHandler(markup=True)
formatter = logging.Formatter("%(message)s", datefmt="[%X]")
handler.setFormatter(formatter)

# avoid adding duplicate handlers on repeated imports
if not logger.handlers:
    logger.addHandler(handler)
