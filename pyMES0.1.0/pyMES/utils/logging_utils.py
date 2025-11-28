import logging
import sys

def setup_logging(
    log_file='model.log',
    level=logging.INFO,
    log_to_console=True,
    log_to_file=True,
    fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
):
    """
    Set up logging for the modeling framework.

    Args:
        log_file (str): Name of the log file.
        level (int): Logging level (e.g., logging.INFO, logging.DEBUG).
        log_to_console (bool): If True, also log to console (stdout).
        log_to_file (bool): If True, log to file.
        fmt (str): Format string for log messages.
    """
    logger = logging.getLogger()
    logger.setLevel(level)
    formatter = logging.Formatter(fmt)

    # Remove existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # Console handler
    if log_to_console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # File handler
    if log_to_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.info("Logging initialized.")

