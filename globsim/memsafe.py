import psutil
import logging
import sys

logger = logging.getLogger('globsim.interpolate.memsafe')


def require_safe_mem_usage(limit_percent: float, level=logging.DEBUG):
    """ Check if memory usage is safe. Kill globsim if not """
    mem = psutil.virtual_memory()
    logger.log(level, f"Memory usage: {mem.used / 1024**3:.2f} GB ({mem.percent}%)")
    safe = mem.percent < limit_percent
    if not safe:
        logger.critical(f"Memory use exceeds safe limit of {limit_percent}%. Exiting safely")
        sys.exit(1)


def require_memory_overhead(min_gb: float, level=logging.DEBUG):
    """
    Ensures a minimum 'breathing room' of available RAM.
    """
    mem = psutil.virtual_memory()
    available_gb = mem.available / (1024**3)

    logger.log(level, f"Memory check: {available_gb:.2f} GB available.")
    
    if available_gb < min_gb:
        logger.critical(
            f"Insufficient memory overhead! "
            f"Required: {min_gb} GB | Available: {available_gb:.2f} GB. "
            "Exiting to prevent system swap-death."
        )
        sys.exit(1)