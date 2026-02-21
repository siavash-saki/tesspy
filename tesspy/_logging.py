"""
Shared logging utilities for tesspy.
"""

import logging
from typing import TextIO

_TESSPY_HANDLER_ATTR = "_tesspy_owned_handler"


def _resolve_level(level: int | str) -> int:
    """Resolve a numeric logging level from int or named string."""
    if isinstance(level, int):
        return level

    if isinstance(level, str):
        name = level.strip().upper()
        mapped = logging.getLevelNamesMapping().get(name)
        if mapped is None:
            raise ValueError(f"Invalid logging level: {level}")
        return mapped

    raise TypeError("level must be an int or string logging level name")


def configure_logging(
    level: int | str = "INFO",
    *,
    stream: TextIO | None = None,
    fmt: str = "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt: str = "%Y-%m-%d %H:%M:%S",
) -> logging.Logger:
    """
    Configure the package logger for tesspy and return it.

    This function is idempotent: repeated calls do not stack duplicate package-
    owned handlers.
    """
    resolved_level = _resolve_level(level)
    package_logger = logging.getLogger("tesspy")

    owned_handlers = [
        h for h in package_logger.handlers if getattr(h, _TESSPY_HANDLER_ATTR, False)
    ]
    primary_handler: logging.StreamHandler | None = None

    if owned_handlers:
        primary_handler = owned_handlers[0]
        for extra in owned_handlers[1:]:
            package_logger.removeHandler(extra)

    if primary_handler is None:
        primary_handler = logging.StreamHandler(stream=stream)
        setattr(primary_handler, _TESSPY_HANDLER_ATTR, True)
        package_logger.addHandler(primary_handler)
    elif stream is not None and primary_handler.stream is not stream:
        primary_handler.setStream(stream)

    formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
    primary_handler.setFormatter(formatter)
    primary_handler.setLevel(resolved_level)

    package_logger.setLevel(resolved_level)
    package_logger.propagate = False

    return package_logger


def log_progress(
    logger: logging.Logger,
    verbose: bool,
    msg: str,
    *args: object,
    level: int = logging.INFO,
    **kwargs: object,
) -> None:
    """Log a progress message only when verbose mode is enabled."""
    if verbose:
        logger.log(level, msg, *args, **kwargs)
