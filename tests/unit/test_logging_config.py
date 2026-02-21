"""
Unit tests for tesspy logging configuration helper.
"""

import logging

import pytest

from tesspy._logging import configure_logging

_OWNED_ATTR = "_tesspy_owned_handler"


def _remove_owned_handlers() -> None:
    logger = logging.getLogger("tesspy")
    for handler in list(logger.handlers):
        if getattr(handler, _OWNED_ATTR, False):
            logger.removeHandler(handler)


def test_configure_logging_adds_one_owned_handler():
    _remove_owned_handlers()

    logger = configure_logging("INFO")
    owned = [h for h in logger.handlers if getattr(h, _OWNED_ATTR, False)]

    assert len(owned) == 1
    assert logger.level == logging.INFO
    assert owned[0].level == logging.INFO

    _remove_owned_handlers()


def test_configure_logging_is_idempotent():
    _remove_owned_handlers()

    logger = configure_logging("INFO")
    logger = configure_logging("DEBUG")
    logger = configure_logging(logging.WARNING)
    owned = [h for h in logger.handlers if getattr(h, _OWNED_ATTR, False)]

    assert len(owned) == 1
    assert logger.level == logging.WARNING
    assert owned[0].level == logging.WARNING

    _remove_owned_handlers()


def test_configure_logging_rejects_invalid_level():
    _remove_owned_handlers()

    with pytest.raises(ValueError, match="Invalid logging level"):
        configure_logging("NOT_A_LEVEL")

