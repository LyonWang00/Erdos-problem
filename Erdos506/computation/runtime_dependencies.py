#!/usr/bin/env python3
"""Check and activate the exact OR-Tools version used by the verification."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version


ORTOOLS_VERSION = "9.15.6755"


def activate_ortools() -> str:
    """Require the pinned OR-Tools release and return its version."""
    try:
        installed = version("ortools")
    except PackageNotFoundError as exc:
        raise RuntimeError(
            "OR-Tools is absent; install computation/requirements.txt first."
        ) from exc
    if installed != ORTOOLS_VERSION:
        raise RuntimeError(
            f"OR-Tools {ORTOOLS_VERSION} is required, but {installed} is installed."
        )
    from ortools.sat.python import cp_model  # noqa: F401

    return ORTOOLS_VERSION
