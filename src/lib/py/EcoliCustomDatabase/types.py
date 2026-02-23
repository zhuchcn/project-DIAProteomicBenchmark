"""Shared types for the E. coli custom database package."""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class SubcommandSpec:
    """Static metadata for a subcommand."""

    name: str
    help: str
