"""Scaffold for CallVariant subcommand."""

from __future__ import annotations
import argparse


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register CallVariant arguments."""
    parser.description = 'Call moPepGen to generate the custom database FASTA.'


def run(args: argparse.Namespace) -> None:
    """Execute CallVariant subcommand."""
    _ = args
    raise NotImplementedError('CallVariant is not implemented yet.')
