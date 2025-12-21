""" Main """
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import dotenv

MODULE_DIR = Path(__file__).parent
WORK_DIR = MODULE_DIR.parent.parent
dotenv.load_dotenv(WORK_DIR/'.env')
LIB_DIR = WORK_DIR/'src/lib/py'

EXTRA_PATHS = [
    MODULE_DIR,
    LIB_DIR
]
for path in EXTRA_PATHS:
    path = str(path.resolve())
    if path not in sys.path:
        sys.path.append(path)

from subcommands import SUBCOMMANDS

def parse_args():
    """ parse args """
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()

    for command in SUBCOMMANDS:
        p_md = sp.add_parser(
            name=command.PROG_NAME,
            help=command.HELP,
            description=command.DESCRIPTION
        )
        for flag, params in command.ARGS.items():
            p_md.add_argument(flag, **params)
        func = command.func
        p_md.set_defaults(func=func)

    return p.parse_args()

def main():
    args = parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
