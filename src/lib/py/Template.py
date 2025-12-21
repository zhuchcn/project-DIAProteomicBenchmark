from __future__ import annotations
from typing import Any
import dataclasses

class SubCommand:
    ARGS:dict[str, dict[str,Any]] = {}
    PROG_NAME = None
    HELP = None
    DESCRIPTION = None

    def get_argument_type(self):
        return dataclasses.make_dataclass(
            'ArgumentType',
            [
                (
                    attrs.get('dest') or arg.lstrip('-').replace('-', '_'),
                    attrs['type']
                )
                for arg, attrs in self.args.items()
            ]
        )
