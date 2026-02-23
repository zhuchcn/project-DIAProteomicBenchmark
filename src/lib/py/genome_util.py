""" Util modules for genome handling. """
from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path
import dataclasses


if TYPE_CHECKING:
    from typing import Iterable

@dataclasses.dataclass
class GtfEntry:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attribute: dict[str, str | list[str]]


def parse_gtf(path:Path) -> Iterable[GtfEntry]:
    """ Parse GTF file to genome annotation. """
    with open(path, 'r') as fi:
        for line in fi:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attr_fields = fields[8].rstrip(';').split('; ')
            attr_dict = {}
            for attr in attr_fields:
                if attr == '':
                    continue
                key, value = attr.split(' ', 1)
                value = value.strip('"').strip("'")
                if key == 'tag':
                    attr_dict.setdefault(key, []).append(value)
                else:
                    attr_dict[key] = value
            yield GtfEntry(
                seqname=fields[0],
                source=fields[1],
                feature=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                score=fields[5],
                strand=fields[6],
                frame=fields[7],
                attribute=attr_dict
            )
