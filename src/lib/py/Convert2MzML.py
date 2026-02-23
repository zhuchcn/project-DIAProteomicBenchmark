""" Convert raw data files to mzML format """
from __future__ import annotations
from typing import TYPE_CHECKING
import os
import csv
from pathlib import Path
import subprocess as sp
from Template import SubCommand
from Common import setup_logger


if TYPE_CHECKING:
    from typing import Iterable

_PROG_NAME = 'Convert2MzML'
_HELP = 'Convert raw data files to mzML format'
_DESCRIPTION = 'Convert raw data files to mzML format using msconvert from ProteoWizard toolkit.'
LOGGER = setup_logger(_PROG_NAME)

class Convert2MzML(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--input-file': {
            'type': Path,
            'required': True,
            'help': 'Path to the input file containing list of raw data files to convert.'
        },
        '--output-dir': {
            'type': str,
            'required': True,
            'help': 'Directory to save the converted mzML files.'
        },
        '--dia': {
            'action': 'store_true',
            'help': 'Whether the input files are DIA data. If set, additional filters will be applied during conversion.'
        }
    }

    @staticmethod
    def func(args):
        input_file = args.input_file
        output_dir = args.output_dir

        LOGGER.info(f"Reading input file: {input_file}")
        for sample_id, ms_file in _collect_ms_files(input_file, args.dia):
            input_dir = str(ms_file.parent)
            LOGGER.info("Converting %s", sample_id)
            cwd = os.getcwd()
            cmd = f"""
                docker run \\
                    --rm \\
                    --platform linux/amd64 \\
                    -w {cwd} \\
                    -v {cwd}:{cwd} \\
                    -v {input_dir}:{input_dir} \\
                    -v {output_dir}:{output_dir} \\
                    chambm/pwiz-skyline-i-agree-to-the-vendor-licenses \\
                        wine msconvert \\
                        {ms_file} \\
                        -o {output_dir} \\
                        --mzML \\
                        --64 \\
                        --filter "peakPicking true 1-" \\
                        --filter "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>"
            """
            LOGGER.info("Running command: %s", cmd)
            sp.run(cmd, check=True, shell=True)

        LOGGER.info("Conversion completed.")

def _collect_ms_files(input_file: Path, is_dia:bool) -> Iterable[tuple[str, Path]]:
    """Collect raw data files from the input file."""
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if is_dia:
                yield row['sample_id'], Path(row['wiff'])
            else:
                yield row['sample_id'], Path(row['raw'])
