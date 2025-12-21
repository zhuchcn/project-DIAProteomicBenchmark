""" Convert raw data files to mzML format """
from __future__ import annotations
import os
import csv
from pathlib import Path
import subprocess as sp
from Template import SubCommand
from Common import setup_logger


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
        }
    }

    @staticmethod
    def func(args):
        input_file = args.input_file
        output_dir = args.output_dir

        LOGGER.info(f"Reading input file: {input_file}")
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row['sample_id']
                wiff = row['wiff']
                input_dir = str(Path(wiff).parent)
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
                            {wiff} \\
                            -o {output_dir} \\
                            --mzML \\
                            --64 \\
                            --filter "peakPicking true 1-" \\
                            --filter "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>"
                """
                LOGGER.info("Running command: %s", cmd)
                sp.run(cmd, check=True, shell=True)

        LOGGER.info("Conversion completed.")
