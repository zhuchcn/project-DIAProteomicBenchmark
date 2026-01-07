""" Run FragPipe with a config and manifest via Docker """
from __future__ import annotations
import csv
import shlex
from pathlib import Path
import subprocess as sp

from Template import SubCommand
from Common import (
    setup_logger,
    FRAGPIPE_DOCKER_IMAGE,
    FRAGPIPE_EXE,
    collect_mount_points,
    build_docker_run_cmd,
    resolve_file,
    resolve_directory,
)


_PROG_NAME = 'RunFragPipe'
_HELP = 'Run FragPipe via the official Docker image'
_DESCRIPTION = (
    'Run FragPipe in command-line mode against a manifest of mzML files using a provided config file.'
)
LOGGER = setup_logger(_PROG_NAME)
MANIFEST_COLUMNS = ('sample_id', 'mzml_path')


class RunFragPipe(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--workflow': {
            'type': Path,
            'required': True,
            'help': 'Path to the FragPipe workflow file to run.',
        },
        '--manifest': {
            'type': Path,
            'required': True,
            'help': (
                'Manifest describing the mzML inputs. Required columns: ' +
                ', '.join(MANIFEST_COLUMNS)
            )
        },
        '--config-tools-folder': {
            'type': Path,
            'required': True,
            'help': 'Path to the directory containing MSFragger, IonQuant and diaTracer.'
        },
        '--output-dir': {
            'type': Path,
            'required': True,
            'help': 'Directory to save FragPipe output files.'
        },
        '--docker-image': {
            'type': str,
            'required': False,
            'default': FRAGPIPE_DOCKER_IMAGE,
            'help': 'Docker image that provides FragPipe and Philosopher tools.'
        },
        '--dry-run': {
            'action': 'store_true',
            'help': 'If set, only run FragPipe in dry-run mode.'
        },
        '--ram': {
            'type': int,
            'default': 0,
            'help': 'Amount of RAM (in GB) to allocate to FragPipe. Default is'
            ' to let FragPipe decide based on available system memory.'
        },
        '--threads': {
            'type': int,
            'default': 0,
            'help': 'Number of CPU threads to allocate to FragPipe. Default is core number - 1.'
        },
    }

    @staticmethod
    def func(args):
        workflow = resolve_file(args.workflow, 'workflow file')
        manifest = resolve_file(args.manifest, 'manifest file')
        mzml_files = _read_manifest(manifest)
        config_tools_folder = resolve_directory(args.config_tools_folder, 'config tools folder')
        output_dir = resolve_directory(args.output_dir, 'output directory')

        workflow_data = _read_workflow(workflow)
        database_fasta = workflow_data['database.db-path']

        mount_points = collect_mount_points(
            workflow,
            manifest,
            config_tools_folder,
            output_dir,
            database_fasta,
            *mzml_files
        )

        docker_image = args.docker_image
        fragpipe_exe = FRAGPIPE_EXE
        cmd = build_docker_run_cmd(
            image=docker_image,
            mount_points=mount_points,
            work_dir=output_dir
        )
        cmd.extend([
            str(fragpipe_exe),
            '--headless',
            '--workflow',
            str(workflow),
            '--manifest',
            str(manifest),
            '--workdir',
            str(output_dir),
            '--config-tools-folder',
            str(config_tools_folder)
        ])
        if args.dry_run:
            cmd.append('--dry-run')

        if args.ram > 0:
            cmd.extend(['--ram', str(args.ram)])

        if args.threads > 0:
            cmd.extend(['--threads', str(args.threads)])

        LOGGER.info('Running FragPipe command: %s', shlex.join(cmd))
        sp.run(cmd, check=True)


def _read_manifest(manifest_path: Path) -> list[Path]:
    """Load manifest and return absolute mzML paths."""
    with manifest_path.open(newline='') as fh:
        paths: list[Path] = []
        for line in fh:
            mzml_value = line.strip().split('\t')[0]
            mzml_path = Path(mzml_value).expanduser()
            if not mzml_path.exists():
                raise FileNotFoundError(f'MZML file not found for sample {sample_id}: {mzml_path}')
            paths.append(mzml_path.resolve())
        if not paths:
            raise ValueError('Manifest does not contain any entries.')
        return paths

def _read_workflow(workflow_path: Path) -> dict:
    """Load FragPipe workflow file and return as a dictionary."""
    database_file = None
    with open(workflow_path, 'r') as f:
        for line in f:
            if line.startswith('database.db-path'):
                _, v = line.strip().split('=', 1)
                v = v.split('#', 1)[0]
                database_file = v.strip().strip('"').strip("'")
                database_file = Path(database_file).expanduser()
                if not database_file.exists():
                    raise FileNotFoundError(
                        f'Database file not found: {database_file}. Please check the workflow file.'
                    )
    if database_file is None:
        raise ValueError('Could not find database.db-path entry in the workflow file.')
    return {
        "database.db-path": database_file
    }

