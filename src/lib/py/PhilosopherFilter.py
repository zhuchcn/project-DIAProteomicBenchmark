"""Run Philosopher filter to finalize identification results."""
from __future__ import annotations
import shlex
from pathlib import Path
import subprocess as sp

from Template import SubCommand
from Common import (
    setup_logger,
    PHILOSOPHER_EXE,
    FRAGPIPE_DOCKER_IMAGE,
    collect_mount_points,
    build_docker_run_cmd,
    resolve_directory,
    resolve_file,
)

_PROG_NAME = 'PhilosopherFilter'
_HELP = 'Run Philosopher filter against search results.'
_DESCRIPTION = (
    'Run the Philosopher filter step inside the FragPipe Docker container using the provided '
    'pepXML and protXML files, applying picked-FDR, razor peptide assignment, and default FDR cutoffs. '
)
LOGGER = setup_logger(_PROG_NAME)


class PhilosopherFilter(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--pepxml': {
            'type': Path,
            'required': True,
            'help': 'Path to the philosopher pepXML output to filter.',
        },
        '--protxml': {
            'type': Path,
            'required': True,
            'help': 'Path to the philosopher protXML output that accompanies the pepXML.',
        },
        '--pep': {
            'type': float,
            'default': 0.01,
            'help': 'Peptide-level q-value threshold to pass to philosopher filter (default 0.01).',
        },
        '--prot': {
            'type': float,
            'default': 0.01,
            'help': 'Protein-level q-value threshold to pass to philosopher filter (default 0.01).',
        },
        '--psm': {
            'type': float,
            'default': 0.01,
            'help': 'PSM-level q-value threshold to pass to philosopher filter (default 0.01).',
        },
        '--ion': {
            'type': float,
            'default': 0.01,
            'help': 'Peptide ion q-value threshold to pass to philosopher filter (default 0.01).',
        },
        '--min-pep-len': {
            'type': int,
            'default': 8,
            'help': 'Minimum peptide length to require when filtering identifications.',
        },
        '--tag': {
            'type': str,
            'default': 'rev_',
            'help': 'Decoy tag prefix to pass to philosopher.',
        },
        '--database': {
            'type': Path,
            'required': True,
            'help': 'Path to the FASTA database used for the search.',
        },
        '--docker-image': {
            'type': str,
            'default': FRAGPIPE_DOCKER_IMAGE,
            'help': 'Docker image that provides FragPipe and Philosopher binaries.',
        },
        '--work-dir': {
            'type': Path,
            'default': None,
            'help': 'Working directory in which to execute Philosopher (defaults to current directory).',
        },
    }

    @staticmethod
    def func(args):
        pepxml_path = resolve_directory(args.pepxml, 'pepXML file')
        protxml_path = resolve_file(args.protxml, 'protXML file')
        work_dir_input = args.work_dir
        work_dir = resolve_directory(work_dir_input, 'working directory')
        database = resolve_file(args.database, 'FASTA database file')
        docker_image = args.docker_image
        mount_points = collect_mount_points(pepxml_path, protxml_path, work_dir, database)

        workspace_clean_cmd = shlex.join([
            str(PHILOSOPHER_EXE),
            'workspace',
            '--clean',
            '--nocheck'
        ])
        workspace_init_cmd = shlex.join([
            str(PHILOSOPHER_EXE),
            'workspace',
            '--init',
            '--nocheck'
        ])
        datatbase_annotate_cmd = shlex.join([
            str(PHILOSOPHER_EXE),
            'database',
            '--annotate', str(args.database),
            '--prefix', args.tag
        ])

        filter_cmd = shlex.join([
            str(PHILOSOPHER_EXE),
            'filter',
            '--picked',
            '--pep', str(args.pep),
            '--prot', str(args.prot),
            '--psm', str(args.psm),
            '--ion', str(args.ion),
            '--minPepLen', str(args.min_pep_len),
            '--group',
            '--tag', args.tag,
            '--pepxml', str(pepxml_path),
            '--protxml', str(protxml_path),
            '--razor',
        ])
        report_cmd = shlex.join([
            str(PHILOSOPHER_EXE),
            'report',
            '--removecontam'
        ])

        docker_cmd = build_docker_run_cmd(
            image=docker_image,
            mount_points=mount_points,
            work_dir=work_dir,
        )

        commands = [
            workspace_clean_cmd,
            workspace_init_cmd,
            datatbase_annotate_cmd,
            filter_cmd,
            report_cmd
        ]
        for cmd in commands:
            full_cmd = docker_cmd + shlex.split(cmd)
            LOGGER.info('Running Philosopher command: %s', shlex.join(full_cmd))
            sp.run(full_cmd, check=True)
