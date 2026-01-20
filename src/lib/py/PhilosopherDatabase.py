"""Generate decoy database via Philosopher inside FragPipe Docker image."""
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
	resolve_file,
)


_PROG_NAME = 'PhilosopherDatabase'
_HELP = 'Create a decoy database using Philosopher in FragPipe Docker'
_DESCRIPTION = (
	'Initialize a Philosopher workspace and generate a decoy database from the provided proteome FASTA.'
)
LOGGER = setup_logger(_PROG_NAME)


class PhilosopherDatabase(SubCommand):
	PROG_NAME = _PROG_NAME
	HELP = _HELP
	DESCRIPTION = _DESCRIPTION
	ARGS = {
		'--fasta': {
			'type': Path,
			'required': True,
			'help': 'Path to the input proteome FASTA file from which to generate the decoy database.',
		},
		'--docker-image': {
			'type': str,
			'required': False,
			'default': FRAGPIPE_DOCKER_IMAGE,
			'help': 'Docker image that contains FragPipe and Philosopher.',
		}
	}

	@staticmethod
	def func(args):
		fasta_path = resolve_file(args.fasta, 'proteome FASTA file')
		fasta_dir = fasta_path.parent
		mount_points = collect_mount_points(fasta_path)
		docker_image = args.docker_image
		commands = [
			[
				str(PHILOSOPHER_EXE),
				'workspace',
				'--init',
				'--nocheck',
			],
			[
				str(PHILOSOPHER_EXE),
				'database',
				'--custom',
				str(fasta_path),
			],
		]

		for subcommand in commands:
			cmd = build_docker_run_cmd(
				image=docker_image,
                mount_points=mount_points,
				work_dir=fasta_dir
            )
			cmd.extend(subcommand)
			LOGGER.info('Running Philosopher command: %s', shlex.join(cmd))
			sp.run(cmd, check=True)