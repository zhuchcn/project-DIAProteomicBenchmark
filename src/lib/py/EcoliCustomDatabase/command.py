"""Top-level E. coli custom database command scaffold."""

from __future__ import annotations
import argparse

from Template import SubCommand
from Common import setup_logger

from EcoliCustomDatabase import (
    cmd_align_bwa_mem2,
    cmd_annotate_vep,
    cmd_build_bwa_mem2_index,
    cmd_build_index,
    cmd_call_variant,
    cmd_call_variant_bcftools,
    cmd_download_sra,
    cmd_format_variant_fasta,
    cmd_install_vep_cache,
    cmd_parse_vep,
)
from EcoliCustomDatabase.constants import (
    DESCRIPTION,
    HELP,
    PROG_NAME,
    SUBCOMMAND_ANNOTATE_VEP,
    SUBCOMMAND_ALIGN_BWA_MEM2,
    SUBCOMMAND_BUILD_BWA_MEM2_INDEX,
    SUBCOMMAND_BUILD_INDEX,
    SUBCOMMAND_CALL_VARIANT,
    SUBCOMMAND_CALL_VARIANT_BCFTOOLS,
    SUBCOMMAND_DOWNLOAD_SRA,
    SUBCOMMAND_FORMAT_VARIANT_FASTA,
    SUBCOMMAND_INSTALL_VEP_CACHE,
    SUBCOMMAND_PARSE_VEP,
    SUBCOMMANDS,
)

LOGGER = setup_logger(PROG_NAME)

SUBCOMMAND_HELP = {
    SUBCOMMAND_ANNOTATE_VEP: 'Annotate using VEP.',
    SUBCOMMAND_BUILD_INDEX: 'Generate moPepGen index using E. coli genome reference.',
    SUBCOMMAND_PARSE_VEP: 'Parse annotated VEP output.',
    SUBCOMMAND_CALL_VARIANT: 'Call moPepGen to generate database FASTA.',
    SUBCOMMAND_DOWNLOAD_SRA: 'Download SRA sample and export FASTQ files.',
    SUBCOMMAND_BUILD_BWA_MEM2_INDEX: 'Build bwa-mem2 index for E. coli reference genome.',
    SUBCOMMAND_ALIGN_BWA_MEM2: 'Align reads to reference genome with bwa-mem2.',
    SUBCOMMAND_CALL_VARIANT_BCFTOOLS: 'Call and filter variants using bcftools.',
    SUBCOMMAND_INSTALL_VEP_CACHE: 'Install VEP cache for E. coli.',
    SUBCOMMAND_FORMAT_VARIANT_FASTA: 'Format variant FASTA headers for proteomics search engines.',
}

SUBCOMMAND_BUILDERS = {
    SUBCOMMAND_ANNOTATE_VEP: cmd_annotate_vep.build_parser,
    SUBCOMMAND_BUILD_INDEX: cmd_build_index.build_parser,
    SUBCOMMAND_PARSE_VEP: cmd_parse_vep.build_parser,
    SUBCOMMAND_CALL_VARIANT: cmd_call_variant.build_parser,
    SUBCOMMAND_DOWNLOAD_SRA: cmd_download_sra.build_parser,
    SUBCOMMAND_BUILD_BWA_MEM2_INDEX: cmd_build_bwa_mem2_index.build_parser,
    SUBCOMMAND_ALIGN_BWA_MEM2: cmd_align_bwa_mem2.build_parser,
    SUBCOMMAND_CALL_VARIANT_BCFTOOLS: cmd_call_variant_bcftools.build_parser,
    SUBCOMMAND_INSTALL_VEP_CACHE: cmd_install_vep_cache.build_parser,
    SUBCOMMAND_FORMAT_VARIANT_FASTA: cmd_format_variant_fasta.build_parser,
}

SUBCOMMAND_RUNNERS = {
    SUBCOMMAND_ANNOTATE_VEP: cmd_annotate_vep.run,
    SUBCOMMAND_BUILD_INDEX: cmd_build_index.run,
    SUBCOMMAND_PARSE_VEP: cmd_parse_vep.run,
    SUBCOMMAND_CALL_VARIANT: cmd_call_variant.run,
    SUBCOMMAND_DOWNLOAD_SRA: cmd_download_sra.run,
    SUBCOMMAND_BUILD_BWA_MEM2_INDEX: cmd_build_bwa_mem2_index.run,
    SUBCOMMAND_ALIGN_BWA_MEM2: cmd_align_bwa_mem2.run,
    SUBCOMMAND_CALL_VARIANT_BCFTOOLS: cmd_call_variant_bcftools.run,
    SUBCOMMAND_INSTALL_VEP_CACHE: cmd_install_vep_cache.run,
    SUBCOMMAND_FORMAT_VARIANT_FASTA: cmd_format_variant_fasta.run,
}


class EcoliCustomDatabase(SubCommand):
    """Entry command that routes to package subcommands."""

    PROG_NAME = PROG_NAME
    HELP = HELP
    DESCRIPTION = DESCRIPTION
    ARGS = {
        'subcommand': {
            'type': str,
            'choices': SUBCOMMANDS,
            'help': (
                'Subcommand to execute. Available options: '
                + ', '.join(SUBCOMMANDS)
            ),
        },
        'subcommand_args': {
            'nargs': argparse.REMAINDER,
            'help': 'Arguments passed through to the selected subcommand.',
        },
    }

    @staticmethod
    def func(args):
        subcommand = args.subcommand
        subcommand_args = list(args.subcommand_args or [])
        parser = argparse.ArgumentParser(
            prog=f'{PROG_NAME} {subcommand}',
            description=SUBCOMMAND_HELP[subcommand],
        )
        SUBCOMMAND_BUILDERS[subcommand](parser)
        parsed = parser.parse_args(subcommand_args)
        LOGGER.info('Dispatching subcommand: %s', subcommand)
        SUBCOMMAND_RUNNERS[subcommand](parsed)
