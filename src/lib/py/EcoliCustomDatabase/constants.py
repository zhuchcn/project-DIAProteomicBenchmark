"""Constants for the E. coli custom database command scaffold."""
import os
from pathlib import Path

from Common import DATA_DIR

PROG_NAME = 'EcoliCustomDatabase'
HELP = 'Build a custom E. coli variant protein database via a multi-step workflow.'
DESCRIPTION = (
    'Workflow entrypoint for creating an E. coli custom database. '
    'Use one of the subcommands to run an individual stage.'
)

SUBCOMMAND_ANNOTATE_VEP = 'AnnotateVEP'
SUBCOMMAND_BUILD_INDEX = 'BuildIndex'
SUBCOMMAND_PARSE_VEP = 'ParseVEP'
SUBCOMMAND_CALL_VARIANT = 'CallVariant'
SUBCOMMAND_DOWNLOAD_SRA = 'DownloadSRA'
SUBCOMMAND_BUILD_BWA_MEM2_INDEX = 'BuildBwaMem2Index'
SUBCOMMAND_ALIGN_BWA_MEM2 = 'AlignBwaMem2'
SUBCOMMAND_CALL_VARIANT_BCFTOOLS = 'CallVariantBcftools'
SUBCOMMAND_INSTALL_VEP_CACHE = 'InstallVepCache'
SUBCOMMAND_FORMAT_VARIANT_FASTA = 'FormatVariantFasta'

# ---- Shared defaults ----
DEFAULT_SRA_ID = 'SRR1770413'
DEFAULT_FASTQ_DIR = DATA_DIR / f'data/van_puyvelde-2022/DNA/Raw/{DEFAULT_SRA_ID}'
DEFAULT_FASTQ_R1 = DEFAULT_FASTQ_DIR / f'{DEFAULT_SRA_ID}_1.fastq'
DEFAULT_FASTQ_R2 = DEFAULT_FASTQ_DIR / f'{DEFAULT_SRA_ID}_2.fastq'

DEFAULT_REFERENCE_GENOME = DATA_DIR / (
    'ref/reference/escherichia_coli/EnsemblBacteria-62/'
    'Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.'
    'ASM584v2.dna.chromosome.Chromosome.fa'
)
DEFAULT_ANNOTATION_GTF = DATA_DIR / (
    'ref/reference/escherichia_coli/EnsemblBacteria-62/'
    'Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.62.gtf'
)
DEFAULT_PROTEOME_FASTA = DATA_DIR / (
    'ref/reference/escherichia_coli/EnsemblBacteria-62/'
    'Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.pep.all.fa'
)
DEFAULT_MOPEPGEN_INDEX_DIR = DATA_DIR / 'ref/index/moPepGen/1.5.1-rc4'
DEFAULT_MOPEPGEN_REPO = Path('/home/chenghaozhu/project/package-moPepGen')
DEFAULT_MOPEPGEN_DOCKER_IMAGE = 'ghcr.io/zhuchcn/mopepgen:1.5.1-rc4'

DEFAULT_BWA_MEM2_INDEX_DIR = DATA_DIR / 'ref/index/BWA-MEM2/2.3/escherichia_coli/EnsemblBacteria-62'
DEFAULT_BWA_MEM2_INDEX_PREFIX = DEFAULT_BWA_MEM2_INDEX_DIR / 'index'
DEFAULT_SORTED_BAM = DATA_DIR / f'data/van_puyvelde-2022/DNA/BAM/{DEFAULT_SRA_ID}/{DEFAULT_SRA_ID}.sorted.bam'
DEFAULT_VCF_DIR = DATA_DIR / f'data/van_puyvelde-2022/DNA/VCF/{DEFAULT_SRA_ID}'
DEFAULT_RAW_VCF_GZ = DEFAULT_VCF_DIR / f'{DEFAULT_SRA_ID}.vcf.gz'
DEFAULT_FILTERED_VCF_GZ = DEFAULT_VCF_DIR / f'{DEFAULT_SRA_ID}.filtered.vcf.gz'
DEFAULT_VEP_DIR = DATA_DIR / f'data/van_puyvelde-2022/DNA/VEP/{DEFAULT_SRA_ID}'
DEFAULT_VEP_OUTPUT_TSV = DEFAULT_VEP_DIR / f'{DEFAULT_SRA_ID}.vep.tsv'
DEFAULT_GVF_DIR = DATA_DIR / f'data/van_puyvelde-2022/DNA/Database/{DEFAULT_SRA_ID}'
DEFAULT_GVF_OUTPUT = DEFAULT_GVF_DIR / f'{DEFAULT_SRA_ID}.gvf'
DEFAULT_VARIANT_PEPTIDE_FASTA = DEFAULT_GVF_DIR / f'{DEFAULT_SRA_ID}.variant_peptide.fasta'
DEFAULT_FORMATTED_VARIANT_FASTA = DEFAULT_GVF_DIR / f'{DEFAULT_SRA_ID}.variant_peptide.formatted.fasta'
DEFAULT_VEP_CACHE_DIR = DATA_DIR / 'ref/index/vep/115.2'
DEFAULT_VEP_SPECIES = 'escherichia_coli_str_k_12_substr_mg1655_gca_000005845'
DEFAULT_VEP_KINGDOM = 'bacteria'
DEFAULT_VEP_ASSEMBLY = 'ASM584v2'
DEFAULT_CODON_TABLE = 'Bacterial'
DEFAULT_CONDA_ENV = 'bio'
DEFAULT_THREADS = max(1, os.cpu_count() or 1)

SUBCOMMANDS = [
    SUBCOMMAND_ANNOTATE_VEP,
    SUBCOMMAND_BUILD_INDEX,
    SUBCOMMAND_PARSE_VEP,
    SUBCOMMAND_CALL_VARIANT,
    SUBCOMMAND_DOWNLOAD_SRA,
    SUBCOMMAND_BUILD_BWA_MEM2_INDEX,
    SUBCOMMAND_ALIGN_BWA_MEM2,
    SUBCOMMAND_CALL_VARIANT_BCFTOOLS,
    SUBCOMMAND_INSTALL_VEP_CACHE,
    SUBCOMMAND_FORMAT_VARIANT_FASTA,
]
