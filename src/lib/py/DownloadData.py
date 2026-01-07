""" Download dataset files from public repositories """
from __future__ import annotations
from pathlib import Path
import ftplib
import hashlib
from tqdm import tqdm
from Template import SubCommand
import Common


_PROG_NAME = 'DownloadData'
_HELP = 'Download dataset files from public repositories'
_DESCRIPTION = 'Download proteomics dataset files from PRIDE Archive and other sources.'
LOGGER = Common.setup_logger(_PROG_NAME)

DATASETS = {
    "van_puyvelde-2022": {
        "server": "ftp.pride.ebi.ac.uk",
        "path": "pride/data/archive/2022/02/PXD028735/",
        "file_names": [
            "LFQ_TTOF5600_SWATH_Condition_A_Sample_Alpha_01.wiff",
            "LFQ_TTOF5600_SWATH_Condition_A_Sample_Alpha_01.wiff.scan",
            "LFQ_TTOF5600_SWATH_Condition_B_Sample_Alpha_01.wiff",
            "LFQ_TTOF5600_SWATH_Condition_B_Sample_Alpha_01.wiff.scan",
            "LFQ_TTOF5600_SWATH_Ecoli_01.wiff",
            "LFQ_TTOF5600_SWATH_Ecoli_01.wiff.scan",
            "LFQ_TTOF5600_SWATH_Yeast_01.wiff",
            "LFQ_TTOF5600_SWATH_Yeast_01.wiff.scan",
            "LFQ_TTOF5600_SWATH_Human_01.wiff",
            "LFQ_TTOF5600_SWATH_Human_01.wiff.scan"
        ]
    }
}

class DownloadData(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--dataset': {
            'type': str,
            'required': True,
            'choices': list(DATASETS.keys()),
            'help': f'Dataset to download. Available options: {", ".join(DATASETS.keys())}'
        }
    }

    @staticmethod
    def func(args):
        dataset = args.dataset
        if dataset == "van_puyvelde-2022":
            _download_van_puyvelde_2022()

def _calculate_sha512(file_path: Path) -> str:
    """Calculate SHA512 checksum of a file."""
    sha512_hash = hashlib.sha512()
    with open(file_path, 'rb') as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha512_hash.update(byte_block)
    return sha512_hash.hexdigest()

def _download_van_puyvelde_2022():
    server = DATASETS["van_puyvelde-2022"]["server"]
    ftp = ftplib.FTP(server)
    ftp.login()
    ftp.cwd(DATASETS["van_puyvelde-2022"]["path"])

    output_dir = Common.get_dataset_dir(
        dataset_name="van_puyvelde-2022",
        molecule='Protein',
        assay='DIA'
    )

    output_dir /= 'Raw'

    # SHA512 checksum file
    checksum_file = output_dir / "SHA512SUMS.txt"

    # Read existing checksums if file exists
    existing_checksums = {}
    if checksum_file.exists():
        LOGGER.info(f"Reading existing checksums from {checksum_file}")
        with open(checksum_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split(maxsplit=1)
                    if len(parts) == 2:
                        existing_checksums[parts[1]] = parts[0]

    # Get list of available files in the FTP directory
    available_files = set(ftp.nlst())
    LOGGER.info(f"Available files in FTP directory: {len(available_files)} files")

    # Check which files exist and which are missing
    requested_files = DATASETS["van_puyvelde-2022"]["file_names"]
    existing_files = [f for f in requested_files if f in available_files]
    missing_files = [f for f in requested_files if f not in available_files]

    # Report missing files
    if missing_files:
        LOGGER.warning(f"Missing files ({len(missing_files)}):")
        for file_name in missing_files:
            LOGGER.warning(f"  - {file_name}")

    # Check which files are already downloaded with matching checksums
    already_finished = []
    files_to_download = []

    for file_name in existing_files:
        local_file_path = output_dir / file_name
        if local_file_path.exists() and file_name in existing_checksums:
            # File exists locally, verify checksum
            LOGGER.info(f"Verifying {file_name}...")
            actual_checksum = _calculate_sha512(local_file_path)
            expected_checksum = existing_checksums[file_name]
            if actual_checksum == expected_checksum:
                already_finished.append(file_name)
            else:
                LOGGER.warning(f"Checksum mismatch for {file_name}, will re-download")
                files_to_download.append(file_name)
        else:
            files_to_download.append(file_name)

    # Report already finished files
    if already_finished:
        LOGGER.info(f"Already downloaded ({len(already_finished)}):")
        for file_name in already_finished:
            LOGGER.info(f"  ✓ {file_name}")

    # Report files to download
    if files_to_download:
        LOGGER.info(f"Files to download ({len(files_to_download)}):")
        for file_name in files_to_download:
            LOGGER.info(f"  - {file_name}")
    else:
        LOGGER.info("All files already downloaded.")
        ftp.quit()
        return

    # Dictionary to store checksums
    file_checksums = {}

    # Keep existing checksums for already finished files
    for file_name in already_finished:
        file_checksums[file_name] = existing_checksums[file_name]

    # Download files that need to be downloaded
    for file_name in tqdm(files_to_download, desc="Downloading files", unit="file"):
        local_file_path = output_dir / file_name
        LOGGER.info(f"Downloading {file_name} to {local_file_path}")
        with open(local_file_path, 'wb') as f:
            ftp.retrbinary(f"RETR {file_name}", f.write)

        # Calculate SHA512 checksum
        LOGGER.info(f"Calculating SHA512 for {file_name}...")
        checksum = _calculate_sha512(local_file_path)
        file_checksums[file_name] = checksum

    ftp.quit()

    # Save all checksums to file
    LOGGER.info(f"Saving checksums to {checksum_file}")
    with open(checksum_file, 'w') as f:
        for file_name in requested_files:
            if file_name in file_checksums:
                f.write(f"{file_checksums[file_name]}  {file_name}\n")

    LOGGER.info("Download completed.")
