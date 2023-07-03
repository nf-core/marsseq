#!/usr/bin/env python
import argparse
import glob
import itertools
import logging
import multiprocessing as mp
import os
import sys
from gzip import open as gzopen
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

config: Dict[str, int] = {
    "LEFT_ADAPTER": 3,
    "RIGHT_ADAPTER": 2,
    "POOL_BARCODE": 4,
    "CELL_BARCODE": 7,
    "UMI": 8,
    "DEFAULT_SAMPLE_NAME": "GGTTTACT",  # arbitrary sequence
}


def check_folder(folder: str, mandatory: Optional[bool] = False):
    """Checks if folder exists, otherwise create it.

    Args:
        folder (str): Output folder
    """
    if not os.path.exists(folder):
        if mandatory:
            logging.error(f"Folder {folder} not found!")
            sys.exit(-1)
        else:
            logging.warning(f"Missing folder {folder}, creating ...")
            os.mkdir(folder)


def trim_cdna(r1: List[str]) -> Tuple[str, str]:
    """Trim cdna by removing adapters and Pool barcode.

    Args:
        r1 (List[str]): full cdna sequence

    Returns:
        Tuple[str, str]: trimmed cdna and quality
    """
    global config

    # [0]: header, [1]: seq, [2]: quality
    _, seq, qa = r1
    cdna: str = seq[(config["LEFT_ADAPTER"] + config["POOL_BARCODE"]) : -config["RIGHT_ADAPTER"]]
    cdna_quality: str = qa[(config["LEFT_ADAPTER"] + config["POOL_BARCODE"]) : -config["RIGHT_ADAPTER"]]

    return cdna, cdna_quality


def create_r2(r1: FastqGeneralIterator, r2: FastqGeneralIterator) -> Tuple[str, str]:
    """Create R2 consisting of
        - pool barcode (4bp)
        - cell barcode (7bp)
        - UMI          (8bp)

    Args:
        r1 ([FastqGeneralIterator]): R1 read
        r2 ([FastqGeneralIterator]): R2 read

    Returns:
        Tuple[str, str]: New barcode sequence and quality
    """
    global config

    _, r1_seq, r1_qa = r1
    _, r2_seq, r2_qa = r2

    # get pool barcode
    pb_seq: str = r1_seq[config["LEFT_ADAPTER"] : config["LEFT_ADAPTER"] + config["POOL_BARCODE"]]
    pb_qa: str = r1_qa[config["LEFT_ADAPTER"] : config["LEFT_ADAPTER"] + config["POOL_BARCODE"]]

    barcode_len: int = config["POOL_BARCODE"] + config["CELL_BARCODE"] + config["UMI"]
    barcode_seq: str = f"{pb_seq}{r2_seq}"[:barcode_len]
    barcode_qa: str = f"{pb_qa}{r2_qa}"[:barcode_len]

    return barcode_seq, barcode_qa


def convert_to_10x(params: Tuple[str, str, str, str], chunk_size: int = 100) -> None:
    """Main function processing the reads.

    Args:
        params (Tuple[str, str, str, str]): R1, R2, input, output folder
    """
    fastq_r1, fastq_r2, input_folder, output_folder = params

    with gzopen(fastq_r1, "rt") as fq_r1, gzopen(fastq_r2, "rt") as fq_r2, gzopen(
        f"{output_folder}/{os.path.basename(fastq_r1)}", "wt"
    ) as fastq_r1_out, gzopen(f"{output_folder}/{os.path.basename(fastq_r2)}", "wt") as fastq_r2_out:
        fastq_r1 = FastqGeneralIterator(fq_r1)
        fastq_r2 = FastqGeneralIterator(fq_r2)

        counter: int = 0
        fastq_r1_content, fastq_r2_content = "", ""
        for r1, r2 in zip(fastq_r1, fastq_r2):
            # [0]: header, [1]: seq, [2]: quality

            cdna, cdna_quality = trim_cdna(r1)
            barcode, barcode_quality = create_r2(r1, r2)

            fastq_r1_content += f'@{r1[0]}:{config["DEFAULT_SAMPLE_NAME"]}\n{barcode}\n+\n{barcode_quality}\n'
            fastq_r2_content += f'@{r2[0]}:{config["DEFAULT_SAMPLE_NAME"]}\n{cdna}\n+\n{cdna_quality}\n'

            if counter % chunk_size == 0:
                fastq_r1_out.write(fastq_r1_content)
                fastq_r2_out.write(fastq_r2_content)
                counter, fastq_r1_content, fastq_r2_content = 0, "", ""

        # write the rest of the file
        if counter != 0:
            fastq_r1_out.write(fastq_r1_content)
            fastq_r2_out.write(fastq_r2_content)


def convert(fastqs_folder: str, output_folder: str, threads: int):
    """Main function for converting MARS-seq raw files to 10X format.

    Args:
        input_folder (str): Input folder
        output_folder (str): Output folder
    """

    check_folder(fastqs_folder, mandatory=True)
    check_folder(output_folder)

    r1_files = sorted(glob.glob(f"{fastqs_folder}/*R1*.fastq.gz"))
    r2_files = sorted(glob.glob(f"{fastqs_folder}/*R2*.fastq.gz"))

    if len(r1_files) != len(r2_files):
        print(f"Something is off, please check you have the same amount of paired-end files!")
        sys.exit(-1)

    data = zip(
        r1_files,
        r2_files,
        itertools.repeat(fastqs_folder, len(r1_files)),
        itertools.repeat(output_folder, len(r1_files)),
    )
    with mp.Pool(threads) as pool:
        _ = pool.map(convert_to_10x, data)


def whitelist(batch: str, amp_batches: str, well_cells: str):
    """Generate whitelist for StarSolo.
    The file should contain cell tags.

    Args:
        batch (str): batch name
        amp_batches (str): Amplification batches [xls]
        well_cells (str): Well cells [xlsx]
    """
    batches = pd.read_excel(amp_batches)
    wells = pd.read_excel(well_cells)

    # preprocess and clean-up
    batches = batches[batches.Seq_batch_ID == batch].astype("category")
    wells = wells[wells.Amp_batch_ID.isin(batches.Amp_batch_ID)]
    wells = wells[["Amp_batch_ID", "Cell_barcode"]]
    wells["Amp_batch_ID"] = wells["Amp_batch_ID"].str.strip().astype("category")
    wells["Cell_barcode"] = wells["Cell_barcode"].str.strip().astype("category")

    # merge and concat barcodes
    wells["Amp_batch_ID"] = wells["Amp_batch_ID"].cat.rename_categories(batches.Pool_barcode.unique())
    wells["whitelist"] = wells[["Amp_batch_ID", "Cell_barcode"]].agg("".join, axis=1)

    # save
    wells["whitelist"].to_csv(f"whitelist.txt", header=None, index=None)


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s  %(message)s", "%d-%m-%Y %H:%M:%S")

    arg_parser = argparse.ArgumentParser(description="Plugin for converting and running velocity on MARS-seq v2.")

    arg_parser.add_argument("--version", "-v", action="version", version=f"velocity 0.1")
    command_parser = arg_parser.add_subparsers(dest="command")

    # convert
    convert_parser = command_parser.add_parser("convert", help="Convert reads to 10X format")
    convert_parser.add_argument("--input", type=str, help="Input folder with fastq files", required=True)
    convert_parser.add_argument("--output", type=str, help="Output folder", required=True)
    convert_parser.add_argument("--threads", type=int, help="Number of threads", required=True, default=4)

    # whitelist
    whitelist_parser = command_parser.add_parser("whitelist", help="Create whitelist for StarSolo.")
    whitelist_parser.add_argument("--batch", type=str, help="Batch name", required=True)
    whitelist_parser.add_argument("--amp_batches", type=str, help="Amplification batches [xls]", required=True)
    whitelist_parser.add_argument("--well_cells", type=str, help="Well cells [xls]", required=True)

    args = arg_parser.parse_args()

    if args.command == "convert":
        convert(args.input, args.output, args.threads)
    elif args.command == "whitelist":
        whitelist(args.batch, args.amp_batches, args.well_cells)
    else:
        logging.error(f"Command {args.command} not recognized!")
