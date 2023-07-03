#!/usr/bin/env python
import argparse
import glob
import logging
import os
import sys
from gzip import open as gzopen
from typing import Optional, Tuple

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm.contrib import tzip


def check_folder(folder: str, mandatory: Optional[bool] = False) -> None:
    """
    Sanity check if input/output exists

    Parameters
    ----------
    folder: str
    mandatory: Optional[bool]
    """
    if not os.path.exists(folder):
        if mandatory:
            logging.error(f"Folder {folder} not found!")
            sys.exit(-1)
        else:
            logging.warning(f"Missing folder {folder}, creating ...")
            os.mkdir(folder)


def filter_reads(
    fastq_r1: str,
    fastq_r2: str,
    amp_batch: str,
    pool_barcode: str,
    output_folder: str,
    chunk_size: int = 100,
) -> None:
    """
    Filter reads based on barcodes from defined batches by user. Save reads into output folder.

    Parameters
    ----------
    fastq_r1: str,
    fastq_r2: str,
    pool_barcode: str,
    output_folder: str,
    chunk_size: int = 100,
    """

    with gzopen(fastq_r1, "rt") as fq_r1, gzopen(fastq_r2, "rt") as fq_r2, gzopen(
        f"{output_folder}/{amp_batch}_R1.fastq.gz", "wt"
    ) as fastq_r1_out, gzopen(f"{output_folder}/{amp_batch}_R2.fastq.gz", "wt") as fastq_r2_out:
        counter: int = 0
        fastq_r1_content, fastq_r2_content = "", ""
        for r1, r2 in tzip(FastqGeneralIterator(fq_r1), FastqGeneralIterator(fq_r2)):
            # [0]: header, [1]: seq, [2]: quality

            if pool_barcode in r1[1]:
                counter += 1
                fastq_r1_content += f"@{r1[0]}\n{r1[1]}\n+\n{r1[2]}\n"
                fastq_r2_content += f"@{r2[0]}\n{r2[1]}\n+\n{r2[2]}\n"

                if counter % chunk_size == 0:
                    fastq_r1_out.write(fastq_r1_content)
                    fastq_r2_out.write(fastq_r2_content)
                    counter, fastq_r1_content, fastq_r2_content = 0, "", ""

        # write the rest of the file
        if counter != 0:
            fastq_r1_out.write(fastq_r1_content)
            fastq_r2_out.write(fastq_r2_content)


def load_files(input_folder: str, output_folder: str, amp_batch: str) -> Tuple[str, str, str, str]:
    """
    Helper function for loading all required input files

    Parameters
    ----------
    input_folder: str
    output_folder: str
    amp_batch: str

    Returns
    -------
    Tuple[str, str, str, str]
        R1, R2, amp_batches (xls)
    """
    check_folder(input_folder, mandatory=True)
    check_folder(output_folder)

    r1 = sorted(glob.glob(f"{input_folder}/*R1*.fastq.gz"))
    r2 = sorted(glob.glob(f"{input_folder}/*R2*.fastq.gz"))

    if len(r1) != len(r2) and len(r1) > 0 and len(r2) > 0:
        raise "Something is off, please check you have the same amount of paired-end files!"

    r1, r2 = r1[0], r2[0]
    amp_batches = glob.glob(f"{input_folder}/amp*.xls*")[0]
    pool_barcode = pd.read_excel(amp_batches).query("Amp_batch_ID == @amp_batch")["Pool_barcode"][0]

    return r1, r2, pool_barcode


def main() -> None:
    arg_parser = argparse.ArgumentParser(description="Filter reads per amplification batch from MARS-seq v2.")
    arg_parser.add_argument("input", type=str, help="Input folder with fastq files")
    arg_parser.add_argument("output", type=str, help="Output folder")
    arg_parser.add_argument("--amp_batch", type=str, help="Batch name", required=True)
    arg_parser.add_argument("--verbose", help="Verbose mode", action="store_true")

    args = arg_parser.parse_args()

    logging.basicConfig(
        stream=sys.stdout,
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S",
    )

    r1, r2, pool_barcode = load_files(args.input, args.output, args.amp_batch)
    filter_reads(r1, r2, args.amp_batch, pool_barcode, args.output)


if __name__ == "__main__":
    main()
