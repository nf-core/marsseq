#!/usr/bin/env python
import argparse
import pandas as pd
import logging
import os
import sys


logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s  %(message)s", "%d-%m-%Y %H:%M:%S")


def is_unique(df: pd.DataFrame, column: str) -> None:
    """Checks if a provided column is unique.

    Args:
        df (pd.DataFrame): Dataset
        column (str): Column to check
    """

    if df[column].duplicated().sum() > 0:
        sys.exit(f"The {column} is not unique, please make sure there are no duplicates.")


def main(args: argparse.Namespace):
    """Main function

    Args:
        args (argparse.Namespace): arguments
    """
    if not os.path.exists(args.input):
        sys.exit(f"Provided {args.input} does no exists!")

    amp = pd.read_csv(f"{args.input}/amp_batches.txt", sep="\t")
    seq = pd.read_csv(f"{args.input}/seq_batches.txt", sep="\t")
    wells = pd.read_csv(f"{args.input}/wells_cells.txt", sep="\t")

    logging.info("Checking Well_ID, Amp_batch_ID, Seq_batch_ID")
    is_unique(wells, "Well_ID")
    is_unique(amp, "Amp_batch_ID")
    is_unique(seq, "Seq_batch_ID")

    logging.info("Checking if barcodes are non-unique in batches")
    for amp_batch in wells["Amp_batch_ID"].unique():
        is_unique(wells.query("Amp_batch_ID == @amp_batch"), "Cell_barcode")

    if sorted(amp["Seq_batch_ID"].unique()) != sorted(seq["Seq_batch_ID"].unique()):
        sys.exit("Some amplification batches are linked to undefined Seq_batch_ID")

    logging.info("Validation: passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate data folder which contains all txt files.")
    parser.add_argument("--input", type=str, help="Input folder")
    args = parser.parse_args()

    main(args)
