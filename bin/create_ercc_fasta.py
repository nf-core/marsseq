#!/usr/bin/env python
import argparse
import pandas as pd
import os
import sys

ROW_LENGTH: int = 70


def create_fasta(args: argparse.Namespace):
    """Generate ERCC fasta file from spike-seq.txt

    Args:
        args (argparse.Namespace): arguments
    """
    if not os.path.exists(args.input):
        print(f"File {args.input} not found")
        sys.exit(-1)

    spikes = pd.read_csv(args.input, sep="\t")
    with open(args.output, "w") as output:
        for _, row in spikes.iterrows():
            comment = f'>{row["ERCC_ID"]}'
            seq = "\n".join(
                [
                    row["Sequence"][i : i + ROW_LENGTH]
                    for i in range(0, len(row["Sequence"]), ROW_LENGTH)
                ]
            )
            output.write(f"{comment}\n")
            output.write(f"{seq}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=str, help="Input file [annotations/spike-seq.txt]"
    )
    parser.add_argument("--output", type=str, help="Output file")
    args = parser.parse_args()

    create_fasta(args)
