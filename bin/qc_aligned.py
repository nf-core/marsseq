#!/usr/bin/env python
import argparse
import os
import pandas as pd
from subprocess import Popen, PIPE


def exec_proc(command: str):
    with Popen([command], shell=True, stdout=PIPE, stderr=PIPE) as p:
        output, _ = p.communicate()

    return output.decode("utf-8").splitlines()


def qc_read(sam_file: str, qc_file: str, output: str, prefix: str = "_"):
    command: str = f"grep -v ^@ {sam_file} | cut -f2"
    lines = exec_proc(command)

    stats = pd.read_csv(qc_file, sep="\t")
    stats["mapped"] = lines.count("0") + lines.count("16")
    stats.to_csv(f"{output}/{prefix}{os.path.basename(qc_file)}", index=None, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", type=str, help="Sam file [str]", required=True)
    parser.add_argument("--qc", type=str, help="QC of sam [str]", required=True)
    parser.add_argument("--output", type=str, help="Output file")
    args = parser.parse_args()

    qc_read(args.sam, args.qc, args.output)
