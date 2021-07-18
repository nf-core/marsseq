#!/usr/bin/env python
# pip install openpyxl
import argparse
import logging
import os
import pandas as pd
import sys

from typing import Any, Tuple

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s  %(message)s", "%d-%m-%Y %H:%M:%S")


def read_input(
    args: argparse.Namespace,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Reads: amp_batches, seq_batches and well_cells

    Args:
        args (argparse.Namespace): arguments

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: Content of above files
    """
    amp = pd.read_excel(args.amp_batches).query("Seq_batch_ID == @args.batch")
    seq = pd.read_excel(args.seq_batches).query("Seq_batch_ID == @args.batch")
    wells = pd.read_excel(args.well_cells).query("Amp_batch_ID in @amp.Amp_batch_ID")

    return amp, seq, wells


def add_column(df: pd.DataFrame, column_name: str, values: Any) -> None:
    """Add column into dataframe.

    Args:
        df (pd.DataFrame): DataFrame
        column_name (str): Column name
        values (List[Any]): List of new values
        unique (bool): Unique value instead of list of values. Defaults to False.
    """
    column_name = column_name.strip()
    if column_name in df.columns:
        logging.warning(f"Replacing existing {column_name}")

    df[column_name] = values


def prepare_batch_metadata(
    amp: pd.DataFrame,
    seq: pd.DataFrame,
    wells: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    # create amp_batches.txt
    amp_out = pd.DataFrame()
    add_column(amp_out, "Amp_batch_ID", amp["Amp_batch_ID"])
    add_column(amp_out, "Seq_batch_ID", amp["Seq_batch_ID"])
    add_column(amp_out, "Pool_barcode", amp["Pool_barcode"])
    add_column(amp_out, "Spike_type", wells["Spike_type"])
    add_column(amp_out, "Spike_dilution", wells["Spike_dilution"])
    add_column(amp_out, "Spike_volume_ul", wells["Spike_volume_ul"])
    add_column(amp_out, "Experiment_ID", amp["Experiment_ID"])
    add_column(amp_out, "Owner", amp["Owner"])
    add_column(amp_out, "Description", amp["Description"])

    # create seq_batches.txt
    seq_out = pd.DataFrame(
        [
            amp["Seq_batch_ID"].unique()[0],
            seq.iloc[0]["Run_name"],
            seq.iloc[0]["Date"],
            seq.iloc[0]["R1_design"],
            "",
            amp["R2_design"].unique()[0],
            seq["Genome_assembly"],
        ],
        index=[
            "Seq_batch_ID",
            "Run_name",
            "Date",
            "R1_design",
            "I5_design",
            "R2_design",
            "Notes",
        ],
    ).T

    # wells_cells.txt
    wells_out = pd.DataFrame()
    add_column(wells_out, "Well_ID", wells["Well_ID"])
    add_column(wells_out, "Well_coordinates", wells["Well_coordinates"])
    add_column(wells_out, "plate_ID", wells["plate_ID"])
    add_column(wells_out, "Subject_ID", wells["Subject_ID"])
    add_column(wells_out, "Amp_batch_ID", wells["Amp_batch_ID"])
    add_column(wells_out, "Cell_barcode", wells["Cell_barcode"])
    add_column(wells_out, "Number_of_cells", wells["Number_of_cells"])

    return amp_out, seq_out, wells_out


def generate_gene_intervals(args: argparse.Namespace) -> None:
    """Generate `gene_intervals.txt`

    Args:
        args (argparse.Namespace): arguments

    Format:
        chrom	start	end	strand	gene_name
        chr1	11873	12227	1	DDX11L1
    """

    if not os.path.exists(args.gtf):
        sys.exit(f"Provided GTF {args.gtf} does not exists!")

    gtf = pd.read_csv(args.gtf, sep="\t", skiprows=5, header=None)
    gtf.columns = [
        "chrom",
        "annot",
        "type",
        "start",
        "end",
        "dot",
        "strand",
        "dot2",
        "gene_id",
    ]
    gtf.dropna(inplace=True)

    # transcripts
    gtf_transcripts = gtf.query('type == "transcript"').copy()
    gtf_transcripts["gene_name"] = gtf_transcripts.gene_id.str.split(";", expand=True)[
        3
    ].str.replace("gene_name ", "")
    gtf_transcripts = gtf_transcripts[["chrom", "start", "end", "strand", "gene_name"]]

    # genes
    gtf_genes = gtf.query('type == "gene"').copy()
    gtf_genes["gene_name"] = (
        gtf_genes.gene_id.str.split(";", expand=True)[2]
        .str[len("gene_name") + 1 :]
        .replace("gene_name ", "")
    )
    gtf_genes = gtf_genes[["chrom", "start", "end", "strand", "gene_name"]]

    # merge both transcript and genes
    gtf = pd.concat([gtf_transcripts, gtf_genes])

    # fix columns types
    gtf.start = gtf.start.astype(int)
    gtf.end = gtf.end.astype(int)
    gtf.strand = gtf.strand.map({"+": 1, "-": -1})
    gtf.gene_name = gtf.gene_name.str.replace('"', "").str.strip()

    # save unique results
    gtf.drop_duplicates(inplace=True)
    gtf.to_csv(f"{args.output}/gene_intervals.txt", index=False, sep="\t")


def args() -> argparse.Namespace:
    """Argument

    Returns:
        argparse.Namespace: Parsed arguments
    """
    arg_parser = argparse.ArgumentParser(
        description="Preprocessing script for MARS-seq pipeline."
    )

    arg_parser.add_argument("--version", "-v", action="version", version=f"v0.1")
    arg_parser.add_argument("--batch", action="store", type=str, required=True)
    arg_parser.add_argument("--amp_batches", action="store", type=str, required=True)
    arg_parser.add_argument("--seq_batches", action="store", type=str, required=True)
    arg_parser.add_argument("--well_cells", action="store", type=str, required=True)
    arg_parser.add_argument("--gtf", action="store", type=str, required=True)
    arg_parser.add_argument(
        "--output",
        action="store",
        type=str,
        required=True,
        help="Output path to store the txt files",
    )

    return arg_parser.parse_args()


def main(args: argparse.Namespace):
    """Preparation process before executing pipeline. Generates all necessary files like:
        - amp_batches_to_process.txt
        - amp_batches.txt
        - seq_batches.txt
        - wells_cells.txt
        - gene_intervals.txt

    Args:
        args (argparse.Namespace): arguments
    """
    logging.info(f"Saving files into {args.output}...")
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    try:
        # read input XLSX files
        logging.info("Generating batch files ...")
        amp_batches, seq_batches, wells = read_input(args)
    except FileNotFoundError as ex:
        sys.exit(ex)

    # processed new txt files
    amp_out, seq_out, wells_out = prepare_batch_metadata(
        amp_batches, seq_batches, wells
    )

    # store new results
    amp_out.to_csv(f"{args.output}/amp_batches.txt", sep="\t", index=False)
    seq_out.to_csv(f"{args.output}/seq_batches.txt", sep="\t", index=False)
    wells_out.to_csv(f"{args.output}/wells_cells.txt", sep="\t", index=False)

    # store which amplification batches to process
    logging.info("Generating amplification batches ...")
    amp_out["Amp_batch_ID"].to_csv(
        f"{args.output}/amp_batches_to_process.txt", sep="\t", index=False, header=False
    )

    # gene intervals
    logging.info("Generating gene intervals ...")
    generate_gene_intervals(args)
    logging.info("Done")


if __name__ == "__main__":
    args = args()
    main(args)
