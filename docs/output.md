# nf-core/marsseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Prepare genome](#prepare-genome) - Build references needer to run the pipeline
- [Prepare pipeline](#prepare-pipeline)
- [Label reads](#label-reads)
- [Align reads](#align-preads)
- [Demultiplexing](#demultiplexing)
- [Quality control](#quality-control)
  - [QC](#qc) - Internal QC report
  - [FastQC](#fastqc) - Raw read QC
  - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Workflow reporting and genomes](#workflow-reporting-and-genomes)
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

The pipeline is executed per `Batch` and therefore the folder structure looks like following:

```console
results/
├── multiqc
├── pipeline_info
├── references
└── SB26
    ├── data
    ├── fastqc
    ├── output
    ├── QC
    ├── SB26.sam
    └── velocity
```

## Prepare genome

The pipeline requires ERCC (spike-ins) to be included in the reference genome. To
accomdate this, the pipeline requires `fasta` and `gtf` reference files. We recommend
using files from [GENCODE](https://www.gencodegenes.org). Reference indexes are built
based on set `--aligner` parameter.

```console
results/references
├── bowtie2
├── gencode.vM32.annotation.gtf
├── GRCm39.primary_assembly.genome_ercc.fa
├── GRCm39.primary_assembly.genome.fa
├── star
└── versions.yml
```

## Prepare pipeline

<details markdown="1">
<summary>Output files</summary>

- `data/`
  - `amp_batches.txt`: Amplification batches
  - `gene_intervals.txt`: Information about gene (chromosome, start, end, strand and symbol)
  - `seq_batches.txt`: Sequencing batches
  - `wells_cells.txt`: Well cells

</details>

At the beginning the pipeline validates the input `*xlsx` files which are further
converted into `.txt` files. The majority of internal tools are written in Perl
by the original authors. The nf-core/marsseq is trying to simplify these steps
while replicating the necessary formats and executions in the background. Because
the pipeline uses `bowtie2` as a default aligner, reads are furher split into
subfiles of `4,000,000` reads per file so the aligner doesn't run out of memory.

## Label reads

In this intermediate step, reads are further processed and barcodes (CB + UMI) is
transferred from R2 to R1. This is important as the CB is further used in the
demultiplexing step.

- `R1` structure: LA (3bp) + PB (4bp) + cDNA (66bp) + RA (2bp)
- `R2` structure: CB (7bp) + UMI (8bp)

where

- LA: left adapter
- PB: pool barcode
- RA: right adapter
- CB: cell barcode
- UMI: Unique Molecular Indentifier

In this step, files are not saved in the overal results, they reside in the `work`
folder.

## Align reads

Split reads are aligned using `bowtie2`. Next, all the aligned reads are merged
into one `SAM` file which is used as an input for demultiplexing.

If `--aligner` flag is set to `bowtie2_star` or `star`, the reads are also aligned
using `StarSolo` to estimated both spliced and unspliced reads which can be used
for RNA velocity estimation. This is an additional plugin which we developed.
In short MARS-seq2.0 reads are converted to `10X v2` format. Additionally, a
whitelist is generated for aligned to perform demultiplexing.

<details markdown="1">
<summary>Output files</summary>

- `<batch>`
  - `<batch>.sam`: Merged aligned reads into one SAM file with `bowtie2`
  - `velocity/`
    - `<batch>_{1,2}.trim.fastq.gz`: Trimmed pair-end converted `10X v2` reads
    - `<batch>.cutadapt.log`: Log file from running `cutadapt`
    - `<batch>.Log.final.out`: STAR alignment report containing the mapping results summary
    - `<batch>.Log.out` and `<batch>.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes
    - `<batch>.Solo.out/*`: Output from StarSolo (Gene, GeneFull, SJ, Velocyto and Barcode.stats)
    - `whitelist.txt`: File containing cell barcodes (combination of pool and cell barcode)

</details>

## Demultiplexing

<details markdown="1">
<summary>Output files</summary>

- `<batch>`
  - `output/`
    - `_debug`: Debugging information
    - `offset.tab/`: A directory that contains tables with the number of IVT products per each(well,gene) pair of all the processed amplification batches where each row is a gene and each column is a well/cell. The filename of a table corresponds to it’s amplification batch (`[amplification batch ID].txt`)
    - `QC`: Intermediate QC stats later used to generate `QC_reports`
    - `singleton_offset.tab`: A directory that contains tables with the number of molecules with a single IVT products per each(well,gene) pair of all the processed amplification batches where each row is a gene and each column is a well/cell. The filename of a table corresponds to it’s amplification batch (`[amplification batch ID].txt)
    - `umi.tab`: A directory that contains the expression tables of all the processed amplification batches. Each table contains he number of detected molecules for each (well,gene) pair where each row is a gene and each column is a well/cell. The filename of a table corresponds to it’s amplification batch (`[amplification batch ID].txt`)
    - `amp_batches_stats.txt`:
    - `amp_batches_summary.txt`:

</details>

Demultiplexing step where pipeline generates a count matrix for invidivual
amplification batch. The count matrix can be found in `<batch>/output/umi.tab`.
The information above have been taken from original documentation for MARS-seq2.0
which can be round [here](https://tanaylab.github.io/old_resources/pages/672.html).

## QC

### QC report

<details markdown="1">
<summary>Output files</summary>

- `<batch>`
  - `output/`
    - `QC_reports`: For each sequencing batch a QC report is generated in `[sequencing batch ID].pdf`. Each contains a detailed QC reports for all the amplification batches associated with it.

</details>

Internal script for generating QC report after `bowtie2` alignment and demultiplexing.

#### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Workflow reporting and genomes

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
