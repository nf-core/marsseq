// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Preparation process before executing pipeline.
 * Generates all necessary files like
 *  - amp_batches_to_process.txt
 *  - amp_batches.txt
 *  - seq_batches.txt
 *  - wells_cells.txt
 *  - gene_intervals.txt
 */
process PREPARE_PIPELINE {
    tag "$meta.id"
    label 'process_tiny'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"$meta.id/data", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    }

    input:
    tuple val(meta), path(reads)
    
    output:
    path "*.txt", emit: txts

    script:
    def ch_gtf          = file(params.gtf, checkIfExists: true)
    def ch_ercc_regions = file("$projectDir/data/ercc-regions.tsv", checkIfExists: true)
    """
    prepare_pipeline.py \\
        --batch ${meta.id} \\
        --amp_batches ${meta.amp_batches} \\
        --seq_batches ${meta.seq_batches} \\
        --well_cells ${meta.well_cells} \\
        --gtf $ch_gtf \\
        --output .
    cat ${ch_ercc_regions} >> gene_intervals.txt
    validate_data.py --input .
    """
}
