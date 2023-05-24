/*
 * Preparation process before executing pipeline.
 * Generates all necessary files like
 *  - amp_batches_to_process.txt
 *  - amp_batches.txt
 *  - seq_batches.txt
 *  - wells_cells.txt
 *  - gene_intervals.txt
 */
process PREPARE {
    tag "$meta.id"
    label 'process_tiny'

    conda "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' :
        'quay.io/biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' }"

    input:
    tuple val(meta), path(reads)
    path(gtf)
    path(ercc_regions)
    
    output:
    path "amp_batches.txt"      , emit: amp_batches
    path "gene_intervals.txt"   , emit: gene_intervals
    path "seq_batches.txt"      , emit: seq_batches
    path "wells_cells.txt"      , emit: wells_cells
    tuple val(meta), path(reads), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_pipeline.py \\
        --batch ${meta.id} \\
        --amp_batches ${meta.amp_batches} \\
        --seq_batches ${meta.seq_batches} \\
        --well_cells ${meta.well_cells} \\
        --gtf $gtf \\
        --output .
    cat $ercc_regions >> gene_intervals.txt
    validate_data.py --input .
    """
}
