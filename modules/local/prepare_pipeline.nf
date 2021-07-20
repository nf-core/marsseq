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
    tag "${meta}"
    label 'process_tiny'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    }

    input:
    tuple val(meta), path(data), path(metadata)
    path gtf
    path ercc
    
    output:
    path "*.txt", emit: txts

    script:
    """
    prepare_pipeline.py \\
        --batch $meta \\
        --amp_batches ${metadata[0]} \\
        --seq_batches ${metadata[1]} \\
        --well_cells ${metadata[2]} \\
        --gtf $gtf \\
        --output .
    cat ${ercc} >> gene_intervals.txt

    validate_data.py --input .
    """
}
