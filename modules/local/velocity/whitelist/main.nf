/*
 * Generate whitelist for StarSolo.
 */
process VELOCITY_WHITELIST {
    tag "$meta.id"
    label 'process_low'
    
    conda "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' :
        'quay.io/biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "whitelist.txt", emit: whitelist

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    velocity.py whitelist \\
        --batch $meta.id \\
        --amp_batches $meta.amp_batches \\
        --well_cells $meta.well_cells
    """

    stub:
    """
    touch whitelist.txt
    """
}
