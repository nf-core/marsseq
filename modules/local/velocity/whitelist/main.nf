/*
 * Generate whitelist for StarSolo.
 */
process VELOCITY_WHITELIST {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' :
        'biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' }"

    input:
    path(amp_batches)
    path(well_cells)
    tuple val(meta), path(reads)

    output:
    path "whitelist.txt", emit: whitelist
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    velocity.py whitelist \\
        --batch $meta.id \\
        --amp_batches $amp_batches \\
        --well_cells $well_cells

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        velocity.py: \$( velocity.py --version )
    END_VERSIONS
    """

    stub:
    """
    touch whitelist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        velocity.py: \$( velocity.py --version )
    END_VERSIONS
    """
}
