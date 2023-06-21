/*
 * Helper: trim SAM read
 */
process CUT_SAM {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.trimmed.sam"), emit: sam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filename = read.baseName + '.trimmed.sam'
    """
    cut -f1-9,12- $read > $filename

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$( cut --version 2>&1 | sed -n 1p | sed 's/cut (GNU coreutils) //g' )
    END_VERSIONS
    """

    stub:
    def filename = read.baseName + '.trimmed.sam'
    """
    touch ${filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$( cut --version 2>&1 | sed -n 1p | sed 's/cut (GNU coreutils) //g' )
    END_VERSIONS
    """
}
