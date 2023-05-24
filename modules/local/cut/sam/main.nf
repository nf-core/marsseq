/*
 * Helper: trim SAM read
 */
process CUT_SAM {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.trimmed.sam"), emit: sam

    when:
    task.ext.when == null || task.ext.when

    script:
    def filename = read.baseName + '.trimmed.sam'
    """
    cut -f1-9,12- $read > $filename
    """
}
