/*
 * Convert MARS-seq raw reads into 10X format.
 */
process VELOCITY_CONVERT {
    tag "$meta.id"
    label 'process_high_cpu'

    conda "anaconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' :
        'biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("_temp/*.fastq.gz"), emit: reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    velocity.py convert_read \\
        --r1 ${reads[0]} \\
        --r2 ${reads[1]} \\
        --output _temp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        velocity.py: \$( velocity.py --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir _temp
    touch _temp/0001.Undetermined_S0_L001_R{1,2}_001.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        velocity.py: \$( velocity.py --version )
    END_VERSIONS
    """
}
