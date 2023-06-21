/*
 * Fastp module which splits the reads
 * into specified number of reads per file.
 * Number of reads specified in `modules.config`
 * with default 4_000_000.
 */
process FASTP_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp=0.23.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.3--h5f740d0_0' :
        'biocontainers/fastp:0.23.3--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('raw_reads/*fastq.gz'), emit: reads
    tuple val(meta), path('*.json')             , emit: json
    path 'raw_reads/*.log'                      , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir raw_reads/
    fastp \
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o raw_reads/${reads[0]} \\
        -O raw_reads/${reads[1]} \\
        --thread $task.cpus \\
        --disable_quality_filtering \\
        --disable_length_filtering \\
        --disable_adapter_trimming \\
        --disable_trim_poly_g \\
        --json ${meta.id}.fastp.json \\
        $args \\
        2> raw_reads/${meta.id}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastp.json
    mkdir raw_reads/
    touch raw_reads/000{1..3}.${reads[0]}
    touch raw_reads/000{1..3}.${reads[1]}
    touch raw_reads/${meta.id}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
