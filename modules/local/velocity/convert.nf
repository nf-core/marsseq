// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Convert MARS-seq raw reads into 10X format.
 */
process VELOCITY_CONVERT {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.id}/velocity", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    }

    input:
    tuple val(meta), path(fastp_reads)

    output:
    path "*.fastq.gz", emit: reads

    script:
    """
    velocity.py convert \\
        --input $fastp_reads \\
        --output _temp/ \\
        --threads $task.cpus

    for f in _temp/*R1*.fastq.gz; do cat \$f >> Undetermined_S0_R1_001.fastq.gz; done
    for f in _temp/*R2*.fastq.gz; do cat \$f >> Undetermined_S0_R2_001.fastq.gz; done
    """
}
