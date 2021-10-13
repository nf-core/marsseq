// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HISAT2_BUILD {
    tag "$fasta"
    label 'process_high'
    label 'process_high_memory'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::hisat2=2.2.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3"
    } else {
        container "quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3"
    }

    input:
    path fasta

    output:
    path 'hisat2'       , emit: index
    path '*.version.txt', emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    mkdir hisat2
    hisat2-build \\
        --threads $task.cpus \\
        $options.args \\
        $fasta \\
        hisat2/${fasta.baseName}

    echo \$(hisat2-build --version 2>&1) | grep version | cut -d' ' -f 3 > ${software}.version.txt
    """
}
