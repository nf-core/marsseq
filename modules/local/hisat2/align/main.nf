// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high_cpu'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::hisat2=2.2.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3"
    } else {
        container "quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3"
    }

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.sam'), emit: sam
    tuple val(meta), path("*.log"), emit: summary
    path  '*.version.txt'         , emit: version

    script:
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def split_cpus = Math.floor(task.cpus / 2)
    def software   = getSoftwareName(task.process)
    """
    INDEX=`find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//'`
    hisat2 \\
        -x \$INDEX \\
        -U ${reads} \\
        --summary-file ${prefix}.hisat2.summary.log \\
        --threads $task.cpus \\
        $options.args \\
        > ${prefix}.sam
    
    if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
        mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
    fi
    if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
        mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
    fi

    echo \$(hisat2 --version 2>&1) | grep version | cut -d' ' -f 3 > ${software}.version.txt
    """
}
