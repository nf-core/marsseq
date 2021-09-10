// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * MERGE provided SAM files
 */
process CAT_SAMS {
    tag "$meta"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"$meta", meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    }

    input:
    tuple val(meta), path(sams)

    output:
    tuple val(meta), path("*_trimmed.sam"), emit: sam

    script:
    // https://www.systutorials.com/how-to-merge-sam-files-on-linux/
    """
    (grep ^@ ${sams[0]}; for f in *R1*.sam; do grep -v ^@ \$f; done) > ${meta}_trimmed.sam
    """
}
