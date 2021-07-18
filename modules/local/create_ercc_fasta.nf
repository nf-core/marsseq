// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Construct spike-ins fasta (ERCC.fasta)
 */
process CREATE_ERCC_FASTA {
    tag "$spikeins"
    label 'process_tiny'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"references/${params.genome}") }

    conda (params.enable_conda ? "conda-forge::pandas==1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {        
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    input:
    path spikeins
    
    output:
    path "ercc.fasta", emit: ercc

    script:
    """
    create_ercc_fasta.py --input ${spikeins} --output ercc.fasta
    """
}
