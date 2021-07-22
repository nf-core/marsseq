// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Construct spike-ins fasta (ERCC.fasta)
 */
process CREATE_ERCC_FASTA {
    tag "$spikeins"
    label 'process_tiny'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0"
    }

    input:
    path spikeins
    
    output:
    path "ercc.fasta", emit: fasta

    script:
    """
    create_ercc_fasta.py --input ${spikeins} --output ercc.fasta
    """
}
