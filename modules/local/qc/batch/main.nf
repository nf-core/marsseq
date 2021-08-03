// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Construct spike-ins fasta (ERCC.fasta)
 */
process QC_BATCH {
    tag "$meta"
    label 'process_tiny'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"$meta", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::r-mass==7.3_54 conda-forge::r-zoo==1.8_9 conda-forge::r-gplots==3.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-520de0c6650803e8b1dab89be11109011a0418b0:464283c471792ace52fbd7b1eea34a42ec86ac81-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-520de0c6650803e8b1dab89be11109011a0418b0:464283c471792ace52fbd7b1eea34a42ec86ac81-0"
    }

    input:
    tuple val(amp_batch), path(qcs), val(meta), path(amp_batches), path(seq_batches), path(wells_cells), path(spike_concentrations)

    output:
    tuple val(meta), path("output/QC/report_per_amp_batch/*.pdf"), emit: pdf
    tuple val(meta), path("output/QC/rd/*.rd")                   , emit: rd

    script:
    """
    mkdir -p output/QC/
    mkdir -p ./QC/report_per_amp_batch/
    mkdir -p ./QC/rd/

    qc_batch.r \\
        $amp_batch \\
        $wells_cells \\
        $amp_batches \\
        $seq_batches \\
        $spike_concentrations \\
        .
    
    mv QC/* output/QC/
    """
}
