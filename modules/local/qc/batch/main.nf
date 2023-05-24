/*
 * QC for individual Batch
 */
process QC_BATCH {
    tag "$meta"
    label 'process_tiny'

    conda "conda-forge::r-mass==7.3_54 conda-forge::r-zoo==1.8_9 conda-forge::r-gplots==3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-520de0c6650803e8b1dab89be11109011a0418b0:464283c471792ace52fbd7b1eea34a42ec86ac81-0' :
        'quay.io/biocontainers/mulled-v2-520de0c6650803e8b1dab89be11109011a0418b0:464283c471792ace52fbd7b1eea34a42ec86ac81-0' }"

    input:
    tuple val(meta), path(folder), path(amp_batches), path(seq_batches), path(wells_cells), path(spike_concentrations)

    output:
    tuple val(meta), path("report_per_amp_batch/*.pdf"), emit: pdf
    tuple val(meta), path("rd/*.rd")                   , emit: rd

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir report_per_amp_batch/ rd/

    qc_batch.r \\
        $meta.amp_batch \\
        $wells_cells \\
        $amp_batches \\
        $seq_batches \\
        $spike_concentrations \\
        $folder
    """
}
