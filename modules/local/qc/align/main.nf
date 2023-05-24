/*
 * QC report after alignment
 */
process QC_ALIGNED {
    tag "$meta.id"
    label 'process_tiny'

    conda "conda-forge::r-mass==7.3_54 conda-forge::r-zoo==1.8_9 conda-forge::r-gplots==3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-520de0c6650803e8b1dab89be11109011a0418b0:464283c471792ace52fbd7b1eea34a42ec86ac81-0' :
        'quay.io/biocontainers/mulled-v2-520de0c6650803e8b1dab89be11109011a0418b0:464283c471792ace52fbd7b1eea34a42ec86ac81-0' }"

    input:
    tuple val(meta), path(sam)
    path(labeled_qc)
    
    output:
    path ("*.txt"), emit: qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    qc_align.r \\
        $sam \\
        $labeled_qc
    """
}
