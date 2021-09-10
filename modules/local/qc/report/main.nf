// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Generate final QC report
 */
process QC_REPORT {
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
    tuple val(meta), path(rds), path(pdfs), path(amp_batches), path(wells_cells)

    output:
    path "amp_batches_summary.txt", emit: amp_batches_summary
    path "amp_batches_stats.txt"  , emit: amp_batches_stats
    path "output/QC_reports*"     , emit: pdf

    script:
    """
    mkdir -p output/QC_reports
    mkdir _temp/

    export TMPDIR=/tmp
    qc_report.r \\
        $wells_cells \\
        $amp_batches \\
        . 
    """
}
