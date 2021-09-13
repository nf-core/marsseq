// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Demultiplex reads using barcodes. At the same time
 * collect qc reports which will be used later to 
 * construct final QC report per batch.
 */
process DEMULTIPLEX {
    tag "$meta [$amp_batch]"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"$meta", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::conda-forge==5.22.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/perl:5.22.2.1"
    } else {
        container "quay.io/biocontainers/perl:5.22.2.1"
    }

    input:
    tuple val(amp_batch), val(pool_barcode), val(meta), path(trimmed_sam), path(wells_cells), path(gene_intervals), path(spike_seq), path(oligos)

    output:
    tuple val(amp_batch), path("_DS/*")                                    , optional: true, emit: downsampling_factor
    tuple val(amp_batch), path ("output/umi.tab/*.txt")                    , emit: umi_tab
    tuple val(amp_batch), path ("output/offset.tab/*.txt")                 , emit: offset_tab
    tuple val(amp_batch), path ("output/singleton_offset.tab/*.txt")       , emit: singleton_offset_tab
    tuple val(amp_batch), path ("output/QC/read_stats/*.txt")              , emit: qc_read_stats
    tuple val(amp_batch), path ("output/QC/read_stats_amp_batch/*.txt")    , emit: qc_read_stats_amp_batch
    tuple val(amp_batch), path ("output/QC/umi_stats/*.txt")               , emit: qc_umi_stats
    tuple val(amp_batch), path ("output/QC/noffsets_per_umi_distrib/*.txt"), emit: qc_noffsets_per_umi_distrib
    tuple val(amp_batch), path ("output/QC/nreads_per_umi_distrib/*.txt")  , emit: qc_nreads_per_umi_distrib
    tuple val(amp_batch), path ("output/QC/umi_nuc_per_pos/*.txt")         , emit: qc_umi_nuc_per_pos
    tuple val(amp_batch), path ("output/_debug/$amp_batch/*.txt")          , emit: debug

    script:
    """
    mkdir -p output/umi.tab/
    mkdir -p output/offset.tab/
    mkdir -p output/singleton_offset.tab/
    mkdir -p output/QC/read_stats/
    mkdir -p output/QC/read_stats_amp_batch/
    mkdir -p output/QC/umi_stats/
    mkdir -p output/QC/noffsets_per_umi_distrib/
    mkdir -p output/QC/nreads_per_umi_distrib/
    mkdir -p output/QC/umi_nuc_per_pos/
    mkdir -p _debug/$amp_batch/

    demultiplex.pl \\
        $amp_batch \\
        $pool_barcode \\
        $wells_cells \\
        $gene_intervals \\
        $spike_seq \\
        $oligos \\
        $trimmed_sam \\
        . \\
        $options.args
    
    mv _debug output/
    """
}
