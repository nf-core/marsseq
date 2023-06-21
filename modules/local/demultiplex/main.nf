/*
 * Demultiplex reads using barcodes. At the same time
 * collect qc reports which will be used later to
 * construct final QC report per batch.
 */
process DEMULTIPLEX {
    tag "$meta.id [$meta.amp_batch]"
    label 'process_medium'

    conda "bioconda::conda-forge==5.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(read), path(wells_cells), path(gene_intervals), path(spike_seq), path(oligos)

    output:
    tuple val(meta), path("output_tmp"), emit: folder
    path "versions.yml"                , emit: versions

    // required so we can append the results into output folder
    path("_DS/*")                      , optional: true
    path ("output/umi.tab/*.txt")
    path ("output/offset.tab/*.txt")
    path ("output/singleton_offset.tab/*.txt")
    path ("output/QC/read_stats/*.txt")
    path ("output/QC/read_stats_amp_batch/*.txt")
    path ("output/QC/umi_stats/*.txt")
    path ("output/QC/noffsets_per_umi_distrib/*.txt")
    path ("output/QC/nreads_per_umi_distrib/*.txt")
    path ("output/QC/umi_nuc_per_pos/*.txt")
    path ("output/_debug/${meta.amp_batch}/*.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
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
    mkdir -p _debug/${meta.amp_batch}/

    demultiplex.pl \\
        ${meta.amp_batch} \\
        ${meta.pool_barcode} \\
        $wells_cells \\
        $gene_intervals \\
        $spike_seq \\
        $oligos \\
        $read \\
        . \\
        $args

    mv _debug output/
    ln -s output output_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        demultiplex.pl: \$( demultiplex.pl --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/umi.tab/
    mkdir -p output/offset.tab/
    mkdir -p output/singleton_offset.tab/
    mkdir -p output/QC/{read_stats,read_stats_amp_batch,umi_stats,noffsets_per_umi_distrib,nreads_per_umi_distrib,umi_nuc_per_pos}

    touch output/umi.tab/${meta.amp_batch}.txt
    touch output/offset.tab/${meta.amp_batch}.txt
    touch output/singleton_offset.tab/${meta.amp_batch}.txt
    touch output/QC/{read_stats,read_stats_amp_batch,umi_stats,noffsets_per_umi_distrib,nreads_per_umi_distrib,umi_nuc_per_pos}/${meta.amp_batch}.txt

    mkdir -p output/_debug/${meta.amp_batch}/
    touch output/_debug/${meta.amp_batch}/{offsets,UMIs}.txt

    ln -s output output_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        demultiplex.pl: \$( demultiplex.pl --version )
    END_VERSIONS
    """
}
