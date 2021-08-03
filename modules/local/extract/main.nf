// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process EXTRACT_LABELS {
    tag "$meta.id"
    label 'process_tiny'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"$meta.id", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::conda-forge==5.22.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/perl:5.22.2.1"
    } else {
        container "quay.io/biocontainers/perl:5.22.2.1"
    }

    input:
    tuple val(meta), path(reads), path(oligos), path(amp_batches), path(seq_batches)
    
    output:
    tuple val(meta), path("labeled_reads/*.fastq"), emit: labeled_read
    path "labeled_reads/*.txt"                    , emit: qc

    script:
    def r1 = reads[0].baseName - '.gz'
    def r2 = reads[1].baseName - '.gz'
    def qc = r1 - '.fastq' + '.txt'
    """
    gunzip $reads

    mkdir labeled_reads labeled_reads/qc

    extract_labels.pl \\
        $r1 \\
        $r2 \\
        $meta.id \\
        $seq_batches \\
        $oligos \\
        $amp_batches \\
        labeled_reads/$r1 \\
        labeled_reads/$qc \\
        . 
    """
}
