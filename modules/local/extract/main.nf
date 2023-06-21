/*
 * Label R1 using labels from R2.
 * Necessary for demultiplexing step anc QC.
 */
process EXTRACT_LABELS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::conda-forge==5.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads), path(oligos), path(amp_batches), path(seq_batches)

    output:
    tuple val(meta), path("labeled_reads/*.fastq"), emit: labeled_read
    path "labeled_reads/*.txt"                    , emit: qc
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def r1 = reads[0].baseName - '.gz'
    def r2 = reads[1].baseName - '.gz'
    def qc = r1 - '.fastq' + '.txt'
    """
    gunzip -f $reads
    mkdir labeled_reads

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_labels.pl: \$( extract_labels.pl --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def r1 = reads[0].baseName - '.gz'
    def qc = r1 - '.fastq' + '.txt'
    """
    mkdir labeled_reads
    touch labeled_reads/${r1}
    touch labeled_reads/${qc}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_labels.pl: \$( extract_labels.pl --version )
    END_VERSIONS
    """
}
