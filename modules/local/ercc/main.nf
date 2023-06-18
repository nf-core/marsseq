process ERCC_CREATE {
    tag "ercc.fa"
    label "process_low"

    conda "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0':
        'biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' }"

    input:
    path spikeins

    output:
    path "ercc.fa",         emit: fasta
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    create_ercc_fasta.py --input $spikeins --output ercc.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ERCC_CREATE: \$( create_ercc_fasta.py --version )
    END_VERSIONS
    """

    stub:
    """
    touch ercc.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ERCC_CREATE: \$( create_ercc_fasta.py --version )
    END_VERSIONS
    """
}
