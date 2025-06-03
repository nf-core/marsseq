process WGET {
    tag "$filename"
    label "process_low"

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val url
    val filename

    output:
    tuple val(meta), path ("$outfile"), emit: file
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    // NOTE: get last extension, will fail for ther cases like tar.gz
    extension = url.split("\\.")[-1]
    outfile = "${filename}.${extension}"
    meta = [ "id": "$outfile" ]
    """
    wget $args $url -O $outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ""
    // NOTE: get last extension, will fail for ther cases like tar.gz
    extension = url.split("\\.")[-1]
    outfile = "${filename}.${extension}"
    meta = [ "id": "$outfile" ]
    """
    touch ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """
}
