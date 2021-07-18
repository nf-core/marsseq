// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process WGET {
    tag "$filename"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"references/${params.genome}", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gnu-wget=1.18" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5"
    } else {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    }

    input:
    tuple val(filename), val(url)

    output:
    path "*.version.txt"    , emit: version
    path "${filename}"      , emit: output_file

    script:
    def software = getSoftwareName(task.process)
    def original_filename = url.split("/")[-1]
    """
    wget $options.args $url -O $original_filename

    case $original_filename in
    *.gz )  
            gunzip -c $original_filename > $filename
        ;;
    * )
            echo "Extension not supported"
        ;;
    esac

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > ${software}.version.txt
    """
}
