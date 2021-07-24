// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"$meta.id/velocity/", meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::star=2.7.3a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/star:2.7.3a--0"
    } else {
        container "quay.io/biocontainers/star:2.7.3a--0"
    }

    input:
    tuple val(meta), path(reads)
    path  whitelist

    output:
    tuple val(meta), path('Solo.out/')         , emit: star_solo
    tuple val(meta), path('Aligned.out.sam')   , emit: sam
    tuple val(meta), path('Log.final.out')     , emit: log_final
    tuple val(meta), path('Log.out')           , emit: log_out
    tuple val(meta), path('Log.progress.out')  , emit: log_progress
    tuple val(meta), path('SJ.out.tab')        , optional:true, emit: tab
    path  '*.version.txt'                      , emit: version

    script:
    def index    = file(WorkflowMain.getGenomeAttribute(params, 'star'))
    def software = getSoftwareName(task.process)
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads[1]} ${reads[0]} \\
        --runThreadN $task.cpus \\
        --readFilesCommand zcat \\
        --soloCBwhitelist $whitelist \\
        $options.args
    
    STAR --version | sed -e "s/STAR_//g" > ${software}.version.txt
    """
}