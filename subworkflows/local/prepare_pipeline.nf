//
// Subworkflow for setting up all necessary files 
// before running the pipeline
//

def modules = params.modules.clone()
params.options = [:]

include { PREPARE_PIPELINE as INIT } from '../../modules/local/prepare_pipeline' addParams( options: [:] )
include { FASTP_SPLIT              } from '../../modules/local/fastp/split/main' addParams( options: modules['fastp'] )

workflow PREPARE_PIPELINE {
    take:
    batches // channel: [ val(meta), [ reads ] ]

    main:

    // convert XLS metadata into txt format
    INIT ( batches )

    // split fastq reads by predefined number of reads per fastq file
    FASTP_SPLIT ( batches )

    emit:
    data        = INIT.out.txts          // channel: [ val(meta), path: *.txts ]
    fastp_reads = FASTP_SPLIT.out.reads  // channel: [ val(meta), path: *.fastq.gz ]
}
