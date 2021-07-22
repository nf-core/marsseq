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
    batches

    main:

    INIT ( batches )

    FASTP_SPLIT ( batches )

    emit:
    FASTP_SPLIT.out.reads
}
