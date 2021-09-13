//
// Subworkflow for setting up all necessary files 
// before running the pipeline
//

def modules = params.modules.clone()
params.options = [:]

include { PREPARE_PIPELINE as INIT } from '../../modules/local/prepare/pipeline/main' addParams( options: [:] )
include { FASTP_SPLIT              } from '../../modules/local/fastp/split/main'      addParams( options: modules['fastp'] )

workflow PREPARE_PIPELINE {
    take:
    batches       // channel: [ val(meta), [ reads ] ]
    gtf           // channel: gtf  
    ercc_regions  // channel: ercc_regions

    main:

    // convert XLS metadata into txt format
    INIT ( batches, gtf, ercc_regions )

    // split fastq reads by predefined number of reads per fastq file
    FASTP_SPLIT ( batches )

    // verify that split was performed correctly
    // R1 and R2 should always have a same pair
    if ( FASTP_SPLIT.out.reads.last().length() % 2 != 0 ) {
        exit 1, 'Error while splitting FASTQ files. Read pairs don\'t match!'
    }

    emit:
    amp_batches            = INIT.out.amp_batches
    gene_intervals         = INIT.out.gene_intervals
    seq_batches            = INIT.out.seq_batches
    wells_cells            = INIT.out.wells_cells
    fastp_reads            = FASTP_SPLIT.out.reads    // channel: [ val(meta), path: *.fastq.gz ]
    fastp_version          = FASTP_SPLIT.out.version  // path: *.version.txt
}
