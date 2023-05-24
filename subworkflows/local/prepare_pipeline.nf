//
// Subworkflow for setting up all necessary files 
// before running the pipeline
//

include { PREPARE     } from '../../modules/local/prepare/main'
include { FASTP_SPLIT } from '../../modules/local/fastp/split/main'

workflow PREPARE_PIPELINE {
    take:
    batches       // channel: [ val(meta), [ reads ] ]
    gtf           // channel: gtf
    ercc_regions  // channel: ercc_regions

    main:
    ch_reads = Channel.empty()

    // convert XLS metadata into txt format
    PREPARE ( batches, gtf, ercc_regions )

    // split fastq reads by predefined number of reads per fastq file
    ch_reads = FASTP_SPLIT ( PREPARE.out.reads ).reads

    // verify that split was performed correctly
    // R1 and R2 should always have a same pair
    if ( ch_reads.last().length() % 2 != 0 ) {
        exit 1, 'Error while splitting FASTQ files. Read pairs don\'t match!'
    }

    emit:
    amp_batches    = PREPARE.out.amp_batches
    seq_batches    = PREPARE.out.seq_batches
    wells_cells    = PREPARE.out.wells_cells
    gene_intervals = PREPARE.out.gene_intervals
    reads          = ch_reads                   // channel: [ val(meta), path: *.fastq.gz ]
    fastp_multiqc  = FASTP_SPLIT.out.json
    versions       = FASTP_SPLIT.out.versions
}
