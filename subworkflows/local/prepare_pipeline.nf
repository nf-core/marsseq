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

    ch_reads = Channel.empty()
    ch_fastp_version = Channel.empty()

    // convert XLS metadata into txt format
    INIT ( batches, gtf, ercc_regions )

    if (params.aligner == "bowtie2") {
        // split fastq reads by predefined number of reads per fastq file
        ch_reads = FASTP_SPLIT ( batches ).out.reads
        ch_fastp_version = FASTP_SPLIT.out.version

        // verify that split was performed correctly
        // R1 and R2 should always have a same pair
        if ( ch_reads.last().length() % 2 != 0 ) {
            exit 1, 'Error while splitting FASTQ files. Read pairs don\'t match!'
        }
    }

    if (params.aligner == "hisat2") {
        // No need to split the files into smaller ones
        ch_reads = batches
    }

    emit:
    amp_batches            = INIT.out.amp_batches
    gene_intervals         = INIT.out.gene_intervals
    seq_batches            = INIT.out.seq_batches
    wells_cells            = INIT.out.wells_cells
    reads                  = ch_reads               // channel: [ val(meta), path: *.fastq.gz ]
    fastp_version          = ch_fastp_version       // path: *.version.txt
}
