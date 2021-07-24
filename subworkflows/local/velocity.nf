//
// Subworkflow for setting up all necessary files 
// before running the pipeline
//

def modules = params.modules.clone()
params.options = [:]

include { VELOCITY_CONVERT                } from '../../modules/local/velocity_convert'     addParams( options: [:] )
include { VELOCITY_WHITELIST              } from '../../modules/local/velocity_whitelist'   addParams( options: [:] )
include { CUTADAPT as VELOCITY_TRIM       } from '../../modules/local/cutadapt/main'        addParams( options: modules['cutadapt'] )
include { STAR_ALIGN as VELOCITY_STARSOLO } from '../../modules/local/star/align/main'      addParams( options: modules['star_align'] )

workflow VELOCITY {
    take:
    fastp_reads

    main:

    // convert fastq files into 10X format
    VELOCITY_CONVERT ( fastp_reads )

    // build whitelist.txt
    VELOCITY_WHITELIST ( fastp_reads )

    // trim poly-T and low quality reads
    VELOCITY_TRIM ( VELOCITY_CONVERT.out.reads )
    
    // alignment using StarSolo
    VELOCITY_STARSOLO ( VELOCITY_TRIM.out.reads, VELOCITY_WHITELIST.out.whitelist )

}
