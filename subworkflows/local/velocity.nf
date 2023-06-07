//
// Subworkflow for setting up all necessary files
// before running the pipeline
//
include { VELOCITY_CONVERT                } from '../../modules/local/velocity/convert/main'
include { VELOCITY_WHITELIST              } from '../../modules/local/velocity/whitelist/main'
include { CUTADAPT as VELOCITY_TRIM       } from '../../modules/nf-core/cutadapt/main'
include { STAR_ALIGN as VELOCITY_STARSOLO } from '../../modules/local/star/align/main'

workflow VELOCITY {
    take:
    reads   // channel [ meta, reads ]
    index   // channel file(star index)

    main:
    ch_versions = Channel.empty()

    ch_folder = reads.map { meta, reads -> [ meta, reads.first().Parent ] }
    ch_index = reads.map { meta, reads -> [ meta, index ] }

    // convert fastq files into 10X format
    VELOCITY_CONVERT ( ch_folder )

    // build whitelist.txt
    VELOCITY_WHITELIST ( reads )

    // trim poly-T and low quality reads
    VELOCITY_TRIM ( VELOCITY_CONVERT.out.reads )
    ch_versions = ch_versions.mix(VELOCITY_TRIM.out.versions)

    // alignment using StarSolo
    VELOCITY_STARSOLO ( VELOCITY_TRIM.out.reads, index, VELOCITY_WHITELIST.out.whitelist )
    ch_versions = ch_versions.mix(VELOCITY_STARSOLO.out.versions)

    emit:
    catadapt_multiqc = VELOCITY_TRIM.out.log
    star_multiqc     = VELOCITY_STARSOLO.out.log_final
    versions         = ch_versions

}
