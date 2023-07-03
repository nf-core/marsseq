//
// Subworkflow for setting up all necessary files
// before running the pipeline
//
include { VELOCITY_CONVERT                } from '../../modules/local/velocity/convert/main'
include { VELOCITY_WHITELIST              } from '../../modules/local/velocity/whitelist/main'
include { CUTADAPT as VELOCITY_TRIM       } from '../../modules/nf-core/cutadapt/main'
include { STAR_ALIGN as VELOCITY_STARSOLO } from '../../modules/nf-core/star/align/main'

workflow VELOCITY {
    take:
    amp_batches   // channel: amp_batch
    well_cells    // channel: well_cells
    reads         // channel [ meta, reads ]
    index         // channel file(star index)
    gtf           // channel file(gtf)

    main:
    ch_versions = Channel.empty()

    ch_folder = reads.map { meta, reads -> [ meta, reads.first().Parent ] }
    ch_index = reads.map { meta, reads -> [ meta, index ] }

    // convert fastq files into 10X format
    VELOCITY_CONVERT ( ch_folder )
    ch_versions = ch_versions.mix(VELOCITY_CONVERT.out.versions)

    // build whitelist.txt
    VELOCITY_WHITELIST ( amp_batches, well_cells, reads )
    ch_versions = ch_versions.mix(VELOCITY_WHITELIST.out.versions)

    // trim poly-T and low quality reads
    VELOCITY_TRIM ( VELOCITY_CONVERT.out.reads )
    ch_versions = ch_versions.mix(VELOCITY_TRIM.out.versions)

    // alignment using StarSolo
    VELOCITY_STARSOLO (
        VELOCITY_TRIM.out.reads,            // reads
        index,                              // star index
        gtf,                                // gtf annotation
        VELOCITY_WHITELIST.out.whitelist,   // whitelist
        true,                               // star_ignore_sjdbgtf
        false,                              // seq_platform
        false                               // seq_center
    )
    ch_versions = ch_versions.mix(VELOCITY_STARSOLO.out.versions)

    emit:
    catadapt_multiqc = VELOCITY_TRIM.out.log
    star_multiqc     = VELOCITY_STARSOLO.out.log_final
    versions         = ch_versions

}
