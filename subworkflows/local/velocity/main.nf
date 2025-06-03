//
// Subworkflow for setting up all necessary files
// before running the pipeline
//
include { CAT_FASTQ                       } from '../../modules/nf-core/cat/fastq/main'
include { CUTADAPT as VELOCITY_TRIM       } from '../../modules/nf-core/cutadapt/main'
include { VELOCITY_CONVERT                } from '../../modules/local/velocity/convert/main'
include { VELOCITY_WHITELIST              } from '../../modules/local/velocity/whitelist/main'
include { STARSOLO                        } from '../../modules/nf-core/star/starsolo/main'

workflow VELOCITY {
    take:
    amp_batches   // channel: amp_batch
    well_cells    // channel: well_cells
    reads         // channel: [ meta, reads ]
    index         // channel: [ meta, index ]

    main:
    ch_versions = Channel.empty()

    // convert fastq files into 10X format
    ch_reads = reads
        .map { meta, reads -> return [ meta, reads.sort().collate(2) ] }
        .transpose()

    VELOCITY_CONVERT ( ch_reads )
    ch_versions = ch_versions.mix(VELOCITY_CONVERT.out.versions)

    // merge fastq files into one sample
    CAT_FASTQ ( VELOCITY_CONVERT.out.reads.groupTuple().map { meta, reads -> [ meta, reads.flatten() ] } )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // build whitelist.txt
    VELOCITY_WHITELIST ( amp_batches, well_cells, reads )
    ch_versions = ch_versions.mix(VELOCITY_WHITELIST.out.versions)

    // trim poly-T and low quality reads
    VELOCITY_TRIM ( CAT_FASTQ.out.reads )
    ch_versions = ch_versions.mix(VELOCITY_TRIM.out.versions)

    // alignment using StarSolo
    ch_reads = VELOCITY_TRIM.out.reads.map {
        meta, reads -> [
            [
                id: meta.id,
                cb_start: 1,
                cb_len: 11,
                umi_start: 12,
                umi_len: 8,
            ], "CB_UMI_Simple", reads
        ]
    }

    STARSOLO ( ch_reads, VELOCITY_WHITELIST.out.whitelist, index )
    ch_versions = ch_versions.mix(STARSOLO.out.versions)

    emit:
    catadapt_multiqc = VELOCITY_TRIM.out.log
    star_multiqc     = STARSOLO.out.log_final
    versions         = ch_versions

}
