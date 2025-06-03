//
// Align reads using bowtie2
//
include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { CUT_SAM       } from '../../modules/local/cut/sam/main'
include { QC_ALIGNED    } from '../../modules/local/qc/align/main'


workflow ALIGN_READS {
    take:
    reads   // channel: [ meta, reads ]
    index   // channel: [ index ]
    fasta   // channel: [ meta, fasta ]
    qc      // channel: [ file(*.txt) ]

    main:
    ch_versions = Channel.empty()
    ch_reads = reads
        .map { meta, reads -> [ [ "id": meta.id, "single_end": true, "filename": reads.baseName ], reads ]}

    BOWTIE2_ALIGN ( ch_reads, index, fasta, false, false )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    QC_ALIGNED ( BOWTIE2_ALIGN.out.sam, qc )
    ch_versions = ch_versions.mix(QC_ALIGNED.out.versions)

    CUT_SAM ( BOWTIE2_ALIGN.out.sam )
    ch_versions = ch_versions.mix(CUT_SAM.out.versions)

    emit:
    reads           = CUT_SAM.out.sam
    bowtie2_multiqc = BOWTIE2_ALIGN.out.log
    versions        = ch_versions
}
