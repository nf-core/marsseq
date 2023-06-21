//
// Align reads using bowtie2
//
include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { CUT_SAM       } from '../../modules/local/cut/sam/main'
include { QC_ALIGNED    } from '../../modules/local/qc/align/main'


workflow ALIGN_READS {
    take:
    reads   // channel [ meta, reads ]
    index   // channel file(bowtie2 index)
    qc      // channel file(*.txt)

    main:
    ch_versions = Channel.empty()
    ch_reads = reads
        .map { meta, reads -> [ [ "id": meta.id, "single_end": true, "filename": reads.baseName ], reads ]}
    ch_index = ch_reads
        .map { meta, reads -> [ meta, index ] }

    BOWTIE2_ALIGN ( ch_reads, ch_index, false, false )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    QC_ALIGNED ( BOWTIE2_ALIGN.out.aligned, qc )
    ch_versions = ch_versions.mix(QC_ALIGNED.out.versions)

    CUT_SAM ( BOWTIE2_ALIGN.out.aligned )
    ch_versions = ch_versions.mix(CUT_SAM.out.versions)

    emit:
    reads           = CUT_SAM.out.sam
    bowtie2_multiqc = BOWTIE2_ALIGN.out.log
    versions        = ch_versions
}
