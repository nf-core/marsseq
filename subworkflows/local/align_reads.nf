//
// Align reads using bowtie2
//

def modules = params.modules.clone()
params.options = [:]

include { BOWTIE2_ALIGN } from '../../modules/local/bowtie2/align/main' addParams( options: modules['bowtie2_align'] )
include { HISAT2_ALIGN  } from '../../modules/local/hisat2/align/main'  addParams( options: modules['hisat2_align'] )
include { CUT_SAM       } from '../../modules/local/cut/sam/main'       addParams( options: modules['cut_sam'] )
include { QC_ALIGNED    } from '../../modules/local/qc/align/main'      addParams( options: modules['qc_aligned'] )


workflow ALIGN_READS {
    take:
    index      // channel [ bowtie2 index ]
    read       // channel [ meta, reads ]
    qc         // channel file(*.txt)

    main:
    ch_reads = read.combine(index)
    ch_sams = Channel.empty()
    ch_aligner_version = Channel.empty()

    BOWTIE2_ALIGN ( ch_reads )

    ch_sams = BOWTIE2_ALIGN.out.sam
    ch_aligner_version = BOWTIE2_ALIGN.out.version

    QC_ALIGNED ( ch_sams, qc )

    CUT_SAM ( ch_sams )
    
    emit:
    sam             = CUT_SAM.out.sam
    aligner_version = ch_aligner_version  // path: *.version.txt
}