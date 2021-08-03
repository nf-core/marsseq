//
// Align reads using bowtie2
//

def modules = params.modules.clone()
params.options = [:]

include { BOWTIE2_ALIGN   } from '../../modules/local/bowtie2/align/main' addParams( options: [:] )
include { QC_ALIGNED      } from '../../modules/local/qc/align/main'      addParams( options: [:] )
include { TRIM_READ       } from '../../modules/local/cut/sam/main'       addParams( options: [:] )


workflow ALIGN_READS {
    take:
    index      // channel [ bowtie2 index ]
    read       // channel [ meta, reads ]
    qc         // channel file(*.txt)

    main:
    ch_reads = read.combine(index)

    BOWTIE2_ALIGN ( ch_reads )

    QC_ALIGNED ( BOWTIE2_ALIGN.out.sam, qc )

    TRIM_READ ( BOWTIE2_ALIGN.out.sam )
    
    emit:
    sam             = TRIM_READ.out.sam
    bowtie2_version = BOWTIE2_ALIGN.out.version  // path: *.version.txt
}