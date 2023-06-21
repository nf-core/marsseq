//
// Extract tags from reads and move them from R2 to R1
//

include { EXTRACT_LABELS } from '../../modules/local/extract/main'

workflow LABEL_READS {
    take:
    oligos      // channel oligos.txt
    amp_batches // channel amp_batches.txt
    seq_batches // channel seq_batches.txt
    fastp_reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_reads = fastp_reads
        .map { meta, reads -> return [ meta, reads.sort().collate(2) ] }
        .transpose()
        .combine(oligos)
        .combine(amp_batches)
        .combine(seq_batches)

    EXTRACT_LABELS ( ch_reads )

    emit:
    read     = EXTRACT_LABELS.out.labeled_read
    qc       = EXTRACT_LABELS.out.qc
    versions = EXTRACT_LABELS.out.versions

}
