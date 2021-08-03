//
// Extract tags from reads and move them from R2 to R1
//

def modules = params.modules.clone()
params.options = [:]

include { EXTRACT_LABELS } from '../../modules/local/extract/main' addParams( options: modules['extract_labels'] )


workflow LABEL_READS {
    take:
    oligos      // channel oligos.txt
    amp_batches // channel amp_batches.txt
    seq_batches // channel seq_batches.txt
    fastp_reads // channel: [ val(meta), [ reads ] ]

    main:

    fastp_reads
        .map { meta, reads -> return [ meta, reads.sort().collate(2) ] }
        .transpose()
        .set { ch_reads }

    ch_reads
        .combine(oligos)
        .combine(amp_batches)
        .combine(seq_batches)
        .set { ch_reads }

    EXTRACT_LABELS ( ch_reads )

    emit:
    read = EXTRACT_LABELS.out.labeled_read
    qc   = EXTRACT_LABELS.out.qc

}