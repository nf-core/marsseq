//
// Demultiplex reads
//
include { DEMULTIPLEX } from '../../modules/local/demultiplex/main'
include { QC_BATCH    } from '../../modules/local/qc/batch/main'


workflow DEMULTIPLEX_READS {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    amp_batches             // channel: amp_batches.txt
    seq_batches             // channel: seq_batches.txt
    wells_cells             // channel: wells_cells.txt
    gene_intervals          // channel: gene_intervals.txt
    spike_seq               // channel: spike-seq.txt
    spike_concentrations    // channel: spike-concentrations.txt
    oligos                  // channel: oligos.txt

    main:
    ch_versions = Channel.empty()

    amp_batch = amp_batches
        .splitCsv( header:true, sep:'\t' )
        .map { row -> [ row['Amp_batch_ID'], row['Pool_barcode'] ] }
        .combine(reads)
        .map { amp_batch, pool_barcode, meta, reads -> [
            [ "id": meta.id, "amp_batch": amp_batch, "pool_barcode": pool_barcode ], reads ]
        }
        .combine(wells_cells)
        .combine(gene_intervals)
        .combine(spike_seq)
        .combine(oligos)

    ch_demultiplex = DEMULTIPLEX ( amp_batch ).folder
        .combine(amp_batches)
        .combine(seq_batches)
        .combine(wells_cells)
        .combine(spike_concentrations)
    ch_versions = ch_versions.mix(DEMULTIPLEX.out.versions)

    QC_BATCH ( ch_demultiplex )
    ch_versions = ch_versions.mix(QC_BATCH.out.versions)

    emit:
    qc_rd    = QC_BATCH.out.rd
    qc_pdf   = QC_BATCH.out.pdf
    versions = ch_versions
}
