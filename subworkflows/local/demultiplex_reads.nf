//
// Demultiplex reads
//

def modules = params.modules.clone()
params.options = [:]

include { DEMULTIPLEX } from '../../modules/local/demultiplex/main' addParams( options: modules['demultiplex'] )
include { QC_BATCH    } from '../../modules/local/qc/batch/main'    addParams( options: modules['demultiplex'] )


workflow DEMULTIPLEX_READS {
    take:
    sam_reads
    amp_batches
    seq_batches
    wells_cells
    gene_intervals
    spike_seq
    spike_concentrations
    oligos

    main:

    amp_batches
        .splitCsv( header:true, sep:'\t' )
        .map { row -> [ row['Amp_batch_ID'], row['Pool_barcode'] ] }
        .combine(sam_reads)
        .combine(wells_cells)
        .combine(gene_intervals)
        .combine(spike_seq)
        .combine(oligos)
        .set { amp_batch }

    DEMULTIPLEX ( amp_batch )

    DEMULTIPLEX.out.umi_tab.map { meta, f -> [ meta, f.Parent ] }
        .join( DEMULTIPLEX.out.offset_tab.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.singleton_offset_tab.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.qc_read_stats.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.qc_read_stats_amp_batch.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.qc_umi_stats.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.qc_umi_nuc_per_pos.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.qc_noffsets_per_umi_distrib.map { meta, f -> [ meta, f.Parent ] } )
        .join( DEMULTIPLEX.out.qc_nreads_per_umi_distrib.map { meta, f -> [ meta, f.Parent ] } )
        .map { meta, ut, ot, sot, qrs, qrsab, qus, qunpp, qnpud, qnpud2 -> [ meta, [ ut, ot, sot, qrs, qrsab, qus, qunpp, qnpud, qnpud2 ] ] }
        .combine(sam_reads.map { meta, sam -> meta })
        .combine(amp_batches)
        .combine(seq_batches)
        .combine(wells_cells)
        .combine(spike_concentrations)
        .set { ch_demultiplex_output }

    QC_BATCH ( ch_demultiplex_output )

    emit:
    qc_rd  = QC_BATCH.out.rd
    qc_pdf = QC_BATCH.out.pdf
}
