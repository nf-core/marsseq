//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    reads // channel: [ val(batch), [reads] ]
}

def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id  = row.batch

    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }

    if (!file(row.amp_batches).exists()) {
        exit 1, "ERROR: amp_batch file does not exist!"
    }

    if (!file(row.seq_batches).exists()) {
        exit 1, "ERROR: seq_batches file does not exist!"
    }

    if (!file(row.well_cells).exists()) {
        exit 1, "ERROR: well_cells file does not exist!"
    }

    meta.amp_batches = file(row.amp_batches)
    meta.seq_batches = file(row.seq_batches)
    meta.well_cells = file(row.well_cells)

    return [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
}
