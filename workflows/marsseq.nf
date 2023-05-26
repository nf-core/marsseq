/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMarsseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.gtf ]
for (param in checkPathParamList) { if (param && !params.build_references) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { CAT_CAT as MERGE_READS } from '../modules/nf-core/cat/cat/main'
include { QC_REPORT              } from '../modules/local/qc/report/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { PREPARE_PIPELINE  } from '../subworkflows/local/prepare_pipeline'
include { LABEL_READS       } from '../subworkflows/local/label_reads'
include { ALIGN_READS       } from '../subworkflows/local/align_reads'
include { DEMULTIPLEX_READS } from '../subworkflows/local/demultiplex_reads'
include { VELOCITY          } from '../subworkflows/local/velocity'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MARSSEQ {

    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_fasta                = file(params.fasta, checkIfExists: true)
    ch_gtf                  = file(params.gtf, checkIfExists: true)
    ch_bowtie_index         = file(params.bowtie2_index, checkIfExists: true)
    ch_star_index           = file(params.star_index, checkIfExists: true)
    ch_ercc_regions         = Channel.fromPath("$projectDir/data/ercc-regions.tsv")
    ch_oligos               = Channel.fromPath("$projectDir/data/oligos.txt")
    ch_spike_seq            = Channel.fromPath("$projectDir/data/spike-seq.txt")
    ch_spike_concentrations = Channel.fromPath("$projectDir/data/spike-concentrations.txt")

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ( ch_input )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    PREPARE_PIPELINE ( INPUT_CHECK.out.reads, ch_gtf, ch_ercc_regions )
    ch_versions = ch_versions.mix(PREPARE_PIPELINE.out.versions)

    LABEL_READS ( 
        ch_oligos, 
        PREPARE_PIPELINE.out.amp_batches, 
        PREPARE_PIPELINE.out.seq_batches, 
        PREPARE_PIPELINE.out.reads
    )

    ALIGN_READS ( LABEL_READS.out.read, ch_bowtie_index, LABEL_READS.out.qc )
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions)

    // merge sam files into one file
    ch_aligned_reads = ALIGN_READS.out.reads
        .map { meta, sam -> [ meta.id, sam ] }
        .groupTuple(by: [0], sort: { it.name })
        .map { batch, sams -> [ [ "id": batch ], sams ] }

    // merged aligned SAM files
    MERGE_READS ( ch_aligned_reads )
    ch_versions = ch_versions.mix(MERGE_READS.out.versions)

    DEMULTIPLEX_READS ( 
        MERGE_READS.out.file_out,
        PREPARE_PIPELINE.out.amp_batches,
        PREPARE_PIPELINE.out.seq_batches,
        PREPARE_PIPELINE.out.wells_cells,
        PREPARE_PIPELINE.out.gene_intervals,
        ch_spike_seq,
        ch_spike_concentrations,
        ch_oligos
    )

    // QC report from demultiplexing
    ch_qcs = DEMULTIPLEX_READS.out.qc_rd
        .join(DEMULTIPLEX_READS.out.qc_pdf)
        .combine(PREPARE_PIPELINE.out.amp_batches)
        .combine(PREPARE_PIPELINE.out.wells_cells)

    QC_REPORT ( ch_qcs )

    //
    // MODULE: Velocity
    //
    if (params.velocity) {
        VELOCITY ( PREPARE_PIPELINE.out.reads, ch_star_index )
        ch_versions = ch_versions.mix(VELOCITY.out.versions)
        
        ch_multiqc_files = ch_multiqc_files.mix(VELOCITY.out.catadapt_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(VELOCITY.out.star_multiqc.collect{it[1]}.ifEmpty([]))
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMarsseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMarsseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_PIPELINE.out.fastp_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_READS.out.bowtie2_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
