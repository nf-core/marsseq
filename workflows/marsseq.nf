/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMarsseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.gtf ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { MERGE_SAMS            } from '../modules/local/cat/sam/main'          addParams( options: [:] )
include { QC_REPORT             } from '../modules/local/qc/report/main'        addParams( options: [:] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'       addParams( options: [:] )
include { PREPARE_PIPELINE  } from '../subworkflows/local/prepare_pipeline'  addParams( options: [:] )
include { LABEL_READS       } from '../subworkflows/local/label_reads'       addParams( options: [:] )
include { ALIGN_READS       } from '../subworkflows/local/align_reads'       addParams( options: [:] )
include { DEMULTIPLEX_READS } from '../subworkflows/local/demultiplex_reads' addParams( options: [:] )
include { VELOCITY          } from '../subworkflows/local/velocity'          addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// MODULE: local modules

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MARSSEQ {

    ch_software_versions    = Channel.empty()
    ch_fasta                = Channel.from(params.fasta)
    ch_gtf                  = Channel.from(params.gtf)
    ch_ercc_regions         = Channel.from("$projectDir/data/ercc-regions.tsv")
    ch_oligos               = Channel.from("$projectDir/data/oligos.txt")
    ch_spike_seq            = Channel.from("$projectDir/data/spike-seq.txt")
    ch_spike_concentrations = Channel.from("$projectDir/data/spike-concentrations.txt")
    ch_bowtie2_index        = Channel.from(WorkflowMain.getGenomeAttribute(params, 'bowtie2'))

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ( ch_input )
        .map { meta, reads -> [ meta, reads ] }
        .set { ch_batches }

    PREPARE_PIPELINE ( ch_batches, ch_gtf, ch_ercc_regions )

    LABEL_READS ( 
        ch_oligos, 
        PREPARE_PIPELINE.out.amp_batches, 
        PREPARE_PIPELINE.out.seq_batches, 
        PREPARE_PIPELINE.out.fastp_reads
    )

    ALIGN_READS ( ch_bowtie2_index, LABEL_READS.out.read, LABEL_READS.out.qc )

    // merge sam files into one file
    ALIGN_READS.out.sam
        .map { meta, sam -> [ meta.id, sam ] }
        .groupTuple(by: [0], sort: { it.name })
        .set { ch_sams }

    MERGE_SAMS ( ch_sams )

    DEMULTIPLEX_READS ( 
        MERGE_SAMS.out.sam,
        PREPARE_PIPELINE.out.amp_batches,
        PREPARE_PIPELINE.out.seq_batches,
        PREPARE_PIPELINE.out.wells_cells,
        PREPARE_PIPELINE.out.gene_intervals,
        ch_spike_seq,
        ch_spike_concentrations,
        ch_oligos
    )

    DEMULTIPLEX_READS.out.qc_rd.groupTuple()
        .join(DEMULTIPLEX_READS.out.qc_pdf.groupTuple())
        .combine(PREPARE_PIPELINE.out.amp_batches)
        .combine(PREPARE_PIPELINE.out.wells_cells)
        .set { ch_qcs }

    QC_REPORT ( ch_qcs )

    //
    // MODULE: Velocity
    //
    VELOCITY ( PREPARE_PIPELINE.out.fastp_reads )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions
        .mix(PREPARE_PIPELINE.out.fastp_version.ifEmpty(null))
        .mix(FASTQC.out.version.first().ifEmpty(null))
        .mix(ALIGN_READS.out.bowtie2_version.ifEmpty(null))
        .mix(VELOCITY.out.star_version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMarsseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
