/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { WGET as DOWNLOAD_FASTA } from '../modules/local/wget/main'                    addParams( options: modules['wget'] )
include { WGET as DOWNLOAD_GTF   } from '../modules/local/wget/main'                    addParams( options: modules['wget'] )
include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/modules/gunzip/main'        addParams( options: modules['gunzip'] )
include { GUNZIP as GUNZIP_GTF   } from '../modules/nf-core/modules/gunzip/main'        addParams( options: modules['gunzip'] )
include { CREATE_ERCC_FASTA      } from '../modules/local/prepare/ercc/main'            addParams( options: modules['create_ercc_fasta'] )
include { MERGE_FASTA            } from '../modules/local/cat/fasta/main'               addParams( options: modules['merge_fasta'] )
include { BOWTIE2_BUILD          } from '../modules/nf-core/modules/bowtie2/build/main' addParams( options: modules['bowtie2_index'] )
include { STAR_GENOMEGENERATE    } from '../modules/local/star/genomegenerate/main'     addParams( options: modules['star_index'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    // download references
    ch_fetched_fasta = DOWNLOAD_FASTA( [ "${params.fasta.split('/')[-1]}.gz", WorkflowMain.getGenomeAttribute(params, 'fasta_url') ] ).file
    ch_fetched_gtf   = DOWNLOAD_GTF( [ "${params.gtf.split('/')[-1]}.gz", WorkflowMain.getGenomeAttribute(params, 'gtf_url') ] ).file

    // uncompress
    ch_fasta = GUNZIP_FASTA ( ch_fetched_fasta ).gunzip
    ch_gtf = GUNZIP_GTF ( ch_fetched_gtf ).gunzip
    
    // merge ERCC and reference genome
    ch_ercc_fasta = CREATE_ERCC_FASTA( Channel.from("$projectDir/data/spike-seq.txt") ).fasta
    ch_fasta      = MERGE_FASTA ( ch_fasta, ch_ercc_fasta ).fasta
    
    // build bowtie2 index for core alignment
    BOWTIE2_BUILD( ch_fasta )
    
    // build STAR index for velocity
    STAR_GENOMEGENERATE( ch_fasta, ch_gtf )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    // NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
