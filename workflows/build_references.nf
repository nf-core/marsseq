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

include { WGET as DOWNLOAD_FASTA } from '../modules/local/wget/main'
include { WGET as DOWNLOAD_GTF   } from '../modules/local/wget/main'
include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../modules/nf-core/modules/gunzip/main'
include { CREATE_ERCC_FASTA      } from '../modules/local/prepare/ercc/main'
include { CAT_FASTA              } from '../modules/local/cat/fasta/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/modules/bowtie2/build/main'
include { STAR_GENOMEGENERATE    } from '../modules/local/star/genomegenerate/main'

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
    ch_fasta      = CAT_FASTA ( ch_fasta, ch_ercc_fasta ).fasta
    
    // build bowtie2 index
    BOWTIE2_BUILD( ch_fasta )

    // build STAR index for velocity
    if (params.velocity) {
        STAR_GENOMEGENERATE( ch_fasta, ch_gtf )
    }
    
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
