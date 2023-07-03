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

include { WGET as DOWNLOAD_FASTA      } from '../modules/local/wget/main'
include { WGET as DOWNLOAD_GTF        } from '../modules/local/wget/main'
include { GUNZIP as GUNZIP_FASTA      } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF        } from '../modules/nf-core/gunzip/main'
include { ERCC_CREATE                 } from '../modules/local/ercc/main'
include { CAT_CAT as CAT_FASTA        } from '../modules/nf-core/cat/cat/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/bowtie2/build/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/star/genomegenerate/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    ch_versions = Channel.empty()

    // download references
    fasta = params.genomes[params.genome].fasta.split('/')[-1]
    DOWNLOAD_FASTA (
        params.genomes[params.genome].fasta_url,
        "_" + fasta
    )
    DOWNLOAD_GTF (
        params.genomes[params.genome].gtf_url,
        params.genomes[params.genome].gtf.split('/')[-1]
    )

    // uncompress
    ch_fasta = GUNZIP_FASTA ( DOWNLOAD_FASTA.out.file )
        .gunzip
        .map { meta, fasta -> fasta }
    ch_gtf = GUNZIP_GTF ( DOWNLOAD_GTF.out.file ).gunzip

    // create ERCC FASTA
    ch_ercc = ERCC_CREATE( Channel.from("$projectDir/data/spike-seq.txt") ).fasta

    ch_fastas = ch_fasta.merge(ch_ercc)
        .map{ it -> [ ["id": "${fasta - '.fa'}"], it ] }

    ch_genome = CAT_FASTA ( ch_fastas ).file_out

    // build bowtie2 index
    BOWTIE2_BUILD( ch_genome )

    // build STAR index for velocity
    if (params.velocity) {
        STAR_GENOMEGENERATE( ch_genome.map{ meta, fasta -> fasta }, ch_gtf.map{ meta, gtf -> gtf } )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    // gather versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions
            .mix(DOWNLOAD_FASTA.out.versions)
            .mix(GUNZIP_FASTA.out.versions)
            .mix(GUNZIP_GTF.out.versions)
            .mix(ERCC_CREATE.out.versions)
            .mix(CAT_FASTA.out.versions)
            .mix(BOWTIE2_BUILD.out.versions)
            .unique()
            .collectFile(name: 'collated_versions.yml')
    )

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
