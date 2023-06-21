//
// Subworkflow for downloading all
// required references
//

include { WGET as DOWNLOAD_FASTA } from '../../modules/local/wget/main'
include { WGET as DOWNLOAD_GTF   } from '../../modules/local/wget/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../../modules/nf-core/gunzip/main'
include { ERCC_CREATE            } from '../../modules/local/ercc/main'
include { CAT_CAT as CAT_FASTA   } from '../../modules/nf-core/cat/cat/main'
include { BOWTIE2_BUILD          } from '../../modules/nf-core/bowtie2/build/main'
include { STAR_GENOMEGENERATE    } from '../../modules/nf-core/star/genomegenerate/main'

workflow BUILD_REFERENCES {

    main:
    ch_star_index = Channel.empty()
    ch_versions   = Channel.empty()
    params.genomes_base = params.genomes_base ? "${params.outdir}/references/" : params.genomes_base

    // download references
    fasta = params.genomes[params.genome].fasta.split('/')[-1]
    DOWNLOAD_FASTA (
        params.genomes[params.genome].fasta_url,
        "_" + fasta
    )
    ch_versions = ch_versions.mix(DOWNLOAD_FASTA.out.versions)

    DOWNLOAD_GTF (
        params.genomes[params.genome].gtf_url,
        params.genomes[params.genome].gtf.split('/')[-1]
    )
    ch_versions = ch_versions.mix(DOWNLOAD_GTF.out.versions)

    // uncompress
    ch_fasta = GUNZIP_FASTA ( DOWNLOAD_FASTA.out.file ).gunzip.map { it[1] }
    ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)

    ch_gtf = GUNZIP_GTF ( DOWNLOAD_GTF.out.file ).gunzip.map { it[1] }
    ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

    // create ERCC FASTA
    ERCC_CREATE( Channel.from("$projectDir/data/spike-seq.txt") )
    ch_versions = ch_versions.mix(ERCC_CREATE.out.versions)

    ch_fastas = ch_fasta.merge(ERCC_CREATE.out.fasta)
        .map{ it -> [ ["id": "${fasta - '.fa'}"], it ] }

    ch_genome = CAT_FASTA ( ch_fastas ).file_out
    ch_versions = ch_versions.mix(CAT_FASTA.out.versions)

    // build bowtie2 index
    BOWTIE2_BUILD( ch_genome )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    // build STAR index for velocity
    if (params.velocity) {
        ch_star_index = STAR_GENOMEGENERATE( ch_genome.map{ meta, fasta -> fasta }, ch_gtf.map{ meta, gtf -> gtf } ).index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
    fasta         = ch_fasta
    gtf           = ch_gtf
    bowtie2_index = BOWTIE2_BUILD.out.index
    star_index    = ch_star_index
    versions      = ch_versions

}
