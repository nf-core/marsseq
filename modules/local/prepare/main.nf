/*
 * Preparation process before executing pipeline.
 * Generates all necessary files like
 *  - amp_batches_to_process.txt
 *  - amp_batches.txt
 *  - seq_batches.txt
 *  - wells_cells.txt
 *  - gene_intervals.txt
 */
process PREPARE {
    tag "$meta.id"
    label 'process_tiny'

    conda "bioconda::openpyxl==2.6.1 conda-forge::pandas==1.2.4 conda-forge::fsspec==2023.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' :
        'biocontainers/mulled-v2-0bcca2890a3ab7be29a83e813a02d340d6f54660:4cb478c6e57df2ef85ea5f8eae6d717c017962cd-0' }"

    input:
    path(amp_batches)
    path(seq_batches)
    path(well_cells)
    path(gtf)
    path(ercc_regions)
    tuple val(meta), path(reads)

    output:
    path "amp_batches.txt"   , emit: amp_batches
    path "gene_intervals.txt", emit: gene_intervals
    path "seq_batches.txt"   , emit: seq_batches
    path "wells_cells.txt"   , emit: wells_cells
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_pipeline.py \\
        --batch ${meta.id} \\
        --amp_batches $amp_batches \\
        --seq_batches $seq_batches \\
        --well_cells $well_cells \\
        --gtf $gtf \\
        --output .
    cat $ercc_regions >> gene_intervals.txt
    validate_data.py --input .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_pipeline.py: \$( prepare_pipeline.py --version )
        validate_data.py: \$( validate_data.py --version )
    END_VERSIONS
    """

    stub:
    """
    cat <<AMP_BATCH > amp_batches.txt
    Amp_batch_ID\tSeq_batch_ID\tPool_barcode\tSpike_type\tSpike_dilution\tSpike_volume_ul\tExperiment_ID\tOwner\tDescription
    AB339\tSB26\tTGAT\tERCC_mix1\t2.5e-05\t0.01\tTECH_ES\tHadas\tES#7_poolA
    AMP_BATCH

    cat <<SEQ_BATCHES > seq_batches.txt
    Seq_batch_ID\tRun_name\tDate\tR1_design\tI5_design\tR2_design\tNotes
    SB26\tsc_v3_Hadas_Diego_05042015\t150405\t5I.4P.51M\t7W.8R\t\tmm10
    SEQ_BATCHES

    cat <<WELLS_CELLS > wells_cells.txt
    Well_ID\tWell_coordinates\tplate_ID\tSubject_ID\tAmp_batch_ID\tCell_barcode\tNumber_of_cells
    TW1\tA1\t154\t35\tAB339\tCTATTCG\t1
    WELLS_CELLS

    cat <<GENE_INTERVALS > gene_intervals.txt
    chrom\tstart\tend\tstrand\tgene_name
    chr1\t3143476\t3144545\t1\t4933401J01Rik
    GENE_INTERVALS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_pipeline.py: \$( prepare_pipeline.py --version )
        validate_data.py: \$( validate_data.py --version )
    END_VERSIONS
    """
}
