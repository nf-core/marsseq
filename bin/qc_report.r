#!/usr/bin/env Rscript
suppressMessages(library(gplots))

get_stats_per_seq_batch = function(seq_batch) {
    n_umis_per_cell = colSums(umitab)
    cell_mask = intersect(
        rownames(sample_list)[sample_list$sequencing_batch %in%seq_batch],
        single_cell_mask
    )
    n_cells = length(cell_mask)
    avg_mapping_stats = colMeans(t(mapping_stats)[cell_mask, ], na.rm = T)
    avg_n_UMI = mean(n_umis_per_cell[cell_mask])
    total_n_UMI = sum(n_umis_per_cell[cell_mask])
    names(avg_mapping_stats) = c("avg_n_mouse_reads", "avg_n_ecoli_reads", "avg_n_ercc_reads",
        "avg_n_raw_reads")
    total_n_raw_reads = as.numeric(avg_mapping_stats["avg_n_raw_reads"] * n_cells)
    avg_n_reads_per_UMI = total_n_raw_reads/total_n_UMI
    return(round(c(n_cells = n_cells, avg_mapping_stats, avg_n_UMIs = avg_n_UMI,
        total_n_UMIs = total_n_UMI, total_n_raw_reads = total_n_raw_reads, avg_n_reads_per_UMI = avg_n_reads_per_UMI)))
}

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
    well_cells_txt = args[1]
    amp_batches_txt = args[2]
    output_dir = args[3]
} else {
    stop("usage: Rscript qc_report.r [well_cells.txt] [amp_batches.txt] [output]")
}

sample_list = read.table(well_cells_txt, stringsAsFactors = F, header = T, sep = "\t")
amp_batches = read.delim(amp_batches_txt, stringsAsFactors = F, header = T, sep = "\t")

rownames(sample_list) = sample_list$Well_ID

main_stats_list = list()

for (batch in amp_batches$Amp_batch_ID) {
    try({
        load(paste(output_dir, "/", batch, ".rd", sep = ""))
        main_stats_list[[batch]] = main_stats
    }, silent = T)
}

m = t(sapply(main_stats_list, function(x) {
    return(x[names(main_stats_list[[1]])])
}))

pdf(paste("_temp/additional_qc.pdf", sep = ""))
m2 = t(sapply(main_stats_list, function(x) {
    return(x[c("gene_umis", "gene_umis_neg_control", "avg_noffsets_per_umi", "avg_noffsets_per_umi_neg_ctrl")])
}))

# UMIs vs. QC1
try({
    plot(amp_batches$QC1..pooled.cDNA., m2[match(amp_batches[, 1], rownames(m2)),
        1], xlab = "QC1", ylab = "average #UMIs", pch = "", log = "y")
    text(amp_batches$QC1..pooled.cDNA., m2[match(amp_batches[, 1], rownames(m2)),
        1], labels = amp_batches[, 1], cex = 0.5)
})

# UMIs vs. QC2
try({
    plot(amp_batches$QC2..RT.2., m2[match(amp_batches[, 1], rownames(m2)), 1], xlab = "QC2",
        ylab = "average #UMIs", pch = "", log = "y")
    text(amp_batches$QC2..RT.2., m2[match(amp_batches[, 1], rownames(m2)), 1], labels = amp_batches[,
        1], cex = 0.5)
})

# UMIs vs. QC3
try({
    plot(amp_batches$QC3..library.Ct., m2[match(amp_batches[, 1], rownames(m2)),
        1], xlab = "QC3", ylab = "average #UMIs", pch = "", log = "y")
    text(amp_batches$QC3..library.Ct., m2[match(amp_batches[, 1], rownames(m2)),
        1], labels = amp_batches[, 1], cex = 0.5)
})
garbage <- dev.off()

for (seq_batch in unique(amp_batches$Seq_batch_ID)) {
    cur_amp_batches = amp_batches$Amp_batch_ID[amp_batches$Seq_batch_ID == seq_batch]
    pdfs = paste(output_dir, "/", cur_amp_batches, ".pdf", sep = "")
    pdfs = pdfs[file.exists(pdfs)]
    cmd = paste(paste("gs -dBATCH -dNOPAUSE -dAutoRotatePages=/None -q -sDEVICE=pdfwrite -sOutputFile=",
        "output/QC_reports/qc_", seq_batch, ".pdf", sep = ""),
        paste(pdfs, collapse = " "))

    system(cmd)
}

df = data.frame(amp_batch = rownames(m), m)
write.table(file = paste(output_dir, "/amp_batches_stats.txt", sep = ""),
    df, row.names = F, col.names = T, sep = "\t", quote = F)

rownames(amp_batches) = amp_batches[, "Amp_batch_ID"]
m4 = cbind(
    amp_batch_ID = rownames(m),
    # amp_batches[rownames(m), c("Experiment_ID", "Description", "QC2..RT.2.")],
    amp_batches[rownames(m), c("Experiment_ID", "Description")],
    m[, match(c("reads", "spike_yield", "gene_umis", "gene_umis_neg_control", "noise_estimation", "avg_noffsets_per_umi", "avg_reads_per_umi"), colnames(m))]
)

write.table(file = paste(output_dir, "/amp_batches_summary.txt", sep = ""), m4, row.names = F, col.names = T, sep = "\t", quote = F)
