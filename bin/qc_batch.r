#!/usr/bin/env Rscript
suppressMessages(library(MASS))
suppressMessages(library(gplots))
suppressMessages(library(zoo))

downsamp_one = function(v, n) {
    hist(sample(rep(1:length(v), times = v), replace = F, size = n), 0.5 + 0:length(v),
        plot = F)$counts
}

plot_cumulative = function(x, y = NA, nbreaks = 50, col = 1, lty = 1, ylab = "",
    xlab = "", ...) {
    if (any(!is.na(y))) {
        minxy = min(c(x, y), na.rm = T)
        maxxy = max(c(x, y), na.rm = T)

        h1 = hist(x, breaks = seq(minxy, maxxy, l = nbreaks), plot = F)
        h2 = hist(y, breaks = seq(minxy, maxxy, l = nbreaks), plot = F)
        matplot(h2$breaks, cbind(c(0, cumsum(h1$counts)/sum(h1$counts)), c(0, cumsum(h2$counts)/sum(h2$counts))),
            type = "l", col = col, lwd = 2, ylab = ylab, lty = lty, xlab = xlab,
            ...)
        grid(lty = 1)
        matplot(h2$breaks, cbind(c(0, cumsum(h1$counts)/sum(h1$counts)), c(0, cumsum(h2$counts)/sum(h2$counts))),
            type = "l", col = col, lwd = 2, ylab = ylab, lty = lty, xlab = xlab,
            add = T, ...)
    } else {
        h1 = hist(x, breaks = nbreaks, plot = F)
        plot(h1$breaks, c(0, cumsum(h1$counts)/sum(h1$counts)), type = "l", panel.first = grid(lty = 1),
            col = col, lwd = 2, ylab = ylab, xlab = xlab, ...)
        return(median(x))
    }
}

make_amp_batch_qc_fig = function(amp_batch_ID, spike_concentrations_file, output_dir = "output") {
    main_stats = list()
    umitab_batch = read.table(paste(output_dir, "/umi.tab/", amp_batch_ID, ".txt",
        sep = ""), sep = "\t", header = T, check.names = F, row.names = 1)
    offsetab_batch = read.table(paste(output_dir, "/offset.tab/", amp_batch_ID, ".txt",
        sep = ""), sep = "\t", header = T, check.names = F, row.names = 1)
    singleoffsetab_batch = read.table(paste(output_dir, "/singleton_offset.tab/",
        amp_batch_ID, ".txt", sep = ""), sep = "\t", header = T, check.names = F,
        row.names = 1)
    spike_concentration <- read.table(as.character(spike_concentrations_file), header = T,
        sep = "\t", stringsAsFactors = F, row.names = 1)

    batch_mask = sample_list$Well_ID[sample_list$Amp_batch_ID == amp_batch_ID]
    batch_details = amp_batches[amp_batches$Amp_batch_ID == amp_batch_ID, ]
    batch_title = paste("batch #", amp_batch_ID, " (", paste(batch_details["Owner"],
        " Exp_ID=", batch_details["Experiment_ID"], sep = ""), ")", sep = "")

    mask1 = batch_mask[pmax(sample_list[batch_mask, "Number_of_cells"], -1) == 1]
    mask0 = batch_mask[pmax(sample_list[batch_mask, "Number_of_cells"], -1) == 0]
    if (length(mask1) <= 1) {
        warning("No single cells in amplification batch ", amp_batch_ID)
    }


    spike_names = rownames(spike_concentration)
    gene_mask = setdiff(rownames(umitab_batch), spike_names)

    read_stats = read.table(paste(output_dir, "/QC/read_stats/", amp_batch_ID, ".txt",
        sep = ""), header = T, row.names = 1)
    read_stats_total = read.table(paste(output_dir, "/QC/read_stats_amp_batch/",
        amp_batch_ID, ".txt", sep = ""), header = T)
    umi_data = read.table(paste(output_dir, "/QC/umi_stats/", amp_batch_ID, ".txt",
        sep = ""), header = T, row.names = 1)

    x = read_stats[mask1, ]

    layout(matrix(c(1, 1, 1, 1, 2:11, 12, 13, 14:17, 18, 18, 19, 19, 20:23), 7, 4,
        byrow = T), heights = c(1, rep(2, 6)))
    par(mar = c(0, 5, 3, 2), mgp = c(2, 1, 0))
    plot.new()

    title(batch_title, cex = 2)
    par(mar = c(4, 5, 3, 2))
    tot = x$total

    med = plot_cumulative(tot, main = "Sequencing", ylab = "Fraction of cells", xlab = "Number of reads")
    mtext(3, text = paste("median =", med, sep = " "), col = 2, cex = 0.7)
    mtext("A", side = 3, adj = -0.6, line = 2)
    message("A")

    main_stats$reads = round(read_stats_total$nreads_amp_batch/length(mask1))
    main_stats$reads_with_known_cell_barcode = round((read_stats_total$nreads_amp_batch -
        read_stats_total$nreads_with_unknown_cell_barcode)/length(mask1))
    main_stats$reads_bad_quality_well_barcode = round(read_stats_total$nreads_bad_quality_well_barcode/length(mask1))
    main_stats$reads.spike_mapped = round(sum(x[mask1, "spike_mapped"])/length(mask1))
    main_stats$reads.gene_mapped = round(sum(x[mask1, "gene_mapped"])/length(mask1))
    main_stats$gene_umis = round(mean(colSums(umitab_batch[gene_mask, mask1])))

    mbar = cbind(gene_mapped = x[, "gene_mapped"], gene_filtered = rowSums(x[, c("gene_UMI_err",
        "gene_barcode_err", "gene_lonely_offset_err", "gene_readpoor_offset_err")]),
        mapped_to_nongenic = x[, "mapped_to_nongenic"], spike_mapped = x[, "spike_mapped"],
        spike_filtered = rowSums(x[, c("spike_UMI_err", "spike_barcode_err", "spike_lonely_offset_err",
            "spike_readpoor_offset_err")]), low_mapq = x[, "low_mapq"], multimap = x[,
            "non_unique_mapping"])/tot
    mbar = mbar[order(mbar[, "gene_mapped"]), ]
    cols = c(colors()[132], "cyan4", "gray", "firebrick", "firebrick1", "green",
        "yellow")

    barplot(t(mbar), col = cols, border = F, names.arg = rep("", dim(mbar)[1]), main = "Mapping per cell",
        ylab = "Reads fraction", xlab = "Single cells", xlim = c(0, round(dim(mbar)[1] *
            2.2)))
    legend("bottomright", legend = c("gene", "gene filtered", "non-exonic", "spike",
        "spike_filtered", "low mapq", "multi-map"), pch = 15, col = cols, cex = 0.5,
        bg = "white")
    mtext("B", side = 3, adj = -0.6, line = 2)

    message("B")
    unmapped_stats = colSums(x[, setdiff(colnames(x), c("well_id", "low_mapq", "non_unique_mapping",
        "mapped_to_nongenic", "spike_UMI_err", "spike_barcode_err", "spike_lonely_offset_err",
        "spike_readpoor_offset_err", "spike_mapped", "gene_UMI_err", "gene_barcode_err",
        "gene_lonely_offset_err", "gene_readpoor_offset_err", "gene_mapped"))])

    unmapped_stats = c(unmapped_stats[-match("unmapped", names(unmapped_stats))],
        other = as.numeric(unmapped_stats["unmapped"]))
    unmapped_stats[-match("total", names(unmapped_stats))]/unmapped_stats["total"]
    unmapped_fracs = unmapped_stats[-match("total", names(unmapped_stats))]/unmapped_stats["total"]
    unmapped_fracs[unmapped_fracs > 0.001]

    barplot(rev(unmapped_fracs), horiz = T, las = 1, border = F, names = rep("",
        length(unmapped_fracs)), xlab = "fraction of total reads", main = "Unmapped Reads")
    mtext(side = 2, text = names(unmapped_fracs), at = seq(length(unmapped_fracs) +
        1.1, 0.7, length = length(unmapped_fracs)), las = T, cex = 0.6)
    mtext("C", side = 3, adj = -0.6, line = 2)
    message("C")

    m = read.table(paste(output_dir, "/QC/umi_nuc_per_pos/", amp_batch_ID, ".txt",
        sep = ""), header = T, row.names = 1)
    freq_per_pos = m/rowSums(m)
    matplot(freq_per_pos, pch = c("A", "C", "G", "T"), col = c(3, 4, "orange", "red"),
        cex = 0.8, ylab = "Freq.", xlab = "position")
    title("UMI nuc/position")
    mtext("D", side = 3, adj = -0.6, line = 2)
    message("D")

    genetot = rowSums(umi_data[mask1, c("gene_UMI_err", "gene_barcode_err", "gene_readpoor_offset_err",
        "gene_lonely_offset_err", "gene_ok")])
    spiketot = rowSums(umi_data[mask1, c("spike_UMI_err", "spike_barcode_err", "spike_readpoor_offset_err",
        "spike_lonely_offset_err", "spike_ok")])

    plot_cumulative(x = umi_data[mask1, "gene_UMI_err"]/genetot, y = umi_data[mask1,
        "spike_UMI_err"]/spiketot, main = "UMI Errors", xlab = "Raw UMI fraction",
        ylab = "Fraction of cells", col = gene_spike_col)
    legend("bottomright", legend = c("Genes", "Spikes"), lty = 1, col = c(colors()[132],
        "firebrick"), cex = 0.6, bg = "white")
    mtext("E", side = 3, adj = -0.6, line = 2)
    message("E")

    plot_cumulative(x = umi_data[mask1, "gene_barcode_err"]/genetot, y = umi_data[mask1,
        "spike_barcode_err"]/spiketot, main = "Barcode Errors", xlab = "Raw UMI fraction",
        ylab = "Fraction of cells", col = gene_spike_col)
    legend("bottomright", legend = c("Genes", "Spikes"), lty = 1, col = c(colors()[132],
        "firebrick"), cex = 0.6, bg = "white")

    mtext("F", side = 3, adj = -0.6, line = 2)
    message("F")

    gene_offset_errors = (umi_data[mask1, "gene_lonely_offset_err"] + umi_data[mask1,
        "gene_readpoor_offset_err"])/genetot
    spike_offset_errors = (umi_data[mask1, "spike_lonely_offset_err"] + umi_data[mask1,
        "spike_readpoor_offset_err"])/spiketot

    if (any(is.na(gene_offset_errors)) | any(is.na(spike_offset_errors))) {
        plot.new()
    } else {
        plot_cumulative(x = gene_offset_errors, y = spike_offset_errors, main = "Switchers",
            xlab = "Raw UMI fraction", ylab = "Fraction of cells", col = gene_spike_col)
        legend("bottomright", legend = c("Genes", "Spikes"), lty = 1, col = c(colors()[132],
            "firebrick"), cex = 0.6, bg = "white")
    }
    mtext("G", side = 3, adj = -0.6, line = 2)
    message("G")

    spike_amount = batch_details$Spike_dilution * batch_details$Spike_volume_ul
    spike_concentration_mat = spike_concentration[, batch_details$Spike_type]
    const = 10^(-18) * (6.02214129 * (10^23))
    spike_molecules_per_ul_mat <- spike_concentration_mat * const

    expected = spike_molecules_per_ul_mat * spike_amount * length(batch_mask)
    observed = rowSums(umitab_batch[names(expected), batch_mask], na.rm = T)
    main_stats[["spike_yield"]] = round(100 * sum(observed, na.rm = T)/sum(expected,
        na.rm = T), digits = 1)

    if (length(mask0) > 0) {
        neg_control_spike_umis = colSums(umitab_batch[spike_names, as.character(mask0),
            drop = F], na.rm = T)
        neg_control_gene_umis = colSums(umitab_batch[setdiff(rownames(umitab_batch),
            spike_names), as.character(mask0), drop = F])
        main_stats$gene_umis_neg_control = round(mean(neg_control_gene_umis))
        main_stats$noise_estimation = round(mean(neg_control_gene_umis)/main_stats$gene_umis,
            digits = 3)
        barplot(rbind(neg_control_gene_umis, neg_control_spike_umis), col = c(colors()[132],
            "firebrick"), border = F, names.arg = sample_list[mask0, "well_coordinates"],
            main = "#UMI (neg control)", ylab = "#UMIs")
        legend("bottomright", legend = c("Genes", "Spikes"), pch = 15, col = c(colors()[132],
            "firebrick"), cex = 0.6, bg = "white")
        mtext(side = 3, text = paste("Est. noise = ", round(100 * main_stats$noise_estimation,
            digits = 1), "%", sep = ""), col = 2, cex = 0.7)

    } else {
        plot.new()
    }
    mtext("H", side = 3, adj = -0.6, line = 2)
    message("H")

    nmols_genes = colSums(umitab_batch[gene_mask, ])
    nmols_spike = colSums(umitab_batch[spike_names, ], na.rm = T)

    plot(nmols_spike, log10(nmols_genes + 1), xlab = "#spike_umis", ylab = "log10(#gene_umis)",
        panel.first = grid(), col = ifelse(colnames(umitab_batch) %in% mask1, 1,
            2), pch = ifelse(colnames(umitab_batch) %in% mask1, 20, 4), cex = ifelse(colnames(umitab_batch) %in%
            mask1, 0.5, 2))
    # text(z[,1],z[,2],labels=z$well_coor,col=ifelse(!z$is_empty,colors()[132],'firebrick'),cex=.6)
    mtext("I", side = 3, adj = -0.6, line = 2)
    message("I")

    single_cell_mask = colnames(umitab_batch) %in% mask1
    plot(colSums(umitab_batch[gene_mask, ]), 100 * colSums(singleoffsetab_batch[gene_mask,
        ])/colSums(umitab_batch[gene_mask, ]), log = "x", pch = ifelse(single_cell_mask,
        20, 4), col = ifelse(single_cell_mask, 1, 2), ylim = c(0, 100), xlab = "#UMIs",
        ylab = "%Singleton IVT prod.", panel.first = grid(lty = 1), cex = ifelse(single_cell_mask,
            0.5, 2))
    mtext("J", side = 3, adj = -0.6, line = 2)
    message("J")

    plate_id = plot_plate(umitab_batch[gene_mask, ], breaks = seq(-3, 3, l = 99))
    title(paste("UMI plate", plate_id))

    mtext("K", side = 3, adj = -0.6, line = 2)
    message("K")

    plate_id = plot_plate(umitab_batch[spike_names, ], breaks = seq(-3, 3, l = 99))
    mtext("L", side = 3, adj = -0.6, line = 2)
    title(paste("spike UMI plate", plate_id))
    message("L")

    nmols = colSums(umitab_batch[gene_mask, mask1])
    main_stats$gene_umis = round(mean(nmols))
    med = plot_cumulative(nmols, main = "#net gene UMIs", ylab = "Fraction of cells")
    main_stats$singleton_gene_umis = round(mean(colSums(singleoffsetab_batch[gene_mask,
        mask1])))
    mtext(side = 3, text = paste("Median =", med, "Mean = ", main_stats$gene_umis,
        sep = " "), col = 2, cex = 0.7)
    mtext("M", side = 3, adj = -0.6, line = 2)
    message("M")

    nmols = colSums(umitab_batch[spike_names, mask1], na.rm = T)
    med = plot_cumulative(nmols, main = "#net spike UMIs", ylab = "Fraction of wells")
    mtext(side = 3, text = paste("Median=", med, paste("Yield=", main_stats[["spike_yield"]],
        "%"), sep = " "), col = 2, cex = 0.7)
    mtext("N", side = 3, adj = -0.6, line = 2)
    message("N")

    main_stats$spike_umis = round(mean(nmols))

    par(mar = c(4, 5, 4, 2))

    noffsets_per_umi = read.table(paste(output_dir, "/QC/noffsets_per_umi_distrib/",
        amp_batch_ID, ".txt", sep = ""), header = T, row.names = 1)
    barplot2(log10(1 + t(noffsets_per_umi[2:20, 3])), col = "gray", main = "#IVT products/UMI",
        names = 1:19, plot.grid = T, border = F, ylab = "log10(counts)", space = 0)
    avg = round(sum(noffsets_per_umi[, 3] * 0:20)/sum(noffsets_per_umi[, 3]), digits = 2)
    mtext(side = 3, text = paste("Average =", avg), col = 2, cex = 0.7)
    mtext("O", side = 3, adj = -0.6, line = 2)
    message("O")
    main_stats$avg_noffsets_per_umi = avg


    nreads_per_umi = read.table(paste(output_dir, "/QC/nreads_per_umi_distrib/",
        amp_batch_ID, ".txt", sep = ""), header = T, row.names = 1)
    barplot2(log10(1 + t(nreads_per_umi[2:51, 3])), col = "gray", main = "#Reads/UMI",
        names = 1:50, plot.grid = T, border = F, ylab = "log10(counts)", space = 0)
    avg = round(sum(nreads_per_umi[, 3] * 0:100)/sum(nreads_per_umi[, 3]), digits = 2)
    mtext(side = 3, text = paste("Average =", avg), col = 2, cex = 0.7)
    mtext("P", side = 3, adj = -0.6, line = 2)
    message("P")
    main_stats$avg_reads_per_umi = avg

    par(mar = c(7, 5, 1, 2))
    mv = rowMeans(umitab_batch[setdiff(rownames(umitab_batch), spike_names), ])
    vv = apply(umitab_batch[setdiff(rownames(umitab_batch), spike_names), ], 1, var)
    gene_names = substr(names(mv), 1, 30)
    names(mv) = gene_names
    names(vv) = gene_names
    barplot(log2(sort(mv, decreasing = T)[1:25]), las = 2, ylab = "log2(#UMI/cell)",
        cex.names = 0.8, border = F, space = 0)
    mtext("Q", side = 3, adj = -0.2, line = 2)
    message("Q")

    barplot(log2(sort((vv/mv)[mv > 0.1], decreasing = T)[1:25]), las = 2, ylab = "log2(var/mean)",
        cex.names = 0.8, border = F, space = 0)
    mtext("R", side = 3, adj = -0.2, line = 2)
    message("R")

    return(unlist(main_stats))
}

plot_plate = function(umitab, breaks = seq(-4, 4, l = 99)) {
    wells = colnames(umitab)
    nrow = 16
    ncol = 24
    letters = LETTERS[1:nrow]
    nums = 1:ncol
    match1 = match(sample_list[wells, "well_coordinates"], as.vector(paste(matrix(letters,
        nrow, ncol), matrix(nums, nrow, ncol, byrow = T), sep = "")))
    plate_id = unique(sample_list[wells, "plate_ID"])
    v = rep(NA, nrow * ncol)
    v[match1] = colSums(umitab)
    m = matrix(v, nrow, ncol)
    m2 = log2((1 + m)/(1 + mean(m, na.rm = T)))

    rownames(m2) = c(letters)
    colnames(m2) = c(nums)
    m2 = m2[dim(m2)[1]:1, ]

    image(t(m2), col = colorRampPalette(c("blue", "white", "red"))(100), axes = F,
        breaks = c(-100, breaks, 100))
    mtext(3, text = colnames(m2), at = seq(0, 1, l = ncol), line = 0, cex = 0.2)
    mtext(2, text = rownames(m2), at = seq(0, 1, l = nrow), las = 1, line = 0.2,
        cex = 0.2)
    grid(lty = 1, col = "gray", nx = ncol, ny = nrow)
    box()

    return(plate_id)
}

old_clustering_MAP = function(p_table_fn = "submission1_likelihood_seed_p.txt", umitab) {
    ptab = read.table(p_table_fn, sep = "\t", check.names = F)[, -1] + 1e-04
    ptab = t(t(ptab)/colSums(ptab))
    get_ll = function(v, ptab) {
        return(colSums(v * log10(ptab), na.rm = T))
    }
    umitab2 = umitab[rownames(ptab), ]
    llv = apply(umitab2 + 0.1, 2, get_ll, ptab)
    llm = matrix(llv, , dim(ptab)[2])
    colnames(llm) = colnames(ptab)
    rownames(llm) = colnames(umitab)
    print(table(apply(llm, 1, which.max)))
}

get_stats_per_seq_batch = function(seq_batch) {
    n_umis_per_cell = colSums(umitab)
    cell_mask = intersect(rownames(sample_list)[sample_list$sequencing_batch %in%
        seq_batch], single_cell_mask)
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
# args=c('AB339','output')
if (length(args) == 6) {
    amp_batch = args[1]
    wells_cells_txt = args[2]
    amp_batches_txt = args[3]
    seq_batches_txt = args[4]
    spike_concentrations_txt = args[5]
    output_dir = args[6] # output where umi.tab and others are stored
} else {
    message("Usage: Rscript qc_batch.r [amp_batch] [well_cells.txt] [amp_batches.txt] [seq_batches.txt] [spike_concentrations.txt] [output_dir]")
    stop()
}

gene_spike_col = c(colors()[132], "firebrick")

sample_list = read.table(wells_cells_txt, stringsAsFactors = F, header = T, sep = "\t")
amp_batches = read.delim(amp_batches_txt, stringsAsFactors = F, header = T, sep = "\t")
seq_batches = read.delim(seq_batches_txt, stringsAsFactors = F, header = T, sep = "\t")

rownames(sample_list) = sample_list$Well_ID

pdf(paste("report_per_amp_batch/", amp_batch, ".pdf", sep = ""), 8.3,
    11.7, paper = "a4")

message(amp_batch)

try({
    main_stats = make_amp_batch_qc_fig(amp_batch, spike_concentrations_txt, output_dir)
    save(main_stats, file = paste("rd/", amp_batch, ".rd", sep = ""))
})
garbage <- dev.off()
