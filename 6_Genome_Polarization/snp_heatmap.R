#!/usr/bin/env Rscript

library(vcfR)
library(pheatmap)

# ---------------------------
# Arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript snp_heatmap.R input.vcf.gz output.pdf [sample_order.txt]")
}

vcf_file <- args[1]
output_pdf <- args[2]
order_file <- if (length(args) == 3) args[3] else NULL

# ---------------------------
# Read VCF (UNCHANGED)
# ---------------------------
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT")

fix <- getFIX(vcf)
chrom <- fix[, "CHROM"]
pos <- as.numeric(fix[, "POS"])

# ---------------------------
# Convert genotypes (0/1/2)
# ---------------------------
gt_num <- matrix(NA_real_, nrow = nrow(gt), ncol = ncol(gt))

gt_num[gt %in% c("0/0","0|0")] <- 0
gt_num[gt %in% c("0/1","1/0","0|1","1|0")] <- 1
gt_num[gt %in% c("1/1","1|1")] <- 2

gt_num <- t(gt_num)
rownames(gt_num) <- colnames(gt)

# ---------------------------
# Species assignment
# ---------------------------
species <- sub("[0-9]+$", "", rownames(gt_num))

# ---------------------------
# Ordering
# ---------------------------
if (!is.null(order_file)) {
  
  desired_order <- scan(order_file, what = character(), quiet = TRUE)
  
  missing_samples <- setdiff(desired_order, rownames(gt_num))
  
  if (length(missing_samples) > 0) {
    stop(paste("Samples not found:", paste(missing_samples, collapse = ", ")))
  }
  
  gt_num <- gt_num[desired_order, , drop = FALSE]
  species <- species[match(desired_order, rownames(gt_num))]
  
} else {
  
  species_levels <- c("Lgra", "Ltim", "Lcor", "Lcas")
  
  species <- factor(species, levels = species_levels)
  ord <- order(species)
  
  gt_num <- gt_num[ord, , drop = FALSE]
  species <- species[ord]
}

# ---------------------------
# Recode relative to Lcas majority
# ---------------------------
lcas_idx <- which(species == "Lcas")

if (length(lcas_idx) == 0) {
  stop("No Lcas samples found in dataset.")
}

for (j in seq_len(ncol(gt_num))) {
  
  vals <- gt_num[lcas_idx, j]
  vals <- vals[!is.na(vals)]
  
  if (length(vals) == 0) next
  
  major <- as.numeric(names(sort(table(vals), decreasing = TRUE)[1]))
  
  gt_num[, j] <- gt_num[, j] - major
  gt_num[, j] <- ifelse(gt_num[, j] < 0, 2, gt_num[, j])
}

# =========================================================
# NEW PART: 50 KB WINDOWING (ALL SNPs USED)
# =========================================================

window_size <- 50000

windowed_list <- list()

for (chr in unique(chrom)) {
  
  cat("Processing chromosome:", chr, "\n")
  
  chr_idx <- which(chrom == chr)
  
  chr_pos <- pos[chr_idx]
  chr_gt  <- gt_num[, chr_idx, drop = FALSE]
  
  window_id <- floor(chr_pos / window_size)
  
  windows <- sort(unique(window_id))
  
  window_mat <- matrix(
    NA_real_,
    nrow = nrow(chr_gt),
    ncol = length(windows)
  )
  
  rownames(window_mat) <- rownames(chr_gt)
  colnames(window_mat) <- paste0(chr, "_", windows)
  
  for (i in seq_along(windows)) {
    
    idx <- which(window_id == windows[i])
    
    window_mat[, i] <- rowMeans(
      chr_gt[, idx, drop = FALSE],
      na.rm = TRUE
    )
  }
  
  windowed_list[[chr]] <- window_mat
}

# ---------------------------
# Counts (still SNP-based, unchanged meaning)
# ---------------------------
count_table <- data.frame(
  Sample = rownames(gt_num),
  Species = as.character(species),
  Blue = rowSums(gt_num == 0, na.rm = TRUE),
  Yellow = rowSums(gt_num == 1, na.rm = TRUE),
  Red = rowSums(gt_num == 2, na.rm = TRUE),
  Missing = rowSums(is.na(gt_num)),
  stringsAsFactors = FALSE
)

counts_file <- sub("\\.pdf$", "_counts.tsv", output_pdf)
if (counts_file == output_pdf) {
  counts_file <- paste0(output_pdf, "_counts.tsv")
}

write.table(count_table,
            file = counts_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# ---------------------------
# Annotation
# ---------------------------
annotation_row <- data.frame(Species = species)
rownames(annotation_row) <- rownames(gt_num)

ann_colors <- list(
  Species = c(
    Lgra = "#c5cca0ff",
    Ltim = "#a3bedcff",
    Lcor = "#fca64dff",
    Lcas = "#f4a6b8ff"
  )
)

# ---------------------------
# PLOTTING (WINDOWED)
# ---------------------------

for (chr in names(windowed_list)) {
  
  mat <- windowed_list[[chr]]
  
  out_pdf <- sub("\\.pdf$", paste0("_", chr, ".pdf"), output_pdf)
  
  pdf(out_pdf,
      width = max(12, ncol(mat) / 50),
      height = max(6, nrow(mat) / 3),
      useDingbats = FALSE)
  
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    border_color = NA,
    fontsize_row = 8,
    color = colorRampPalette(c("blue", "purple", "red"))(100),
    breaks = seq(0, 2, length.out = 101),
    na_col = "grey90",
    main = chr
  )
  
  dev.off()
}

cat("Saved counts:", counts_file, "\n")
cat("Saved per-chromosome heatmaps.\n")
