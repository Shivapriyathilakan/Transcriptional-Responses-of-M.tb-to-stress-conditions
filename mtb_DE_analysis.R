# =============================================================================
# Mtb pathogen-side transcriptomics: per-condition differential expression 
#
# Three independent studies (hypoxia, rifampicin, iron starvation)
#3 replicates per group each = 18 samples total. condition is fully confounded with study,
# so every DE comparison is done WITHIN a single study/condition, never pooled.
#
# Covers: per-condition DE, volcano plots, cross-condition overlap with
# direction concordance, and KEGG pathway enrichment per condition.
# =============================================================================


cran_pkgs <- c("dplyr", "ggplot2", "ggVennDiagram", "pheatmap", "rlang")
bioc_pkgs <- c("DESeq2", "apeglm", "clusterProfiler", "enrichplot")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

# Load featureCounts matrix
counts_raw <- read.delim("mtb_gene_counts.txt", comment.char = "#", stringsAsFactors = FALSE)

gene_ids   <- counts_raw$Geneid
counts_mat <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(counts_mat) <- gene_ids

colnames(counts_mat) <- sub(".*?(SRR\\d+).*", "\\1", colnames(counts_mat))

# Load metadata 
sample_metadata <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)

stopifnot(
  "Some sample_metadata$sample_id are missing from counts_mat" =
    all(sample_metadata$sample_id %in% colnames(counts_mat))
)

counts_mat <- counts_mat[, sample_metadata$sample_id]

stopifnot(
  "Column order mismatch between counts_mat and sample_metadata after reorder" =
    identical(colnames(counts_mat), sample_metadata$sample_id)
)

# Combined QC object across all 18 samples for PCA 
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sample_metadata,
  design    = ~condition
)
dds_all <- dds_all[rowSums(counts(dds_all)) > 10, ]
vsd_all <- vst(dds_all, blind = TRUE)

pca_by_condition <- plotPCA(vsd_all, intgroup = "condition") +
  ggtitle("QC PCA - colored by condition")
pca_by_study <- plotPCA(vsd_all, intgroup = "study_id") +
  ggtitle("QC PCA - colored by study_id (batch)")

ggsave("QC_PCA_by_condition.png", pca_by_condition, width = 6, height = 5)
ggsave("QC_PCA_by_study.png", pca_by_study, width = 6, height = 5)
# since condition == study_id here, the two plot will look identical


condition_config <- list(
  hypoxia = list(
    condition   = "hypoxia",
    treat_level = "treatment",
    ctrl_level  = "control"
  ),
  rifampicin = list(
    condition   = "rifampicin",
    treat_level = "treatment",  
    ctrl_level  = "control"     
  ),
  iron_starvation = list(
    condition   = "iron_starvation",
    treat_level = "treatment",  
    ctrl_level  = "control"     
  )
)

# STEP 6: Per-condition DE
run_de <- function(meta, counts, cfg) {
  sub_meta <- meta[meta$condition == cfg$condition, ]
  sub_counts <- counts[, sub_meta$sample_id]
  
  sub_meta$group <- factor(sub_meta$group)
  stopifnot(
    "treat_level/ctrl_level not found in this condition's group factor" =
      all(c(cfg$treat_level, cfg$ctrl_level) %in% levels(sub_meta$group))
  )
  sub_meta$group <- relevel(sub_meta$group, ref = cfg$ctrl_level)
  
  dds <- DESeqDataSetFromMatrix(sub_counts, sub_meta, design = ~group)
  
  min_group_n <- min(table(sub_meta$group))
  keep <- rowSums(counts(dds) >= 10) >= min_group_n
  dds <- dds[keep, ]
  
  DESeq(dds)
}

get_shrunk_results <- function(dds) {
  coef_name <- resultsNames(dds)[grep("group", resultsNames(dds))]
  lfcShrink(dds, coef = coef_name, type = "apeglm")
}

dds_by_condition <- lapply(condition_config, run_de, meta = sample_metadata, counts = counts_mat)
res_by_condition <- lapply(dds_by_condition, get_shrunk_results)

# Save full and significant-only results per condition
sig_by_condition <- list()

for (cond in names(res_by_condition)) {
  res_df <- as.data.frame(res_by_condition[[cond]])
  write.csv(res_df, sprintf("DE_%s_full.csv", cond))
  
  sig_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
  sig_df <- sig_df[order(-abs(sig_df$log2FoldChange)), ]
  write.csv(sig_df, sprintf("DE_%s_significant.csv", cond), row.names = FALSE)
  sig_by_condition[[cond]] <- sig_df
  
  message(sprintf("%s: %d significant DEGs (padj<0.05 & |log2FC|>1) out of %d genes tested",
                  cond, nrow(sig_df), nrow(res_df)))
}

# Volcano plot per condition
plot_volcano <- function(res_df, title, outfile) {
  df <- as.data.frame(res_df)
  df$neglog10padj <- -log10(df$padj)
  df$sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1
  
  p <- ggplot(df, aes(log2FoldChange, neglog10padj, color = sig)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = title, x = "log2 fold change", y = "-log10(padj)", color = "Significant") +
    theme_minimal()
  
  ggsave(outfile, p, width = 6, height = 5)
}

for (cond in names(res_by_condition)) {
  plot_volcano(
    res_by_condition[[cond]],
    title   = sprintf("DE volcano - %s", cond),
    outfile = sprintf("volcano_%s.png", cond)
  )
}

# Cross-condition overlap - simple significance-only intersection
sig_gene_sets <- lapply(sig_by_condition, rownames)

shared_all3 <- Reduce(intersect, sig_gene_sets)
unique_by_condition <- lapply(names(sig_gene_sets), function(cond) {
  others <- setdiff(names(sig_gene_sets), cond)
  setdiff(sig_gene_sets[[cond]], Reduce(union, sig_gene_sets[others]))
})
names(unique_by_condition) <- names(sig_gene_sets)

message(sprintf("Genes significant in ALL THREE conditions: %d", length(shared_all3)))
for (cond in names(unique_by_condition)) {
  message(sprintf("Genes significant ONLY in %s: %d", cond, length(unique_by_condition[[cond]])))
}

writeLines(shared_all3, "DEGs_shared_all3_conditions.txt")
for (cond in names(unique_by_condition)) {
  writeLines(unique_by_condition[[cond]], sprintf("DEGs_unique_to_%s.txt", cond))
}

# Cross-condition Venn of significant DEGs 
venn_plot <- ggVennDiagram(sig_gene_sets, label = "count") +
  ggtitle("Significant DEGs by condition (direction not shown)")
ggsave("DEG_venn_overlap.png", venn_plot, width = 6, height = 6, bg = "white")

# Direction concordance among the 3-way shared genes
direction_table <- sapply(names(res_by_condition), function(cond) {
  df <- as.data.frame(res_by_condition[[cond]])
  sign(df[shared_all3, "log2FoldChange"])
})
rownames(direction_table) <- shared_all3

concordant_genes <- shared_all3[apply(direction_table, 1, function(x) length(unique(x)) == 1)]
discordant_genes <- setdiff(shared_all3, concordant_genes)

write.csv(
  data.frame(gene = shared_all3, direction_table, check.names = FALSE),
  "DEGs_shared_all3_direction_table.csv", row.names = FALSE
)
message(sprintf(
  "Of %d genes significant in all 3 conditions: %d concordant (same direction everywhere), %d discordant",
  length(shared_all3), length(concordant_genes), length(discordant_genes)
))

# Heatmap of the shared, concordant DEGs across all samples------
vsd_by_condition <- lapply(dds_by_condition, function(dds) vst(dds, blind = FALSE))
expr_by_condition <- lapply(vsd_by_condition, assay)

common_genes <- concordant_genes
combined_expr <- do.call(cbind, lapply(expr_by_condition, function(m) m[common_genes, , drop = FALSE]))

z_combined <- t(scale(t(combined_expr)))

sample_annotation <- sample_metadata %>%
  filter(sample_id %in% colnames(z_combined)) %>%
  select(sample_id, condition, group) %>%
  tibble::column_to_rownames("sample_id")

pheatmap(
  z_combined,
  annotation_col = sample_annotation[, c("condition", "group")],
  show_rownames = nrow(z_combined) <= 50,  # avoid unreadable labels if the list is large
  filename = "heatmap_shared_concordant_DEGs.png",
  width = 8, height = 10
)

# KEGG pathway enrichment per condition
run_kegg <- function(sig_genes, background_genes, title, outfile) {
  ekegg <- enrichKEGG(
    gene         = sig_genes,
    universe     = background_genes,
    organism     = "mtu",
    keyType      = "kegg",
    pvalueCutoff = 0.05
  )
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    message(sprintf("%s: no significant KEGG pathways returned (check gene ID format if unexpected)", title))
    return(NULL)
  }
  write.csv(as.data.frame(ekegg), sprintf("KEGG_%s.csv", tolower(gsub(" ", "_", title))))
  p <- barplot(ekegg, showCategory = 15, title = sprintf("KEGG enrichment - %s", title))
  ggsave(outfile, p, width = 8, height = 6)
  ekegg
}

kegg_by_condition <- lapply(names(sig_by_condition), function(cond) {
  run_kegg(
    sig_genes       = rownames(sig_by_condition[[cond]]),
    background_genes = rownames(res_by_condition[[cond]]),
    title           = cond,
    outfile         = sprintf("KEGG_barplot_%s.png", cond)
  )
})
names(kegg_by_condition) <- names(sig_by_condition)