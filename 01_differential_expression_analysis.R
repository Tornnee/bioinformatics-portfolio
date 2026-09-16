# =============================================================================
# SCRIPT 01: Differential Gene Expression Analysis
# Author:    gh (Postgraduate Researcher, Bioinformatics & Molecular Biology)
# Dataset:   GEO GSE120103 — Stage IV Ovarian Endometriosis vs Disease-Free Endometrium
# Platform:  Agilent 014850 Whole Human Genome Microarray 4x44K (GPL6480)
# Tools:     GEOquery, limma, ggplot2, clusterProfiler, enrichplot
# Date:      2026
# =============================================================================
# DESCRIPTION:
#   This script performs a full differential expression analysis pipeline:
#   1. Downloads and loads GEO microarray data
#   2. Quality control and quantile normalisation
#   3. Identifies differentially expressed genes (DEGs) using limma
#   4. Generates a volcano plot
#   5. Performs GO and KEGG pathway enrichment using clusterProfiler
#   6. Exports results tables
# =============================================================================

# ── 0. INSTALL & LOAD PACKAGES ─────────────────────────────────────────────────

# Install Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("GEOquery", "limma", "org.Hs.eg.db",
               "clusterProfiler", "enrichplot", "AnnotationDbi")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

cran_pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "pheatmap", "RColorBrewer")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# ── 1. LOAD GEO DATASET ────────────────────────────────────────────────────────

cat("Downloading GEO dataset GSE120103...\n")

gse <- getGEO("GSE120103", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]  # Extract the ExpressionSet

# Inspect the data
cat("Dataset dimensions:", dim(exprs(gse)), "\n")
cat("Sample metadata columns:", colnames(pData(gse)), "\n")

# ── 2. DEFINE SAMPLE GROUPS ────────────────────────────────────────────────────

# Extract the phenotype/sample metadata
pheno <- pData(gse)

# Identify control vs endometriosis samples
# Adjust the column name if needed — check pheno column names first
cat("\nUnique values in 'characteristics_ch1':\n")
print(unique(pheno$characteristics_ch1))

# Create a group factor based on disease status
# Modify the string match to match your dataset exactly
group <- ifelse(grepl("disease-free|control|normal", pheno$characteristics_ch1,
                       ignore.case = TRUE), "Control", "Endometriosis")
group <- factor(group, levels = c("Control", "Endometriosis"))
cat("\nSample group assignments:\n")
print(table(group))

# ── 3. SELECT QUALITY-CONTROLLED SAMPLES ───────────────────────────────────────

# Based on QC (box plot review), retain 4 control + 4 endometriosis samples
# These are the 8 accessions confirmed in the analysis
control_ids      <- c("GSM3393509", "GSM3393510", "GSM3393511", "GSM3393512")
endometriosis_ids <- c("GSM3393518", "GSM3393519", "GSM3393520", "GSM3393521")
keep_ids         <- c(control_ids, endometriosis_ids)

# Subset ExpressionSet to selected samples
gse_filtered <- gse[, sampleNames(gse) %in% keep_ids]
group_filtered <- group[sampleNames(gse) %in% keep_ids]
cat("\nSamples retained after QC:", ncol(gse_filtered), "\n")

# ── 4. EXTRACT AND NORMALISE EXPRESSION DATA ───────────────────────────────────

expr_raw <- exprs(gse_filtered)

# Log2 transform if not already done (check data range first)
if (max(expr_raw, na.rm = TRUE) > 50) {
  cat("Log2 transforming expression data...\n")
  expr_log2 <- log2(expr_raw + 1)
} else {
  cat("Data appears already log2-transformed.\n")
  expr_log2 <- expr_raw
}

# Quantile normalisation using limma
expr_norm <- normalizeBetweenArrays(expr_log2, method = "quantile")
cat("Quantile normalisation complete.\n")

# ── 5. QUALITY CONTROL PLOTS ───────────────────────────────────────────────────

# 5a. Box plot of normalised expression
pdf("QC_boxplot_normalised.pdf", width = 10, height = 6)
boxplot(expr_norm,
        col      = c(rep("#2C6B2F", 4), rep("#8B3A00", 4)),
        las      = 2,
        cex.axis = 0.7,
        main     = "Quantile-Normalised Expression — GSE120103",
        ylab     = "log2 Expression",
        names    = colnames(expr_norm))
legend("topright",
       legend = c("Control (n=4)", "Stage IV Endometriosis (n=4)"),
       fill   = c("#2C6B2F", "#8B3A00"),
       bty    = "n")
dev.off()
cat("Box plot saved: QC_boxplot_normalised.pdf\n")

# 5b. PCA plot
pca <- prcomp(t(expr_norm), scale. = TRUE)
pca_df <- data.frame(
  PC1   = pca$x[, 1],
  PC2   = pca$x[, 2],
  Group = group_filtered,
  Label = colnames(expr_norm)
)
var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

pdf("QC_PCA_plot.pdf", width = 7, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, label = Label)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c("Control" = "#2C6B2F",
                                 "Endometriosis" = "#8B3A00")) +
  labs(
    title    = "PCA of Normalised Microarray Expression",
    subtitle = "GSE120103 — Stage IV Ovarian Endometriosis",
    x        = paste0("PC1 (", var_explained[1], "% variance)"),
    y        = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")
dev.off()
cat("PCA plot saved: QC_PCA_plot.pdf\n")

# ── 6. DIFFERENTIAL EXPRESSION ANALYSIS (limma) ────────────────────────────────

# Build design matrix (no intercept, group-coded)
design <- model.matrix(~0 + group_filtered)
colnames(design) <- levels(group_filtered)

# Fit linear model
fit <- lmFit(expr_norm, design)

# Define contrast: Endometriosis vs Control
contrast_matrix <- makeContrasts(
  Endometriosis_vs_Control = Endometriosis - Control,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract all results with gene annotation
results_all <- topTable(fit2,
                        coef       = "Endometriosis_vs_Control",
                        number     = Inf,
                        adjust.method = "BH",
                        sort.by    = "P")

# Add gene symbols from feature data
fdata <- fData(gse_filtered)
gene_col <- intersect(c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment"),
                      colnames(fdata))[1]
results_all$GeneSymbol <- fdata[rownames(results_all), gene_col]

cat("\nTop 10 DEGs:\n")
print(head(results_all[, c("GeneSymbol", "logFC", "AveExpr",
                             "t", "P.Value", "adj.P.Val")], 10))

# ── 7. FILTER DEGs ─────────────────────────────────────────────────────────────

# Apply significance thresholds: adj.P.Val < 0.05 and |logFC| > 1
# NOTE: |logFC| > 0 was used in the original analysis but > 1 is more standard
degs <- results_all %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  mutate(
    Direction = ifelse(logFC > 0, "Upregulated", "Downregulated")
  )

cat("\nTotal DEGs (adj.P.Val < 0.05, |logFC| > 1):", nrow(degs), "\n")
cat("  Upregulated:  ", sum(degs$Direction == "Upregulated"), "\n")
cat("  Downregulated:", sum(degs$Direction == "Downregulated"), "\n")

# ── 8. VOLCANO PLOT ────────────────────────────────────────────────────────────

# Prepare plot data
plot_data <- results_all %>%
  mutate(
    neg_log10_p = -log10(adj.P.Val),
    Significance = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE                           ~ "Not Significant"
    ),
    Label = ifelse(Significance != "Not Significant" &
                     abs(logFC) > 3 & neg_log10_p > 3,
                   GeneSymbol, NA)
  )

pdf("DEG_volcano_plot.pdf", width = 8, height = 7)
ggplot(plot_data, aes(x = logFC, y = neg_log10_p, colour = Significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_text_repel(aes(label = Label),
                  size = 3, max.overlaps = 30,
                  box.padding = 0.4, segment.colour = "grey50") +
  scale_colour_manual(values = c(
    "Upregulated"     = "#B22222",
    "Downregulated"   = "#1A6B8A",
    "Not Significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  labs(
    title    = "Volcano Plot — Stage IV Ovarian Endometriosis vs Disease-Free Endometrium",
    subtitle = paste0("GSE120103 | ", nrow(degs), " DEGs (adj.P < 0.05, |logFC| > 1)"),
    x        = "log2 Fold Change (Endometriosis / Control)",
    y        = "-log10 Adjusted P-Value",
    colour   = "Direction"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")
dev.off()
cat("Volcano plot saved: DEG_volcano_plot.pdf\n")

# ── 9. PATHWAY ENRICHMENT ANALYSIS ─────────────────────────────────────────────

# Convert gene symbols to Entrez IDs
deg_symbols <- na.omit(unique(degs$GeneSymbol))

entrez_ids <- bitr(deg_symbols,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

cat("\nGenes mapped to Entrez IDs:", nrow(entrez_ids), "of", length(deg_symbols), "\n")

# 9a. GO Biological Process enrichment
go_bp <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("\nSignificant GO-BP terms:", nrow(as.data.frame(go_bp)), "\n")

if (nrow(as.data.frame(go_bp)) > 0) {
  pdf("Enrichment_GO_BP_dotplot.pdf", width = 10, height = 7)
  print(dotplot(go_bp,
                showCategory = 15,
                title        = "GO Biological Process Enrichment (adj.P < 0.05)",
                font.size    = 10))
  dev.off()
  cat("GO-BP dot plot saved: Enrichment_GO_BP_dotplot.pdf\n")
}

# 9b. KEGG Pathway enrichment
kegg <- enrichKEGG(
  gene          = entrez_ids$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

cat("Significant KEGG pathways:", nrow(as.data.frame(kegg)), "\n")

if (nrow(as.data.frame(kegg)) > 0) {
  pdf("Enrichment_KEGG_dotplot.pdf", width = 10, height = 6)
  print(dotplot(kegg,
                showCategory = 10,
                title        = "KEGG Pathway Enrichment (adj.P < 0.05)",
                font.size    = 10))
  dev.off()
  cat("KEGG dot plot saved: Enrichment_KEGG_dotplot.pdf\n")
}

# ── 10. EXPORT RESULTS ─────────────────────────────────────────────────────────

# All DEGs table
write.csv(degs, "DEG_results_filtered.csv", row.names = TRUE)
cat("DEG results saved: DEG_results_filtered.csv\n")

# All results (unfiltered)
write.csv(results_all, "DEG_results_all.csv", row.names = TRUE)
cat("Full results saved: DEG_results_all.csv\n")

# GO-BP results
if (exists("go_bp") && nrow(as.data.frame(go_bp)) > 0) {
  write.csv(as.data.frame(go_bp), "GO_BP_enrichment_results.csv", row.names = FALSE)
  cat("GO-BP results saved: GO_BP_enrichment_results.csv\n")
}

# KEGG results
if (exists("kegg") && nrow(as.data.frame(kegg)) > 0) {
  write.csv(as.data.frame(kegg), "KEGG_enrichment_results.csv", row.names = FALSE)
  cat("KEGG results saved: KEGG_enrichment_results.csv\n")
}

cat("\n=== Analysis complete. All output files saved. ===\n")
