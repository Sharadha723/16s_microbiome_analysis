# ==============================================================================
# PCOS 16S rRNA Microbiome Analysis Pipeline
# ==============================================================================
#
# Developed and substantially expanded by:
#   Sarada Giridharan
#
# Based on an original DADA2/phyloseq workflow provided by:
#   Donald Porchia and Zhigang Li
#
# Version:
#   v1.0.0
#
# Description:
#   An R-based workflow for processing paired-end 16S rRNA gene sequencing
#   data and performing downstream microbiome analyses.
#
# Major components include:
#   - Quality filtering and trimming
#   - DADA2 ASV inference
#   - Chimera removal
#   - Taxonomic assignment using SILVA
#   - Phylogenetic tree construction
#   - Phyloseq object generation
#   - Alpha and beta diversity analyses
#   - Differential abundance analyses
#   - Microbiome visualization and result export
#
# Attribution:
#   This pipeline substantially expands and modifies an original
#   DADA2/phyloseq workflow provided by Donald Porchia and Zhigang Li.
#
# License:
#   PolyForm Noncommercial License 1.0.0
#   See the LICENSE file for full terms.
#
# ==============================================================================


# ==============================================================================
# 0. SETUP — REQUIRED PACKAGES
# ==============================================================================

required_packages <- c(
  "dada2", "phyloseq", "Biostrings", "data.table", "dplyr",
  "ggplot2", "vegan", "DESeq2", "ANCOMBC", 
  "pheatmap", "RColorBrewer", "scales",
  "EnhancedVolcano", "ComplexHeatmap", "microbiomeMarker",
  "DECIPHER", "phangorn", "circlize"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install the missing package(s) before running the pipeline.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(vegan)
  library(DESeq2)
  library(ANCOMBC)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(microbiomeMarker)
  library(DECIPHER)
  library(phangorn)
})

set.seed(123)


# ==============================================================================
# 1. INPUT AND OUTPUT PATHS
# ==============================================================================

raw_data_dir <- file.path("data", "raw_data")
silva_train <- file.path("data", "reference","silva_nr99_v138.1_train_set.fa.gz")
silva_species <- file.path("data", "reference","silva_species_assignment_v138.1.fa.gz")
metadata_file <- file.path("data", "Metadata.csv")


# ==============================================================================
# OUTPUT DIRECTORIES
# ==============================================================================

results_dir <- "results"

alpha_dir <- file.path(results_dir, "alpha_diversity")
beta_dir  <- file.path(results_dir, "beta_diversity")
ord_dir   <- file.path(results_dir, "ordination")
rare_dir  <- file.path(results_dir, "rarefaction")
fig_dir   <- file.path(results_dir, "figures")
da_dir    <- file.path(results_dir, "differential_abundance")

for (d in c(results_dir,alpha_dir,beta_dir,ord_dir,rare_dir,fig_dir,da_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# 2. Locate and unzip FASTQ files
# ==============================================================================
zip_files <- list.files(raw_data_dir, pattern = "\\.zip$", full.names = TRUE)

for (zf in zip_files) {
  sample_base  <- tools::file_path_sans_ext(basename(zf))
  expected_fq  <- file.path(raw_data_dir, paste0(sample_base, "_1.fastq"))
  if (!file.exists(expected_fq)) unzip(zf, exdir = raw_data_dir)
}

fnFs <- sort(list.files(raw_data_dir, pattern = "_1\\.(fastq|fq)(\\.(gz|bz2))?$",
                        full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(raw_data_dir, pattern = "_2\\.(fastq|fq)(\\.(gz|bz2))?$",
                        full.names = TRUE, recursive = TRUE))

stopifnot(length(fnFs) == length(fnRs), length(fnFs) > 0)
sample.names <- sub("_1\\.(fastq|fq)(\\.(gz|bz2))?$", "", basename(fnFs))
stopifnot(identical(sample.names,
                    sub("_2\\.(fastq|fq)(\\.(gz|bz2))?$", "", basename(fnRs))))

# ==============================================================================
# 3. Quality profiles
# ==============================================================================
plotQualityProfile(fnFs, aggregate = TRUE)
plotQualityProfile(fnRs, aggregate = TRUE)

# ==============================================================================
# 4. Filter and trim
# ==============================================================================
filt_path <- file.path(raw_data_dir, "filtered")
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)
filtFs <- file.path(filt_path, paste0(sample.names, "_R1.filt.fastq.bz2"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R2.filt.fastq.bz2"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     trimLeft  = c(5, 5),
                     truncLen  = c(240, 200), #decided after looking at quality profiles
                     truncQ    = 2,
                     maxN      = 0,
                     maxEE     = c(2, 2),
                     rm.phix   = TRUE,
                     multithread = TRUE,
                     compress  = TRUE)

print(head(out, 10))

plotQualityProfile(filtFs, aggregate = TRUE)
plotQualityProfile(filtRs, aggregate = TRUE)

# ==============================================================================
# 5. Error models, dereplication, ASV inference
# ==============================================================================
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF)

derepFs <- derepFastq(filtFs); names(derepFs) <- sample.names
derepRs <- derepFastq(filtRs); names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# ==============================================================================
# 6. Merge, build sequence table, remove chimeras
# ==============================================================================
mergers    <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                         minOverlap = 12, verbose = TRUE)
seqtab     <- makeSequenceTable(mergers)

hist(nchar(colnames(seqtab)), breaks = 50,
     main = "ASV length distribution", xlab = "Length (bp)")

seq_lengths <- nchar(colnames(seqtab))
seqtab      <- seqtab[, seq_lengths >= 240 & seq_lengths <= 245, drop = FALSE]

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE)

cat(sprintf("Non-chimeric reads retained: %.1f%%\n",100 * sum(seqtab.nochim) / sum(seqtab)))

# ==============================================================================
# 7. Taxonomy assignment
# ==============================================================================
tax <- assignTaxonomy(seqtab.nochim, silva_train, minBoot = 80, tryRC = TRUE, multithread = TRUE)
tax <- addSpecies(tax, silva_species)

genus_counts_raw <- sort(table(tax[, "Genus"]), decreasing = TRUE)
barplot(genus_counts_raw[1:10], las = 2, main = "Top 10 Genera", ylab = "Number of ASVs")

# ==============================================================================
# 8. Build phyloseq object with phylogenetic tree
# ==============================================================================
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_len(ncol(seqtab.nochim)))
asv_map  <- data.frame(ASV = asv_ids, Sequence = asv_seqs, stringsAsFactors = FALSE)

colnames(seqtab.nochim) <- asv_ids
rownames(tax)            <- asv_ids

dna       <- DNAStringSet(asv_seqs); names(dna) <- asv_ids
alignment <- AlignSeqs(dna, anchor = NA)

phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm          <- phangorn::dist.ml(phang.align)
treeNJ      <- phangorn::NJ(dm)
fit         <- phangorn::pml(treeNJ, data = phang.align)
fit         <- phangorn::optim.pml(fit, model = "GTR",
                                   optInv = TRUE, optGamma = TRUE,
                                   rearrangement = "stochastic",
                                   control = phangorn::pml.control(trace = 0))

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax)),
  phy_tree(fit$tree),
  dna
)

# Attach metadata
metadata_df  <- fread(metadata_file)
metadata_df2 <- as.data.frame(metadata_df)
metadata_df2 <- metadata_df2[match(sample_names(ps), metadata_df2$run), ]
rownames(metadata_df2) <- metadata_df2$run
stopifnot(all(rownames(metadata_df2) == sample_names(ps)))
sample_data(ps) <- sample_data(metadata_df2)

# Fill missing genus labels
tt <- tax_table(ps)

if ("Genus" %in% colnames(tt)) {
  tt[, "Genus"][is.na(tt[, "Genus"]) | tt[, "Genus"] == ""] <- "Unassigned"
  tax_table(ps) <- tt
} else {
  txm <- cbind(as.matrix(tt), Genus = "Unassigned")
  tax_table(ps) <- tax_table(txm)
}

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

genus_counts <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) genus_counts <- t(genus_counts)
genus_rel    <- sweep(genus_counts, 1, pmax(rowSums(genus_counts), 1), "/")
genus_names  <- as.vector(tax_table(ps_genus)[, "Genus"])
colnames(genus_counts) <- genus_names
colnames(genus_rel)    <- genus_names

# ==============================================================================
# 9. ALPHA DIVERSITY
# ==============================================================================
# --  Estimate richness metrics  -----------------------------------------------
alpha_df        <- estimate_richness(ps, measures = c("Observed", "Shannon","Simpson", "Chao1"))
alpha_df$group  <- sample_data(ps)$group
alpha_df$Sample <- rownames(alpha_df)

# --  Wilcoxon tests for each metric  -----------------------------------------
alpha_tests <- lapply(c("Shannon", "Observed", "Simpson", "Chao1"), function(m) {
  wt <- wilcox.test(as.formula(paste(m, "~ group")), data = alpha_df)
  data.frame(Metric = m, W = wt$statistic, p.value = wt$p.value)
})

alpha_test_df <- do.call(rbind, alpha_tests)
alpha_test_df$p.adjusted <- p.adjust(alpha_test_df$p.value, method = "BH")
print(alpha_test_df)
fwrite(alpha_test_df, file = file.path(alpha_dir, "alpha_wilcoxon_tests.csv"))

# --  Plot all four metrics  ---------------------------------------------------
p_alpha <- plot_richness(ps,
                         x       = "group",
                         measures = c("Shannon", "Observed", "Simpson", "Chao1"),
                         color   = "group") +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(title = "Alpha Diversity by Group") +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "alpha_diversity.png"),
       p_alpha, width = 10, height = 6)

fwrite(alpha_df, file = file.path(alpha_dir, "alpha_diversity.csv"))

# ==============================================================================
# 10. BETA DIVERSITY — distances, ordination, statistical tests
# ==============================================================================

# --  Distance matrices  -------------------------------------------------------
bray_dist     <- phyloseq::distance(ps, method = "bray")
wunifrac_dist <- phyloseq::distance(ps, method = "wunifrac")

write.csv(as.matrix(bray_dist),
          file = file.path(beta_dir, "bray_curtis_distance_matrix.csv"))
write.csv(as.matrix(wunifrac_dist),
          file = file.path(beta_dir, "weighted_unifrac_distance_matrix.csv"))

# --  Hierarchical clustering dendrogram  -------------------------------------
png(file.path(fig_dir, "bray_dendrogram.png"), width = 2000, height = 1200, res = 150)
plot(hclust(bray_dist), main = "Bray-Curtis Hierarchical Clustering",
     xlab = "", sub = "")
dev.off()

# --  PCoA: Bray-Curtis  -------------------------------------------------------
ordu      <- ordinate(ps, method = "PCoA", distance = "bray")
p_bray    <- plot_ordination(ps, ordu, color = "group") +
  geom_point(size = 4) +
  stat_ellipse(type = "t", linetype = 2) +
  theme_bw() +
  labs(title = "PCoA – Bray-Curtis")

ggsave(file.path(fig_dir, "PCoA_bray.png"), p_bray, width = 8, height = 6)

bray_coords        <- as.data.frame(ordu$vectors)
bray_coords$Sample <- rownames(bray_coords)
fwrite(bray_coords, file = file.path(ord_dir, "PCoA_bray_coordinates.csv"))

# --  PCoA: Weighted UniFrac  --------------------------------------------------
ordu_unifrac    <- ordinate(ps, method = "PCoA", distance = "wunifrac")
p_unifrac       <- plot_ordination(ps, ordu_unifrac, color = "group") +
  geom_point(size = 4) +
  stat_ellipse(type = "t", linetype = 2) +
  theme_bw() +
  labs(title = "PCoA – Weighted UniFrac")

ggsave(file.path(fig_dir, "PCoA_unifrac.png"), p_unifrac, width = 8, height = 6)

unifrac_coords        <- as.data.frame(ordu_unifrac$vectors)
unifrac_coords$Sample <- rownames(unifrac_coords)
fwrite(unifrac_coords, file = file.path(ord_dir, "PCoA_unifrac_coordinates.csv"))

# --  NMDS: Bray-Curtis  -------------------------------------------------------
ordu_nmds <- ordinate(ps, method = "NMDS", distance = "bray")
p_nmds    <- plot_ordination(ps, ordu_nmds, color = "group") +
  geom_point(size = 4) +
  stat_ellipse(type = "t", linetype = 2) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Stress = ", round(ordu_nmds$stress, 3)),
           hjust = -0.1, vjust = 1.5, size = 4) +
  theme_bw() +
  labs(title = "NMDS – Bray-Curtis")

ggsave(file.path(fig_dir, "NMDS_bray.png"), p_nmds, width = 8, height = 6)

nmds_coords        <- as.data.frame(ordu_nmds$points)
nmds_coords$Sample <- rownames(nmds_coords)
fwrite(nmds_coords, file = file.path(ord_dir, "NMDS_bray_coordinates.csv"))

# --  PERMANOVA (adonis2)  -----------------------------------------------------
sample_df   <- as(sample_data(ps), "data.frame")

perm_bray   <- adonis2(bray_dist ~ group, data = sample_df, permutations = 999)
perm_unifrac <- adonis2(wunifrac_dist ~ group, data = sample_df, permutations = 999)

cat("=== PERMANOVA Bray-Curtis ===\n"); print(perm_bray)
cat("=== PERMANOVA Weighted UniFrac ===\n"); print(perm_unifrac)

write.csv(as.data.frame(perm_bray),
          file = file.path(beta_dir, "PERMANOVA_bray.csv"))
write.csv(as.data.frame(perm_unifrac),
          file = file.path(beta_dir, "PERMANOVA_unifrac.csv"))

# --  ANOSIM  ------------------------------------------------------------------
anosim_bray    <- anosim(bray_dist,    sample_df$group, permutations = 999)
anosim_unifrac <- anosim(wunifrac_dist, sample_df$group, permutations = 999)

cat("=== ANOSIM Bray-Curtis ===\n");       print(summary(anosim_bray))
cat("=== ANOSIM Weighted UniFrac ===\n");  print(summary(anosim_unifrac))

anosim_results <- data.frame(
  Distance  = c("Bray-Curtis", "Weighted UniFrac"),
  R_statistic = c(anosim_bray$statistic,    anosim_unifrac$statistic),
  p_value   = c(anosim_bray$signif,         anosim_unifrac$signif)
)
write.csv(anosim_results, file = file.path(beta_dir, "ANOSIM_results.csv"),
          row.names = FALSE)

# --  Betadisper (homogeneity of within-group dispersion)  --------------------
betadisp_bray    <- betadisper(bray_dist,    sample_df$group)
betadisp_unifrac <- betadisper(wunifrac_dist, sample_df$group)

betadisp_test_bray    <- permutest(betadisp_bray,    permutations = 999)
betadisp_test_unifrac <- permutest(betadisp_unifrac, permutations = 999)

cat("=== Betadisper Bray-Curtis ===\n");       print(betadisp_test_bray)
cat("=== Betadisper Weighted UniFrac ===\n");  print(betadisp_test_unifrac)

betadisp_results <- data.frame(
  Distance    = c("Bray-Curtis", "Weighted UniFrac"),
  F_statistic = c(betadisp_test_bray$tab$F[1],    betadisp_test_unifrac$tab$F[1]),
  p_value     = c(betadisp_test_bray$tab$`Pr(>F)`[1],
                  betadisp_test_unifrac$tab$`Pr(>F)`[1])
)
write.csv(betadisp_results, file = file.path(beta_dir, "Betadisper_results.csv"),
          row.names = FALSE)

png(file.path(fig_dir, "betadisper_bray.png"), width = 1600, height = 1200, res = 150)
plot(betadisp_bray, main = "Betadisper – Bray-Curtis")
dev.off()

# ==============================================================================
# 11. DIFFERENTIAL ABUNDANCE
# ==============================================================================

# Helper: make group a factor with reference level (change "Control" if needed)
ref_level <- "Control"     # <-- edit to match your metadata

# ==============================================================================
# 11.a. DESeq2
# ==============================================================================
ps_deseq <- ps_genus

# DESeq2 requires integer counts — round if needed
otu_table(ps_deseq) <- otu_table(
  round(as(otu_table(ps_deseq), "matrix")),
  taxa_are_rows = taxa_are_rows(ps_deseq)
)

sample_data(ps_deseq)$group <- relevel(
  factor(sample_data(ps_deseq)$group), ref = ref_level
)

dds <- phyloseq_to_deseq2(ps_deseq, ~ group)

# Geometric mean with zero-robust estimator
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, test = "Wald", fitType = "local")

res_deseq <- results(dds, alpha = 0.05)
res_deseq_df <- as.data.frame(res_deseq)
res_deseq_df$Taxon     <- rownames(res_deseq_df)
res_deseq_df$Genus     <- as.vector(tax_table(ps_genus)[rownames(res_deseq_df), "Genus"])
res_deseq_df$Direction <- ifelse(res_deseq_df$log2FoldChange > 0,
                                 paste0("Up in non-", ref_level),
                                 paste0("Up in ", ref_level))
res_deseq_df           <- res_deseq_df[order(res_deseq_df$padj, na.last = TRUE), ]

fwrite(res_deseq_df,
       file = file.path(da_dir, "DESeq2_results.csv"))

cat("\n=== Top DESeq2 significant taxa (padj < 0.05) ===\n")
print(head(res_deseq_df[!is.na(res_deseq_df$padj) & res_deseq_df$padj < 0.05, ], 20))

# --  Volcano plot (DESeq2)  ---------------------------------------------------
sig_deseq <- res_deseq_df[!is.na(res_deseq_df$padj), ]

p_volcano <- EnhancedVolcano::EnhancedVolcano(sig_deseq,
                                              lab          = sig_deseq$Genus,
                                              x            = "log2FoldChange",
                                              y            = "padj",
                                              pCutoff      = 0.05,
                                              FCcutoff     = 1,
                                              title        = "DESeq2 – Volcano Plot",
                                              subtitle     = paste0("Reference: ", ref_level),
                                              legendPosition = "right",
                                              drawConnectors = TRUE,
                                              widthConnectors = 0.3
)
ggsave(file.path(fig_dir, "DESeq2_volcano.png"),
       p_volcano, width = 12, height = 8)

# ==============================================================================
# 11.b. ANCOM-BC
# ==============================================================================
ancombc_out <- ancombc2(
  data         = ps_genus,
  fix_formula  = "group",
  p_adj_method = "BH",
  lib_cut      = 1000,
  group        = "group",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  global       = TRUE
)

res_ancombc        <- ancombc_out$res
res_ancombc$Genus  <- as.vector(
  tax_table(ps_genus)[res_ancombc$taxon, "Genus"]
)

fwrite(res_ancombc,
       file = file.path(da_dir, "ANCOMBC_results.csv"))

cat("\n=== ANCOM-BC significant taxa ===\n")
ancombc_sig_col <- grep("^diff_", colnames(res_ancombc), value = TRUE)[1]
if (!is.null(ancombc_sig_col)) {
  print(res_ancombc[res_ancombc[[ancombc_sig_col]], ])
}

# ==============================================================================
# 11.c. Kruskal-Wallis + FDR correction (per taxon)
# ==============================================================================
genus_rel_df        <- as.data.frame(genus_rel)
genus_rel_df$group  <- sample_data(ps_genus)$group

kw_results <- lapply(genus_names, function(g) {
  vals  <- genus_rel_df[[g]]
  grps  <- genus_rel_df$group
  # skip taxa with zero variance
  if (length(unique(vals)) <= 1) return(NULL)
  kt    <- kruskal.test(vals ~ grps)
  data.frame(Genus   = g,
             chi_sq  = kt$statistic,
             df      = kt$parameter,
             p.value = kt$p.value,
             stringsAsFactors = FALSE)
})
kw_df           <- do.call(rbind, Filter(Negate(is.null), kw_results))
kw_df$p.adjusted <- p.adjust(kw_df$p.value, method = "BH")
kw_df           <- kw_df[order(kw_df$p.adjusted), ]

fwrite(kw_df, file = file.path(da_dir, "KruskalWallis_FDR_results.csv"))

cat("\n=== Kruskal-Wallis significant taxa (FDR < 0.05) ===\n")
print(head(kw_df[kw_df$p.adjusted < 0.05, ], 20))

# ==============================================================================
# 11.d. LEfSe-style LDA score plot
# ==============================================================================
# lefser requires a SummarizedExperiment object built from a phyloseq
# We use the genus-level relative abundance


# How many samples per group?
table(sample_data(ps_genus)$group)

mm_res <- run_lefse(
  ps_genus,
  group    = "group",
  taxa_rank = "Genus",
  kw_cutoff  = 0.05,
  wilcoxon_cutoff = 0.05,
  lda_cutoff = 1.0  # changed from 2.0 to 1.0 
)
marker_table(mm_res)

p_lefse <- plot_ef_bar(mm_res) +
  labs(title = "LEfSe – LDA Score Plot") +
  theme_bw()

# Export
write.csv(
  data.frame(mm_res@marker_table),
  file = file.path(da_dir, "LEfSe_results.csv"),
  row.names = FALSE
)

ggsave(file.path(fig_dir, "LEfSe_LDA_scores.png"),
       p_lefse, width = 10, height = 8)

# ==============================================================================
# 12. HEATMAP OF SIGNIFICANT TAXA
# ==============================================================================

# Use top significant taxa from DESeq2 (adjust source as needed)
n_top_taxa <- 30   # number of top taxa to show

sig_taxa <- res_deseq_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = n_top_taxa) %>%
  pull(Genus)

# If fewer than n_top_taxa are significant, fall back to top by p-value
if (length(sig_taxa) < 5) {
  sig_taxa <- res_deseq_df %>%
    arrange(pvalue) %>%
    slice_head(n = n_top_taxa) %>%
    pull(Genus)
  cat("Note: fewer than 5 significant taxa at padj<0.05; using top by p-value.\n")
}

# Build matrix: rows = taxa, cols = samples
heat_mat <- t(genus_rel)[sig_taxa, , drop = FALSE]

# Z-score per taxon (row-wise)
heat_z <- t(scale(t(heat_mat)))

# Annotation for columns (samples)
col_annot_df <- data.frame(
  Group = sample_df$group,
  row.names = rownames(sample_df)
)

# Get groups in a consistent order
group_levels <- unique(sample_df$group)
n_groups <- length(group_levels)

group_colours <- setNames(
  brewer.pal(max(n_groups, 3), "Set1")[1:n_groups],  # subset to exact n
  group_levels
)


annot_colours <- list(Group = group_colours)

# --  pheatmap  ----------------------------------------------------------------
png(file.path(fig_dir, "heatmap_sig_taxa_pheatmap.png"),
    width = 2400, height = 2000, res = 200)
pheatmap(
  heat_z,
  annotation_col     = col_annot_df,
  annotation_colors  = annot_colours,
  color              = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  cluster_rows       = TRUE,
  cluster_cols       = TRUE,
  clustering_method  = "ward.D2",
  show_colnames      = FALSE,
  fontsize_row       = 9,
  main               = paste0("Top ", length(sig_taxa), " Significant Taxa\n(Z-scored Relative Abundance)")
)
dev.off()

# --  ComplexHeatmap (publication-quality alternative)  -----------------------
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#D6604D")
)
group_vec    <- col_annot_df$Group
col_annot_ch <- HeatmapAnnotation(
  Group = group_vec,
  col   = list(Group = group_colours),
  annotation_legend_param = list(Group = list(title = "Group"))
)

png(file.path(fig_dir, "heatmap_sig_taxa_ComplexHeatmap.png"),
    width = 2400, height = 2000, res = 200)
Heatmap(
  heat_z,
  name               = "Z-score",
  col                = col_fun,
  top_annotation     = col_annot_ch,
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_column_names  = FALSE,
  row_names_gp       = gpar(fontsize = 9),
  column_title       = paste0("Top ", length(sig_taxa),
                              " Significant Taxa (Z-scored Rel. Abundance)")
)
dev.off()

# ==============================================================================
# 13. STACKED BAR CHARTS
# ==============================================================================

# Top N genera by mean relative abundance
n_bar   <- 15
top_gen <- names(sort(colMeans(genus_rel), decreasing = TRUE))[1:n_bar]

bar_long <- as.data.table(genus_rel[, top_gen, drop = FALSE],
                          keep.rownames = "Sample")
bar_long <- melt(bar_long, id.vars = "Sample",
                 variable.name = "Genus", value.name = "RelAbund")
bar_long$Group <- sample_df[bar_long$Sample, "group"]

# Per-group mean
bar_mean <- bar_long[, .(MeanRelAbund = mean(RelAbund)),
                     by = .(Group, Genus)]

col_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Accent"))[1:n_bar]

p_bar_sample <- ggplot(bar_long, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  labs(title = paste0("Top ", n_bar, " Genera – Per Sample"),
       x = "Sample", y = "Relative Abundance")

ggsave(file.path(fig_dir, "stacked_bar_per_sample.png"),
       p_bar_sample, width = 14, height = 6)

p_bar_group <- ggplot(bar_mean, aes(x = Group, y = MeanRelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  labs(title = paste0("Top ", n_bar, " Genera – Group Mean"),
       x = "Group", y = "Mean Relative Abundance")

ggsave(file.path(fig_dir, "stacked_bar_group_mean.png"),
       p_bar_group, width = 7, height = 6)

# ==============================================================================
# 14. BOX / VIOLIN PLOTS 
# ==============================================================================

# Pick top 9 taxa from Kruskal-Wallis for individual plots
top_kw_taxa <- head(kw_df$Genus[kw_df$p.adjusted < 0.05], 9)
if (length(top_kw_taxa) == 0) top_kw_taxa <- head(kw_df$Genus, 9)

violin_long <- as.data.table(genus_rel[, top_kw_taxa, drop = FALSE],
                             keep.rownames = "Sample")
violin_long$Group <- sample_df[violin_long$Sample, "group"]
violin_long <- melt(violin_long, id.vars = c("Sample", "Group"),
                    variable.name = "Genus", value.name = "RelAbund")

p_violin <- ggplot(violin_long,
                   aes(x = Group, y = RelAbund, fill = Group)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 1, fill = "white") +
  geom_jitter(width = 0.07, size = 1, alpha = 0.5) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Top Differential Taxa – Relative Abundance",
       x = "Group", y = "Relative Abundance")

ggsave(file.path(fig_dir, "violin_top_taxa.png"),
       p_violin, width = 12, height = 10)

# ==============================================================================
# 15. PHYLOGENETIC TREE PLOTS
# ==============================================================================

top50 <- names(sort(taxa_sums(ps), TRUE))[1:50]

ggsave(
  file.path(fig_dir, "genus_tree.png"),
  plot  = plot_tree(ps_genus, color = "Phylum", label.tips = "Genus"),
  width = 2200/250, height = 1800/250, dpi = 250
)

ggsave(
  file.path(fig_dir, "top50_ASV_tree.png"),
  plot  = plot_tree(prune_taxa(top50, ps)),
  width = 2200/250, height = 1800/250, dpi = 250
)

# ==============================================================================
# 16. BASIC PHYLOSEQ HEATMAP
# ==============================================================================
png(file.path(fig_dir, "genus_heatmap_phyloseq.png"),
    width = 2200, height = 1800, res = 250)
plot_heatmap(ps_genus)
dev.off()

# ==============================================================================
# 17. EXPORT RESULTS
# ==============================================================================
fwrite(asv_map, file = file.path(results_dir, "asv_sequences.csv"))

# Clean up sample names for export - adjust to your sample names if needed
rn  <- rownames(genus_counts)
pat <- ".*(MP\\d+)_(BL|Y[12]|T[012]).*"
new <- ifelse(grepl(pat, rn), sub(pat, "\\1_\\2", rn), rn)
new <- sub("T([12])$", "Y\\1", new)
new <- sub("T0$",      "BL",   new)
rownames(genus_counts) <- new
rownames(genus_rel)    <- new

fwrite(cbind(Sample = rownames(genus_counts),
             as.data.frame(genus_counts, check.names = FALSE)),
       file = file.path(results_dir, "counts_by_Genus.csv"))
fwrite(cbind(Sample = rownames(genus_rel),
             as.data.frame(genus_rel, check.names = FALSE)),
       file = file.path(results_dir, "relative_abundance_by_Genus.csv"))

genus_long <- as.data.table(genus_counts, keep.rownames = "Sample")
genus_long <- melt(genus_long, id.vars = "Sample",
                   variable.name = "Genus", value.name = "Count")
fwrite(genus_long[Count > 0],
       file = file.path(results_dir, "counts_by_Genus_long.csv"))

saveRDS(ps,       file = file.path(results_dir, "phyloseq_ASV.rds"))
saveRDS(ps_genus, file = file.path(results_dir, "phyloseq_Genus.rds"))

fwrite(alpha_df, file = file.path(alpha_dir, "alpha_diversity.csv"))

# ==============================================================================
# 18. READ TRACKING
# ==============================================================================
getN <- function(x) sum(getUniques(x))
merged_counts <- sapply(mergers, function(x) {
  if (is.null(x) || nrow(x) == 0) 0 else sum(x$abundance)
})

track <- data.frame(
  Sample    = sample.names,
  input     = out[, "reads.in"],
  filtered  = out[, "reads.out"],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged    = merged_counts[sample.names],
  nonchim   = rowSums(seqtab.nochim)[sample.names],
  row.names = NULL
)
fwrite(track, file = file.path(results_dir, "read_tracking.csv"))

# ==============================================================================
# End of pipeline
# ==============================================================================
cat("\n=== Pipeline complete. All outputs saved to:", results_dir, "===\n")
cat("\nKey output files:\n")
cat("  figures/              – all plots (PNG)\n")
cat("  alpha_diversity/      – Shannon/Observed/Chao1/Simpson + Wilcoxon tests\n")
cat("  beta_diversity/       – PERMANOVA, ANOSIM, Betadisper, distance matrices\n")
cat("  ordination/           – PCoA + NMDS coordinates\n")
cat("  differential_abundance/ – DESeq2, ANCOM-BC, KW+FDR, LEfSe\n")
cat("  results/              – count/relabund tables, phyloseq .rds objects\n")

