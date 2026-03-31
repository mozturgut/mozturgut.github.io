#!/usr/bin/env Rscript
# CELLCALL NEUTROPHIL SIGNALING ANALYSIS - QIN DATASET ONLY
# Research question: Signals received by Arg1+ vs Arg1- neutrophils (Arg1+ ~ anti-inflammatory bias;
#   Arg1- ~ pro-inflammatory bias as a working hypothesis). Baseline + global Arg1-labeled analyses.
#
# RUN ORDER: after upstream Seurat objects exist; independent of Lee/Wang CellCall but often run after CellChat Qin.
# Required in getwd(): QinDat.rds (or qindat.rds — master normalizes name)
# Rerun from: line 15 (suppressPackageStartupMessages) for full run.
# Optional snapshot: # load("CELLCALL_MAIN_RESULTS/CELLCALL_QIN_ONLY_Final_Environment.RData")
#
# ---- Working directory ----
CELLCOMM_PROJECT_ROOT <- Sys.getenv("CELLCOMM_PROJECT_ROOT", unset = "")
if (nzchar(CELLCOMM_PROJECT_ROOT) && dir.exists(CELLCOMM_PROJECT_ROOT)) setwd(CELLCOMM_PROJECT_ROOT)

suppressPackageStartupMessages({
  library(Seurat)
  library(cellcall)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(circlize)
  library(grid)
  library(pheatmap)
})

options(repr.plot.width = 8, repr.plot.height = 6)
if (capabilities("png")) {
  png_width <- 800
  png_height <- 600
} else {
  png_width <- 800
  png_height <- 600
}

OUTPUT_DIR <- "CELLCALL_MAIN_RESULTS"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory created:", OUTPUT_DIR, "\n")
cat("All results will be saved in:", OUTPUT_DIR, "\n\n")

# =============================================================================
# LOAD QIN DATASET ONLY
# =============================================================================
cat("=== LOAD QIN DATASET ===\n")
qin_data <- readRDS("QinDat.rds")
DefaultAssay(qin_data) <- "RNA"

# CRITICAL: Convert Qin Seurat v5 to v4 format for CellCall compatibility
if (inherits(qin_data[["RNA"]], "Assay5")) {
  cat("Detected Seurat v5 Assay5, converting to v4 format...\n")
  qin_data[["RNA"]] <- JoinLayers(qin_data[["RNA"]])
  options(Seurat.object.assay.version = "v3")
  qin_data[["RNA"]] <- as(qin_data[["RNA"]], Class = "Assay")
  cat("Conversion complete. Assay class:", class(qin_data[["RNA"]]), "\n")
} else {
  cat("Already in v4 format. Assay class:", class(qin_data[["RNA"]]), "\n")
}

cat("Qin cells:", ncol(qin_data), " | genes:", nrow(qin_data), "\n")
gc()

# =============================================================================
# STANDARDIZE TIMEPOINTS AND DEFINE ARG1 STATUS
# =============================================================================
cat("\n=== STANDARDIZE TIMEPOINTS (QIN) ===\n")
qin_data$time <- as.character(qin_data$time)
qin_data$time <- gsub("^Uninjured$", "0", qin_data$time)
qin_data$time <- gsub("^1dpi$", "1", qin_data$time)
qin_data$time <- gsub("^3dpi$", "3", qin_data$time)
qin_data$time <- gsub("^7dpi$", "7", qin_data$time)
qin_data$time <- as.numeric(qin_data$time)
cat("QIN timepoints:", paste(sort(unique(qin_data$time)), collapse = ", "), "\n")

cat("\n=== DEFINE ARG1 STATUS BASED ON ACTUAL EXPRESSION ===\n")
qin_arg1_expr <- GetAssayData(qin_data, assay = "RNA", slot = "data")["Arg1", ]
qin_data$Arg1_status <- "Arg1neg"
qin_data$Arg1_status[qin_arg1_expr > 0] <- "Arg1pos"
cat("QIN Arg1+:", sum(qin_data$Arg1_status == "Arg1pos"), " | Arg1-:", sum(qin_data$Arg1_status == "Arg1neg"), "\n")
gc()

# =============================================================================
# ANALYSIS 1: QIN DAY 1 UNMODIFIED (BASELINE)
# =============================================================================
cat("\n=== ANALYSIS 1: QIN DAY 1 UNMODIFIED (BASELINE) ===\n")
qin_day1_unmodified <- subset(qin_data, time == 1)
cat("Qin Day 1 cells:", ncol(qin_day1_unmodified), "\n")
new_cell_ids_qin_d1 <- gsub("-", "_", colnames(qin_day1_unmodified))
qin_day1_unmodified <- RenameCells(qin_day1_unmodified, new.names = new_cell_ids_qin_d1)
Idents(qin_day1_unmodified) <- as.character(qin_day1_unmodified$pruned_labels)
gc()
cc_qin_d1_unmod <- CreateObject_fromSeurat(Seurat.object = qin_day1_unmodified,
                                           slot = "counts",
                                           cell_type = "pruned_labels",
                                           data_source = "UMI",
                                           scale.factor = 10^6,
                                           Org = "Mus musculus")
cc_qin_d1_unmod <- TransCommuProfile(object = cc_qin_d1_unmod,
                                     pValueCor = 0.1,
                                     CorValue = 0.05,
                                     topTargetCor = 1,
                                     p.adjust = 0.1,
                                     use.type = "mean",
                                     probs = 0.1,
                                     method = "weighted",
                                     IS_core = TRUE,
                                     Org = "Mus musculus")
saveRDS(cc_qin_d1_unmod, file.path(OUTPUT_DIR, "CellCall_QIN_Day1_UNMODIFIED.rds"))
save.image(file.path(OUTPUT_DIR, "Analysis1_Qin_Day1_Unmodified_After_TransCommu.RData"))
cat("Analysis 1 (Qin Day 1 Unmodified) TransCommuProfile complete\n")

# Visualizations Analysis 1
neutrophil_cols_unmod1 <- grepl("-Neutrophil", colnames(cc_qin_d1_unmod@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_sending_cols_unmod1 <- grepl("^Neutrophil.*-", colnames(cc_qin_d1_unmod@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_cols_unmod1 <- neutrophil_cols_unmod1 & !neutrophil_sending_cols_unmod1
neutrophil_matrix_unmod1 <- cc_qin_d1_unmod@data$expr_l_r_log2_scale[, neutrophil_cols_unmod1, drop = FALSE]
top_n <- min(20, nrow(neutrophil_matrix_unmod1))
mean_signals <- rowMeans(neutrophil_matrix_unmod1, na.rm = TRUE)
top_lr_pairs <- names(sort(mean_signals, decreasing = TRUE))[1:top_n]
neutrophil_matrix_unmod1_top <- neutrophil_matrix_unmod1[top_lr_pairs, , drop = FALSE]
set.seed(123)
noise_matrix <- matrix(rnorm(nrow(neutrophil_matrix_unmod1_top) * ncol(neutrophil_matrix_unmod1_top), mean = 0, sd = 1e-10),
                       nrow = nrow(neutrophil_matrix_unmod1_top), ncol = ncol(neutrophil_matrix_unmod1_top))
neutrophil_matrix_unmod1_noise <- neutrophil_matrix_unmod1_top + noise_matrix
p_hm_unmod1 <- pheatmap::pheatmap(neutrophil_matrix_unmod1_noise,
  color = colorRampPalette(c("blue", "white", "yellow", "orange", "red"))(100),
  breaks = seq(0, 1, length.out = 101), show_rownames = TRUE, show_colnames = TRUE,
  treeheight_row = 0, treeheight_col = 10, cluster_rows = TRUE, cluster_cols = FALSE,
  fontsize = 12, angle_col = 90, fontsize_row = 12, fontsize_col = 12,
  main = "Analysis 1: Qin Day1 Unmodified (Top 20 Neutrophil-Receiving Signals)")
print(p_hm_unmod1)
circlize::circos.clear()
unique_cell_types1 <- unique(gsub("-.*", "", colnames(cc_qin_d1_unmod@data$expr_l_r_log2_scale)))
cell_color_unmod1 <- data.frame(color = rainbow(length(unique_cell_types1)), stringsAsFactors = FALSE)
rownames(cell_color_unmod1) <- unique_cell_types1
ViewInterCircos(object = cc_qin_d1_unmod, font = 2, cellColor = cell_color_unmod1, lrColor = c("#F16B6F", "#84B1ED"),
  arr.type = "big.arrow", arr.length = 0.04, trackhight1 = 0.05, slot = "expr_l_r_log2_scale",
  linkcolor.from.sender = TRUE, linkcolor = NULL, gap.degree = 0.1, trackhight2 = 0.032, track.margin2 = c(0.01, 0.12), DIY = FALSE)
grid.text("Analysis 1: Qin Day 1 Unmodified - Neutrophil-Receiving Cell-Cell Communication", y = 0.95, gp = gpar(fontsize = 16))
n_unmod1_neutrophil <- cc_qin_d1_unmod@data$expr_l_r_log2_scale[, neutrophil_cols_unmod1, drop = FALSE]
pathway.hyper.list_unmod1 <- list()
for(i in colnames(n_unmod1_neutrophil)) { tmp <- getHyperPathway(data = n_unmod1_neutrophil, object = cc_qin_d1_unmod, cella_cellb = i, Org = "Mus musculus", IS_core = TRUE); pathway.hyper.list_unmod1[[i]] <- tmp }
pathway.hyper.list_unmod1 <- pathway.hyper.list_unmod1[!sapply(pathway.hyper.list_unmod1, is.null)]
myPub.df_unmod1 <- getForBubble(pathway.hyper.list_unmod1, cella_cellb = colnames(n_unmod1_neutrophil))
myPub.df_unmod1$label_text <- ""; myPub.df_unmod1$label_text[myPub.df_unmod1$p.adjust < 0.05] <- format(myPub.df_unmod1$p.adjust[myPub.df_unmod1$p.adjust < 0.05], digits = 2, scientific = TRUE)
p_bubble_unmod1 <- plotBubble(myPub.df_unmod1) + geom_text_repel(aes(label = label_text), size = 3) + scale_size_continuous(limits = c(0, 100), range = c(1, 10)) + scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.05)) + ggtitle("Qin Day 1 UNMODIFIED: Neutrophil-Receiving Pathway Analysis (with p.adjust)") + theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_bubble_unmod1)
rm(cc_qin_d1_unmod, qin_day1_unmodified, p_hm_unmod1, cell_color_unmod1, unique_cell_types1, neutrophil_matrix_unmod1, neutrophil_matrix_unmod1_top, neutrophil_matrix_unmod1_noise, pathway.hyper.list_unmod1, myPub.df_unmod1, p_bubble_unmod1, neutrophil_cols_unmod1, neutrophil_sending_cols_unmod1, n_unmod1_neutrophil)
gc()
cat("Memory cleanup completed after Analysis 1\n")

# =============================================================================
# ANALYSIS 2: QIN DAY 3 UNMODIFIED (BASELINE)
# =============================================================================
cat("\n=== ANALYSIS 2: QIN DAY 3 UNMODIFIED (BASELINE) ===\n")
qin_day3_unmodified <- subset(qin_data, time == 3)
cat("Qin Day 3 cells:", ncol(qin_day3_unmodified), "\n")
new_cell_ids_qin_d3 <- gsub("-", "_", colnames(qin_day3_unmodified))
qin_day3_unmodified <- RenameCells(qin_day3_unmodified, new.names = new_cell_ids_qin_d3)
Idents(qin_day3_unmodified) <- as.character(qin_day3_unmodified$pruned_labels)
gc()
cc_qin_d3_unmod <- CreateObject_fromSeurat(Seurat.object = qin_day3_unmodified,
                                           slot = "counts",
                                           cell_type = "pruned_labels",
                                           data_source = "UMI",
                                           scale.factor = 10^6,
                                           Org = "Mus musculus")
cc_qin_d3_unmod <- TransCommuProfile(object = cc_qin_d3_unmod,
                                     pValueCor = 0.1,
                                     CorValue = 0.05,
                                     topTargetCor = 1,
                                     p.adjust = 0.1,
                                     use.type = "mean",
                                     probs = 0.1,
                                     method = "weighted",
                                     IS_core = TRUE,
                                     Org = "Mus musculus")
saveRDS(cc_qin_d3_unmod, file.path(OUTPUT_DIR, "CellCall_QIN_Day3_UNMODIFIED.rds"))
save.image(file.path(OUTPUT_DIR, "Analysis2_Qin_Day3_Unmodified_After_TransCommu.RData"))
cat("Analysis 2 (Qin Day 3 Unmodified) TransCommuProfile complete\n")

# Visualizations Analysis 2
neutrophil_cols_unmod2 <- grepl("-Neutrophil", colnames(cc_qin_d3_unmod@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_sending_cols_unmod2 <- grepl("^Neutrophil.*-", colnames(cc_qin_d3_unmod@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_cols_unmod2 <- neutrophil_cols_unmod2 & !neutrophil_sending_cols_unmod2
neutrophil_matrix_unmod2 <- cc_qin_d3_unmod@data$expr_l_r_log2_scale[, neutrophil_cols_unmod2, drop = FALSE]
top_n <- min(20, nrow(neutrophil_matrix_unmod2))
mean_signals <- rowMeans(neutrophil_matrix_unmod2, na.rm = TRUE)
top_lr_pairs <- names(sort(mean_signals, decreasing = TRUE))[1:top_n]
neutrophil_matrix_unmod2_top <- neutrophil_matrix_unmod2[top_lr_pairs, , drop = FALSE]
set.seed(123)
noise_matrix <- matrix(rnorm(nrow(neutrophil_matrix_unmod2_top) * ncol(neutrophil_matrix_unmod2_top), mean = 0, sd = 1e-10),
                       nrow = nrow(neutrophil_matrix_unmod2_top), ncol = ncol(neutrophil_matrix_unmod2_top))
neutrophil_matrix_unmod2_noise <- neutrophil_matrix_unmod2_top + noise_matrix
p_hm_unmod2 <- pheatmap::pheatmap(neutrophil_matrix_unmod2_noise,
  color = colorRampPalette(c("blue", "white", "yellow", "orange", "red"))(100),
  breaks = seq(0, 1, length.out = 101), show_rownames = TRUE, show_colnames = TRUE,
  treeheight_row = 0, treeheight_col = 10, cluster_rows = TRUE, cluster_cols = FALSE,
  fontsize = 12, angle_col = 90, fontsize_row = 12, fontsize_col = 12,
  main = "Analysis 2: Qin Day3 Unmodified (Top 20 Neutrophil-Receiving Signals)")
print(p_hm_unmod2)
circlize::circos.clear()
unique_cell_types2 <- unique(gsub("-.*", "", colnames(cc_qin_d3_unmod@data$expr_l_r_log2_scale)))
cell_color_unmod2 <- data.frame(color = rainbow(length(unique_cell_types2)), stringsAsFactors = FALSE)
rownames(cell_color_unmod2) <- unique_cell_types2
ViewInterCircos(object = cc_qin_d3_unmod, font = 2, cellColor = cell_color_unmod2, lrColor = c("#F16B6F", "#84B1ED"),
  arr.type = "big.arrow", arr.length = 0.04, trackhight1 = 0.05, slot = "expr_l_r_log2_scale",
  linkcolor.from.sender = TRUE, linkcolor = NULL, gap.degree = 0.1, trackhight2 = 0.032, track.margin2 = c(0.01, 0.12), DIY = FALSE)
grid.text("Analysis 2: Qin Day 3 Unmodified - Neutrophil-Receiving Cell-Cell Communication", y = 0.95, gp = gpar(fontsize = 16))
n_unmod2_neutrophil <- cc_qin_d3_unmod@data$expr_l_r_log2_scale[, neutrophil_cols_unmod2, drop = FALSE]
pathway.hyper.list_unmod2 <- list()
for(i in colnames(n_unmod2_neutrophil)) { tmp <- getHyperPathway(data = n_unmod2_neutrophil, object = cc_qin_d3_unmod, cella_cellb = i, Org = "Mus musculus", IS_core = TRUE); pathway.hyper.list_unmod2[[i]] <- tmp }
pathway.hyper.list_unmod2 <- pathway.hyper.list_unmod2[!sapply(pathway.hyper.list_unmod2, is.null)]
myPub.df_unmod2 <- getForBubble(pathway.hyper.list_unmod2, cella_cellb = colnames(n_unmod2_neutrophil))
myPub.df_unmod2$label_text <- ""; myPub.df_unmod2$label_text[myPub.df_unmod2$p.adjust < 0.05] <- format(myPub.df_unmod2$p.adjust[myPub.df_unmod2$p.adjust < 0.05], digits = 2, scientific = TRUE)
p_bubble_unmod2 <- plotBubble(myPub.df_unmod2) + geom_text_repel(aes(label = label_text), size = 3) + scale_size_continuous(limits = c(0, 100), range = c(1, 10)) + scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.05)) + ggtitle("Qin Day 3 UNMODIFIED: Neutrophil-Receiving Pathway Analysis (with p.adjust)") + theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_bubble_unmod2)
rm(cc_qin_d3_unmod, qin_day3_unmodified, p_hm_unmod2, cell_color_unmod2, unique_cell_types2, neutrophil_matrix_unmod2, neutrophil_matrix_unmod2_top, neutrophil_matrix_unmod2_noise, pathway.hyper.list_unmod2, myPub.df_unmod2, p_bubble_unmod2, neutrophil_cols_unmod2, neutrophil_sending_cols_unmod2, n_unmod2_neutrophil)
gc()
cat("Memory cleanup completed after Analysis 2\n")

# =============================================================================
# GLOBAL ANALYSIS 1: QIN DAY 1 (Arg1-labeled)
# =============================================================================
cat("\n=== GLOBAL ANALYSIS 1: QIN DAY 1 (Arg1+ vs Arg1- Comparison) ===\n")
qin_day1_global <- subset(qin_data, time == 1)
qin_day1_global$cellcall_label <- gsub("-", "", as.character(qin_day1_global$pruned_labels))
neut_idx_q1 <- grepl("Neutrophil", qin_day1_global$cellcall_label, ignore.case = TRUE)
qin_day1_global$cellcall_label[neut_idx_q1] <- paste0("Neutrophil", qin_day1_global$Arg1_status[neut_idx_q1])
cat("Qin Day 1 cells:", ncol(qin_day1_global), "\n"); print(table(qin_day1_global$cellcall_label))
qin_day1_global <- RenameCells(qin_day1_global, new.names = gsub("-", "_", colnames(qin_day1_global)))
Idents(qin_day1_global) <- as.character(qin_day1_global$cellcall_label)
cc_qin_d1_global <- CreateObject_fromSeurat(Seurat.object = qin_day1_global, slot = "counts", cell_type = "cellcall_label", data_source = "UMI", scale.factor = 10^6, Org = "Mus musculus")
cc_qin_d1_global <- TransCommuProfile(object = cc_qin_d1_global, pValueCor = 0.1, CorValue = 0.05, topTargetCor = 1, p.adjust = 0.1, use.type = "mean", probs = 0.1, method = "weighted", IS_core = TRUE, Org = "Mus musculus")
saveRDS(cc_qin_d1_global, file.path(OUTPUT_DIR, "CellCall_QIN_Day1_GLOBAL.rds"))
save.image(file.path(OUTPUT_DIR, "GlobalAnalysis1_Qin_Day1_After_TransCommu.RData"))

neutrophil_cols_global1 <- grepl("-Neutrophil", colnames(cc_qin_d1_global@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_sending_global1 <- grepl("^Neutrophil.*-", colnames(cc_qin_d1_global@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_cols_global1 <- neutrophil_cols_global1 & !neutrophil_sending_global1
neutrophil_receiving_global1 <- cc_qin_d1_global@data$expr_l_r_log2_scale[, neutrophil_cols_global1, drop = FALSE]
arg1pos_receiving_global1 <- neutrophil_receiving_global1[, grepl("NeutrophilArg1pos$", colnames(neutrophil_receiving_global1)), drop = FALSE]
arg1neg_receiving_global1 <- neutrophil_receiving_global1[, grepl("NeutrophilArg1neg$", colnames(neutrophil_receiving_global1)), drop = FALSE]
signals_to_arg1pos_global1 <- rowMeans(arg1pos_receiving_global1, na.rm = TRUE)
signals_to_arg1neg_global1 <- rowMeans(arg1neg_receiving_global1, na.rm = TRUE)
global1_comparison <- data.frame(LR_Pair = names(signals_to_arg1pos_global1), Signals_to_Arg1pos = signals_to_arg1pos_global1, Signals_to_Arg1neg = signals_to_arg1neg_global1[names(signals_to_arg1pos_global1)], Difference_Pos_minus_Neg = signals_to_arg1pos_global1 - signals_to_arg1neg_global1[names(signals_to_arg1pos_global1)], Log2FC_Pos_vs_Neg = log2((signals_to_arg1pos_global1 + 1e-6)/(signals_to_arg1neg_global1[names(signals_to_arg1pos_global1)] + 1e-6)))
global1_comparison <- global1_comparison[order(global1_comparison$Difference_Pos_minus_Neg, decreasing = TRUE), ]
global1_comparison$pvalue <- NA_real_
for (lr in global1_comparison$LR_Pair) {
  arg1pos_scores <- as.numeric(arg1pos_receiving_global1[lr, ]); arg1neg_scores <- as.numeric(arg1neg_receiving_global1[lr, ])
  arg1pos_scores <- arg1pos_scores[!is.na(arg1pos_scores)]; arg1neg_scores <- arg1neg_scores[!is.na(arg1neg_scores)]
  if (length(arg1pos_scores) > 0 && length(arg1neg_scores) > 0) tryCatch({ wt <- wilcox.test(arg1pos_scores, arg1neg_scores); global1_comparison$pvalue[global1_comparison$LR_Pair == lr] <- wt$p.value }, error = function(e) {})
}
global1_comparison$padj <- p.adjust(global1_comparison$pvalue, method = "BH")
write_tsv(global1_comparison, file.path(OUTPUT_DIR, "QIN_Day1_Arg1pos_vs_Arg1neg_GLOBAL_Comparison.tsv"))
cat("Qin Day 1 Global Analysis complete. Top signals to Arg1+ neutrophils:\n"); print(head(global1_comparison, 10))

interaction_data_global1 <- cc_qin_d1_global@data$expr_l_r_log2_scale
receiver_cells_global1 <- sapply(strsplit(colnames(interaction_data_global1), "-"), `[`, 2)
neutrophil_receiving_routes_global1 <- grepl("NeutrophilArg1(pos|neg)$", receiver_cells_global1)
neutrophil_matrix_global1 <- interaction_data_global1[, neutrophil_receiving_routes_global1, drop = FALSE]
top_signals_global1 <- names(sort(rowMeans(neutrophil_matrix_global1, na.rm = TRUE), decreasing = TRUE))[1:20]
neutrophil_matrix_global1_top <- neutrophil_matrix_global1[top_signals_global1, , drop = FALSE]
p_hm_global1 <- pheatmap::pheatmap(neutrophil_matrix_global1_top, color = colorRampPalette(c("blue", "white", "yellow", "orange", "red"))(100), breaks = seq(0, 1, length.out = 101), show_rownames = TRUE, show_colnames = TRUE, treeheight_row = 0, treeheight_col = 10, cluster_rows = TRUE, cluster_cols = FALSE, fontsize = 12, angle_col = 90, fontsize_row = 12, fontsize_col = 12, main = "Global 1: Qin Day1 (Top 20 Neutrophil-Receiving Signals)")
print(p_hm_global1)
circlize::circos.clear()
unique_cell_types_global1 <- unique(gsub("-.*", "", colnames(cc_qin_d1_global@data$expr_l_r_log2_scale)))
cell_color_global1 <- data.frame(color = rainbow(length(unique_cell_types_global1)), stringsAsFactors = FALSE)
rownames(cell_color_global1) <- unique_cell_types_global1
ViewInterCircos(object = cc_qin_d1_global, font = 2, cellColor = cell_color_global1, lrColor = c("#F16B6F", "#84B1ED"), arr.type = "big.arrow", arr.length = 0.04, trackhight1 = 0.05, slot = "expr_l_r_log2_scale", linkcolor.from.sender = TRUE, linkcolor = NULL, gap.degree = 0.1, trackhight2 = 0.032, track.margin2 = c(0.01, 0.12), DIY = FALSE)
grid.text("Global 1: Qin Day 1 - Neutrophil-Receiving Cell-Cell Communication", y = 0.95, gp = gpar(fontsize = 16))
pathway.hyper.list_global1 <- list()
for(i in colnames(neutrophil_matrix_global1)) { tmp <- getHyperPathway(data = neutrophil_matrix_global1, object = cc_qin_d1_global, cella_cellb = i, Org = "Mus musculus", IS_core = TRUE); pathway.hyper.list_global1[[i]] <- tmp }
pathway.hyper.list_global1 <- pathway.hyper.list_global1[!sapply(pathway.hyper.list_global1, is.null)]
myPub.df_global1 <- getForBubble(pathway.hyper.list_global1, cella_cellb = colnames(neutrophil_matrix_global1))
myPub.df_global1$label_text <- ""; myPub.df_global1$label_text[myPub.df_global1$p.adjust < 0.05] <- format(myPub.df_global1$p.adjust[myPub.df_global1$p.adjust < 0.05], digits = 2, scientific = TRUE)
p_bubble_global1 <- plotBubble(myPub.df_global1) + geom_text_repel(aes(label = label_text), size = 3) + scale_size_continuous(limits = c(0, 100), range = c(1, 10)) + scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.05)) + ggtitle("Global 1: Qin Day 1 Neutrophil-Receiving Pathway Analysis (with p.adjust)") + theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_bubble_global1)
rm(qin_day1_global, cc_qin_d1_global, neutrophil_receiving_global1, arg1pos_receiving_global1, arg1neg_receiving_global1, neutrophil_matrix_global1, neutrophil_matrix_global1_top, pathway.hyper.list_global1, myPub.df_global1, p_hm_global1, p_bubble_global1, cell_color_global1, unique_cell_types_global1, neutrophil_cols_global1, neutrophil_sending_global1)
gc()
cat("Memory cleanup completed after Global Analysis 1\n")

# =============================================================================
# GLOBAL ANALYSIS 2: QIN DAY 3 (Arg1-labeled)
# =============================================================================
cat("\n=== GLOBAL ANALYSIS 2: QIN DAY 3 (Arg1+ vs Arg1- Comparison) ===\n")
qin_day3_global <- subset(qin_data, time == 3)
qin_day3_global$cellcall_label <- gsub("-", "", as.character(qin_day3_global$pruned_labels))
neut_idx_q3 <- grepl("Neutrophil", qin_day3_global$cellcall_label, ignore.case = TRUE)
qin_day3_global$cellcall_label[neut_idx_q3] <- paste0("Neutrophil", qin_day3_global$Arg1_status[neut_idx_q3])
cat("Qin Day 3 cells:", ncol(qin_day3_global), "\n"); print(table(qin_day3_global$cellcall_label))
qin_day3_global <- RenameCells(qin_day3_global, new.names = gsub("-", "_", colnames(qin_day3_global)))
Idents(qin_day3_global) <- as.character(qin_day3_global$cellcall_label)
cc_qin_d3_global <- CreateObject_fromSeurat(Seurat.object = qin_day3_global, slot = "counts", cell_type = "cellcall_label", data_source = "UMI", scale.factor = 10^6, Org = "Mus musculus")
cc_qin_d3_global <- TransCommuProfile(object = cc_qin_d3_global, pValueCor = 0.1, CorValue = 0.05, topTargetCor = 1, p.adjust = 0.1, use.type = "mean", probs = 0.1, method = "weighted", IS_core = TRUE, Org = "Mus musculus")
saveRDS(cc_qin_d3_global, file.path(OUTPUT_DIR, "CellCall_QIN_Day3_GLOBAL.rds"))
save.image(file.path(OUTPUT_DIR, "GlobalAnalysis2_Qin_Day3_After_TransCommu.RData"))

neutrophil_cols_global2 <- grepl("-Neutrophil", colnames(cc_qin_d3_global@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_sending_global2 <- grepl("^Neutrophil.*-", colnames(cc_qin_d3_global@data$expr_l_r_log2_scale), ignore.case = TRUE)
neutrophil_cols_global2 <- neutrophil_cols_global2 & !neutrophil_sending_global2
neutrophil_receiving_global2 <- cc_qin_d3_global@data$expr_l_r_log2_scale[, neutrophil_cols_global2, drop = FALSE]
arg1pos_receiving_global2 <- neutrophil_receiving_global2[, grepl("NeutrophilArg1pos$", colnames(neutrophil_receiving_global2)), drop = FALSE]
arg1neg_receiving_global2 <- neutrophil_receiving_global2[, grepl("NeutrophilArg1neg$", colnames(neutrophil_receiving_global2)), drop = FALSE]
signals_to_arg1pos_global2 <- rowMeans(arg1pos_receiving_global2, na.rm = TRUE)
signals_to_arg1neg_global2 <- rowMeans(arg1neg_receiving_global2, na.rm = TRUE)
global2_comparison <- data.frame(LR_Pair = names(signals_to_arg1pos_global2), Signals_to_Arg1pos = signals_to_arg1pos_global2, Signals_to_Arg1neg = signals_to_arg1neg_global2[names(signals_to_arg1pos_global2)], Difference_Pos_minus_Neg = signals_to_arg1pos_global2 - signals_to_arg1neg_global2[names(signals_to_arg1pos_global2)], Log2FC_Pos_vs_Neg = log2((signals_to_arg1pos_global2 + 1e-6)/(signals_to_arg1neg_global2[names(signals_to_arg1pos_global2)] + 1e-6)))
global2_comparison <- global2_comparison[order(global2_comparison$Difference_Pos_minus_Neg, decreasing = TRUE), ]
global2_comparison$pvalue <- NA_real_
for (lr in global2_comparison$LR_Pair) {
  arg1pos_scores <- as.numeric(arg1pos_receiving_global2[lr, ]); arg1neg_scores <- as.numeric(arg1neg_receiving_global2[lr, ])
  arg1pos_scores <- arg1pos_scores[!is.na(arg1pos_scores)]; arg1neg_scores <- arg1neg_scores[!is.na(arg1neg_scores)]
  if (length(arg1pos_scores) > 0 && length(arg1neg_scores) > 0) tryCatch({ wt <- wilcox.test(arg1pos_scores, arg1neg_scores); global2_comparison$pvalue[global2_comparison$LR_Pair == lr] <- wt$p.value }, error = function(e) {})
}
global2_comparison$padj <- p.adjust(global2_comparison$pvalue, method = "BH")
write_tsv(global2_comparison, file.path(OUTPUT_DIR, "QIN_Day3_Arg1pos_vs_Arg1neg_GLOBAL_Comparison.tsv"))
cat("Qin Day 3 Global Analysis complete. Top signals to Arg1+ neutrophils:\n"); print(head(global2_comparison, 10))

interaction_data_global2 <- cc_qin_d3_global@data$expr_l_r_log2_scale
receiver_cells_global2 <- sapply(strsplit(colnames(interaction_data_global2), "-"), `[`, 2)
neutrophil_receiving_routes_global2 <- grepl("NeutrophilArg1(pos|neg)$", receiver_cells_global2)
neutrophil_matrix_global2 <- interaction_data_global2[, neutrophil_receiving_routes_global2, drop = FALSE]
top_signals_global2 <- names(sort(rowMeans(neutrophil_matrix_global2, na.rm = TRUE), decreasing = TRUE))[1:20]
neutrophil_matrix_global2_top <- neutrophil_matrix_global2[top_signals_global2, , drop = FALSE]
p_hm_global2 <- pheatmap::pheatmap(neutrophil_matrix_global2_top, color = colorRampPalette(c("blue", "white", "yellow", "orange", "red"))(100), breaks = seq(0, 1, length.out = 101), show_rownames = TRUE, show_colnames = TRUE, treeheight_row = 0, treeheight_col = 10, cluster_rows = TRUE, cluster_cols = FALSE, fontsize = 12, angle_col = 90, fontsize_row = 12, fontsize_col = 12, main = "Global 2: Qin Day3 (Top 20 Neutrophil-Receiving Signals)")
print(p_hm_global2)
circlize::circos.clear()
unique_cell_types_global2 <- unique(gsub("-.*", "", colnames(cc_qin_d3_global@data$expr_l_r_log2_scale)))
cell_color_global2 <- data.frame(color = rainbow(length(unique_cell_types_global2)), stringsAsFactors = FALSE)
rownames(cell_color_global2) <- unique_cell_types_global2
ViewInterCircos(object = cc_qin_d3_global, font = 2, cellColor = cell_color_global2, lrColor = c("#F16B6F", "#84B1ED"), arr.type = "big.arrow", arr.length = 0.04, trackhight1 = 0.05, slot = "expr_l_r_log2_scale", linkcolor.from.sender = TRUE, linkcolor = NULL, gap.degree = 0.1, trackhight2 = 0.032, track.margin2 = c(0.01, 0.12), DIY = FALSE)
grid.text("Global 2: Qin Day 3 - Neutrophil-Receiving Cell-Cell Communication", y = 0.95, gp = gpar(fontsize = 16))
pathway.hyper.list_global2 <- list()
for(i in colnames(neutrophil_matrix_global2)) { tmp <- getHyperPathway(data = neutrophil_matrix_global2, object = cc_qin_d3_global, cella_cellb = i, Org = "Mus musculus", IS_core = TRUE); pathway.hyper.list_global2[[i]] <- tmp }
pathway.hyper.list_global2 <- pathway.hyper.list_global2[!sapply(pathway.hyper.list_global2, is.null)]
myPub.df_global2 <- getForBubble(pathway.hyper.list_global2, cella_cellb = colnames(neutrophil_matrix_global2))
myPub.df_global2$label_text <- ""; myPub.df_global2$label_text[myPub.df_global2$p.adjust < 0.05] <- format(myPub.df_global2$p.adjust[myPub.df_global2$p.adjust < 0.05], digits = 2, scientific = TRUE)
p_bubble_global2 <- plotBubble(myPub.df_global2) + geom_text_repel(aes(label = label_text), size = 3) + scale_size_continuous(limits = c(0, 100), range = c(1, 10)) + scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.05)) + ggtitle("Global 2: Qin Day 3 Neutrophil-Receiving Pathway Analysis (with p.adjust)") + theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_bubble_global2)
rm(qin_day3_global, cc_qin_d3_global, neutrophil_receiving_global2, arg1pos_receiving_global2, arg1neg_receiving_global2, neutrophil_matrix_global2, neutrophil_matrix_global2_top, pathway.hyper.list_global2, myPub.df_global2, p_hm_global2, p_bubble_global2, cell_color_global2, unique_cell_types_global2, neutrophil_cols_global2, neutrophil_sending_global2)
gc()
cat("Memory cleanup completed after Global Analysis 2\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("CELLCALL QIN-ONLY ANALYSES COMPLETE!\n")
cat("Total analyses: 4 (2 unmodified + 2 global)\n")
cat("  - Unmodified: Qin Day 1, Qin Day 3\n")
cat("  - Global (Arg1-labeled): Qin Day 1, Qin Day 3\n")
cat("Files: CellCall_QIN_Day1_UNMODIFIED.rds, CellCall_QIN_Day3_UNMODIFIED.rds,\n")
cat("       CellCall_QIN_Day1_GLOBAL.rds, CellCall_QIN_Day3_GLOBAL.rds,\n")
cat("       QIN_Day1_Arg1pos_vs_Arg1neg_GLOBAL_Comparison.tsv, QIN_Day3_Arg1pos_vs_Arg1neg_GLOBAL_Comparison.tsv\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

save.image(file.path(OUTPUT_DIR, "CELLCALL_QIN_ONLY_Final_Environment.RData"))
