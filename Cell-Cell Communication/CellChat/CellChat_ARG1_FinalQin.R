# CellChat ARG1 Neutrophil Analysis (Qin)
# Research question: Which incoming ligand–receptor signals differ between Arg1+ and Arg1- neutrophils?
#   Arg1+ as candidate anti-inflammatory–biased; Arg1- as candidate pro-inflammatory–biased (hypothesis).
# Note: Correlation from communication inference, not causation.
#
# RUN ORDER: run CellChat_ARG1_Final.R first (Lee/Wang RDS), then this script for Qin CellChat objects.
# Required in getwd(): LeeDat.rds, WangDat.rds, QinDat.rds
# Rerun from: line 16 (library(CellChat)) for full pipeline including Lee/Wang baselines + Qin; or line 419
#   (# ANALYSIS 4: QIN DAY 1 BASELINE) only if LeeDat/WangDat RDS outputs already exist and you only refresh Qin sections
#   (not recommended unless you know the script’s dependencies).
#
# ---- Working directory ----
CELLCOMM_PROJECT_ROOT <- Sys.getenv("CELLCOMM_PROJECT_ROOT", unset = "")
if (nzchar(CELLCOMM_PROJECT_ROOT) && dir.exists(CELLCOMM_PROJECT_ROOT)) setwd(CELLCOMM_PROJECT_ROOT)

library(CellChat)
library(NMF)
library(patchwork)
library(readr)
library(Seurat)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(ggalluvial)

# =============================================================================
# SECTION 1: SHARED VARIABLES AND FUNCTIONS
# =============================================================================

# Global variables used across all analyses
neutrophil_base_label <- "Neutrophil"
neutrophil_aliases <- c("Neutrophil", "Neutrophils")
gene_candidates <- c("Arg1", "ARG1")
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
lee_neutrophil_label <- "Neutrophil"
wang_neutrophil_label <- "Neutrophil"
qin_neutrophil_label <- "Neutrophil"

# CellChat computeCommunProb: nboot=100 for stable p-values; raw.use=FALSE uses normalized data after merge
CELLCHAT_NBOOT <- 100L
CELLCHAT_MARGIN_CIRCLE <- 0.35
N_PATTERNS_INCOMING <- 3L

# Seurat v5: layer="data"; v4: slot="data"
cellchat_norm_matrix <- function(seu, genes, cells) {
  if (utils::packageVersion("Seurat") >= "5.0.0") {
    Seurat::GetAssayData(seu, layer = "data")[genes, cells, drop = FALSE]
  } else {
    Seurat::GetAssayData(seu, slot = "data")[genes, cells, drop = FALSE]
  }
}

print_signaling_role_network <- function(cellchat_obj, column_title) {
  ht <- netAnalysis_signalingRole_network(cellchat_obj, width = 8, height = 2.5, font.size = 10)
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    ComplexHeatmap::draw(ht, column_title = column_title)
  } else {
    print(ht)
  }
}

netVisual_heatmap_tnf_if_present <- function(cellchat_obj, title.name) {
  paths <- tryCatch(cellchat_obj@netP$pathways, error = function(e) character(0))
  if (length(paths) > 0 && "TNF" %in% paths) {
    netVisual_heatmap(cellchat_obj, signaling = "TNF", title.name = title.name)
  } else {
    message("Skipping TNF heatmap (TNF not in @netP$pathways): ", title.name)
  }
}

# Fixed scaling parameters for consistent comparisons
GLOBAL_MIN_COUNT <- NA
GLOBAL_MAX_COUNT <- NA
GLOBAL_MIN_PROB <- NA
GLOBAL_MAX_PROB <- NA

# Function to create standardized circle plots with auto-scaling
create_standard_circle_plot <- function(net_count, vertex_weight, targets_use, title_name, 
                                        vertex_label_cex = 2, margin = CELLCHAT_MARGIN_CIRCLE) {
  netVisual_circle(net_count, 
                   vertex.weight = vertex_weight, 
                   weight.scale = TRUE,  # Auto-scaling for better visibility
                   label.edge = FALSE, 
                   targets.use = targets_use, 
                   vertex.label.cex = vertex_label_cex, 
                   margin = margin, 
                   title.name = title_name)
}

# Function to create standardized bubble plots with fixed scaling
create_standard_bubble_plot <- function(cellchat_obj, targets_use, title_name) {
  netVisual_bubble(cellchat_obj, 
                   targets.use = targets_use, 
                   remove.isolate = FALSE, 
                   title.name = title_name)
}

# =============================================================================
# ANALYSIS 1: LEE DAY 1 BASELINE - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 1: LEE DAY 1 BASELINE")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Lee Day 1 Baseline
print("Loading and preparing Lee Day 1 data...")
leeDat <- readRDS(file = "LeeDat.rds")
wangDat <- readRDS(file = "WangDat.rds")

# Process labels and merge datasets
leeDat$pruned_labels <- as.character(leeDat$celltype)
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
leeDat$dataset <- "Lee"
wangDat$dataset <- "Wang"

seuObj.integrated <- merge(leeDat, wangDat)
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Lee Day 1
oneDPI <- subset(seuObj.integrated, time == "1")
DefaultAssay(oneDPI) <- "RNA"

# Standardize neutrophil labels
oneDPI$pruned_labels <- as.character(oneDPI$pruned_labels)
oneDPI$pruned_labels[oneDPI$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
oneDPI$pruned_labels[is.na(oneDPI$pruned_labels)] <- "NA"

print("✅ Lee Day 1 data preparation completed")
print("ANALYSIS 1: Lee Day 1 Baseline - Signals TO Original Neutrophils")
print(paste(rep("=", 80), collapse=""))

cellchatOne_baseline <- createCellChat(object = oneDPI, group.by = "celltype")
CellChatDB <- CellChatDB.mouse
cellchatOne_baseline@DB <- CellChatDB

cellchatOne_baseline <- subsetData(cellchatOne_baseline)
cellchatOne_baseline <- identifyOverExpressedGenes(cellchatOne_baseline)
cellchatOne_baseline <- identifyOverExpressedInteractions(cellchatOne_baseline)
cellchatOne_baseline <- computeCommunProb(cellchatOne_baseline, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatOne_baseline <- filterCommunication(cellchatOne_baseline, min.cells = 10)
cellchatOne_baseline <- computeCommunProbPathway(cellchatOne_baseline)
cellchatOne_baseline <- aggregateNet(cellchatOne_baseline)

groupSizeOne_baseline <- as.numeric(table(cellchatOne_baseline@idents)[rownames(cellchatOne_baseline@net$count)])
groupSizeOne_baseline[is.na(groupSizeOne_baseline)] <- 0

cellchatOne_baseline <- netAnalysis_computeCentrality(cellchatOne_baseline, slot.name = "netP")

selectK(cellchatOne_baseline, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatOne_baseline <- identifyCommunicationPatterns(cellchatOne_baseline, pattern = "incoming", k = nPatterns)

# Save RDS file immediately after processing
saveRDS(cellchatOne_baseline, "cellchatOne_baseline_complete.rds")
print("✅ Saved: cellchatOne_baseline_complete.rds (Lee Day 1 Baseline - Complete)")

print("Generating visualizations for Lee Day 1 Baseline...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
lee_neutrophil_label <- "Neutrophil"
cellchatOne_baseline <- readRDS("cellchatOne_baseline_complete.rds")
groupSizeOne_baseline <- as.numeric(table(cellchatOne_baseline@idents)[rownames(cellchatOne_baseline@net$count)])
groupSizeOne_baseline[is.na(groupSizeOne_baseline)] <- 0

pdf("1dpiLee_Baseline_Target_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatOne_baseline@net$count, vertex.weight = groupSizeOne_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = lee_neutrophil_label, vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Lee Day 1 Baseline - Signals TO Neutrophils",
)
dev.off()

tiff("1dpiLee_Baseline_Target_Neutrophils.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatOne_baseline@net$count, vertex.weight = groupSizeOne_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = lee_neutrophil_label, vertex.label.cex = 1.8, margin = 0, title.name = "Lee Day 1 Baseline - Signals TO Neutrophils",
)
dev.off()

png("Neu1dpiBaseline_TargetBubbleLee.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatOne_baseline, targets.use = lee_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Lee Day 1 Baseline - Signals TO Neutrophils")
dev.off()

print("Inline plots for Lee Day 1 Baseline:")
netVisual_circle(cellchatOne_baseline@net$count, vertex.weight = groupSizeOne_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = lee_neutrophil_label, title.name = "Lee Day 1 Baseline - Signals TO Neutrophils")
netVisual_bubble(cellchatOne_baseline, targets.use = lee_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Lee Day 1 Baseline - Signals TO Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatOne_baseline, measure = "weight", targets.use = lee_neutrophil_label, 
                    title.name = "Lee Day 1 Baseline - Signals TO Neutrophils")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Lee Day 1 Baseline (insufficient color range)\n")
})
netVisual_heatmap(cellchatOne_baseline, signaling = "TNF", title.name = "Lee Day 1 Baseline - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatOne_baseline) + ggtitle("Lee Day 1 Baseline - Communication Strength"))
print(netAnalysis_signalingRole_network(cellchatOne_baseline, width = 8, height = 2.5, font.size = 10) + ggtitle("Lee Day 1 Baseline - Signaling Role Network"))

ht2_day1_baseline <- netAnalysis_signalingRole_heatmap(cellchatOne_baseline, pattern = "incoming", font.size = 5)
print(ht2_day1_baseline + plot_annotation(title = "Lee Day 1 Baseline - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatOne_baseline, pattern = "incoming") + ggtitle("Lee Day 1 Baseline - Communication Patterns"))
print(netAnalysis_dot(cellchatOne_baseline, pattern = "incoming") + ggtitle("Lee Day 1 Baseline - Pattern Strength"))

save.image("environment_analysis1_baseline.RData")
print("✅ Saved: environment_analysis1_baseline.RData (Analysis 1 Environment)")

# =============================================================================
# ANALYSIS 2: LEE DAY 3 BASELINE - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 2: LEE DAY 3 BASELINE")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Lee Day 3 Baseline
print("Loading and preparing Lee Day 3 data...")
leeDat <- readRDS(file = "LeeDat.rds")
wangDat <- readRDS(file = "WangDat.rds")

# Process labels and merge datasets
leeDat$pruned_labels <- as.character(leeDat$celltype)
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
leeDat$dataset <- "Lee"
wangDat$dataset <- "Wang"

seuObj.integrated <- merge(leeDat, wangDat)
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Lee Day 3
threeDPILee <- subset(seuObj.integrated, time == "3" & dataset == "Lee")
DefaultAssay(threeDPILee) <- "RNA"

# Standardize neutrophil labels
threeDPILee$pruned_labels <- as.character(threeDPILee$pruned_labels)
threeDPILee$pruned_labels[threeDPILee$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
threeDPILee$pruned_labels[is.na(threeDPILee$pruned_labels)] <- "NA"

print("✅ Lee Day 3 data preparation completed")
print("ANALYSIS 2: Lee Day 3 Baseline - Signals TO Original Neutrophils")
print(paste(rep("=", 80), collapse=""))

cellchatThreeLee_baseline <- createCellChat(object = threeDPILee, group.by = "celltype")
CellChatDB <- CellChatDB.mouse
cellchatThreeLee_baseline@DB <- CellChatDB

cellchatThreeLee_baseline <- subsetData(cellchatThreeLee_baseline)
cellchatThreeLee_baseline <- identifyOverExpressedGenes(cellchatThreeLee_baseline)
cellchatThreeLee_baseline <- identifyOverExpressedInteractions(cellchatThreeLee_baseline)
cellchatThreeLee_baseline <- computeCommunProb(cellchatThreeLee_baseline, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatThreeLee_baseline <- filterCommunication(cellchatThreeLee_baseline, min.cells = 10)
cellchatThreeLee_baseline <- computeCommunProbPathway(cellchatThreeLee_baseline)
cellchatThreeLee_baseline <- aggregateNet(cellchatThreeLee_baseline)

groupSizeThreeLee_baseline <- as.numeric(table(cellchatThreeLee_baseline@idents)[rownames(cellchatThreeLee_baseline@net$count)])
groupSizeThreeLee_baseline[is.na(groupSizeThreeLee_baseline)] <- 0

cellchatThreeLee_baseline <- netAnalysis_computeCentrality(cellchatThreeLee_baseline, slot.name = "netP")

cellchatThreeLee_baseline <- identifyCommunicationPatterns(cellchatThreeLee_baseline, pattern = "incoming", k = nPatterns)

# Save RDS file immediately after processing
saveRDS(cellchatThreeLee_baseline, "cellchatThreeLee_baseline_complete.rds")
print("✅ Saved: cellchatThreeLee_baseline_complete.rds (Lee Day 3 Baseline - Complete)")

print("Generating visualizations for Lee Day 3 Baseline...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
lee_neutrophil_label <- "Neutrophil"
cellchatThreeLee_baseline <- readRDS("cellchatThreeLee_baseline_complete.rds")
groupSizeThreeLee_baseline <- as.numeric(table(cellchatThreeLee_baseline@idents)[rownames(cellchatThreeLee_baseline@net$count)])
groupSizeThreeLee_baseline[is.na(groupSizeThreeLee_baseline)] <- 0

tiff("3dpiLee_Baseline_Target_Neutrophils.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeLee_baseline@net$count, vertex.weight = groupSizeThreeLee_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = lee_neutrophil_label, vertex.label.cex = 1.8, margin = 0, title.name = "Lee Day 3 Baseline - Signals TO Neutrophils")
dev.off()

png("Neu3dpiBaseline_TargetCircleLee.png", width = 6, height = 6, units = "in", res = 1200)
netVisual_circle(cellchatThreeLee_baseline@net$count, vertex.weight = groupSizeThreeLee_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = lee_neutrophil_label, title.name = "Lee Day 3 Baseline - Signals TO Neutrophils")
dev.off()

print("Inline plots for Lee Day 3 Baseline:")
netVisual_circle(cellchatThreeLee_baseline@net$count, vertex.weight = groupSizeThreeLee_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = lee_neutrophil_label, title.name = "Lee Day 3 Baseline - Signals TO Neutrophils")
netVisual_bubble(cellchatThreeLee_baseline, targets.use = lee_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Lee Day 3 Baseline - Signals TO Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatThreeLee_baseline, measure = "weight", targets.use = lee_neutrophil_label, 
                    title.name = "Lee Day 3 Baseline - Signals TO Neutrophils")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Lee Day 3 Baseline (insufficient color range)\n")
})
netVisual_heatmap_tnf_if_present(cellchatThreeLee_baseline, "Lee Day 3 Baseline - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatThreeLee_baseline) + ggtitle("Lee Day 3 Baseline - Communication Strength"))
print_signaling_role_network(cellchatThreeLee_baseline, "Lee Day 3 Baseline - Signaling Role Network")

ht2_lee_day3_baseline <- netAnalysis_signalingRole_heatmap(cellchatThreeLee_baseline, pattern = "incoming", font.size = 5)
print(ht2_lee_day3_baseline + plot_annotation(title = "Lee Day 3 Baseline - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatThreeLee_baseline, pattern = "incoming") + ggtitle("Lee Day 3 Baseline - Communication Patterns"))
print(netAnalysis_dot(cellchatThreeLee_baseline, pattern = "incoming") + ggtitle("Lee Day 3 Baseline - Pattern Strength"))

save.image("environment_analysis2_baseline.RData")
print("✅ Saved: environment_analysis2_baseline.RData (Analysis 2 Environment)")

# =============================================================================
# ANALYSIS 3: WANG DAY 3 BASELINE - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 3: WANG DAY 3 BASELINE")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Wang Day 3 Baseline
print("Loading and preparing Wang Day 3 data...")
leeDat <- readRDS(file = "LeeDat.rds")
wangDat <- readRDS(file = "WangDat.rds")

# Process labels and merge datasets
leeDat$pruned_labels <- as.character(leeDat$celltype)
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
leeDat$dataset <- "Lee"
wangDat$dataset <- "Wang"

seuObj.integrated <- merge(leeDat, wangDat)
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Wang Day 3
threeDPIWang <- subset(seuObj.integrated, time == "3" & dataset == "Wang")
DefaultAssay(threeDPIWang) <- "RNA"

# Standardize neutrophil labels
threeDPIWang$pruned_labels <- as.character(threeDPIWang$pruned_labels)
threeDPIWang$pruned_labels[threeDPIWang$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
threeDPIWang$pruned_labels[is.na(threeDPIWang$pruned_labels)] <- "NA"

print("✅ Wang Day 3 data preparation completed")
print("ANALYSIS 3: Wang Day 3 Baseline - Signals TO Original Neutrophils")
print(paste(rep("=", 80), collapse=""))

cellchatThreeWang_baseline <- createCellChat(object = threeDPIWang, group.by = "pruned_labels")
CellChatDB <- CellChatDB.mouse
cellchatThreeWang_baseline@DB <- CellChatDB

cellchatThreeWang_baseline <- subsetData(cellchatThreeWang_baseline)
cellchatThreeWang_baseline <- identifyOverExpressedGenes(cellchatThreeWang_baseline)
cellchatThreeWang_baseline <- identifyOverExpressedInteractions(cellchatThreeWang_baseline)
cellchatThreeWang_baseline <- computeCommunProb(cellchatThreeWang_baseline, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatThreeWang_baseline <- filterCommunication(cellchatThreeWang_baseline, min.cells = 10)
cellchatThreeWang_baseline <- computeCommunProbPathway(cellchatThreeWang_baseline)
cellchatThreeWang_baseline <- aggregateNet(cellchatThreeWang_baseline)

groupSizeThreeWang_baseline <- as.numeric(table(cellchatThreeWang_baseline@idents)[rownames(cellchatThreeWang_baseline@net$count)])
groupSizeThreeWang_baseline[is.na(groupSizeThreeWang_baseline)] <- 0

cellchatThreeWang_baseline <- netAnalysis_computeCentrality(cellchatThreeWang_baseline, slot.name = "netP")

cellchatThreeWang_baseline <- identifyCommunicationPatterns(cellchatThreeWang_baseline, pattern = "incoming", k = nPatterns)

# Save RDS file immediately after processing
saveRDS(cellchatThreeWang_baseline, "cellchatThreeWang_baseline_complete.rds")
print("✅ Saved: cellchatThreeWang_baseline_complete.rds (Wang Day 3 Baseline - Complete)")

print("Generating visualizations for Wang Day 3 Baseline...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
wang_neutrophil_label <- "Neutrophil"
cellchatThreeWang_baseline <- readRDS("cellchatThreeWang_baseline_complete.rds")
groupSizeThreeWang_baseline <- as.numeric(table(cellchatThreeWang_baseline@idents)[rownames(cellchatThreeWang_baseline@net$count)])
groupSizeThreeWang_baseline[is.na(groupSizeThreeWang_baseline)] <- 0

tiff("3dpiWang_Baseline_Target_Neutrophils.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeWang_baseline@net$count, vertex.weight = groupSizeThreeWang_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = wang_neutrophil_label, vertex.label.cex = 1.8, margin = 0, title.name = "Wang Day 3 Baseline - Signals TO Neutrophils")
dev.off()

png("Neu3dpiBaseline_TargetBubbleWang.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatThreeWang_baseline, targets.use = wang_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Wang Day 3 Baseline - Signals TO Neutrophils")
dev.off()

png("Neu3dpiBaseline_TargetCircleWang.png", width = 6, height = 6, units = "in", res = 1200)
netVisual_circle(cellchatThreeWang_baseline@net$count, vertex.weight = groupSizeThreeWang_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = wang_neutrophil_label, title.name = "Wang Day 3 Baseline - Signals TO Neutrophils")
dev.off()

print("Inline plots for Wang Day 3 Baseline:")
netVisual_circle(cellchatThreeWang_baseline@net$count, vertex.weight = groupSizeThreeWang_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = wang_neutrophil_label, title.name = "Wang Day 3 Baseline - Signals TO Neutrophils")
netVisual_bubble(cellchatThreeWang_baseline, targets.use = wang_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Wang Day 3 Baseline - Signals TO Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatThreeWang_baseline, measure = "weight", targets.use = wang_neutrophil_label, 
                    title.name = "Wang Day 3 Baseline - Signals TO Neutrophils")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Wang Day 3 Baseline (insufficient color range)\n")
})
netVisual_heatmap(cellchatThreeWang_baseline, signaling = "TNF", title.name = "Wang Day 3 Baseline - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatThreeWang_baseline) + ggtitle("Wang Day 3 Baseline - Communication Strength"))
print(netAnalysis_signalingRole_network(cellchatThreeWang_baseline, width = 8, height = 2.5, font.size = 10) + ggtitle("Wang Day 3 Baseline - Signaling Role Network"))

ht2_wang_day3_baseline <- netAnalysis_signalingRole_heatmap(cellchatThreeWang_baseline, pattern = "incoming", font.size = 5)
print(ht2_wang_day3_baseline + plot_annotation(title = "Wang Day 3 Baseline - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatThreeWang_baseline, pattern = "incoming") + ggtitle("Wang Day 3 Baseline - Communication Patterns"))
print(netAnalysis_dot(cellchatThreeWang_baseline, pattern = "incoming") + ggtitle("Wang Day 3 Baseline - Pattern Strength"))

save.image("environment_analysis3_baseline.RData")
print("✅ Saved: environment_analysis3_baseline.RData (Analysis 3 Environment)")

# =============================================================================
# ANALYSIS 4: QIN DAY 1 BASELINE - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 4: QIN DAY 1 BASELINE")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Qin Day 1 Baseline
print("Loading and preparing Qin Day 1 data...")
qinDat <- readRDS(file = "QinDat.rds")
wangDat <- readRDS(file = "WangDat.rds")
leeDat <- readRDS(file = "LeeDat.rds")

# Process labels and merge datasets
qinDat$pruned_labels <- as.character(qinDat$celltype)
qinDat$dataset <- "Qin"
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
wangDat$dataset <- "Wang"
leeDat$pruned_labels <- as.character(leeDat$celltype)
leeDat$dataset <- "Lee"

seuObj.integrated <- merge(leeDat, y = c(wangDat, qinDat))
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Qin Day 1
oneDPIQin <- subset(seuObj.integrated, time == "1" & dataset == "Qin")
DefaultAssay(oneDPIQin) <- "RNA"

# Standardize neutrophil labels
oneDPIQin$pruned_labels <- as.character(oneDPIQin$pruned_labels)
oneDPIQin$pruned_labels[oneDPIQin$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
oneDPIQin$pruned_labels[is.na(oneDPIQin$pruned_labels)] <- "NA"

print("✅ Qin Day 1 data preparation completed")
print("ANALYSIS 4: Qin Day 1 Baseline - Signals TO Original Neutrophils")
print(paste(rep("=", 80), collapse=""))

cellchatOneQin_baseline <- createCellChat(object = oneDPIQin, group.by = "celltype")
CellChatDB <- CellChatDB.mouse
cellchatOneQin_baseline@DB <- CellChatDB

cellchatOneQin_baseline <- subsetData(cellchatOneQin_baseline)
cellchatOneQin_baseline <- identifyOverExpressedGenes(cellchatOneQin_baseline)
cellchatOneQin_baseline <- identifyOverExpressedInteractions(cellchatOneQin_baseline)
cellchatOneQin_baseline <- computeCommunProb(cellchatOneQin_baseline, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatOneQin_baseline <- filterCommunication(cellchatOneQin_baseline, min.cells = 10)
cellchatOneQin_baseline <- computeCommunProbPathway(cellchatOneQin_baseline)
cellchatOneQin_baseline <- aggregateNet(cellchatOneQin_baseline)

groupSizeOneQin_baseline <- as.numeric(table(cellchatOneQin_baseline@idents)[rownames(cellchatOneQin_baseline@net$count)])
groupSizeOneQin_baseline[is.na(groupSizeOneQin_baseline)] <- 0

cellchatOneQin_baseline <- netAnalysis_computeCentrality(cellchatOneQin_baseline, slot.name = "netP")

selectK(cellchatOneQin_baseline, pattern = "incoming")
nPatterns <- 3
cellchatOneQin_baseline <- identifyCommunicationPatterns(cellchatOneQin_baseline, pattern = "incoming", k = nPatterns)

# Save RDS file immediately after processing
saveRDS(cellchatOneQin_baseline, "cellchatOneQin_baseline_complete.rds")
print("✅ Saved: cellchatOneQin_baseline_complete.rds (Qin Day 1 Baseline - Complete)")

print("Generating visualizations for Qin Day 1 Baseline...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
qin_neutrophil_label <- "Neutrophil"
cellchatOneQin_baseline <- readRDS("cellchatOneQin_baseline_complete.rds")
groupSizeOneQin_baseline <- as.numeric(table(cellchatOneQin_baseline@idents)[rownames(cellchatOneQin_baseline@net$count)])
groupSizeOneQin_baseline[is.na(groupSizeOneQin_baseline)] <- 0

pdf("1dpiQin_Baseline_Target_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatOneQin_baseline@net$count, vertex.weight = groupSizeOneQin_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = qin_neutrophil_label, vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Qin Day 1 Baseline - Signals TO Neutrophils")
dev.off()

tiff("1dpiQin_Baseline_Target_Neutrophils.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatOneQin_baseline@net$count, vertex.weight = groupSizeOneQin_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = qin_neutrophil_label, vertex.label.cex = 1.8, margin = 0, title.name = "Qin Day 1 Baseline - Signals TO Neutrophils")
dev.off()

png("Neu1dpiBaseline_TargetBubbleQin.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatOneQin_baseline, targets.use = qin_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Qin Day 1 Baseline - Signals TO Neutrophils")
dev.off()

print("Inline plots for Qin Day 1 Baseline:")
netVisual_circle(cellchatOneQin_baseline@net$count, vertex.weight = groupSizeOneQin_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = qin_neutrophil_label, title.name = "Qin Day 1 Baseline - Signals TO Neutrophils")
netVisual_bubble(cellchatOneQin_baseline, targets.use = qin_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Qin Day 1 Baseline - Signals TO Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatOneQin_baseline, measure = "weight", targets.use = qin_neutrophil_label, 
                    title.name = "Qin Day 1 Baseline - Signals TO Neutrophils")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Qin Day 1 Baseline (insufficient color range)\n")
})
netVisual_heatmap(cellchatOneQin_baseline, signaling = "TNF", title.name = "Qin Day 1 Baseline - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatOneQin_baseline) + ggtitle("Qin Day 1 Baseline - Communication Strength"))
print(netAnalysis_signalingRole_network(cellchatOneQin_baseline, width = 8, height = 2.5, font.size = 10) + ggtitle("Qin Day 1 Baseline - Signaling Role Network"))

ht2_qin_day1_baseline <- netAnalysis_signalingRole_heatmap(cellchatOneQin_baseline, pattern = "incoming", font.size = 5)
print(ht2_qin_day1_baseline + plot_annotation(title = "Qin Day 1 Baseline - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatOneQin_baseline, pattern = "incoming") + ggtitle("Qin Day 1 Baseline - Communication Patterns"))
print(netAnalysis_dot(cellchatOneQin_baseline, pattern = "incoming") + ggtitle("Qin Day 1 Baseline - Pattern Strength"))

save.image("environment_analysis4_baseline.RData")
print("✅ Saved: environment_analysis4_baseline.RData (Analysis 4 Environment)")

# =============================================================================
# ANALYSIS 5: QIN DAY 3 BASELINE - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 5: QIN DAY 3 BASELINE")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Qin Day 3 Baseline
print("Loading and preparing Qin Day 3 data...")
qinDat <- readRDS(file = "QinDat.rds")
wangDat <- readRDS(file = "WangDat.rds")
leeDat <- readRDS(file = "LeeDat.rds")

# Process labels and merge datasets
qinDat$pruned_labels <- as.character(qinDat$celltype)
qinDat$dataset <- "Qin"
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
wangDat$dataset <- "Wang"
leeDat$pruned_labels <- as.character(leeDat$celltype)
leeDat$dataset <- "Lee"

seuObj.integrated <- merge(leeDat, y = c(wangDat, qinDat))
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Qin Day 3
threeDPIQin <- subset(seuObj.integrated, time == "3" & dataset == "Qin")
DefaultAssay(threeDPIQin) <- "RNA"

# Standardize neutrophil labels
threeDPIQin$pruned_labels <- as.character(threeDPIQin$pruned_labels)
threeDPIQin$pruned_labels[threeDPIQin$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
threeDPIQin$pruned_labels[is.na(threeDPIQin$pruned_labels)] <- "NA"

print("✅ Qin Day 3 data preparation completed")
print("ANALYSIS 5: Qin Day 3 Baseline - Signals TO Original Neutrophils")
print(paste(rep("=", 80), collapse=""))

cellchatThreeQin_baseline <- createCellChat(object = threeDPIQin, group.by = "celltype")
CellChatDB <- CellChatDB.mouse
cellchatThreeQin_baseline@DB <- CellChatDB

cellchatThreeQin_baseline <- subsetData(cellchatThreeQin_baseline)
cellchatThreeQin_baseline <- identifyOverExpressedGenes(cellchatThreeQin_baseline)
cellchatThreeQin_baseline <- identifyOverExpressedInteractions(cellchatThreeQin_baseline)
cellchatThreeQin_baseline <- computeCommunProb(cellchatThreeQin_baseline, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatThreeQin_baseline <- filterCommunication(cellchatThreeQin_baseline, min.cells = 10)
cellchatThreeQin_baseline <- computeCommunProbPathway(cellchatThreeQin_baseline)
cellchatThreeQin_baseline <- aggregateNet(cellchatThreeQin_baseline)

groupSizeThreeQin_baseline <- as.numeric(table(cellchatThreeQin_baseline@idents)[rownames(cellchatThreeQin_baseline@net$count)])
groupSizeThreeQin_baseline[is.na(groupSizeThreeQin_baseline)] <- 0

cellchatThreeQin_baseline <- netAnalysis_computeCentrality(cellchatThreeQin_baseline, slot.name = "netP")

selectK(cellchatThreeQin_baseline, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatThreeQin_baseline <- identifyCommunicationPatterns(cellchatThreeQin_baseline, pattern = "incoming", k = nPatterns)

# Save RDS file immediately after processing
saveRDS(cellchatThreeQin_baseline, "cellchatThreeQin_baseline_complete.rds")
print("✅ Saved: cellchatThreeQin_baseline_complete.rds (Qin Day 3 Baseline - Complete)")

print("Generating visualizations for Qin Day 3 Baseline...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
qin_neutrophil_label <- "Neutrophil"
cellchatThreeQin_baseline <- readRDS("cellchatThreeQin_baseline_complete.rds")
groupSizeThreeQin_baseline <- as.numeric(table(cellchatThreeQin_baseline@idents)[rownames(cellchatThreeQin_baseline@net$count)])
groupSizeThreeQin_baseline[is.na(groupSizeThreeQin_baseline)] <- 0

tiff("3dpiQin_Baseline_Target_Neutrophils.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeQin_baseline@net$count, vertex.weight = groupSizeThreeQin_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = qin_neutrophil_label, vertex.label.cex = 1.8, margin = 0, title.name = "Qin Day 3 Baseline - Signals TO Neutrophils")
dev.off()

png("Neu3dpiBaseline_TargetCircleQin.png", width = 6, height = 6, units = "in", res = 1200)
netVisual_circle(cellchatThreeQin_baseline@net$count, vertex.weight = groupSizeThreeQin_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = qin_neutrophil_label, title.name = "Qin Day 3 Baseline - Signals TO Neutrophils")
dev.off()

print("Inline plots for Qin Day 3 Baseline:")
netVisual_circle(cellchatThreeQin_baseline@net$count, vertex.weight = groupSizeThreeQin_baseline, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = qin_neutrophil_label, title.name = "Qin Day 3 Baseline - Signals TO Neutrophils")
netVisual_bubble(cellchatThreeQin_baseline, targets.use = qin_neutrophil_label, remove.isolate = FALSE, 
                 title.name = "Qin Day 3 Baseline - Signals TO Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatThreeQin_baseline, measure = "weight", targets.use = qin_neutrophil_label, 
                    title.name = "Qin Day 3 Baseline - Signals TO Neutrophils")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Qin Day 3 Baseline (insufficient color range)\n")
})
netVisual_heatmap(cellchatThreeQin_baseline, signaling = "TNF", title.name = "Qin Day 3 Baseline - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatThreeQin_baseline) + ggtitle("Qin Day 3 Baseline - Communication Strength"))
print(netAnalysis_signalingRole_network(cellchatThreeQin_baseline, width = 8, height = 2.5, font.size = 10) + ggtitle("Qin Day 3 Baseline - Signaling Role Network"))

ht2_qin_day3_baseline <- netAnalysis_signalingRole_heatmap(cellchatThreeQin_baseline, pattern = "incoming", font.size = 5)
print(ht2_qin_day3_baseline + plot_annotation(title = "Qin Day 3 Baseline - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatThreeQin_baseline, pattern = "incoming") + ggtitle("Qin Day 3 Baseline - Communication Patterns"))
print(netAnalysis_dot(cellchatThreeQin_baseline, pattern = "incoming") + ggtitle("Qin Day 3 Baseline - Pattern Strength"))

save.image("environment_analysis5_baseline.RData")
print("✅ Saved: environment_analysis5_baseline.RData (Analysis 5 Environment)")

# =============================================================================
# ANALYSIS 6: LEE DAY 1 ARG1-SPLIT - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 6: LEE DAY 1 ARG1-SPLIT")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Lee Day 1 ARG1-split
print("Loading and preparing Lee Day 1 data with ARG1 annotation...")
leeDat <- readRDS(file = "LeeDat.rds")
wangDat <- readRDS(file = "WangDat.rds")

# Process labels and merge datasets
leeDat$pruned_labels <- as.character(leeDat$celltype)
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
leeDat$dataset <- "Lee"
wangDat$dataset <- "Wang"

seuObj.integrated <- merge(leeDat, wangDat)
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Lee Day 1
oneDPI <- subset(seuObj.integrated, time == "1")
DefaultAssay(oneDPI) <- "RNA"

# Standardize neutrophil labels
oneDPI$pruned_labels <- as.character(oneDPI$pruned_labels)
oneDPI$pruned_labels[oneDPI$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
oneDPI$pruned_labels[is.na(oneDPI$pruned_labels)] <- "NA"

# ARG1 annotation for neutrophils
print("Annotating ARG1 status for neutrophils...")
neut_one <- which(oneDPI$pruned_labels == neutrophil_base_label)
gene_one <- head(gene_candidates[gene_candidates %in% rownames(oneDPI)], 1)
gene_one_value <- c(gene_one, NA_character_)[1]

expr_one <- Seurat::GetAssayData(oneDPI, slot = "data")[gene_one, neut_one, drop = FALSE]
expr_one_vec <- as.numeric(expr_one)
expr_one_pos <- expr_one_vec[expr_one_vec > 0]
threshold_one <- stats::median(expr_one_pos)
threshold_one[is.na(threshold_one)] <- 0

indicator_one <- expr_one_vec >= threshold_one
indicator_one <- indicator_one & rep(length(expr_one_pos) > 0, length(indicator_one))

pos_one <- neut_one[indicator_one]
neg_one <- setdiff(neut_one, pos_one)

oneDPI$arg1_status <- rep("Other", ncol(oneDPI))
oneDPI$arg1_status[neut_one] <- "ARG1neg"
oneDPI$arg1_status[pos_one] <- "ARG1pos"

oneDPI$cellchat_labels <- oneDPI$pruned_labels
oneDPI$cellchat_labels[pos_one] <- paste0(oneDPI$cellchat_labels[pos_one], "Arg1pos")
oneDPI$cellchat_labels[neg_one] <- paste0(oneDPI$cellchat_labels[neg_one], "Arg1neg")

print("✅ Lee Day 1 ARG1-split data preparation completed")
print("ANALYSIS 6: Lee Day 1 ARG1-split - Signals TO ARG1 Neutrophil States")
print(paste(rep("=", 80), collapse=""))

cellchatOne <- createCellChat(object = oneDPI, group.by = "cellchat_labels")
CellChatDB <- CellChatDB.mouse
cellchatOne@DB <- CellChatDB

cellchatOne <- subsetData(cellchatOne)
cellchatOne <- identifyOverExpressedGenes(cellchatOne)
cellchatOne <- identifyOverExpressedInteractions(cellchatOne)
cellchatOne <- computeCommunProb(cellchatOne, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatOne <- filterCommunication(cellchatOne, min.cells = 10)
cellchatOne <- computeCommunProbPathway(cellchatOne)
cellchatOne <- aggregateNet(cellchatOne)

cellchatOne@idents <- factor(cellchatOne@idents, levels = unique(c(levels(cellchatOne@idents), neutrophil_states)))
groupSizeOne <- as.numeric(table(cellchatOne@idents)[rownames(cellchatOne@net$count)])
groupSizeOne[is.na(groupSizeOne)] <- 0

cellchatOne <- netAnalysis_computeCentrality(cellchatOne, slot.name = "netP")

selectK(cellchatOne, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatOne <- identifyCommunicationPatterns(cellchatOne, pattern = "incoming", k = nPatterns)

comm_day1 <- subsetCommunication(cellchatOne)
comm_day1_pos <- dplyr::select(dplyr::filter(comm_day1, target == neutrophil_states[1]), 
                               source, ligand, receptor, pathway_name, interaction_name_2, prob.pos = prob)
comm_day1_neg <- dplyr::select(dplyr::filter(comm_day1, target == neutrophil_states[2]), 
                               source, ligand, receptor, pathway_name, interaction_name_2, prob.neg = prob)
arg1_diff_day1 <- dplyr::full_join(comm_day1_pos, comm_day1_neg, by = c("source", "ligand", "receptor", "pathway_name", "interaction_name_2"))
arg1_diff_day1 <- dplyr::mutate(arg1_diff_day1, dataset = "Lee_Day1_ARG1split")
arg1_diff_day1$prob.pos[is.na(arg1_diff_day1$prob.pos)] <- 0
arg1_diff_day1$prob.neg[is.na(arg1_diff_day1$prob.neg)] <- 0
arg1_diff_day1$prob.diff <- arg1_diff_day1$prob.pos - arg1_diff_day1$prob.neg
arg1_diff_day1 <- dplyr::arrange(arg1_diff_day1, dplyr::desc(prob.diff))
readr::write_csv(arg1_diff_day1, "differential_signals_Lee_Day1_ARG1split.csv")

# Save RDS file immediately after processing
saveRDS(cellchatOne, "cellchatOne_ARG1_complete.rds")
print("✅ Saved: cellchatOne_ARG1_complete.rds (Lee Day 1 ARG1-split - Complete)")

print("Generating visualizations for Lee Day 1 ARG1-split...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
neutrophil_base_label <- "Neutrophil"
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
cellchatOne <- readRDS("cellchatOne_ARG1_complete.rds")
groupSizeOne <- as.numeric(table(cellchatOne@idents)[rownames(cellchatOne@net$count)])
groupSizeOne[is.na(groupSizeOne)] <- 0

pdf("1dpiLee_Target_Arg1pos.pdf", width = 25, height = 15)
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

pdf("1dpiLee_Target_Arg1neg.pdf", width = 25, height = 15)
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

pdf("Lee_Day1_Signals_TO_ARG1_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states, vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Lee Day 1 ARG1-split - Signals TO ARG1 Neutrophil States")
dev.off()

png("Neu1dpiTargetBubbleLee_Arg1pos.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

png("Neu1dpiTargetBubbleLee_Arg1neg.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

print("Inline plots for Lee Day 1 ARG1-split:")
print("ARG1+ Neutrophils:")
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1pos")
print("ARG1- Neutrophils:")
netVisual_circle(cellchatOne@net$count, vertex.weight = groupSizeOne, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1neg")
print("Combined ARG1 Neutrophils:")
netVisual_bubble(cellchatOne, targets.use = neutrophil_states, remove.isolate = FALSE, 
                 title.name = "Lee Day 1 ARG1-split - Signals TO ARG1 Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatOne, measure = "weight", targets.use = neutrophil_states, 
                    title.name = "Lee Day 1 ARG1-split - Signals TO ARG1 States")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Lee Day 1 ARG1-split (insufficient color range)\n")
})
netVisual_heatmap_tnf_if_present(cellchatOne, "Lee Day 1 ARG1-split - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatOne) + ggtitle("Lee Day 1 ARG1-split - Communication Strength"))
print_signaling_role_network(cellchatOne, "Lee Day 1 ARG1-split - Signaling Role Network")

ht2_day1 <- netAnalysis_signalingRole_heatmap(cellchatOne, pattern = "incoming", font.size = 5)
print(ht2_day1 + plot_annotation(title = "Lee Day 1 ARG1-split - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatOne, pattern = "incoming") + ggtitle("Lee Day 1 ARG1-split - Communication Patterns"))
print(netAnalysis_dot(cellchatOne, pattern = "incoming") + ggtitle("Lee Day 1 ARG1-split - Pattern Strength"))

save.image("environment_analysis6_ARG1.RData")
print("✅ Saved: environment_analysis6_ARG1.RData (Analysis 6 Environment)")

# =============================================================================
# ANALYSIS 7: LEE DAY 3 ARG1-SPLIT - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 7: LEE DAY 3 ARG1-SPLIT")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Lee Day 3 ARG1-split
print("Loading and preparing Lee Day 3 data with ARG1 annotation...")
leeDat <- readRDS(file = "LeeDat.rds")
wangDat <- readRDS(file = "WangDat.rds")

# Process labels and merge datasets
leeDat$pruned_labels <- as.character(leeDat$celltype)
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
leeDat$dataset <- "Lee"
wangDat$dataset <- "Wang"

seuObj.integrated <- merge(leeDat, wangDat)
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Lee Day 3
threeDPILee <- subset(seuObj.integrated, time == "3" & dataset == "Lee")
DefaultAssay(threeDPILee) <- "RNA"

# Standardize neutrophil labels
threeDPILee$pruned_labels <- as.character(threeDPILee$pruned_labels)
threeDPILee$pruned_labels[threeDPILee$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
threeDPILee$pruned_labels[is.na(threeDPILee$pruned_labels)] <- "NA"

# ARG1 annotation for neutrophils
print("Annotating ARG1 status for neutrophils...")
neut_three_lee <- which(threeDPILee$pruned_labels == neutrophil_base_label)
gene_three_lee <- head(gene_candidates[gene_candidates %in% rownames(threeDPILee)], 1)
gene_three_lee_value <- c(gene_three_lee, NA_character_)[1]

expr_three_lee <- cellchat_norm_matrix(threeDPILee, gene_three_lee, neut_three_lee)
expr_three_lee_vec <- as.numeric(expr_three_lee)
expr_three_lee_pos <- expr_three_lee_vec[expr_three_lee_vec > 0]
if (length(expr_three_lee_pos) == 0) {
  warning("No Arg1-expressing neutrophils (Lee Day 3 ARG1-split); all neutrophils labeled Arg1neg.")
  threshold_three_lee <- NA_real_
  indicator_three_lee <- rep(FALSE, length(expr_three_lee_vec))
} else {
  threshold_three_lee <- stats::median(expr_three_lee_pos)
  indicator_three_lee <- expr_three_lee_vec >= threshold_three_lee
}
print(paste0("Lee Day 3 ARG1-split: Arg1 threshold=", threshold_three_lee,
             " | fraction Arg1pos among neutrophils=", round(sum(indicator_three_lee) / max(length(expr_three_lee_vec), 1L), 4)))

pos_three_lee <- neut_three_lee[indicator_three_lee]
neg_three_lee <- setdiff(neut_three_lee, pos_three_lee)

threeDPILee$arg1_status <- rep("Other", ncol(threeDPILee))
threeDPILee$arg1_status[neut_three_lee] <- "ARG1neg"
threeDPILee$arg1_status[pos_three_lee] <- "ARG1pos"

threeDPILee$cellchat_labels <- threeDPILee$pruned_labels
threeDPILee$cellchat_labels[pos_three_lee] <- neutrophil_states[1]
threeDPILee$cellchat_labels[neg_three_lee] <- neutrophil_states[2]

print("✅ Lee Day 3 ARG1-split data preparation completed")
print("ANALYSIS 7: Lee Day 3 ARG1-split - Signals TO ARG1 Neutrophil States")
print(paste(rep("=", 80), collapse=""))

cellchatThreeLee <- createCellChat(object = threeDPILee, group.by = "cellchat_labels")
CellChatDB <- CellChatDB.mouse
cellchatThreeLee@DB <- CellChatDB

cellchatThreeLee <- subsetData(cellchatThreeLee)
cellchatThreeLee <- identifyOverExpressedGenes(cellchatThreeLee)
cellchatThreeLee <- identifyOverExpressedInteractions(cellchatThreeLee)
cellchatThreeLee <- computeCommunProb(cellchatThreeLee, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatThreeLee <- filterCommunication(cellchatThreeLee, min.cells = 10)
cellchatThreeLee <- computeCommunProbPathway(cellchatThreeLee)
cellchatThreeLee <- aggregateNet(cellchatThreeLee)

cellchatThreeLee@idents <- factor(cellchatThreeLee@idents, levels = unique(c(levels(cellchatThreeLee@idents), neutrophil_states)))
groupSizeThreeLee <- as.numeric(table(cellchatThreeLee@idents)[rownames(cellchatThreeLee@net$count)])
groupSizeThreeLee[is.na(groupSizeThreeLee)] <- 0

cellchatThreeLee <- netAnalysis_computeCentrality(cellchatThreeLee, slot.name = "netP")

selectK(cellchatThreeLee, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatThreeLee <- identifyCommunicationPatterns(cellchatThreeLee, pattern = "incoming", k = nPatterns)

comm_day3_lee <- subsetCommunication(cellchatThreeLee)
comm_day3_lee_pos <- dplyr::select(dplyr::filter(comm_day3_lee, target == neutrophil_states[1]), 
                                   source, ligand, receptor, pathway_name, interaction_name_2, prob.pos = prob)
comm_day3_lee_neg <- dplyr::select(dplyr::filter(comm_day3_lee, target == neutrophil_states[2]), 
                                   source, ligand, receptor, pathway_name, interaction_name_2, prob.neg = prob)
arg1_diff_day3_lee <- dplyr::full_join(comm_day3_lee_pos, comm_day3_lee_neg, by = c("source", "ligand", "receptor", "pathway_name", "interaction_name_2"))
arg1_diff_day3_lee <- dplyr::mutate(arg1_diff_day3_lee, prob.pos = prob.pos, prob.neg = prob.neg, prob.diff = prob.pos - prob.neg, dataset = "Lee_Day3_ARG1split")
arg1_diff_day3_lee$prob.pos[is.na(arg1_diff_day3_lee$prob.pos)] <- 0
arg1_diff_day3_lee$prob.neg[is.na(arg1_diff_day3_lee$prob.neg)] <- 0
arg1_diff_day3_lee$prob.diff <- arg1_diff_day3_lee$prob.pos - arg1_diff_day3_lee$prob.neg
arg1_diff_day3_lee <- dplyr::arrange(arg1_diff_day3_lee, dplyr::desc(prob.diff))
readr::write_csv(arg1_diff_day3_lee, "differential_signals_Lee_Day3_ARG1split.csv")

# Save RDS file immediately after processing
saveRDS(cellchatThreeLee, "cellchatThreeLee_ARG1_complete.rds")
print("✅ Saved: cellchatThreeLee_ARG1_complete.rds (Lee Day 3 ARG1-split - Complete)")

print("Generating visualizations for Lee Day 3 ARG1-split...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
neutrophil_base_label <- "Neutrophil"
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
cellchatThreeLee <- readRDS("cellchatThreeLee_ARG1_complete.rds")
groupSizeThreeLee <- as.numeric(table(cellchatThreeLee@idents)[rownames(cellchatThreeLee@net$count)])
groupSizeThreeLee[is.na(groupSizeThreeLee)] <- 0

tiff("3dpiLee_Target_Arg1pos.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], vertex.label.cex = 1.8, 
                 title.name = "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1pos", margin = 0)
dev.off()

tiff("3dpiLee_Target_Arg1neg.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], vertex.label.cex = 1.8, 
                 title.name = "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1neg", margin = 0)
dev.off()

pdf("Lee_Day3_Signals_TO_ARG1_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states, vertex.label.cex = 2, 
                 title.name = "Lee Day 3 ARG1-split - Signals TO ARG1 Neutrophil States", margin = 1)
dev.off()

png("Neu3dpiTargetCircleLee_Arg1pos.png", width = 6, height = 6, units = "in", res = 1200)
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

png("Neu3dpiTargetCircleLee_Arg1neg.png", width = 6, height = 6, units = "in", res = 1200)
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

print("Inline plots for Lee Day 3 ARG1-split:")
print("ARG1+ Neutrophils:")
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1pos")
print("ARG1- Neutrophils:")
netVisual_circle(cellchatThreeLee@net$count, vertex.weight = groupSizeThreeLee, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1neg")
print("Combined ARG1 Neutrophils:")
netVisual_bubble(cellchatThreeLee, targets.use = neutrophil_states, remove.isolate = FALSE, 
                 title.name = "Lee Day 3 ARG1-split - Signals TO ARG1 Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatThreeLee, measure = "weight", targets.use = neutrophil_states, 
                    title.name = "Lee Day 3 ARG1-split - Signals TO ARG1 States")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Lee Day 3 ARG1-split (insufficient color range)\n")
})
netVisual_heatmap(cellchatThreeLee, signaling = "TNF", title.name = "Lee Day 3 ARG1-split - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatThreeLee) + ggtitle("Lee Day 3 ARG1-split - Communication Strength"))
print(netAnalysis_signalingRole_network(cellchatThreeLee, width = 8, height = 2.5, font.size = 10) + ggtitle("Lee Day 3 ARG1-split - Signaling Role Network"))

ht2_lee_day3 <- netAnalysis_signalingRole_heatmap(cellchatThreeLee, pattern = "incoming", font.size = 5)
print(ht2_lee_day3 + plot_annotation(title = "Lee Day 3 ARG1-split - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatThreeLee, pattern = "incoming") + ggtitle("Lee Day 3 ARG1-split - Communication Patterns"))
print(netAnalysis_dot(cellchatThreeLee, pattern = "incoming") + ggtitle("Lee Day 3 ARG1-split - Pattern Strength"))

save.image("environment_analysis7_ARG1.RData")
print("✅ Saved: environment_analysis7_ARG1.RData (Analysis 7 Environment)")

# =============================================================================
# ANALYSIS 8: WANG DAY 3 ARG1-SPLIT - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 8: WANG DAY 3 ARG1-SPLIT")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Wang Day 3 ARG1-split
print("Loading and preparing Wang Day 3 data with ARG1 annotation...")
leeDat <- readRDS(file = "LeeDat.rds")
wangDat <- readRDS(file = "WangDat.rds")

# Process labels and merge datasets
leeDat$pruned_labels <- as.character(leeDat$celltype)
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
leeDat$dataset <- "Lee"
wangDat$dataset <- "Wang"

seuObj.integrated <- merge(leeDat, wangDat)
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Wang Day 3
threeDPIWang <- subset(seuObj.integrated, time == "3" & dataset == "Wang")
DefaultAssay(threeDPIWang) <- "RNA"

# Standardize neutrophil labels
threeDPIWang$pruned_labels <- as.character(threeDPIWang$pruned_labels)
threeDPIWang$pruned_labels[threeDPIWang$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
threeDPIWang$pruned_labels[is.na(threeDPIWang$pruned_labels)] <- "NA"

# ARG1 annotation for neutrophils
print("Annotating ARG1 status for neutrophils...")
neut_three_wang <- which(threeDPIWang$pruned_labels == neutrophil_base_label)
gene_three_wang <- head(gene_candidates[gene_candidates %in% rownames(threeDPIWang)], 1)
gene_three_wang_value <- c(gene_three_wang, NA_character_)[1]

expr_three_wang <- Seurat::GetAssayData(threeDPIWang, slot = "data")[gene_three_wang, neut_three_wang, drop = FALSE]
expr_three_wang_vec <- as.numeric(expr_three_wang)
expr_three_wang_pos <- expr_three_wang_vec[expr_three_wang_vec > 0]
threshold_three_wang <- stats::median(expr_three_wang_pos)
threshold_three_wang[is.na(threshold_three_wang)] <- 0

indicator_three_wang <- expr_three_wang_vec >= threshold_three_wang
indicator_three_wang <- indicator_three_wang & rep(length(expr_three_wang_pos) > 0, length(indicator_three_wang))

pos_three_wang <- neut_three_wang[indicator_three_wang]
neg_three_wang <- setdiff(neut_three_wang, pos_three_wang)

threeDPIWang$arg1_status <- rep("Other", ncol(threeDPIWang))
threeDPIWang$arg1_status[neut_three_wang] <- "ARG1neg"
threeDPIWang$arg1_status[pos_three_wang] <- "ARG1pos"

threeDPIWang$cellchat_labels <- threeDPIWang$pruned_labels
threeDPIWang$cellchat_labels[pos_three_wang] <- paste0(threeDPIWang$cellchat_labels[pos_three_wang], "Arg1pos")
threeDPIWang$cellchat_labels[neg_three_wang] <- paste0(threeDPIWang$cellchat_labels[neg_three_wang], "Arg1neg")

print("✅ Wang Day 3 ARG1-split data preparation completed")
print("ANALYSIS 8: Wang Day 3 ARG1-split - Signals TO ARG1 Neutrophil States")
print(paste(rep("=", 80), collapse=""))

cellchatThreeWang <- createCellChat(object = threeDPIWang, group.by = "cellchat_labels")
CellChatDB <- CellChatDB.mouse
cellchatThreeWang@DB <- CellChatDB

cellchatThreeWang <- subsetData(cellchatThreeWang)
cellchatThreeWang <- identifyOverExpressedGenes(cellchatThreeWang)
cellchatThreeWang <- identifyOverExpressedInteractions(cellchatThreeWang)
cellchatThreeWang <- computeCommunProb(cellchatThreeWang, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatThreeWang <- filterCommunication(cellchatThreeWang, min.cells = 10)
cellchatThreeWang <- computeCommunProbPathway(cellchatThreeWang)
cellchatThreeWang <- aggregateNet(cellchatThreeWang)

cellchatThreeWang@idents <- factor(cellchatThreeWang@idents, levels = unique(c(levels(cellchatThreeWang@idents), neutrophil_states)))
groupSizeThreeWang <- as.numeric(table(cellchatThreeWang@idents)[rownames(cellchatThreeWang@net$count)])
groupSizeThreeWang[is.na(groupSizeThreeWang)] <- 0

cellchatThreeWang <- netAnalysis_computeCentrality(cellchatThreeWang, slot.name = "netP")

selectK(cellchatThreeWang, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatThreeWang <- identifyCommunicationPatterns(cellchatThreeWang, pattern = "incoming", k = nPatterns)

comm_day3_wang <- subsetCommunication(cellchatThreeWang)
comm_day3_wang_pos <- dplyr::select(dplyr::filter(comm_day3_wang, target == neutrophil_states[1]), 
                                    source, ligand, receptor, pathway_name, interaction_name_2, prob.pos = prob)
comm_day3_wang_neg <- dplyr::select(dplyr::filter(comm_day3_wang, target == neutrophil_states[2]), 
                                    source, ligand, receptor, pathway_name, interaction_name_2, prob.neg = prob)
arg1_diff_day3_wang <- dplyr::full_join(comm_day3_wang_pos, comm_day3_wang_neg, by = c("source", "ligand", "receptor", "pathway_name", "interaction_name_2"))
arg1_diff_day3_wang <- dplyr::mutate(arg1_diff_day3_wang, dataset = "Wang_Day3_ARG1split")
arg1_diff_day3_wang$prob.pos[is.na(arg1_diff_day3_wang$prob.pos)] <- 0
arg1_diff_day3_wang$prob.neg[is.na(arg1_diff_day3_wang$prob.neg)] <- 0
arg1_diff_day3_wang$prob.diff <- arg1_diff_day3_wang$prob.pos - arg1_diff_day3_wang$prob.neg
arg1_diff_day3_wang <- dplyr::arrange(arg1_diff_day3_wang, dplyr::desc(prob.diff))
readr::write_csv(arg1_diff_day3_wang, "differential_signals_Wang_Day3_ARG1split.csv")

# Save RDS file immediately after processing
saveRDS(cellchatThreeWang, "cellchatThreeWang_ARG1_complete.rds")
print("✅ Saved: cellchatThreeWang_ARG1_complete.rds (Wang Day 3 ARG1-split - Complete)")

print("Generating visualizations for Wang Day 3 ARG1-split...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
neutrophil_base_label <- "Neutrophil"
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
cellchatThreeWang <- readRDS("cellchatThreeWang_ARG1_complete.rds")
groupSizeThreeWang <- as.numeric(table(cellchatThreeWang@idents)[rownames(cellchatThreeWang@net$count)])
groupSizeThreeWang[is.na(groupSizeThreeWang)] <- 0

tiff("3dpiWang_Target_Arg1pos.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeWang@net$count, vertex.weight = groupSizeThreeWang, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], vertex.label.cex = 1.8, 
                 title.name = "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1pos", margin = 0)
dev.off()

tiff("3dpiWang_Target_Arg1neg.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeWang@net$count, vertex.weight = groupSizeThreeWang, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], vertex.label.cex = 1.8, 
                 title.name = "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1neg", margin = 0)
dev.off()

pdf("Wang_Day3_Signals_TO_ARG1_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatThreeWang@net$count, vertex.weight = groupSizeThreeWang, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states, vertex.label.cex = 2, 
                 title.name = "Wang Day 3 ARG1-split - Signals TO ARG1 Neutrophil States", margin = CELLCHAT_MARGIN_CIRCLE)
dev.off()

png("Neu3dpiTargetBubbleWang_Arg1pos.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatThreeWang, targets.use = neutrophil_states[1], remove.isolate = FALSE, 
                 title.name = "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

png("Neu3dpiTargetBubbleWang_Arg1neg.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatThreeWang, targets.use = neutrophil_states[2], remove.isolate = FALSE, 
                 title.name = "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

print("Inline plots for Wang Day 3 ARG1-split:")
print("ARG1+ Neutrophils:")
netVisual_circle(cellchatThreeWang@net$count, vertex.weight = groupSizeThreeWang, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1pos")
print("ARG1- Neutrophils:")
netVisual_circle(cellchatThreeWang@net$count, vertex.weight = groupSizeThreeWang, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1neg")
print("Combined ARG1 Neutrophils:")
netVisual_bubble(cellchatThreeWang, targets.use = neutrophil_states, remove.isolate = FALSE, 
                 title.name = "Wang Day 3 ARG1-split - Signals TO ARG1 Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatThreeWang, measure = "weight", targets.use = neutrophil_states, 
                    title.name = "Wang Day 3 ARG1-split - Signals TO ARG1 States")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Wang Day 3 ARG1-split (insufficient color range)\n")
})
netVisual_heatmap(cellchatThreeWang, signaling = "TNF", title.name = "Wang Day 3 ARG1-split - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatThreeWang) + ggtitle("Wang Day 3 ARG1-split - Communication Strength"))
print(netAnalysis_signalingRole_network(cellchatThreeWang, width = 8, height = 2.5, font.size = 10) + ggtitle("Wang Day 3 ARG1-split - Signaling Role Network"))

ht2_wang_day3 <- netAnalysis_signalingRole_heatmap(cellchatThreeWang, pattern = "incoming", font.size = 5)
print(ht2_wang_day3 + plot_annotation(title = "Wang Day 3 ARG1-split - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatThreeWang, pattern = "incoming") + ggtitle("Wang Day 3 ARG1-split - Communication Patterns"))
print(netAnalysis_dot(cellchatThreeWang, pattern = "incoming") + ggtitle("Wang Day 3 ARG1-split - Pattern Strength"))

save.image("environment_analysis8_ARG1.RData")
print("✅ Saved: environment_analysis8_ARG1.RData (Analysis 8 Environment)")

# =============================================================================
# ANALYSIS 9: QIN DAY 1 ARG1-SPLIT - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 9: QIN DAY 1 ARG1-SPLIT")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Qin Day 1 ARG1-split
print("Loading and preparing Qin Day 1 data with ARG1 annotation...")
qinDat <- readRDS(file = "QinDat.rds")
wangDat <- readRDS(file = "WangDat.rds")
leeDat <- readRDS(file = "LeeDat.rds")

# Process labels and merge datasets
qinDat$pruned_labels <- as.character(qinDat$celltype)
qinDat$dataset <- "Qin"
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
wangDat$dataset <- "Wang"
leeDat$pruned_labels <- as.character(leeDat$celltype)
leeDat$dataset <- "Lee"

seuObj.integrated <- merge(leeDat, y = c(wangDat, qinDat))
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Qin Day 1
oneDPIQin <- subset(seuObj.integrated, time == "1" & dataset == "Qin")
DefaultAssay(oneDPIQin) <- "RNA"

# Standardize neutrophil labels
oneDPIQin$pruned_labels <- as.character(oneDPIQin$pruned_labels)
oneDPIQin$pruned_labels[oneDPIQin$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
oneDPIQin$pruned_labels[is.na(oneDPIQin$pruned_labels)] <- "NA"

# ARG1 annotation for neutrophils
print("Annotating ARG1 status for neutrophils...")
neut_one_qin <- which(oneDPIQin$pruned_labels == neutrophil_base_label)
gene_one_qin <- head(gene_candidates[gene_candidates %in% rownames(oneDPIQin)], 1)
gene_one_qin_value <- c(gene_one_qin, NA_character_)[1]

expr_one_qin <- cellchat_norm_matrix(oneDPIQin, gene_one_qin, neut_one_qin)
expr_one_qin_vec <- as.numeric(expr_one_qin)
expr_one_qin_pos <- expr_one_qin_vec[expr_one_qin_vec > 0]
if (length(expr_one_qin_pos) == 0) {
  warning("No Arg1-expressing neutrophils (Qin Day 1 ARG1-split); all neutrophils labeled Arg1neg.")
  threshold_one_qin <- NA_real_
  indicator_one_qin <- rep(FALSE, length(expr_one_qin_vec))
} else {
  threshold_one_qin <- stats::median(expr_one_qin_pos)
  indicator_one_qin <- expr_one_qin_vec >= threshold_one_qin
}
print(paste0("Qin Day 1 ARG1-split: Arg1 threshold=", threshold_one_qin,
             " | fraction Arg1pos among neutrophils=", round(sum(indicator_one_qin) / max(length(expr_one_qin_vec), 1L), 4)))

pos_one_qin <- neut_one_qin[indicator_one_qin]
neg_one_qin <- setdiff(neut_one_qin, pos_one_qin)

oneDPIQin$arg1_status <- rep("Other", ncol(oneDPIQin))
oneDPIQin$arg1_status[neut_one_qin] <- "ARG1neg"
oneDPIQin$arg1_status[pos_one_qin] <- "ARG1pos"

oneDPIQin$cellchat_labels <- oneDPIQin$pruned_labels
oneDPIQin$cellchat_labels[pos_one_qin] <- neutrophil_states[1]
oneDPIQin$cellchat_labels[neg_one_qin] <- neutrophil_states[2]

print("✅ Qin Day 1 ARG1-split data preparation completed")
print("ANALYSIS 9: Qin Day 1 ARG1-split - Signals TO ARG1 Neutrophil States")
print(paste(rep("=", 80), collapse=""))

cellchatOneQin <- createCellChat(object = oneDPIQin, group.by = "cellchat_labels")
CellChatDB <- CellChatDB.mouse
cellchatOneQin@DB <- CellChatDB

cellchatOneQin <- subsetData(cellchatOneQin)
cellchatOneQin <- identifyOverExpressedGenes(cellchatOneQin)
cellchatOneQin <- identifyOverExpressedInteractions(cellchatOneQin)
cellchatOneQin <- computeCommunProb(cellchatOneQin, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatOneQin <- filterCommunication(cellchatOneQin, min.cells = 10)
cellchatOneQin <- computeCommunProbPathway(cellchatOneQin)
cellchatOneQin <- aggregateNet(cellchatOneQin)

cellchatOneQin@idents <- factor(cellchatOneQin@idents, levels = unique(c(levels(cellchatOneQin@idents), neutrophil_states)))
groupSizeOneQin <- as.numeric(table(cellchatOneQin@idents)[rownames(cellchatOneQin@net$count)])
groupSizeOneQin[is.na(groupSizeOneQin)] <- 0

cellchatOneQin <- netAnalysis_computeCentrality(cellchatOneQin, slot.name = "netP")

selectK(cellchatOneQin, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatOneQin <- identifyCommunicationPatterns(cellchatOneQin, pattern = "incoming", k = nPatterns)

comm_qin_day1 <- subsetCommunication(cellchatOneQin)
comm_qin_day1_pos <- dplyr::select(dplyr::filter(comm_qin_day1, target == neutrophil_states[1]), 
                                   source, ligand, receptor, pathway_name, interaction_name_2, prob.pos = prob)
comm_qin_day1_neg <- dplyr::select(dplyr::filter(comm_qin_day1, target == neutrophil_states[2]), 
                                   source, ligand, receptor, pathway_name, interaction_name_2, prob.neg = prob)
arg1_diff_qin_day1 <- dplyr::full_join(comm_qin_day1_pos, comm_qin_day1_neg, by = c("source", "ligand", "receptor", "pathway_name", "interaction_name_2"))
arg1_diff_qin_day1 <- dplyr::mutate(arg1_diff_qin_day1, prob.pos = prob.pos, prob.neg = prob.neg, prob.diff = prob.pos - prob.neg, dataset = "Qin_Day1_ARG1split")
arg1_diff_qin_day1$prob.pos[is.na(arg1_diff_qin_day1$prob.pos)] <- 0
arg1_diff_qin_day1$prob.neg[is.na(arg1_diff_qin_day1$prob.neg)] <- 0
arg1_diff_qin_day1$prob.diff <- arg1_diff_qin_day1$prob.pos - arg1_diff_qin_day1$prob.neg
arg1_diff_qin_day1 <- dplyr::arrange(arg1_diff_qin_day1, dplyr::desc(prob.diff))
readr::write_csv(arg1_diff_qin_day1, "differential_signals_Qin_Day1_ARG1split.csv")

# Save RDS file immediately after processing
saveRDS(cellchatOneQin, "cellchatOneQin_ARG1_complete.rds")
print("✅ Saved: cellchatOneQin_ARG1_complete.rds (Qin Day 1 ARG1-split - Complete)")

print("Generating visualizations for Qin Day 1 ARG1-split...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
neutrophil_base_label <- "Neutrophil"
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
cellchatOneQin <- readRDS("cellchatOneQin_ARG1_complete.rds")
groupSizeOneQin <- as.numeric(table(cellchatOneQin@idents)[rownames(cellchatOneQin@net$count)])
groupSizeOneQin[is.na(groupSizeOneQin)] <- 0

pdf("1dpiQin_Target_Arg1pos.pdf", width = 25, height = 15)
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

pdf("1dpiQin_Target_Arg1neg.pdf", width = 25, height = 15)
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

pdf("Qin_Day1_Signals_TO_ARG1_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states, vertex.label.cex = 2, margin = CELLCHAT_MARGIN_CIRCLE, title.name = "Qin Day 1 ARG1-split - Signals TO ARG1 Neutrophil States")
dev.off()

png("Neu1dpiTargetBubbleQin_Arg1pos.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

png("Neu1dpiTargetBubbleQin_Arg1neg.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

print("Inline plots for Qin Day 1 ARG1-split:")
print("ARG1+ Neutrophils:")
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1pos")
print("ARG1- Neutrophils:")
netVisual_circle(cellchatOneQin@net$count, vertex.weight = groupSizeOneQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1neg")
print("Combined ARG1 Neutrophils:")
netVisual_bubble(cellchatOneQin, targets.use = neutrophil_states, remove.isolate = FALSE, 
                 title.name = "Qin Day 1 ARG1-split - Signals TO ARG1 Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatOneQin, measure = "weight", targets.use = neutrophil_states, 
                    title.name = "Qin Day 1 ARG1-split - Signals TO ARG1 States")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Qin Day 1 ARG1-split (insufficient color range)\n")
})
netVisual_heatmap_tnf_if_present(cellchatOneQin, "Qin Day 1 ARG1-split - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatOneQin) + ggtitle("Qin Day 1 ARG1-split - Communication Strength"))
print_signaling_role_network(cellchatOneQin, "Qin Day 1 ARG1-split - Signaling Role Network")

ht2_qin_day1 <- netAnalysis_signalingRole_heatmap(cellchatOneQin, pattern = "incoming", font.size = 5)
print(ht2_qin_day1 + plot_annotation(title = "Qin Day 1 ARG1-split - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatOneQin, pattern = "incoming") + ggtitle("Qin Day 1 ARG1-split - Communication Patterns"))
print(netAnalysis_dot(cellchatOneQin, pattern = "incoming") + ggtitle("Qin Day 1 ARG1-split - Pattern Strength"))

save.image("environment_analysis9_ARG1.RData")
print("✅ Saved: environment_analysis9_ARG1.RData (Analysis 9 Environment)")

# =============================================================================
# ANALYSIS 10: QIN DAY 3 ARG1-SPLIT - COMPLETE ANALYSIS
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("ANALYSIS 10: QIN DAY 3 ARG1-SPLIT")
print(paste(rep("=", 80), collapse=""))

# Data preparation for Qin Day 3 ARG1-split
print("Loading and preparing Qin Day 3 data with ARG1 annotation...")
qinDat <- readRDS(file = "QinDat.rds")
wangDat <- readRDS(file = "WangDat.rds")
leeDat <- readRDS(file = "LeeDat.rds")

# Process labels and merge datasets
qinDat$pruned_labels <- as.character(qinDat$celltype)
qinDat$dataset <- "Qin"
wangDat$pruned_labels <- as.character(wangDat$pruned_labels)
wangDat$dataset <- "Wang"
leeDat$pruned_labels <- as.character(leeDat$celltype)
leeDat$dataset <- "Lee"

seuObj.integrated <- merge(leeDat, y = c(wangDat, qinDat))
seuObj.integrated$pruned_labels <- as.character(seuObj.integrated$pruned_labels)
seuObj.integrated$pruned_labels[is.na(seuObj.integrated$pruned_labels)] <- "NA"

# Process time information
integrated_time <- as.character(seuObj.integrated$time)
integrated_time[is.na(integrated_time)] <- ""
lee_mask <- seuObj.integrated$dataset == "Lee"
integrated_time[lee_mask] <- as.character(seuObj.integrated$orig.ident)[lee_mask]
integrated_time <- gsub("_.*", "", integrated_time)
integrated_time <- gsub("dpi", "", integrated_time)
integrated_time <- gsub("uninj", "0", integrated_time)
integrated_time[integrated_time == ""] <- NA_character_
seuObj.integrated$time <- integrated_time

# Subset for Qin Day 3
threeDPIQin <- subset(seuObj.integrated, time == "3" & dataset == "Qin")
DefaultAssay(threeDPIQin) <- "RNA"

# Standardize neutrophil labels
threeDPIQin$pruned_labels <- as.character(threeDPIQin$pruned_labels)
threeDPIQin$pruned_labels[threeDPIQin$pruned_labels %in% neutrophil_aliases] <- neutrophil_base_label
threeDPIQin$pruned_labels[is.na(threeDPIQin$pruned_labels)] <- "NA"

# ARG1 annotation for neutrophils
print("Annotating ARG1 status for neutrophils...")
neut_three_qin <- which(threeDPIQin$pruned_labels == neutrophil_base_label)
gene_three_qin <- head(gene_candidates[gene_candidates %in% rownames(threeDPIQin)], 1)
gene_three_qin_value <- c(gene_three_qin, NA_character_)[1]

expr_three_qin <- Seurat::GetAssayData(threeDPIQin, slot = "data")[gene_three_qin, neut_three_qin, drop = FALSE]
expr_three_qin_vec <- as.numeric(expr_three_qin)
expr_three_qin_pos <- expr_three_qin_vec[expr_three_qin_vec > 0]
threshold_three_qin <- stats::median(expr_three_qin_pos)
threshold_three_qin[is.na(threshold_three_qin)] <- 0

indicator_three_qin <- expr_three_qin_vec >= threshold_three_qin
indicator_three_qin <- indicator_three_qin & rep(length(expr_three_qin_pos) > 0, length(indicator_three_qin))

pos_three_qin <- neut_three_qin[indicator_three_qin]
neg_three_qin <- setdiff(neut_three_qin, pos_three_qin)

threeDPIQin$arg1_status <- rep("Other", ncol(threeDPIQin))
threeDPIQin$arg1_status[neut_three_qin] <- "ARG1neg"
threeDPIQin$arg1_status[pos_three_qin] <- "ARG1pos"

threeDPIQin$cellchat_labels <- threeDPIQin$pruned_labels
threeDPIQin$cellchat_labels[pos_three_qin] <- paste0(threeDPIQin$cellchat_labels[pos_three_qin], "Arg1pos")
threeDPIQin$cellchat_labels[neg_three_qin] <- paste0(threeDPIQin$cellchat_labels[neg_three_qin], "Arg1neg")

print("✅ Qin Day 3 ARG1-split data preparation completed")
print("ANALYSIS 10: Qin Day 3 ARG1-split - Signals TO ARG1 Neutrophil States")
print(paste(rep("=", 80), collapse=""))

cellchatThreeQin <- createCellChat(object = threeDPIQin, group.by = "cellchat_labels")
CellChatDB <- CellChatDB.mouse
cellchatThreeQin@DB <- CellChatDB

cellchatThreeQin <- subsetData(cellchatThreeQin)
cellchatThreeQin <- identifyOverExpressedGenes(cellchatThreeQin)
cellchatThreeQin <- identifyOverExpressedInteractions(cellchatThreeQin)
cellchatThreeQin <- computeCommunProb(cellchatThreeQin, type = "triMean", nboot = CELLCHAT_NBOOT, raw.use = FALSE)
cellchatThreeQin <- filterCommunication(cellchatThreeQin, min.cells = 10)
cellchatThreeQin <- computeCommunProbPathway(cellchatThreeQin)
cellchatThreeQin <- aggregateNet(cellchatThreeQin)

cellchatThreeQin@idents <- factor(cellchatThreeQin@idents, levels = unique(c(levels(cellchatThreeQin@idents), neutrophil_states)))
groupSizeThreeQin <- as.numeric(table(cellchatThreeQin@idents)[rownames(cellchatThreeQin@net$count)])
groupSizeThreeQin[is.na(groupSizeThreeQin)] <- 0

cellchatThreeQin <- netAnalysis_computeCentrality(cellchatThreeQin, slot.name = "netP")

selectK(cellchatThreeQin, pattern = "incoming")
nPatterns <- N_PATTERNS_INCOMING
cellchatThreeQin <- identifyCommunicationPatterns(cellchatThreeQin, pattern = "incoming", k = nPatterns)

comm_qin_day3 <- subsetCommunication(cellchatThreeQin)
comm_qin_day3_pos <- dplyr::select(dplyr::filter(comm_qin_day3, target == neutrophil_states[1]), 
                                   source, ligand, receptor, pathway_name, interaction_name_2, prob.pos = prob)
comm_qin_day3_neg <- dplyr::select(dplyr::filter(comm_qin_day3, target == neutrophil_states[2]), 
                                   source, ligand, receptor, pathway_name, interaction_name_2, prob.neg = prob)
arg1_diff_qin_day3 <- dplyr::full_join(comm_qin_day3_pos, comm_qin_day3_neg, by = c("source", "ligand", "receptor", "pathway_name", "interaction_name_2"))
arg1_diff_qin_day3 <- dplyr::mutate(arg1_diff_qin_day3, prob.pos = prob.pos, prob.neg = prob.neg, prob.diff = prob.pos - prob.neg, dataset = "Qin_Day3_ARG1split")
arg1_diff_qin_day3$prob.pos[is.na(arg1_diff_qin_day3$prob.pos)] <- 0
arg1_diff_qin_day3$prob.neg[is.na(arg1_diff_qin_day3$prob.neg)] <- 0
arg1_diff_qin_day3$prob.diff <- arg1_diff_qin_day3$prob.pos - arg1_diff_qin_day3$prob.neg
arg1_diff_qin_day3 <- dplyr::arrange(arg1_diff_qin_day3, dplyr::desc(prob.diff))
readr::write_csv(arg1_diff_qin_day3, "differential_signals_Qin_Day3_ARG1split.csv")

# Save RDS file immediately after processing
saveRDS(cellchatThreeQin, "cellchatThreeQin_ARG1_complete.rds")
print("✅ Saved: cellchatThreeQin_ARG1_complete.rds (Qin Day 3 ARG1-split - Complete)")

print("Generating visualizations for Qin Day 3 ARG1-split...")
print("Loading pre-processed CellChat object...")
# Load the RDS file we just saved
neutrophil_base_label <- "Neutrophil"
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
cellchatThreeQin <- readRDS("cellchatThreeQin_ARG1_complete.rds")
groupSizeThreeQin <- as.numeric(table(cellchatThreeQin@idents)[rownames(cellchatThreeQin@net$count)])
groupSizeThreeQin[is.na(groupSizeThreeQin)] <- 0

tiff("3dpiQin_Target_Arg1pos.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeQin@net$count, vertex.weight = groupSizeThreeQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], vertex.label.cex = 1.8, 
                 title.name = "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1pos", margin = 0)
dev.off()

tiff("3dpiQin_Target_Arg1neg.tiff", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchatThreeQin@net$count, vertex.weight = groupSizeThreeQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], vertex.label.cex = 1.8, 
                 title.name = "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1neg", margin = 0)
dev.off()

pdf("Qin_Day3_Signals_TO_ARG1_Neutrophils.pdf", width = 25, height = 15)
netVisual_circle(cellchatThreeQin@net$count, vertex.weight = groupSizeThreeQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states, vertex.label.cex = 2, 
                 title.name = "Qin Day 3 ARG1-split - Signals TO ARG1 Neutrophil States", margin = CELLCHAT_MARGIN_CIRCLE)
dev.off()

png("Neu3dpiTargetBubbleQin_Arg1pos.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatThreeQin, targets.use = neutrophil_states[1], remove.isolate = FALSE, 
                 title.name = "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1pos")
dev.off()

png("Neu3dpiTargetBubbleQin_Arg1neg.png", width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
netVisual_bubble(cellchatThreeQin, targets.use = neutrophil_states[2], remove.isolate = FALSE, 
                 title.name = "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1neg")
dev.off()

print("Inline plots for Qin Day 3 ARG1-split:")
print("ARG1+ Neutrophils:")
netVisual_circle(cellchatThreeQin@net$count, vertex.weight = groupSizeThreeQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[1], title.name = "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1pos")
print("ARG1- Neutrophils:")
netVisual_circle(cellchatThreeQin@net$count, vertex.weight = groupSizeThreeQin, weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = neutrophil_states[2], title.name = "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1neg")
print("Combined ARG1 Neutrophils:")
netVisual_bubble(cellchatThreeQin, targets.use = neutrophil_states, remove.isolate = FALSE, 
                 title.name = "Qin Day 3 ARG1-split - Signals TO ARG1 Neutrophils")
tryCatch({
  netVisual_heatmap(cellchatThreeQin, measure = "weight", targets.use = neutrophil_states, 
                    title.name = "Qin Day 3 ARG1-split - Signals TO ARG1 States")
}, error = function(e) {
  cat("Note: Weight heatmap skipped for Qin Day 3 ARG1-split (insufficient color range)\n")
})
netVisual_heatmap_tnf_if_present(cellchatThreeQin, "Qin Day 3 ARG1-split - TNF Pathway")

print(netAnalysis_signalingRole_scatter(cellchatThreeQin) + ggtitle("Qin Day 3 ARG1-split - Communication Strength"))
print_signaling_role_network(cellchatThreeQin, "Qin Day 3 ARG1-split - Signaling Role Network")

ht2_qin_day3 <- netAnalysis_signalingRole_heatmap(cellchatThreeQin, pattern = "incoming", font.size = 5)
print(ht2_qin_day3 + plot_annotation(title = "Qin Day 3 ARG1-split - Incoming Signals TO Neutrophils"))

print(netAnalysis_river(cellchatThreeQin, pattern = "incoming") + ggtitle("Qin Day 3 ARG1-split - Communication Patterns"))
print(netAnalysis_dot(cellchatThreeQin, pattern = "incoming") + ggtitle("Qin Day 3 ARG1-split - Pattern Strength"))

save.image("environment_analysis10_ARG1.RData")
print("✅ Saved: environment_analysis10_ARG1.RData (Analysis 10 Environment)")

# =============================================================================
# SAVE FINAL ENVIRONMENT
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("🎉 All 10 analyses complete!")
print(paste(rep("=", 80), collapse=""))

# =============================================================================
# CALCULATE GLOBAL SCALING PARAMETERS FOR CONSISTENT COMPARISONS
# (CAN BE RUN INDEPENDENTLY - All variables defined here)
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("CALCULATING GLOBAL SCALING PARAMETERS")
print(paste(rep("=", 80), collapse=""))

# Define required variables for standalone execution
neutrophil_base_label <- "Neutrophil"
neutrophil_states <- c(paste0(neutrophil_base_label, "sArg1pos"), paste0(neutrophil_base_label, "sArg1neg"))
lee_neutrophil_label <- "Neutrophil"
wang_neutrophil_label <- "Neutrophil"
qin_neutrophil_label <- "Neutrophil"

# Load all CellChat objects to calculate global scales
cellchat_objects <- list(
  "Lee_Day1_Baseline" = readRDS("cellchatOne_baseline_complete.rds"),
  "Lee_Day3_Baseline" = readRDS("cellchatThreeLee_baseline_complete.rds"),
  "Wang_Day3_Baseline" = readRDS("cellchatThreeWang_baseline_complete.rds"),
  "Qin_Day1_Baseline" = readRDS("cellchatOneQin_baseline_complete.rds"),
  "Qin_Day3_Baseline" = readRDS("cellchatThreeQin_baseline_complete.rds"),
  "Lee_Day1_ARG1" = readRDS("cellchatOne_ARG1_complete.rds"),
  "Lee_Day3_ARG1" = readRDS("cellchatThreeLee_ARG1_complete.rds"),
  "Wang_Day3_ARG1" = readRDS("cellchatThreeWang_ARG1_complete.rds"),
  "Qin_Day1_ARG1" = readRDS("cellchatOneQin_ARG1_complete.rds"),
  "Qin_Day3_ARG1" = readRDS("cellchatThreeQin_ARG1_complete.rds")
)

# Calculate global min/max values across all analyses
all_count_values <- c()
all_prob_values <- c()

for (obj_name in names(cellchat_objects)) {
  obj <- cellchat_objects[[obj_name]]
  all_count_values <- c(all_count_values, as.vector(obj@net$count))
  all_prob_values <- c(all_prob_values, as.vector(obj@net$prob))
}

# Set global scaling parameters
GLOBAL_MIN_COUNT <<- min(all_count_values[all_count_values > 0], na.rm = TRUE)
GLOBAL_MAX_COUNT <<- max(all_count_values, na.rm = TRUE)
GLOBAL_MIN_PROB <<- min(all_prob_values[all_prob_values > 0], na.rm = TRUE)
GLOBAL_MAX_PROB <<- max(all_prob_values, na.rm = TRUE)

# Create fixed color scale
FIXED_COLOR_SCALE <<- colorRampPalette(c("white", "red"))(100)

print(paste("✅ Global scaling parameters calculated:"))
print(paste("  Count range:", round(GLOBAL_MIN_COUNT, 3), "to", round(GLOBAL_MAX_COUNT, 3)))
print(paste("  Probability range:", round(GLOBAL_MIN_PROB, 3), "to", round(GLOBAL_MAX_PROB, 3)))

# =============================================================================
# REGENERATE ALL PLOTS WITH FIXED SCALING
# (This entire section can be run independently after defining variables above)
# =============================================================================
print(paste(rep("=", 80), collapse=""))
print("REGENERATING ALL PLOTS WITH FIXED SCALING")
print(paste(rep("=", 80), collapse=""))

# Function to create auto-scaled circle plot
create_fixed_circle_plot <- function(cellchat_obj, group_sizes, target_label, title_name, filename) {
  pdf(filename, width = 25, height = 15)
  netVisual_circle(cellchat_obj@net$count, 
                   vertex.weight = group_sizes, 
                   weight.scale = TRUE,  # Auto-scaling for better visibility
                   label.edge = FALSE, 
                   targets.use = target_label, 
                   vertex.label.cex = 2, 
                   margin = CELLCHAT_MARGIN_CIRCLE, 
                   title.name = title_name)
  dev.off()
}

# Function to create fixed-scale bubble plot
create_fixed_bubble_plot <- function(cellchat_obj, target_label, title_name, filename) {
  png(filename, width = 6, height = 8, units = "in", res = 1200, pointsize = 3.5)
  netVisual_bubble(cellchat_obj, 
                   targets.use = target_label, 
                   remove.isolate = FALSE, 
                   title.name = title_name,
                   color.heatmap = "Spectral",
                   direction = 1)
  dev.off()
}

# Regenerate all plots with fixed scaling
print("Regenerating Lee Day 1 Baseline plots...")
obj <- cellchat_objects[["Lee_Day1_Baseline"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, lee_neutrophil_label, 
                         "Lee Day 1 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "1dpiLee_Baseline_Target_Neutrophils_FixedScale.pdf")
create_fixed_bubble_plot(obj, lee_neutrophil_label,
                         "Lee Day 1 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "Neu1dpiBaseline_TargetBubbleLee_FixedScale.png")

print("Regenerating Lee Day 3 Baseline plots...")
obj <- cellchat_objects[["Lee_Day3_Baseline"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, lee_neutrophil_label,
                         "Lee Day 3 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "3dpiLee_Baseline_Target_Neutrophils_FixedScale.pdf")

print("Regenerating Wang Day 3 Baseline plots...")
obj <- cellchat_objects[["Wang_Day3_Baseline"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, wang_neutrophil_label,
                         "Wang Day 3 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "3dpiWang_Baseline_Target_Neutrophils_FixedScale.pdf")
create_fixed_bubble_plot(obj, wang_neutrophil_label,
                         "Wang Day 3 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "Neu3dpiBaseline_TargetBubbleWang_FixedScale.png")

print("Regenerating Lee Day 1 ARG1-split plots...")
obj <- cellchat_objects[["Lee_Day1_ARG1"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, neutrophil_states[1],
                         "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1pos (Fixed Scale)",
                         "1dpiLee_Target_Arg1pos_FixedScale.pdf")
create_fixed_circle_plot(obj, group_sizes, neutrophil_states[2],
                         "Lee Day 1 ARG1-split - Signals TO NeutrophilsArg1neg (Fixed Scale)",
                         "1dpiLee_Target_Arg1neg_FixedScale.pdf")
create_fixed_bubble_plot(obj, neutrophil_states,
                         "Lee Day 1 ARG1-split - Signals TO ARG1 Neutrophils (Fixed Scale)",
                         "Neu1dpiTargetBubbleLee_ARG1_FixedScale.png")

print("Regenerating Lee Day 3 ARG1-split plots...")
obj <- cellchat_objects[["Lee_Day3_ARG1"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, neutrophil_states[1],
                         "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1pos (Fixed Scale)",
                         "3dpiLee_Target_Arg1pos_FixedScale.pdf")
create_fixed_circle_plot(obj, group_sizes, neutrophil_states[2],
                         "Lee Day 3 ARG1-split - Signals TO NeutrophilsArg1neg (Fixed Scale)",
                         "3dpiLee_Target_Arg1neg_FixedScale.pdf")
create_fixed_bubble_plot(obj, neutrophil_states,
                         "Lee Day 3 ARG1-split - Signals TO ARG1 Neutrophils (Fixed Scale)",
                         "Neu3dpiTargetBubbleLee_ARG1_FixedScale.png")

print("Regenerating Wang Day 3 ARG1-split plots...")
obj <- cellchat_objects[["Wang_Day3_ARG1"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, neutrophil_states[1],
                         "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1pos (Fixed Scale)",
                         "3dpiWang_Target_Arg1pos_FixedScale.pdf")
create_fixed_circle_plot(obj, group_sizes, neutrophil_states[2],
                         "Wang Day 3 ARG1-split - Signals TO NeutrophilsArg1neg (Fixed Scale)",
                         "3dpiWang_Target_Arg1neg_FixedScale.pdf")
create_fixed_bubble_plot(obj, neutrophil_states,
                         "Wang Day 3 ARG1-split - Signals TO ARG1 Neutrophils (Fixed Scale)",
                         "Neu3dpiTargetBubbleWang_ARG1_FixedScale.png")

print("Regenerating Qin Day 1 Baseline plots...")
obj <- cellchat_objects[["Qin_Day1_Baseline"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, qin_neutrophil_label,
                         "Qin Day 1 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "1dpiQin_Baseline_Target_Neutrophils_FixedScale.pdf")
create_fixed_bubble_plot(obj, qin_neutrophil_label,
                         "Qin Day 1 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "Neu1dpiBaseline_TargetBubbleQin_FixedScale.png")

print("Regenerating Qin Day 3 Baseline plots...")
obj <- cellchat_objects[["Qin_Day3_Baseline"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, qin_neutrophil_label,
                         "Qin Day 3 Baseline - Signals TO Neutrophil (Fixed Scale)",
                         "3dpiQin_Baseline_Target_Neutrophils_FixedScale.pdf")

print("Regenerating Qin Day 1 ARG1-split plots...")
obj <- cellchat_objects[["Qin_Day1_ARG1"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, neutrophil_states[1],
                         "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1pos (Fixed Scale)",
                         "1dpiQin_Target_Arg1pos_FixedScale.pdf")
create_fixed_circle_plot(obj, group_sizes, neutrophil_states[2],
                         "Qin Day 1 ARG1-split - Signals TO NeutrophilsArg1neg (Fixed Scale)",
                         "1dpiQin_Target_Arg1neg_FixedScale.pdf")
create_fixed_bubble_plot(obj, neutrophil_states,
                         "Qin Day 1 ARG1-split - Signals TO ARG1 Neutrophils (Fixed Scale)",
                         "Neu1dpiTargetBubbleQin_ARG1_FixedScale.png")

print("Regenerating Qin Day 3 ARG1-split plots...")
obj <- cellchat_objects[["Qin_Day3_ARG1"]]
group_sizes <- as.numeric(table(obj@idents)[rownames(obj@net$count)])
group_sizes[is.na(group_sizes)] <- 0
names(group_sizes) <- rownames(obj@net$count)

create_fixed_circle_plot(obj, group_sizes, neutrophil_states[1],
                         "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1pos (Fixed Scale)",
                         "3dpiQin_Target_Arg1pos_FixedScale.pdf")
create_fixed_circle_plot(obj, group_sizes, neutrophil_states[2],
                         "Qin Day 3 ARG1-split - Signals TO NeutrophilsArg1neg (Fixed Scale)",
                         "3dpiQin_Target_Arg1neg_FixedScale.pdf")
create_fixed_bubble_plot(obj, neutrophil_states,
                         "Qin Day 3 ARG1-split - Signals TO ARG1 Neutrophils (Fixed Scale)",
                         "Neu3dpiTargetBubbleQin_ARG1_FixedScale.png")

save.image("environment_final_ARG1.RData")
print("✅ Saved: environment_final_ARG1.RData (Complete environment)")

print("")
print("📊 SUMMARY:")
print("✅ Analysis 1: Lee Day 1 Baseline → cellchatOne_baseline_complete.rds")
print("✅ Analysis 2: Lee Day 3 Baseline → cellchatThreeLee_baseline_complete.rds")
print("✅ Analysis 3: Wang Day 3 Baseline → cellchatThreeWang_baseline_complete.rds")
print("✅ Analysis 4: Qin Day 1 Baseline → cellchatOneQin_baseline_complete.rds")
print("✅ Analysis 5: Qin Day 3 Baseline → cellchatThreeQin_baseline_complete.rds")
print("✅ Analysis 6: Lee Day 1 ARG1-split → cellchatOne_ARG1_complete.rds")
print("✅ Analysis 7: Lee Day 3 ARG1-split → cellchatThreeLee_ARG1_complete.rds")
print("✅ Analysis 8: Wang Day 3 ARG1-split → cellchatThreeWang_ARG1_complete.rds")
print("✅ Analysis 9: Qin Day 1 ARG1-split → cellchatOneQin_ARG1_complete.rds")
print("✅ Analysis 10: Qin Day 3 ARG1-split → cellchatThreeQin_ARG1_complete.rds")
print("")
print("🎯 MIXED SCALING STRATEGY APPLIED:")
print("✅ Circle plots: Auto-scaling for optimal visibility")
print("✅ Bubble plots: Fixed color scale (white-red) for comparison")
print("✅ Heatmaps: Fixed color scale for comparison")
print("✅ Files saved with '_FixedScale' suffix for non-circle plots")
print("")
print("All 10 RDS files contain complete processing (centrality + pattern analysis)")
print("Ready for visualization without preprocessing!")
print("")
print("📈 DATASETS INCLUDED:")
print("  • Lee Day 1 & Day 3")
print("  • Wang Day 3")
print("  • Qin Day 1 & Day 3 (NEW!)")
print("")
print("🧬 ANALYSIS TYPES:")
print("  • 5 Baseline analyses (original neutrophils)")
print("  • 5 ARG1-split analyses (ARG1+ vs ARG1- neutrophils)")
