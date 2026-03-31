#!/usr/bin/env Rscript
# ============================================================================
# MASTER COMPREHENSIVE POST-SCI CELL COMMUNICATION PIPELINE
# ============================================================================
# Research question: What signals do Arg1+ neutrophils receive that may associate with
#   anti-inflammatory activation, and what signals do Arg1- neutrophils receive that may associate
#   with pro-inflammatory activation? (Working hypotheses; not causal claims.)
# Datasets: Lee (Day 1, Day 3) + Wang (Day 3) + Qin (Day 1, Day 3)
# Analysis Tools: LIANA, BARTsc, decoupleR+PROGENy
#
# VISUALIZATION AUDIT (aligned with official tool docs):
#   LIANA (saezlab.github.io/liana): liana_dotplot, heat_freq, chord_freq, liana_heatmap (mat from trunc source x target)
#     Custom LR prior: Lee D1/D3, Wang D3, Qin D1/D3 all use resource=custom + external_resource when LIANA_USE_CUSTOM_LR and cohort CellChat/CellCall RDS load; else MouseConsensus with an explicit console note.
#     - Filter by aggregate_rank <= 0.01 only (no target filter) to avoid heat_freq annotation error
#     - chord_freq: source/target converted to character so labels show cell-type names (not factor 1,2,3); source_groups/target_groups = unique(source)/unique(target) after conversion; PNG 1400x1400
#   PROGENy+decoupleR: pathway activity via decoupleR::run_wmean + model_mouse_full (not progeny() wrapper); bar + RdBu heatmap; ligand→pathway dot map is exploratory (ligand genes vs footprint targets), not canonical receptor-centric inference.
#   BARTsc: analysis blocks save *BARTsc_Object.rds + CSVs; separate per-cohort "BARTsc visualizations" sections reload RDS (optional) and call BARTsc::dot_plot, deviation_heatmap, key_regulator_scatter_unified into BARTsc_* folders (no re-run of bartsc/crossCT).
#   Integration (per cohort after BARTsc): LIANA x PROGENy footprint co-membership (exploratory CSV + sender x pathway pheatmap) + OmniPath receptor→TF (primary mechanistic chain); Lee D1/D3, Wang D3, Qin D1/D3
#
# OFFICIAL PIPELINES (follow these for correctness):
#   LIANA:   https://saezlab.github.io/liana/  (liana_wrap, liana_aggregate, resource, external_resource)
#   PROGENy: https://saezlab.github.io/progeny/ + decoupleR for pathway activity (progeny(), model_mouse_full)
#   BARTsc:  https://github.com/hongpan-uva/BARTsc (install + initialize/load_bart2, then bartsc pipeline)
#
# STRUCTURE: Linear/flattened per cohort — LIANA/PROGENy ggplot objects built in each section; each comparison figure is printed immediately after it is built using UNIFY_* limits pooled from every cohort object that exists() at that line (so scales tighten as later cohorts run; full cross-cohort limits once the pipeline reaches Qin D3). Then CellChat/CellCall (Lee D3), VlnPlots, BARTsc, integration.
# FILE INTEGRITY: The saved .r on disk is authoritative. Chat or email excerpts can break mid-line (e.g. lee_day1_arg1pos_cells_in_mat) without the real file being truncated.
# LINE-BY-LINE RUN: Use source("MASTER_POST_COMMUNICATION_PIPELINE.r", echo = TRUE) so R prints each line as it executes (RStudio: run with echo, or step through with Cursor line run).
# STANDALONE SECTIONS: Lee Day 3 LIANA begins with a preamble (checkpoint and/or LeeDat.rds + defaults + libraries). Other days/sections still expect prior blocks unless extended the same way.
# Qin Day 3: After QinDay3_LIANA_Results.rds the script continues in cohort order: LIANA plots (dotplot, heat_freq, chord_freq PNG, liana_heatmap, …) → top-receptor VlnPlot → PROGENy → BARTsc analysis (CSVs + *BARTsc_Object.rds) → BARTsc visualizations section → LIANA–PROGENy–BARTsc integration (same pattern as Lee Day 1).
# BARTsc figures only: each cohort has a "BARTsc visualizations" block after its analysis block; load *BARTsc_Object.rds + BARTsc::load_bart2() (see comment at each block).
#
# RUN ORDER (recommended):
#   1) CellChat_ARG1_Final.R
#   2) CellChat_ARG1_FinalQin.R
#   3) CELLCALL_NEUTROPHIL_6_ANALYSES Linear.R
#   4) CELLCALL_QIN_ONLY_ANALYSES.R
#   5) Upstream CellChat/CellCall .rds paths default to D:/Active Analysis/Mustafa/post-communication/ (overridable via CELLCOMM_POST_RDS_DIR).
#   6) Working directory for LeeDat/WangDat/QinDat defaults to D:/Mustafa/cellcall (overridable via CELLCOMM_PROJECT_ROOT).
# Rerun from: line 1 (full pipeline) or line 42 (SECTION 0) / line 70 (library(Seurat)). Resume blocks are documented inside.
# ============================================================================

# ============================================================================
# SECTION 0: INSTALL OPTIONAL LIANA METHOD PACKAGES (if missing)
# ============================================================================
# LIANA methods call_cellchat, call_italk, cytotalk require these packages.
# Install once; then load is optional (LIANA calls them when method is used).
# Official install commands (run in R if a package is missing):
#   CellChat:  devtools::install_github("sqjin/CellChat")   # or jinworks/CellChat
#   iTALK:    devtools::install_github("Coolgenome/iTALK")
#   CytoTalk: devtools::install_github("huBioinfo/CytoTalk")  # or tanlabcode/CytoTalk
if (!requireNamespace("CellChat", quietly = TRUE)) {
  message("CellChat not installed. LIANA method 'call_cellchat' will be skipped. Install: devtools::install_github('sqjin/CellChat')")
}
if (!requireNamespace("iTALK", quietly = TRUE)) {
  message("iTALK not installed. LIANA method 'call_italk' will be skipped. Install: devtools::install_github('Coolgenome/iTALK')")
}
if (!requireNamespace("CytoTalk", quietly = TRUE)) {
  message("CytoTalk not installed. LIANA method 'cytotalk' will be skipped. Install: devtools::install_github('huBioinfo/CytoTalk')")
}
if (!requireNamespace("OmnipathR", quietly = TRUE)) {
  message("OmnipathR not installed. Receptor-TF integration will be limited. Install: BiocManager::install('OmnipathR')")
}

# ============================================================================
# SECTION 1: LOAD LIBRARIES
# ============================================================================
# NOTE: Do NOT load dplyr/tidyr explicitly - they are imported by Seurat, liana,
# etc. Explicit library(dplyr) triggers unload/reload and fails when dplyr is
# already in use. They will be available after loading Seurat.
library(Seurat)
library(ggplot2)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(grid)
library(ggrepel)
library(liana)
if (requireNamespace("OmnipathR", quietly = TRUE)) library(OmnipathR)
library(decoupleR)
library(nichenetr)
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(progeny)
library(SingleCellExperiment)
if (requireNamespace("scuttle", quietly = TRUE)) library(scuttle)
library(GENIE3)
library(AUCell)
library(msigdbr)
library(igraph)
library(ggalluvial)
library(viridis)
library(circlize)
if (requireNamespace("BARTsc", quietly = TRUE)) {
  library(BARTsc)
}
# scCustomize: optional - used to convert Assay5 to V3 for LIANA compatibility
if (requireNamespace("scCustomize", quietly = TRUE)) {
  library(scCustomize)
}

# ============================================================================
# SECTION 2: CONFIGURE GLOBAL VARIABLES
# ============================================================================

# Define neutrophil identifiers (aliases = alternate labels in datasets so Arg1+/- split applies to all neutrophils)
NEUTROPHIL_BASE_LABEL <- "Neutrophil"
NEUTROPHIL_ALIASES <- c("cNeutrophil", "Neutrophils")
GENE_CANDIDATES <- c("Arg1", "ARG1")
# LIANA methods: use internal methods only. call_cellchat/call_italk cause duplicate row.names in liana_aggregate when mixed with internal methods (LIANA warns "Using internal and external methods should be done with caution!"). Your CellChat/CellCall L-R prior is already used via external_resource when LIANA_USE_CUSTOM_LR = TRUE.
LIANA_METHODS <- c("natmi", "connectome", "logfc", "sca", "cellphonedb")
if (requireNamespace("CytoTalk", quietly = TRUE)) {
  LIANA_METHODS <- c(LIANA_METHODS, "cytotalk")
}
NEUTROPHIL_STATES <- c(
  paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos"),
  paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
)
# LIANA liana_wrap: min cells per group (official default often 5; 3 is safer when Arg1+ neutrophil n is small).
LIANA_MIN_CELLS <- 3L
# Top single-gene LIANA receptors to plot with Seurat::VlnPlot per cohort (after LIANA).
LIANA_TOP_RECEPTOR_VLN <- 5L

# Windows defaults (your analysis machine). Override: Sys.setenv(CELLCOMM_PROJECT_ROOT = "D:/path") for Seurat .rds folder;
# Sys.setenv(CELLCOMM_POST_RDS_DIR = "D:/path") for CellChat/CellCall upstream .rds folder.
POST_COMM_RDS_DIR <- Sys.getenv("CELLCOMM_POST_RDS_DIR", unset = "D:/Active Analysis/Mustafa/post-communication")
# Create output directory (relative to getwd() when you source/run; setwd() to project root first)
OUTPUT_DIR <- "MASTER_PIPELINE_RESULTS"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
print(paste0("OUTPUT_DIR: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), " | getwd(): ", normalizePath(getwd(), winslash = "/", mustWork = FALSE), " | POST_COMM_RDS_DIR: ", POST_COMM_RDS_DIR))
# Lee Day 3 LIANA only: set TRUE to load lee_day3 + labels from a prior run (needs Checkpoint_AfterLeeDay3_LIANA.rds under OUTPUT_DIR). First run: keep FALSE and execute Lee Day 3 block (~252–270) after lee_dat exists.
LEE_DAY3_LOAD_LIANA_CHECKPOINT <- FALSE

# Custom LR resource from CellChat/CellCall - use MIF-CD74, TNF-TNFRSF1B style names instead of F7-F10-F3
# DESIGN: We merge only LR *identities* (ligand-receptor pairs) from CellChat and CellCall into a union,
#   deduplicated list. We do NOT merge their scores/probabilities. LIANA then scores this custom set
#   uniformly (logfc, magnitude_rank, etc.). CellChat/CellCall define the candidate LR universe;
#   LIANA provides the harmonized scoring and Arg1+ vs Arg1- comparison (Methods: "biologically informed LR prior").
# CellChat: cellchatOne_ARG1_complete.rds (Day 1), cellchatThreeLee_ARG1_complete.rds (Day 3) - from CellChat_ARG1_Final.R
# CellCall: CellCall_LEE_Day1_GLOBAL.rds, CellCall_LEE_Day3_GLOBAL.rds - from CELLCALL_NEUTROPHIL_6_ANALYSES Linear.R
# TRUE = use CellChat+CellCall L-R pairs as external_resource; FALSE = MouseConsensus only
LIANA_USE_CUSTOM_LR <- TRUE
# Upstream outputs: CellChat_ARG1_Final.R / CellChat_ARG1_FinalQin.R; CellCall: CELLCALL_NEUTROPHIL_6_ANALYSES / CELLCALL_QIN_ONLY_ANALYSES.R (GLOBAL or merged LR universe).
CELLCHAT_LEE_DAY1_RDS <- file.path(POST_COMM_RDS_DIR, "cellchatOne_ARG1_complete.rds")
CELLCALL_LEE_DAY1_RDS <- file.path(POST_COMM_RDS_DIR, "CellCall_LEE_Day1_GLOBAL.rds")
CELLCHAT_LEE_DAY3_RDS <- file.path(POST_COMM_RDS_DIR, "cellchatThreeLee_ARG1_complete.rds")
CELLCALL_LEE_DAY3_RDS <- file.path(POST_COMM_RDS_DIR, "CellCall_LEE_Day3_GLOBAL.rds")
CELLCHAT_WANG_DAY3_RDS <- file.path(POST_COMM_RDS_DIR, "cellchatThreeWang_ARG1_complete.rds")
CELLCALL_WANG_DAY3_RDS <- file.path(POST_COMM_RDS_DIR, "CellCall_WANG_Day3_GLOBAL_FULL.rds")
CELLCHAT_QIN_DAY1_RDS <- file.path(POST_COMM_RDS_DIR, "cellchatOneQin_ARG1_complete.rds")
CELLCALL_QIN_DAY1_RDS <- file.path(POST_COMM_RDS_DIR, "CellCall_QIN_Day1_GLOBAL.rds")
CELLCHAT_QIN_DAY3_RDS <- file.path(POST_COMM_RDS_DIR, "cellchatThreeQin_ARG1_complete.rds")
CELLCALL_QIN_DAY3_RDS <- file.path(POST_COMM_RDS_DIR, "CellCall_QIN_Day3_GLOBAL.rds")

# Set plot options. Plots draw inline in RStudio (no png/pdf/ggsave for most ggplots). Do not run dev.off() in the console.
options(repr.plot.width = 8, repr.plot.height = 6)
if (interactive() && grDevices::dev.cur() == 1L) { grDevices::dev.new(); message("Graphics device opened for plots (current device ", grDevices::dev.cur(), "). Check Plots pane or any new window.") }
PLOT_TITLE_THEME <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 12, margin = margin(t = 10)),
  axis.title.y = element_text(size = 12, margin = margin(r = 10)),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  plot.margin = margin(10, 10, 20, 10)
)

print("✓ Global variables configured")

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# SECTION 3: LOAD DATA
# ============================================================================

print(">>> LOADING DATA OBJECTS <<<")

# Default working directory for LeeDat.rds / WangDat.rds / QinDat.rds (Windows). Override with Sys.setenv(CELLCOMM_PROJECT_ROOT = "D:/your/path").
CELLCOMM_PROJECT_ROOT <- Sys.getenv("CELLCOMM_PROJECT_ROOT", unset = "")
if (nzchar(CELLCOMM_PROJECT_ROOT) && dir.exists(CELLCOMM_PROJECT_ROOT)) {
  setwd(CELLCOMM_PROJECT_ROOT)
} else {
  setwd("D:/Mustafa/cellcall")
}

# Fix OmnipathR logging directory issue (required for PROGENy)
current_wd_omnipath <- tryCatch(getwd(), error = function(e) "D:/Mustafa/cellcall")
log_dir_path_omnipath <- tryCatch(file.path(current_wd_omnipath, "omnipathr-log"), error = function(e) "D:/Mustafa/cellcall/omnipathr-log")
dir_exists_omnipath <- tryCatch(dir.exists(log_dir_path_omnipath), error = function(e) FALSE)
dir_created_omnipath <- tryCatch(if (!dir_exists_omnipath) dir.create(log_dir_path_omnipath, recursive = TRUE) else TRUE, error = function(e) FALSE)
Sys.setenv(OMNIPATHR_LOG_DIR = log_dir_path_omnipath)
# Set OmnipathR option to use current directory (ignore old F:/ path)
tryCatch(options(omnipathr.log_dir = log_dir_path_omnipath), error = function(e) NULL)
tryCatch(options(OmnipathR.log_dir = log_dir_path_omnipath), error = function(e) NULL)
tryCatch(options(omnipath.log_dir = log_dir_path_omnipath), error = function(e) NULL)

# Validate input files exist (linear execution)
lee_dat_file_exists <- file.exists("LeeDat.rds")
wang_dat_file_exists <- file.exists("WangDat.rds")
qin_dat_file_exists <- file.exists("QinDat.rds") | file.exists("qindat.rds")
stopifnot(lee_dat_file_exists)
stopifnot(wang_dat_file_exists)
stopifnot(qin_dat_file_exists)
print("✓ All required input files found")

# Load Lee dataset
lee_dat <- readRDS("LeeDat.rds")
DefaultAssay(lee_dat) <- "RNA"

# Load Wang dataset
wang_dat <- readRDS("WangDat.rds")
DefaultAssay(wang_dat) <- "RNA"

# Load Qin dataset (QinDat.rds or qindat.rds)
qin_dat <- readRDS(if (file.exists("QinDat.rds")) "QinDat.rds" else "qindat.rds")
DefaultAssay(qin_dat) <- "RNA"

# Prepare time variables
lee_dat$time <- as.character(lee_dat$time)
lee_dat$time <- gsub("Uninjured", "0", lee_dat$time)
lee_dat$time <- gsub("1dpi", "1", lee_dat$time)
lee_dat$time <- gsub("3dpi", "3", lee_dat$time)
lee_dat$time <- gsub("7dpi", "7", lee_dat$time)
lee_dat$time <- as.numeric(lee_dat$time)

wang_dat$time <- as.character(wang_dat$time)
wang_dat$time <- gsub("Uninjured", "0", wang_dat$time)
wang_dat$time <- gsub("1dpi", "1", wang_dat$time)
wang_dat$time <- gsub("3dpi", "3", wang_dat$time)
wang_dat$time <- gsub("7dpi", "7", wang_dat$time)
wang_dat$time <- as.numeric(wang_dat$time)

qin_dat$time <- as.character(qin_dat$time)
qin_dat$time <- gsub("Uninjured", "0", qin_dat$time)
qin_dat$time <- gsub("1dpi|Day1", "1", qin_dat$time)
qin_dat$time <- gsub("3dpi|Day3", "3", qin_dat$time)
qin_dat$time <- gsub("7dpi|Day7", "7", qin_dat$time)
qin_dat$time <- gsub("14dpi|Day14", "14", qin_dat$time)
qin_dat$time <- gsub("28dpi|Day28", "28", qin_dat$time)
qin_dat$time <- as.numeric(qin_dat$time)

print("✓ Data loaded successfully")

# ============================================================================
# SECTION 4: STEP 0 - LIANA CONSENSUS ANALYSIS
# ============================================================================
# Official LIANA: https://saezlab.github.io/liana/
#   - liana_wrap(sce, method, resource, idents_col, external_resource) then liana_aggregate() for consensus ranks.
#   - Custom L-R: use resource = "custom" and external_resource in OmniPath format (source_genesymbol, target_genesymbol).
#   - We use LR pairs from your pre-run CellChat and CellCall results (loaded from RDS) as the custom resource when LIANA_USE_CUSTOM_LR = TRUE.
#   - Visuals: liana_dotplot(source_groups, target_groups, ntop); heat_freq/chord_freq on filter(aggregate_rank <= 0.01).
# ============================================================================

print(">>> STEP 0: LIANA CONSENSUS CELL-CELL COMMUNICATION ANALYSIS <<<")

# -------- Define Arg1 Status for All Datasets (ONCE, REUSED EVERYWHERE) --------
print("--- Defining Arg1 Status for All Datasets ---")

# Lee Day 1: Create object and define Arg1 status (before any LIANA or PROGENy)
# Arg1+ vs Arg1- defined within neutrophils only (biologically correct)
lee_day1 <- subset(lee_dat, subset = time == 1)
DefaultAssay(lee_day1) <- "RNA"
lee_day1_pruned_labels <- as.character(lee_day1$celltype)
lee_day1_pruned_labels <- gsub("-", "", lee_day1_pruned_labels)
lee_day1_pruned_labels[lee_day1_pruned_labels %in% NEUTROPHIL_ALIASES] <- NEUTROPHIL_BASE_LABEL
lee_day1_neut <- which(lee_day1_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL))
lee_day1_gene_arg1 <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(lee_day1)][1]
stopifnot(!is.na(lee_day1_gene_arg1) && nchar(lee_day1_gene_arg1) > 0)
lee_day1_expr_arg1 <- as.numeric(Seurat::GetAssayData(lee_day1, layer = "data")[lee_day1_gene_arg1, lee_day1_neut, drop = FALSE])
lee_day1_indicator <- lee_day1_expr_arg1 > 0
lee_day1_indicator[is.na(lee_day1_indicator)] <- FALSE
lee_day1_pos <- lee_day1_neut[lee_day1_indicator]
lee_day1_neg <- setdiff(lee_day1_neut, lee_day1_pos)
lee_day1_arg1_status <- rep("Arg1neg", ncol(lee_day1))
lee_day1_arg1_status[lee_day1_pos] <- "Arg1pos"
lee_day1$arg1_status <- lee_day1_arg1_status
Idents(lee_day1) <- factor(lee_day1_pruned_labels)

# Lee Day 3: Create object and define Arg1 status (same pattern as Lee Day 1)
lee_day3 <- subset(lee_dat, subset = time == 3)
DefaultAssay(lee_day3) <- "RNA"
lee_day3_pruned_labels <- as.character(lee_day3$celltype)
lee_day3_pruned_labels <- gsub("-", "", lee_day3_pruned_labels)
lee_day3_pruned_labels[lee_day3_pruned_labels %in% NEUTROPHIL_ALIASES] <- NEUTROPHIL_BASE_LABEL
lee_day3_neut <- which(lee_day3_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL))
lee_day3_gene_arg1 <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(lee_day3)][1]
stopifnot(!is.na(lee_day3_gene_arg1) && nchar(lee_day3_gene_arg1) > 0)
lee_day3_expr_arg1 <- as.numeric(Seurat::GetAssayData(lee_day3, layer = "data")[lee_day3_gene_arg1, lee_day3_neut, drop = FALSE])
lee_day3_indicator <- lee_day3_expr_arg1 > 0
lee_day3_indicator[is.na(lee_day3_indicator)] <- FALSE
lee_day3_pos <- lee_day3_neut[lee_day3_indicator]
lee_day3_neg <- setdiff(lee_day3_neut, lee_day3_pos)
lee_day3_arg1_status <- rep("Arg1neg", ncol(lee_day3))
lee_day3_arg1_status[lee_day3_pos] <- "Arg1pos"
lee_day3$arg1_status <- lee_day3_arg1_status
Idents(lee_day3) <- factor(lee_day3_pruned_labels)
lee_day3_neut_cells <- colnames(lee_day3)[c(lee_day3_pos, lee_day3_neg)]

# Wang Day 3: Create object and define Arg1 status (same label source pattern as Lee Day 1)
wang_day3 <- subset(wang_dat, time == 3)
wang_day3_label_col <- if ("celltype" %in% colnames(wang_day3@meta.data)) wang_day3$celltype else wang_day3$pruned_labels
wang_day3_pruned_labels <- as.character(wang_day3_label_col)
wang_day3_pruned_labels <- gsub("-", "", wang_day3_pruned_labels)
wang_day3_pruned_labels[wang_day3_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL)] <- NEUTROPHIL_BASE_LABEL
wang_day3_pruned_labels[is.na(wang_day3_pruned_labels) | wang_day3_pruned_labels == ""] <- "Other"
wang_day3_neut <- which(wang_day3_pruned_labels == NEUTROPHIL_BASE_LABEL)
wang_day3_gene_candidates <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(wang_day3)]
wang_day3_gene <- head(wang_day3_gene_candidates, 1)
wang_day3_expr <- tryCatch(Seurat::GetAssayData(wang_day3, layer = "data")[wang_day3_gene, wang_day3_neut, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
wang_day3_expr_vec <- as.numeric(wang_day3_expr)
wang_day3_indicator <- wang_day3_expr_vec > 0
wang_day3_indicator[is.na(wang_day3_indicator)] <- FALSE
wang_day3_pos <- wang_day3_neut[wang_day3_indicator]
wang_day3_neg <- setdiff(wang_day3_neut, wang_day3_pos)
wang_day3_arg1_status <- rep("Arg1neg", ncol(wang_day3))
wang_day3_arg1_status[wang_day3_pos] <- "Arg1pos"
wang_day3$arg1_status <- wang_day3_arg1_status
Idents(wang_day3) <- factor(wang_day3_pruned_labels)

# Qin Day 1: Create object and define Arg1 status (same label pattern as Lee Day 1)
qin_day1 <- subset(qin_dat, subset = time == 1)
DefaultAssay(qin_day1) <- "RNA"
qin_day1_label_col <- if ("celltype" %in% colnames(qin_day1@meta.data)) qin_day1$celltype else qin_day1$pruned_labels
qin_day1_pruned_labels <- as.character(qin_day1_label_col)
qin_day1_pruned_labels <- gsub("-", "", qin_day1_pruned_labels)
qin_day1_pruned_labels[qin_day1_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL)] <- NEUTROPHIL_BASE_LABEL
qin_day1_pruned_labels[is.na(qin_day1_pruned_labels) | qin_day1_pruned_labels == ""] <- "Other"
qin_day1_neut <- which(qin_day1_pruned_labels == NEUTROPHIL_BASE_LABEL)
qin_day1_gene_arg1 <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(qin_day1)][1]
stopifnot(!is.na(qin_day1_gene_arg1) && nchar(qin_day1_gene_arg1) > 0)
qin_day1_expr_arg1 <- as.numeric(Seurat::GetAssayData(qin_day1, layer = "data")[qin_day1_gene_arg1, qin_day1_neut, drop = FALSE])
qin_day1_indicator <- qin_day1_expr_arg1 > 0
qin_day1_indicator[is.na(qin_day1_indicator)] <- FALSE
qin_day1_pos <- qin_day1_neut[qin_day1_indicator]
qin_day1_neg <- setdiff(qin_day1_neut, qin_day1_pos)
qin_day1_arg1_status <- rep("Arg1neg", ncol(qin_day1))
qin_day1_arg1_status[qin_day1_pos] <- "Arg1pos"
qin_day1$arg1_status <- qin_day1_arg1_status
Idents(qin_day1) <- factor(qin_day1_pruned_labels)

# Qin Day 3: Create object and define Arg1 status (same label pattern as Lee Day 1)
qin_day3 <- subset(qin_dat, subset = time == 3)
DefaultAssay(qin_day3) <- "RNA"
qin_day3_label_col <- if ("celltype" %in% colnames(qin_day3@meta.data)) qin_day3$celltype else qin_day3$pruned_labels
qin_day3_pruned_labels <- as.character(qin_day3_label_col)
qin_day3_pruned_labels <- gsub("-", "", qin_day3_pruned_labels)
qin_day3_pruned_labels[qin_day3_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL)] <- NEUTROPHIL_BASE_LABEL
qin_day3_pruned_labels[is.na(qin_day3_pruned_labels) | qin_day3_pruned_labels == ""] <- "Other"
qin_day3_neut <- which(qin_day3_pruned_labels == NEUTROPHIL_BASE_LABEL)
qin_day3_gene_arg1 <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(qin_day3)][1]
stopifnot(!is.na(qin_day3_gene_arg1) && nchar(qin_day3_gene_arg1) > 0)
qin_day3_expr_arg1 <- as.numeric(Seurat::GetAssayData(qin_day3, layer = "data")[qin_day3_gene_arg1, qin_day3_neut, drop = FALSE])
qin_day3_indicator <- qin_day3_expr_arg1 > 0
qin_day3_indicator[is.na(qin_day3_indicator)] <- FALSE
qin_day3_pos <- qin_day3_neut[qin_day3_indicator]
qin_day3_neg <- setdiff(qin_day3_neut, qin_day3_pos)
qin_day3_arg1_status <- rep("Arg1neg", ncol(qin_day3))
qin_day3_arg1_status[qin_day3_pos] <- "Arg1pos"
qin_day3$arg1_status <- qin_day3_arg1_status
Idents(qin_day3) <- factor(qin_day3_pruned_labels)

print("✓ Arg1 status defined for all datasets")

# -------- Lee Day 1: LIANA → LIANA plots → CellChat/CellCall overlap → top-receptor Vln → PROGENy (+ plots) → BARTsc → integration (one cohort, top-to-bottom) --------
# LIANA: official liana_wrap -> liana_aggregate (https://saezlab.github.io/liana/). Uses CellChat + CellCall RDS as custom resource when LIANA_USE_CUSTOM_LR = TRUE.
print("--- LIANA Analysis: Lee Day 1 ---")

## 1) Build LIANA labels (NeutrophilArg1pos / NeutrophilArg1neg / Other)
lee_day1_liana_labels <- as.character(lee_day1_pruned_labels)
lee_day1_liana_labels[lee_day1_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
lee_day1_liana_labels[lee_day1_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
lee_day1_liana_labels[is.na(lee_day1_liana_labels) | lee_day1_liana_labels == ""] <- "Other"

lee_day1$liana_label <- factor(lee_day1_liana_labels)
Seurat::Idents(lee_day1) <- lee_day1$liana_label

## 2) Clean Seurat for LIANA: RNA only, drop other assays/reductions/graphs (per official workflow)
SeuratObject::DefaultAssay(lee_day1) <- "RNA"
for (assay_name in names(lee_day1@assays)) {
  if (assay_name != "RNA") lee_day1[[assay_name]] <- NULL
}
for (red_name in names(lee_day1@reductions)) lee_day1[[red_name]] <- NULL
for (graph_name in names(lee_day1@graphs)) lee_day1[[graph_name]] <- NULL
## 2b) Build custom L-R resource from CellChat and CellCall (Lee Day 1). Linear: one check per step; errors on exact line.
##     LIANA expects "_" as complex subunit separator. CellCall: split on first "-"; receptor subunits join with "_".
lee_day1_liana_resource <- "MouseConsensus"
lee_day1_liana_external <- NULL
custom_lr_list_d1 <- list()
path_cc_d1 <- CELLCHAT_LEE_DAY1_RDS
path_ccall_d1 <- CELLCALL_LEE_DAY1_RDS
cc_d1_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_cc_d1) && nzchar(trimws(path_cc_d1)) && file.exists(path_cc_d1) && requireNamespace("CellChat", quietly = TRUE)
cc_d1 <- NULL
if (cc_d1_ok) cc_d1 <- readRDS(path_cc_d1)
comm_d1 <- NULL
if (!is.null(cc_d1) && inherits(cc_d1, "CellChat")) comm_d1 <- CellChat::subsetCommunication(cc_d1)
has_comm_d1 <- !is.null(comm_d1) && nrow(comm_d1) > 0 && "ligand" %in% colnames(comm_d1) && "receptor" %in% colnames(comm_d1)
if (has_comm_d1) {
  cc_lr_d1 <- data.frame(source_genesymbol = comm_d1$ligand, target_genesymbol = comm_d1$receptor, stringsAsFactors = FALSE)
  cc_lr_d1 <- cc_lr_d1[!duplicated(cc_lr_d1), ]
  custom_lr_list_d1[["CellChat"]] <- cc_lr_d1
}
ccall_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_ccall_d1) && nzchar(trimws(path_ccall_d1)) && file.exists(path_ccall_d1)
ccall_d1 <- NULL
if (ccall_ok) ccall_d1 <- readRDS(path_ccall_d1)
has_ccall_lr <- !is.null(ccall_d1) && !is.null(ccall_d1@data$expr_l_r_log2_scale)
if (has_ccall_lr) {
  lr_rownames_d1 <- rownames(ccall_d1@data$expr_l_r_log2_scale)
  ccall_lig_d1 <- sub("-.*", "", lr_rownames_d1)
  ccall_rec_d1 <- sub("^[^-]+-", "", lr_rownames_d1)
  ccall_rec_d1 <- gsub("-", "_", ccall_rec_d1)
  ccall_lr_d1 <- data.frame(source_genesymbol = ccall_lig_d1, target_genesymbol = ccall_rec_d1, stringsAsFactors = FALSE)
  ccall_lr_d1 <- ccall_lr_d1[nzchar(ccall_lr_d1$target_genesymbol), ]
  ccall_lr_d1 <- unique(ccall_lr_d1)
  custom_lr_list_d1[["CellCall"]] <- ccall_lr_d1
}
has_custom_lr <- length(custom_lr_list_d1) > 0
if (has_custom_lr) {
  custom_lr_d1 <- do.call(rbind, custom_lr_list_d1)
  custom_lr_d1 <- custom_lr_d1[!duplicated(custom_lr_d1[, c("source_genesymbol", "target_genesymbol")]), ]
  lee_day1_liana_external <- data.frame(source_genesymbol = custom_lr_d1$source_genesymbol, target_genesymbol = custom_lr_d1$target_genesymbol, stringsAsFactors = FALSE)
  lee_day1_liana_resource <- "custom"
  print(paste("Loaded", nrow(lee_day1_liana_external), "custom L-R pairs from CellChat + CellCall (Lee Day 1)"))
}
## 3) Build clean SCE to bypass validObject() bug; use Seurat data slot as logcounts + base = exp(1) for parity with Seurat path.
lee_day1_counts <- tryCatch(Seurat::GetAssayData(lee_day1, slot = "counts"), error = function(e) Seurat::GetAssayData(lee_day1, layer = "counts"))
lee_day1_logcounts <- tryCatch(Seurat::GetAssayData(lee_day1, slot = "data"), error = function(e) Seurat::GetAssayData(lee_day1, layer = "data"))
lee_day1_sce <- SingleCellExperiment(assays = list(counts = lee_day1_counts, logcounts = lee_day1_logcounts))
lee_day1_sce$liana_label <- lee_day1$liana_label
SingleCellExperiment::colLabels(lee_day1_sce) <- lee_day1_sce$liana_label
liana_args_ld1 <- list(sce = lee_day1_sce, method = LIANA_METHODS, resource = lee_day1_liana_resource, idents_col = "liana_label", expr_prop = 0.05, verbose = TRUE, min_cells = LIANA_MIN_CELLS, base = exp(1))
if (!is.null(lee_day1_liana_external)) liana_args_ld1$external_resource <- lee_day1_liana_external
lee_day1_liana_result <- do.call(liana::liana_wrap, liana_args_ld1)
lee_day1_liana_result_df <- liana::liana_aggregate(lee_day1_liana_result)

## 4) Downstream processing (unchanged logic)
lee_day1_arg1pos_received <- dplyr::filter(lee_day1_liana_result_df, target == "NeutrophilArg1pos")
lee_day1_arg1pos_received <- dplyr::arrange(lee_day1_arg1pos_received, aggregate_rank)
lee_day1_arg1pos_top_signals <- head(lee_day1_arg1pos_received, 50)

lee_day1_arg1neg_received <- dplyr::filter(lee_day1_liana_result_df, target == "NeutrophilArg1neg")
lee_day1_arg1neg_received <- dplyr::arrange(lee_day1_arg1neg_received, aggregate_rank)

topN <- 50
lee_day1_arg1pos_top <- head(lee_day1_arg1pos_received, topN)
lee_day1_arg1neg_top <- head(lee_day1_arg1neg_received, topN)
pos_pairs <- paste0(lee_day1_arg1pos_top$ligand.complex, "_", lee_day1_arg1pos_top$receptor.complex)
neg_pairs_all <- paste0(lee_day1_arg1neg_received$ligand.complex, "_", lee_day1_arg1neg_received$receptor.complex)
neg_rank_lookup <- setNames(lee_day1_arg1neg_received$aggregate_rank, neg_pairs_all)
neg_ranks_matched <- neg_rank_lookup[pos_pairs]
neg_ranks_matched[is.na(neg_ranks_matched)] <- 1.0
spec_index <- (neg_ranks_matched - lee_day1_arg1pos_top$aggregate_rank) / (neg_ranks_matched + lee_day1_arg1pos_top$aggregate_rank + 1e-10)
lee_day1_arg1pos_top$specificity_index <- spec_index
lee_day1_arg1pos_specific <- dplyr::filter(lee_day1_arg1pos_top, specificity_index > 0.3 | !(pos_pairs %in% neg_pairs_all))

lee_day1_liana_arg1_pos <- lee_day1_arg1pos_received
lee_day1_liana_arg1_neg <- lee_day1_arg1neg_received
lee_day1_liana_consensus <- head(lee_day1_arg1pos_received, 30)
lee_day1_liana_arg1pos_topranked <- head(lee_day1_arg1pos_received, 20)
lee_day1_liana_arg1neg_topranked <- head(lee_day1_arg1neg_received, 20)

print(paste("Arg1pos received signals (all):", nrow(lee_day1_arg1pos_received)))
print(paste("Arg1pos top signals (top 50):", nrow(lee_day1_arg1pos_top_signals)))
print(paste("Arg1pos-specific signals:", nrow(lee_day1_arg1pos_specific)))

lr_pairs_lee_d1 <- dplyr::distinct(lee_day1_liana_result_df, ligand.complex, receptor.complex, .keep_all = FALSE)
print("LIANA ligand.complex and receptor.complex columns - unique L-R pairs in full result:")
print(head(lr_pairs_lee_d1, 30))
all_sym <- unique(c(lr_pairs_lee_d1$ligand.complex, lr_pairs_lee_d1$receptor.complex))
all_sym_single <- unique(trimws(unlist(strsplit(all_sym[grepl("[_+]", all_sym)], "[_+]"))))
all_sym_single <- unique(c(all_sym[!grepl("[_+]", all_sym)], all_sym_single))
tbl_lr <- data.frame(SYMBOL = character(0), GENENAME = character(0), stringsAsFactors = FALSE)
if (length(all_sym_single) > 0) {
  tbl_lr <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = all_sym_single, columns = "GENENAME", keytype = "SYMBOL"))
  tbl_lr <- dplyr::distinct(tbl_lr, SYMBOL, .keep_all = TRUE)
}
sym_to_name <- if (nrow(tbl_lr) > 0) stats::setNames(tbl_lr$GENENAME, tbl_lr$SYMBOL) else character(0)
lr_ref <- lr_pairs_lee_d1
lr_ref$ligand_genename <- rep(NA_character_, nrow(lr_ref))
lr_ref$receptor_genename <- rep(NA_character_, nrow(lr_ref))
has_sym_to_name <- length(sym_to_name) > 0
if (!has_sym_to_name) {
  lr_ref$ligand_genename <- lr_ref$ligand.complex
  lr_ref$receptor_genename <- lr_ref$receptor.complex
}
if (has_sym_to_name) {
  ligand_is_complex <- grepl("[_+]", lr_ref$ligand.complex)
  lr_ref$ligand_genename[!ligand_is_complex] <- ifelse(lr_ref$ligand.complex[!ligand_is_complex] %in% names(sym_to_name), sym_to_name[lr_ref$ligand.complex[!ligand_is_complex]], lr_ref$ligand.complex[!ligand_is_complex])
  complex_ligands <- lr_ref$ligand.complex[ligand_is_complex]
  mapped_complex_lig <- vapply(complex_ligands, function(s) {
    parts <- trimws(strsplit(s, "[_+]")[[1]])
    mapped_parts <- ifelse(parts %in% names(sym_to_name), sym_to_name[parts], parts)
    paste(mapped_parts, collapse = "_")
  }, character(1))
  lr_ref$ligand_genename[ligand_is_complex] <- mapped_complex_lig
  receptor_is_complex <- grepl("[_+]", lr_ref$receptor.complex)
  lr_ref$receptor_genename[!receptor_is_complex] <- ifelse(lr_ref$receptor.complex[!receptor_is_complex] %in% names(sym_to_name), sym_to_name[lr_ref$receptor.complex[!receptor_is_complex]], lr_ref$receptor.complex[!receptor_is_complex])
  complex_receptors <- lr_ref$receptor.complex[receptor_is_complex]
  mapped_complex_rec <- vapply(complex_receptors, function(s) {
    parts <- trimws(strsplit(s, "[_+]")[[1]])
    mapped_parts <- ifelse(parts %in% names(sym_to_name), sym_to_name[parts], parts)
    paste(mapped_parts, collapse = "_")
  }, character(1))
  lr_ref$receptor_genename[receptor_is_complex] <- mapped_complex_rec
}
write.csv(lr_ref, file.path(OUTPUT_DIR, "LeeDay1_LIANA_LigandReceptor_Reference.csv"), row.names = FALSE)
print(paste("Saved L-R reference with gene names to LeeDay1_LIANA_LigandReceptor_Reference.csv"))
print_gene_mapping <- has_sym_to_name
if (print_gene_mapping) {
  gene_map_df <- data.frame(symbol = names(sym_to_name), full_name = unname(sym_to_name), stringsAsFactors = FALSE)
  print("LIANA gene symbols in Lee Day 1 -> full names (org.Mm.eg.db):")
  print(gene_map_df)
}

# Save Lee Day 1 LIANA
write.csv(lee_day1_liana_result_df, file.path(OUTPUT_DIR, "LeeDay1_LIANA_AllResults.csv"), row.names = FALSE)
write.csv(lee_day1_arg1pos_received, file.path(OUTPUT_DIR, "LeeDay1_Arg1pos_ReceivedSignals.csv"), row.names = FALSE)
write.csv(lee_day1_arg1neg_received, file.path(OUTPUT_DIR, "LeeDay1_Arg1neg_ReceivedSignals.csv"), row.names = FALSE)
write.csv(lee_day1_arg1pos_top_signals, file.path(OUTPUT_DIR, "LeeDay1_LIANA_Arg1pos_Top50Signals.csv"), row.names = FALSE)
write.csv(lee_day1_arg1pos_specific, file.path(OUTPUT_DIR, "LeeDay1_LIANA_Arg1pos_Specific.csv"), row.names = FALSE)
write.csv(lee_day1_liana_consensus, file.path(OUTPUT_DIR, "LeeDay1_LIANA_NeutrophilConsensus.csv"), row.names = FALSE)
write.csv(lee_day1_liana_arg1pos_topranked, file.path(OUTPUT_DIR, "LeeDay1_LIANA_Arg1pos_TopRanked.csv"), row.names = FALSE)
write.csv(lee_day1_liana_arg1neg_topranked, file.path(OUTPUT_DIR, "LeeDay1_LIANA_Arg1neg_TopRanked.csv"), row.names = FALSE)
saveRDS(list(
  liana_result = lee_day1_liana_result,
  liana_aggregated = lee_day1_liana_result_df,
  arg1pos_received = lee_day1_arg1pos_received,
  arg1pos_top_signals = lee_day1_arg1pos_top_signals,
  arg1pos_specific = lee_day1_arg1pos_specific,
  neutrophil_consensus = lee_day1_liana_consensus,
  arg1pos_topranked = lee_day1_liana_arg1pos_topranked,
  arg1neg_topranked = lee_day1_liana_arg1neg_topranked
), file.path(OUTPUT_DIR, "LeeDay1_LIANA_Results.rds"))
# Uncomment block below to load and skip re-running Lee Day 1 LIANA:
# lee_day1_liana_loaded <- readRDS(file.path(OUTPUT_DIR, "LeeDay1_LIANA_Results.rds"))
# lee_day1_liana_result <- lee_day1_liana_loaded$liana_result
# lee_day1_liana_result_df <- lee_day1_liana_loaded$liana_aggregated
# lee_day1_arg1pos_received <- lee_day1_liana_loaded$arg1pos_received
# lee_day1_arg1pos_top_signals <- lee_day1_liana_loaded$arg1pos_top_signals
# lee_day1_arg1pos_specific <- lee_day1_liana_loaded$arg1pos_specific
# lee_day1_liana_consensus <- lee_day1_liana_loaded$neutrophil_consensus
# lee_day1_liana_arg1pos_topranked <- lee_day1_liana_loaded$arg1pos_topranked
# lee_day1_liana_arg1neg_topranked <- lee_day1_liana_loaded$arg1neg_topranked

# -------- Lee Day 1: LIANA visualizations (official: liana_dotplot, heat_freq, chord_freq; https://saezlab.github.io/liana/) --------
# Plot tweaks (change these to adjust appearance)
USE_FULL_GENE_NAMES_P01F <- FALSE
LIANA_P01F_LEGEND_DOT_SIZE <- 3
LIANA_P01F_LEGEND_KEY_CM <- 0.6
LIANA_DOTPLOT_SIZE_RANGE <- c(2, 10)
LIANA_TOP_SIGNAL_POINT_SIZE_RANGE <- c(2, 8)
LIANA_P01F_WIDTH <- 12
LIANA_P01F_HEIGHT_PER_ROW <- 0.35
LIANA_P01G_WIDTH <- 10
LIANA_P01G_HEIGHT <- 8
LIANA_P01G2_WIDTH <- 12
LIANA_P01G2_HEIGHT <- 10
liana_network_lee_d1_both <- lee_day1_liana_result_df |>
  dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |>
  dplyr::filter(!(source %in% NEUTROPHIL_STATES)) |>
  dplyr::group_by(target) |>
  dplyr::arrange(aggregate_rank) |>
  dplyr::slice_head(n = 20) |>
  dplyr::ungroup()
symbols_lr <- unique(c(liana_network_lee_d1_both$ligand.complex, liana_network_lee_d1_both$receptor.complex))
symbols_complex <- symbols_lr[grepl("[_+]", symbols_lr)]
symbols_from_complex <- unique(trimws(unlist(strsplit(symbols_complex, "[_+]"))))
symbols_single <- unique(c(symbols_lr[!grepl("[_+]", symbols_lr)], symbols_from_complex))
gene_info <- character(0)
if (length(symbols_single) > 0) {
  tbl <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = symbols_single, columns = "GENENAME", keytype = "SYMBOL"))
  tbl <- dplyr::distinct(tbl, SYMBOL, .keep_all = TRUE)
  gene_info <- stats::setNames(tbl$GENENAME, tbl$SYMBOL)
}
liana_network_lee_d1_both$ligand_pretty <- liana_network_lee_d1_both$ligand.complex
liana_network_lee_d1_both$receptor_pretty <- liana_network_lee_d1_both$receptor.complex
n_p01f_rows <- nrow(liana_network_lee_d1_both)
if (length(gene_info) > 0 && n_p01f_rows > 0) {
  i_p01f <- 1
  while (i_p01f <= n_p01f_rows) {
    sym_lig <- liana_network_lee_d1_both$ligand.complex[i_p01f]
    if (!grepl("[_+]", sym_lig)) {
      nm <- gene_info[sym_lig]
      if (!is.na(nm) && nzchar(nm)) liana_network_lee_d1_both$ligand_pretty[i_p01f] <- paste0(nm, " (", sym_lig, ")")
    }
    if (grepl("[_+]", sym_lig)) {
      parts_lig <- trimws(strsplit(sym_lig, "[_+]")[[1]])
      out_lig <- character(length(parts_lig))
      j_lig <- 1
      while (j_lig <= length(parts_lig)) {
        p <- parts_lig[j_lig]
        nm <- gene_info[p]
        out_lig[j_lig] <- if (!is.na(nm) && nzchar(nm)) paste0(nm, " (", p, ")") else p
        j_lig <- j_lig + 1
      }
      liana_network_lee_d1_both$ligand_pretty[i_p01f] <- paste(out_lig, collapse = "_")
    }
    i_p01f <- i_p01f + 1
  }
  i_p01f <- 1
  while (i_p01f <= n_p01f_rows) {
    sym_rec <- liana_network_lee_d1_both$receptor.complex[i_p01f]
    if (!grepl("[_+]", sym_rec)) {
      nm <- gene_info[sym_rec]
      if (!is.na(nm) && nzchar(nm)) liana_network_lee_d1_both$receptor_pretty[i_p01f] <- paste0(nm, " (", sym_rec, ")")
    }
    if (grepl("[_+]", sym_rec)) {
      parts_rec <- trimws(strsplit(sym_rec, "[_+]")[[1]])
      out_rec <- character(length(parts_rec))
      j_rec <- 1
      while (j_rec <= length(parts_rec)) {
        p <- parts_rec[j_rec]
        nm <- gene_info[p]
        out_rec[j_rec] <- if (!is.na(nm) && nzchar(nm)) paste0(nm, " (", p, ")") else p
        j_rec <- j_rec + 1
      }
      liana_network_lee_d1_both$receptor_pretty[i_p01f] <- paste(out_rec, collapse = "_")
    }
    i_p01f <- i_p01f + 1
  }
}
liana_network_lee_d1_both$interaction_label <- if (USE_FULL_GENE_NAMES_P01F) paste0(liana_network_lee_d1_both$ligand_pretty, "\u2013(", liana_network_lee_d1_both$receptor_pretty, ")") else paste0(liana_network_lee_d1_both$ligand.complex, "\u2013(", liana_network_lee_d1_both$receptor.complex, ")")
p_01f <- ggplot(liana_network_lee_d1_both, aes(x = target, y = interaction_label, size = -log10(aggregate_rank + 1e-10), color = source)) + geom_point(alpha = 0.8) + scale_x_discrete(limits = NEUTROPHIL_STATES, drop = FALSE) + theme_minimal() + labs(title = "Top Signals Received: Arg1+ vs Arg1- (Lee Day 1)", subtitle = "Y = Ligand\u2013(Receptor); X = receiver; color = sender (neutrophil\u2192neutrophil excluded from data)", x = "Receiver (neutrophil state)", y = "Interaction", color = "Sender cell type", size = "Consensus support\n(-log10 aggregate rank)") + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), plot.margin = margin(10, 80, 10, 10), legend.key.size = unit(LIANA_P01F_LEGEND_KEY_CM, "cm"), legend.text = element_text(size = 11)) + guides(color = guide_legend(override.aes = list(size = LIANA_P01F_LEGEND_DOT_SIZE)))
# Unified-scale print: pool values from every cohort object that exists() so far, then print (same pattern for all LIANA/PROGENy comparison ggplots below).
vals_nlr <- numeric(0)
if (exists("liana_network_lee_d1_both") && is.data.frame(liana_network_lee_d1_both) && nrow(liana_network_lee_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d1_both$aggregate_rank + 1e-10))
if (exists("liana_network_lee_d3_both") && is.data.frame(liana_network_lee_d3_both) && nrow(liana_network_lee_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d3_both$aggregate_rank + 1e-10))
if (exists("liana_network_wang_d3_both") && is.data.frame(liana_network_wang_d3_both) && nrow(liana_network_wang_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_wang_d3_both$aggregate_rank + 1e-10))
if (exists("liana_network_qin_d1_both") && is.data.frame(liana_network_qin_d1_both) && nrow(liana_network_qin_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d1_both$aggregate_rank + 1e-10))
if (exists("liana_network_qin_d3_both") && is.data.frame(liana_network_qin_d3_both) && nrow(liana_network_qin_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d3_both$aggregate_rank + 1e-10))
vals_nlr <- vals_nlr[is.finite(vals_nlr)]
UNIFY_LIANA_NEGLOG10 <- if (length(vals_nlr) > 0) range(vals_nlr) else c(0, 10)
if (length(UNIFY_LIANA_NEGLOG10) != 2 || !all(is.finite(UNIFY_LIANA_NEGLOG10))) UNIFY_LIANA_NEGLOG10 <- c(0, 10)
if (UNIFY_LIANA_NEGLOG10[1] == UNIFY_LIANA_NEGLOG10[2]) UNIFY_LIANA_NEGLOG10[2] <- UNIFY_LIANA_NEGLOG10[1] + 1e-6
print(p_01f + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = (if (exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE else c(2, 8))))
# Exclude neutrophil states as sources so plot shows only external signals TO neutrophils
lee_day1_external_to_neutrophils <- lee_day1_liana_result_df |>
  dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |>
  dplyr::filter(!(source %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))) |>
  dplyr::arrange(aggregate_rank)
liana_top_lee_d1 <- dplyr::slice_head(lee_day1_external_to_neutrophils, n = 20)
liana_top_lee_d1$lr_label <- paste0(liana_top_lee_d1$source, " -> ", liana_top_lee_d1$target, "  ", liana_top_lee_d1$ligand.complex, "\u2013(", liana_top_lee_d1$receptor.complex, ")")
p_01g <- ggplot(liana_top_lee_d1, aes(x = reorder(lr_label, aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 20 L-R Pairs - LIANA (Lee Day 1)\nExternal signals to Arg1+ / Arg1- neutrophils", x = "Source -> Target  Ligand\u2013(Receptor)", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME
vals_ext <- numeric(0)
if (exists("lee_day1_external_to_neutrophils") && is.data.frame(lee_day1_external_to_neutrophils) && nrow(lee_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day1_external_to_neutrophils$aggregate_rank)
if (exists("lee_day3_external_to_neutrophils") && is.data.frame(lee_day3_external_to_neutrophils) && nrow(lee_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day3_external_to_neutrophils$aggregate_rank)
if (exists("wang_day3_external_to_neutrophils") && is.data.frame(wang_day3_external_to_neutrophils) && nrow(wang_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(wang_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, wang_day3_external_to_neutrophils$aggregate_rank)
if (exists("qin_day1_external_to_neutrophils") && is.data.frame(qin_day1_external_to_neutrophils) && nrow(qin_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day1_external_to_neutrophils$aggregate_rank)
if (exists("qin_day3_external_to_neutrophils") && is.data.frame(qin_day3_external_to_neutrophils) && nrow(qin_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day3_external_to_neutrophils$aggregate_rank)
if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0 && "aggregate_rank" %in% names(liana_top_lee_d1)) vals_ext <- c(vals_ext, liana_top_lee_d1$aggregate_rank)
if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_lee_d3_external)) vals_ext <- c(vals_ext, liana_top_lee_d3_external$aggregate_rank)
if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_wang_d3_external)) vals_ext <- c(vals_ext, liana_top_wang_d3_external$aggregate_rank)
if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d1_external)) vals_ext <- c(vals_ext, liana_top_qin_d1_external$aggregate_rank)
if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d3_external)) vals_ext <- c(vals_ext, liana_top_qin_d3_external$aggregate_rank)
vals_ext <- suppressWarnings(as.numeric(vals_ext))
vals_ext <- vals_ext[is.finite(vals_ext)]
plot_ext_ar <- suppressWarnings(as.numeric(liana_top_lee_d1[["aggregate_rank"]]))
plot_ext_ar <- plot_ext_ar[is.finite(plot_ext_ar)]
UNIFY_EXT_AR <- range(c(vals_ext, plot_ext_ar), na.rm = TRUE)
if (length(plot_ext_ar) == 0 && length(vals_ext) == 0) UNIFY_EXT_AR <- c(0, 1)
if (!all(is.finite(UNIFY_EXT_AR))) UNIFY_EXT_AR <- c(0, 1)
if (UNIFY_EXT_AR[1] == UNIFY_EXT_AR[2]) UNIFY_EXT_AR[2] <- UNIFY_EXT_AR[1] + max(abs(UNIFY_EXT_AR[1]) * 1e-6, 1e-12)
print(p_01g + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
unique_sources_lee_d1 <- unique(lee_day1_liana_result_df$source)
neutrophil_targets_ld1 <- NEUTROPHIL_STATES[NEUTROPHIL_STATES %in% unique(lee_day1_liana_result_df$target)]
dot_sources_lee_d1 <- setdiff(unique_sources_lee_d1, neutrophil_targets_ld1)
if (length(dot_sources_lee_d1) == 0) dot_sources_lee_d1 <- unique_sources_lee_d1
do_ld1_dotplot <- length(neutrophil_targets_ld1) > 0 && length(dot_sources_lee_d1) > 0
p_01g2 <- NULL
if (do_ld1_dotplot) p_01g2 <- liana::liana_dotplot(lee_day1_liana_result_df, source_groups = dot_sources_lee_d1, target_groups = neutrophil_targets_ld1, ntop = 20, size_range = LIANA_DOTPLOT_SIZE_RANGE)
if (!is.null(p_01g2)) p_01g2 <- p_01g2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
if (!is.null(p_01g2)) print(p_01g2)
# Official LIANA: filter(aggregate_rank <= 0.01) only. Chord: source/target as character so labels show cell-type names.
liana_trunc_lee_d1 <- dplyr::filter(lee_day1_liana_result_df, aggregate_rank <= 0.01)
has_trunc_ld1 <- nrow(liana_trunc_lee_d1) > 0
if (has_trunc_ld1) liana_trunc_lee_d1$source <- as.character(liana_trunc_lee_d1$source)
if (has_trunc_ld1) liana_trunc_lee_d1$target <- as.character(liana_trunc_lee_d1$target)
if (!has_trunc_ld1) print("No interactions with aggregate_rank <= 0.01; skip heat_freq and chord_freq (Lee D1)")
p_heat_freq_ld1 <- NULL
if (has_trunc_ld1) p_heat_freq_ld1 <- liana::heat_freq(liana_trunc_lee_d1)
if (!is.null(p_heat_freq_ld1)) print(p_heat_freq_ld1)
unique_sources_chord_ld1 <- unique(liana_trunc_lee_d1$source)
unique_targets_chord_ld1 <- unique(liana_trunc_lee_d1$target)
p_chord_freq_ld1 <- NULL
if (has_trunc_ld1) grDevices::png(file.path(OUTPUT_DIR, "LeeDay1_LIANA_ChordFreq.png"), width = 1400, height = 1400, res = 150)
if (has_trunc_ld1) tryCatch(liana::chord_freq(liana_trunc_lee_d1, source_groups = unique_sources_chord_ld1, target_groups = unique_targets_chord_ld1), error = function(e) message("chord_freq Lee D1 (PNG): ", conditionMessage(e)))
if (has_trunc_ld1) grDevices::dev.off()
if (has_trunc_ld1) p_chord_freq_ld1 <- tryCatch(liana::chord_freq(liana_trunc_lee_d1, source_groups = unique_sources_chord_ld1, target_groups = unique_targets_chord_ld1), error = function(e) { message("chord_freq Lee D1: ", conditionMessage(e)); NULL })
if (!is.null(p_chord_freq_ld1)) print(p_chord_freq_ld1)
liana_mat_ld1 <- if (has_trunc_ld1) as.matrix(table(liana_trunc_lee_d1$source, liana_trunc_lee_d1$target)) else matrix(0, 0, 0)
p_liana_heatmap_ld1 <- NULL
if (has_trunc_ld1 && nrow(liana_mat_ld1) > 0 && ncol(liana_mat_ld1) > 0) p_liana_heatmap_ld1 <- liana::liana_heatmap(liana_mat_ld1)
if (!is.null(p_liana_heatmap_ld1)) ComplexHeatmap::draw(p_liana_heatmap_ld1)
receptor_freq_lee_d1 <- dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(lee_day1_liana_consensus, receptor.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15)
p_01h <- ggplot(receptor_freq_lee_d1, aes(x = reorder(receptor.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Receptors - Lee Day 1", x = "Receptor", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_01h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
ligand_freq_lee_d1 <- dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(lee_day1_liana_consensus, ligand.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15)
p_01i <- ggplot(ligand_freq_lee_d1, aes(x = reorder(ligand.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Ligands - Lee Day 1", x = "Ligand", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_01i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
interaction_matrix_lee_d1 <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(lee_day1_liana_consensus, source, target), interaction_count = dplyr::n(), .groups = "drop"), names_from = target, values_from = interaction_count, values_fill = 0)
cols_num_ld1 <- setdiff(colnames(interaction_matrix_lee_d1), "source")
interaction_matrix_lee_d1_mat <- as.matrix(interaction_matrix_lee_d1[, cols_num_ld1, drop = FALSE])
rownames(interaction_matrix_lee_d1_mat) <- interaction_matrix_lee_d1$source
p_01j <- NULL
if (nrow(interaction_matrix_lee_d1_mat) >= 2 && ncol(interaction_matrix_lee_d1_mat) >= 2) p_01j <- pheatmap::pheatmap(interaction_matrix_lee_d1_mat, color = colorRampPalette(c("white", "yellow", "orange", "red"))(100), main = "Cell Type Interaction Frequency - Lee Day 1", display_numbers = TRUE)
if (!is.null(p_01j)) { print(p_01j); grid::grid.newpage(); grid::grid.draw(p_01j$gtable) }
source_importance_lee_d1 <- dplyr::arrange(dplyr::summarise(dplyr::group_by(lee_day1_liana_consensus, source), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), importance_score = dplyr::n() * (1 - mean(aggregate_rank))), desc(importance_score))
p_01k <- ggplot(source_importance_lee_d1, aes(x = reorder(source, importance_score), y = importance_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkred") + labs(title = "Source Cell Importance - Lee Day 1", x = "Cell Type", y = "Importance Score") + theme_minimal() + PLOT_TITLE_THEME
vals_iy <- numeric(0)
vals_imr <- numeric(0)
imp_tabnames <- c("source_importance_lee_d1", "source_importance_lee_d3", "source_importance_wang_d3", "source_importance_qin_d1", "source_importance_qin_d3")
for (fn in imp_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("importance_score" %in% names(d)) vals_iy <- c(vals_iy, d$importance_score)
  if ("mean_rank" %in% names(d)) vals_imr <- c(vals_imr, d$mean_rank)
}
vals_iy <- vals_iy[is.finite(vals_iy)]
vals_imr <- vals_imr[is.finite(vals_imr)]
UNIFY_SRC_IMP_Y <- if (length(vals_iy) > 0) range(vals_iy) else c(0, 1)
UNIFY_SRC_IMP_MR <- if (length(vals_imr) > 0) range(vals_imr) else c(0, 1)
if (UNIFY_SRC_IMP_Y[1] == UNIFY_SRC_IMP_Y[2]) UNIFY_SRC_IMP_Y[2] <- UNIFY_SRC_IMP_Y[1] + 1e-6
if (UNIFY_SRC_IMP_MR[1] == UNIFY_SRC_IMP_MR[2]) UNIFY_SRC_IMP_MR[2] <- UNIFY_SRC_IMP_MR[1] + 1e-6
print(p_01k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
p_01n <- ggplot(head(lee_day1_liana_arg1pos_topranked, 15), aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), aggregate_rank), y = -aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top Arg1pos Interactions - Lee Day 1", x = "L-R Pair", y = "Rank Score") + theme_minimal() + PLOT_TITLE_THEME
vals_a1 <- numeric(0)
if (exists("lee_day1_liana_arg1pos_topranked") && is.data.frame(lee_day1_liana_arg1pos_topranked) && nrow(lee_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("lee_day3_liana_arg1pos_topranked") && is.data.frame(lee_day3_liana_arg1pos_topranked) && nrow(lee_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("wang_day3_liana_arg1pos_topranked") && is.data.frame(wang_day3_liana_arg1pos_topranked) && nrow(wang_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(wang_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day1_liana_arg1pos_topranked") && is.data.frame(qin_day1_liana_arg1pos_topranked) && nrow(qin_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day3_liana_arg1pos_topranked") && is.data.frame(qin_day3_liana_arg1pos_topranked) && nrow(qin_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day3_liana_arg1pos_topranked$aggregate_rank, 15))
vals_a1 <- vals_a1[is.finite(vals_a1)]
UNIFY_ARG1_AR <- if (length(vals_a1) > 0) range(vals_a1) else c(0, 1)
if (UNIFY_ARG1_AR[1] == UNIFY_ARG1_AR[2]) UNIFY_ARG1_AR[2] <- UNIFY_ARG1_AR[1] + 1e-6
UNIFY_ARG1_NEGY <- c(-UNIFY_ARG1_AR[2], -UNIFY_ARG1_AR[1])
print(p_01n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
lee_day1_arg1pos_specific_plot <- head(lee_day1_arg1pos_specific, 15)
p_01n2 <- ggplot(lee_day1_arg1pos_specific_plot, aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), specificity_index), y = specificity_index, fill = source)) + geom_bar(stat = "identity") + coord_flip() + labs(title = "Arg1+-Enriched Signals (Lee Day 1, spec>0.1 or unique)", x = "L-R Pair", y = "Specificity Index") + theme_minimal() + PLOT_TITLE_THEME
print(p_01n2)
target_specificity_lee_d1 <- dplyr::arrange(dplyr::summarise(dplyr::group_by(dplyr::filter(lee_day1_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")), target), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), specificity_score = dplyr::n() * (1 - mean(aggregate_rank)), .groups = "drop"), desc(specificity_score))
p_01l <- ggplot(target_specificity_lee_d1, aes(x = reorder(target, specificity_score), y = specificity_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkgreen") + labs(title = "Target Cell Specificity - Lee Day 1", x = "Cell Type", y = "Specificity Score") + theme_minimal() + PLOT_TITLE_THEME
vals_sy <- numeric(0)
vals_smr <- numeric(0)
spec_tabnames <- c("target_specificity_lee_d1", "target_specificity_lee_d3", "target_specificity_wang_d3", "target_specificity_qin_d1", "target_specificity_qin_d3")
for (fn in spec_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("specificity_score" %in% names(d)) vals_sy <- c(vals_sy, d$specificity_score)
  if ("mean_rank" %in% names(d)) vals_smr <- c(vals_smr, d$mean_rank)
}
vals_sy <- vals_sy[is.finite(vals_sy)]
vals_smr <- vals_smr[is.finite(vals_smr)]
UNIFY_TGT_SPEC_Y <- if (length(vals_sy) > 0) range(vals_sy) else c(0, 1)
UNIFY_TGT_SPEC_MR <- if (length(vals_smr) > 0) range(vals_smr) else c(0, 1)
if (UNIFY_TGT_SPEC_Y[1] == UNIFY_TGT_SPEC_Y[2]) UNIFY_TGT_SPEC_Y[2] <- UNIFY_TGT_SPEC_Y[1] + 1e-6
if (UNIFY_TGT_SPEC_MR[1] == UNIFY_TGT_SPEC_MR[2]) UNIFY_TGT_SPEC_MR[2] <- UNIFY_TGT_SPEC_MR[1] + 1e-6
print(p_01l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
p_01m <- ggplot(lee_day1_liana_consensus, aes(x = aggregate_rank)) + geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) + labs(title = "Distribution of L-R Pair Aggregate Ranks - Lee Day 1", x = "Aggregate Rank Score", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME
vals_hist <- numeric(0)
if (exists("lee_day1_liana_consensus") && is.data.frame(lee_day1_liana_consensus) && nrow(lee_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day1_liana_consensus)) vals_hist <- c(vals_hist, lee_day1_liana_consensus$aggregate_rank)
if (exists("lee_day3_liana_consensus") && is.data.frame(lee_day3_liana_consensus) && nrow(lee_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day3_liana_consensus)) vals_hist <- c(vals_hist, lee_day3_liana_consensus$aggregate_rank)
if (exists("wang_day3_liana_consensus") && is.data.frame(wang_day3_liana_consensus) && nrow(wang_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(wang_day3_liana_consensus)) vals_hist <- c(vals_hist, wang_day3_liana_consensus$aggregate_rank)
if (exists("qin_day1_liana_consensus") && is.data.frame(qin_day1_liana_consensus) && nrow(qin_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day1_liana_consensus)) vals_hist <- c(vals_hist, qin_day1_liana_consensus$aggregate_rank)
if (exists("qin_day3_liana_consensus") && is.data.frame(qin_day3_liana_consensus) && nrow(qin_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day3_liana_consensus)) vals_hist <- c(vals_hist, qin_day3_liana_consensus$aggregate_rank)
vals_hist <- vals_hist[is.finite(vals_hist)]
UNIFY_HIST_AR <- if (length(vals_hist) > 0) range(vals_hist) else c(0, 1)
if (UNIFY_HIST_AR[1] == UNIFY_HIST_AR[2]) UNIFY_HIST_AR[2] <- UNIFY_HIST_AR[1] + 1e-6
print(p_01m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
method_cols_lee_d1 <- colnames(lee_day1_liana_result_df)[grep("pval_", colnames(lee_day1_liana_result_df))]
consensus_lee_d1_step1 <- dplyr::slice_head(lee_day1_liana_consensus, n = 15)
consensus_lee_d1 <- NULL
consensus_lee_d1_has_cols <- length(method_cols_lee_d1) > 0 && nrow(consensus_lee_d1_step1) > 0
if (consensus_lee_d1_has_cols) consensus_lee_d1 <- as.data.frame(dplyr::select(consensus_lee_d1_step1, dplyr::all_of(method_cols_lee_d1)))
if (consensus_lee_d1_has_cols) consensus_lee_d1_rownames <- paste0(consensus_lee_d1_step1$source, " | ", consensus_lee_d1_step1$ligand.complex, " -> ", consensus_lee_d1_step1$receptor.complex)
if (consensus_lee_d1_has_cols) rownames(consensus_lee_d1) <- make.unique(as.character(consensus_lee_d1_rownames))
consensus_lee_d1_ready <- !is.null(consensus_lee_d1) && nrow(consensus_lee_d1) > 0
if (consensus_lee_d1_ready) p_01o <- pheatmap::pheatmap(consensus_lee_d1, color = colorRampPalette(c("red", "white", "blue"))(100), main = "Method Consensus (p-values) - Lee Day 1")
if (consensus_lee_d1_ready) print(p_01o)

print("✓ Lee Day 1 LIANA analysis complete")
# Checkpoint: resume from here if later section crashes. load(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay1_LIANA.RData"))
save.image(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay1_LIANA.RData"))
saveRDS(list(lee_day1 = lee_day1, lee_day1_liana_result_df = lee_day1_liana_result_df, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay1_LIANA.rds"))
# checkpoint_ld1 <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay1_LIANA.rds")); list2env(checkpoint_ld1, envir = .GlobalEnv)
print("✓ Checkpoint saved: Workspace_AfterLeeDay1_LIANA.RData, Checkpoint_AfterLeeDay1_LIANA.rds")

# -------- Lee Day 1: CellChat/CellCall downstream (overlap with LIANA). Linear: one step per line. --------
neutrophil_labels_d1 <- c("NeutrophilArg1pos", "NeutrophilArg1neg")
path_cc_d1 <- CELLCHAT_LEE_DAY1_RDS
path_ccall_d1 <- CELLCALL_LEE_DAY1_RDS
cc_d1_down <- NULL
cc_d1_file_ok <- !is.na(path_cc_d1) && nzchar(trimws(path_cc_d1)) && file.exists(path_cc_d1)
if (cc_d1_file_ok) cc_d1_down <- readRDS(path_cc_d1)
lee_day1_cellchat_neut <- NULL
if (!is.null(cc_d1_down)) comm_d1 <- CellChat::subsetCommunication(cc_d1_down)
if (!is.null(cc_d1_down)) lee_day1_cellchat_neut <- dplyr::filter(comm_d1, source %in% neutrophil_labels_d1 | target %in% neutrophil_labels_d1)
if (!is.null(lee_day1_cellchat_neut) && nrow(lee_day1_cellchat_neut) > 0) write.csv(lee_day1_cellchat_neut, file.path(OUTPUT_DIR, "LeeDay1_CellChat_NeutrophilCommunications.csv"), row.names = FALSE)
if (!is.null(lee_day1_cellchat_neut)) print(paste("Lee Day 1: CellChat neutrophil communications:", nrow(lee_day1_cellchat_neut)))
ccall_d1_file_ok <- !is.na(path_ccall_d1) && nzchar(trimws(path_ccall_d1)) && file.exists(path_ccall_d1)
ccall_d1 <- NULL
if (ccall_d1_file_ok) ccall_d1 <- readRDS(path_ccall_d1)
lee_day1_cellcall_lr <- NULL
if (!is.null(ccall_d1) && !is.null(ccall_d1@data$expr_l_r_log2_scale)) {
  lr_rownames_d1 <- rownames(ccall_d1@data$expr_l_r_log2_scale)
  ccall_lig_d1 <- sub("-.*", "", lr_rownames_d1)
  ccall_rec_d1 <- sub("^[^-]+-", "", lr_rownames_d1)
  ccall_rec_d1 <- gsub("-", "_", ccall_rec_d1)
  lee_day1_cellcall_lr <- data.frame(ligand.complex = ccall_lig_d1, receptor.complex = ccall_rec_d1, stringsAsFactors = FALSE)
  lee_day1_cellcall_lr <- lee_day1_cellcall_lr[nzchar(lee_day1_cellcall_lr$receptor.complex), ]
  lee_day1_cellcall_lr <- dplyr::distinct(lee_day1_cellcall_lr, ligand.complex, receptor.complex)
}
if (!is.null(lee_day1_cellcall_lr)) write.csv(lee_day1_cellcall_lr, file.path(OUTPUT_DIR, "LeeDay1_CellCall_LRPairs.csv"), row.names = FALSE)
if (!is.null(lee_day1_cellcall_lr)) print(paste("Lee Day 1: CellCall LR pairs:", nrow(lee_day1_cellcall_lr)))
lee_day1_liana_lr <- dplyr::distinct(lee_day1_liana_result_df, ligand.complex, receptor.complex)
if (nrow(lee_day1_liana_lr) > 0 && !is.null(lee_day1_cellcall_lr)) {
  lee_day1_lr_overlap_ccall <- dplyr::semi_join(lee_day1_liana_lr, lee_day1_cellcall_lr, by = c("ligand.complex", "receptor.complex"))
  write.csv(lee_day1_lr_overlap_ccall, file.path(OUTPUT_DIR, "LeeDay1_LIANA_CellCall_Overlap.csv"), row.names = FALSE)
  print(paste("Lee Day 1: LIANA-CellCall overlap:", nrow(lee_day1_lr_overlap_ccall), "LR pairs"))
}
has_cc_neut <- !is.null(lee_day1_cellchat_neut) && nrow(lee_day1_cellchat_neut) > 0 && "ligand" %in% names(lee_day1_cellchat_neut) && "receptor" %in% names(lee_day1_cellchat_neut)
if (nrow(lee_day1_liana_lr) > 0 && has_cc_neut) {
  lee_day1_cellchat_lr <- dplyr::distinct(lee_day1_cellchat_neut, ligand, receptor)
  lee_day1_cellchat_lr$ligand.complex <- lee_day1_cellchat_lr$ligand
  lee_day1_cellchat_lr$receptor.complex <- lee_day1_cellchat_lr$receptor
  lee_day1_lr_overlap_cc <- dplyr::semi_join(lee_day1_liana_lr, lee_day1_cellchat_lr, by = c("ligand.complex", "receptor.complex"))
  write.csv(lee_day1_lr_overlap_cc, file.path(OUTPUT_DIR, "LeeDay1_LIANA_CellChat_Overlap.csv"), row.names = FALSE)
  print(paste("Lee Day 1: LIANA-CellChat overlap:", nrow(lee_day1_lr_overlap_cc), "LR pairs"))
}

# Neutrophil cell list for PROGENy and BARTsc (Arg1+ and Arg1- only). Same target population as LIANA.
# Biological order: Ligand (LIANA) -> Pathway (PROGENy) -> TF (BARTsc) -> Gene expression.
lee_day1_neut_cells <- colnames(lee_day1)[c(lee_day1_pos, lee_day1_neg)]

# Top LIANA-predicted receptors (single-symbol only): expression by Arg1 status on neutrophils
top_receptors_ld1 <- head(unique(lee_day1_arg1pos_received$receptor.complex), LIANA_TOP_RECEPTOR_VLN)
top_receptors_ld1_single <- top_receptors_ld1[!grepl("[_+]", top_receptors_ld1)]
top_receptors_ld1_in_data <- top_receptors_ld1_single[top_receptors_ld1_single %in% rownames(lee_day1)]
if (length(top_receptors_ld1_in_data) > 0) {
  lee_day1_neut_obj_vln <- subset(lee_day1, cells = lee_day1_neut_cells)
  p_receptor_vln_ld1 <- Seurat::VlnPlot(lee_day1_neut_obj_vln, features = top_receptors_ld1_in_data, group.by = "arg1_status", pt.size = 0.1, ncol = min(3L, length(top_receptors_ld1_in_data)))
  print(p_receptor_vln_ld1)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "LeeDay1_TopLIANA_Receptors_VlnPlot.png"), p_receptor_vln_ld1, width = 10, height = 6, dpi = 150)
}

# -------- Lee Day 1 PROGENy Pathway Analysis (Exploratory; downstream pathway support only) --------
# Single block: model_mouse_full → decoupleR::run_wmean → stats → CSV/RDS → bar + violin + heatmap + ligand–pathway dot (all plots below stay in this section).
# Runs immediately after LIANA: Ligand (LIANA) -> Pathway (PROGENy) -> TF (BARTsc) -> Gene expression.
# Target population: NeutrophilArg1pos and NeutrophilArg1neg only; pathway activity compared Arg1pos vs Arg1neg.
# Official: https://saezlab.github.io/progeny/ ; use progeny::model_mouse_full and decoupleR for activity.
print("--- PROGENy Pathway Analysis: Lee Day 1 (Exploratory) ---")
# Secondary support only: PROGENy summarizes downstream pathway state in Arg1pos vs Arg1neg neutrophils.

# Get PROGENy pathway gene sets for mouse (official: model_mouse_full). Column names differ by progeny version (Gene vs gene).
progeny_model_mouse_lee_d1 <- progeny::model_mouse_full
colnames(progeny_model_mouse_lee_d1) <- tolower(colnames(progeny_model_mouse_lee_d1))
colnames(progeny_model_mouse_lee_d1)[colnames(progeny_model_mouse_lee_d1) == "p.value"] <- "p_value"
stopifnot(all(c("gene", "pathway", "weight") %in% colnames(progeny_model_mouse_lee_d1)))
progeny_model_mouse_lee_d1_dim <- dim(progeny_model_mouse_lee_d1)
progeny_model_mouse_lee_d1_colnames <- colnames(progeny_model_mouse_lee_d1)
print(paste("PROGENy model dimensions:", paste(progeny_model_mouse_lee_d1_dim, collapse = " x ")))
print(paste("PROGENy model columns:", paste(progeny_model_mouse_lee_d1_colnames, collapse = ", ")))
progeny_network_lee_d1 <- data.frame(
  source = progeny_model_mouse_lee_d1$pathway,
  target = progeny_model_mouse_lee_d1$gene,
  weight = progeny_model_mouse_lee_d1$weight,
  stringsAsFactors = FALSE
)
progeny_network_lee_d1 <- progeny_network_lee_d1[progeny_network_lee_d1$weight != 0, ]
progeny_network_lee_d1_nrow <- nrow(progeny_network_lee_d1)
progeny_network_lee_d1_pathways_unique <- unique(progeny_network_lee_d1$source)
print(paste("PROGENy network: nrow =", progeny_network_lee_d1_nrow, ", unique pathways =", length(progeny_network_lee_d1_pathways_unique)))
print(paste("PROGENy pathways:", paste(progeny_network_lee_d1_pathways_unique, collapse = ", ")))
stopifnot(progeny_network_lee_d1_nrow > 0)
progeny_network_lee_d1_has_pvalue <- "p_value" %in% colnames(progeny_network_lee_d1)
progeny_network_lee_d1_cols <- colnames(progeny_network_lee_d1)
progeny_network_lee_d1_cols_no_pvalue <- progeny_network_lee_d1_cols[progeny_network_lee_d1_cols != "p_value"]
progeny_network_lee_d1 <- progeny_network_lee_d1[, progeny_network_lee_d1_cols_no_pvalue, drop = FALSE]
progeny_network_lee_d1_has_source <- "source" %in% colnames(progeny_network_lee_d1)
progeny_network_lee_d1_has_target <- "target" %in% colnames(progeny_network_lee_d1)
stopifnot(progeny_network_lee_d1_has_source)
stopifnot(progeny_network_lee_d1_has_target)

# Subset to neutrophils only for pathway analysis (lee_day1_neut_cells defined above)
lee_day1_neut_expr_mat <- tryCatch(Seurat::GetAssayData(lee_day1, layer = "data")[, lee_day1_neut_cells, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
lee_day1_expr_mat <- lee_day1_neut_expr_mat

# Calculate pathway activity scores using weighted mean (WMEAN)
lee_day1_progeny_result <- tryCatch(decoupleR::run_wmean(mat = lee_day1_expr_mat, network = progeny_network_lee_d1, .source = "source", .target = "target", .mor = "weight", minsize = 5), error = function(e) data.frame())

# Convert result to wide format (pathways x cells). Use norm_wmean statistic only (not wmean).
lee_day1_progeny_for_wide <- if (nrow(lee_day1_progeny_result) > 0 && "statistic" %in% colnames(lee_day1_progeny_result)) dplyr::filter(lee_day1_progeny_result, statistic == "norm_wmean") else lee_day1_progeny_result
if (nrow(lee_day1_progeny_for_wide) == 0 && nrow(lee_day1_progeny_result) > 0 && "statistic" %in% colnames(lee_day1_progeny_result)) lee_day1_progeny_for_wide <- dplyr::filter(lee_day1_progeny_result, statistic == "wmean")
lee_day1_pw_wide <- tryCatch(tidyr::pivot_wider(lee_day1_progeny_for_wide, names_from = "condition", values_from = "score", id_cols = "source"), error = function(e) data.frame())
lee_day1_pw_cols_num <- tryCatch(sapply(lee_day1_pw_wide[, -1, drop = FALSE], function(x) as.numeric(unlist(x))), error = function(e) matrix(0, nrow = 0, ncol = 0))
lee_day1_progeny_scores_mat <- tryCatch(as.matrix(lee_day1_pw_cols_num), error = function(e) matrix(0, nrow = 0, ncol = 0))
rownames(lee_day1_progeny_scores_mat) <- tryCatch(as.character(lee_day1_pw_wide$source), error = function(e) character(0))
colnames(lee_day1_progeny_scores_mat) <- tryCatch(colnames(lee_day1_pw_wide)[-1], error = function(e) character(0))

# Get cell names for Arg1pos and Arg1neg (neutrophils only)
lee_day1_arg1pos_cells <- colnames(lee_day1)[lee_day1_pos]
lee_day1_arg1neg_cells <- colnames(lee_day1)[lee_day1_neg]

# Extract pathway scores for each group (only if cells exist in matrix)
lee_day1_arg1pos_cells_in_mat <- tryCatch(lee_day1_arg1pos_cells[lee_day1_arg1pos_cells %in% colnames(lee_day1_progeny_scores_mat)], error = function(e) character(0))
lee_day1_arg1neg_cells_in_mat <- tryCatch(lee_day1_arg1neg_cells[lee_day1_arg1neg_cells %in% colnames(lee_day1_progeny_scores_mat)], error = function(e) character(0))

lee_day1_progeny_arg1pos_scores <- tryCatch(lee_day1_progeny_scores_mat[, lee_day1_arg1pos_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
lee_day1_progeny_arg1neg_scores <- tryCatch(lee_day1_progeny_scores_mat[, lee_day1_arg1neg_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))

# Calculate mean and median pathway scores per group (effect sizes)
lee_day1_pathway_means_arg1pos <- tryCatch(rowMeans(lee_day1_progeny_arg1pos_scores, na.rm = TRUE), error = function(e) numeric(0))
lee_day1_pathway_means_arg1neg <- tryCatch(rowMeans(lee_day1_progeny_arg1neg_scores, na.rm = TRUE), error = function(e) numeric(0))
lee_day1_pathway_medians_arg1pos <- tryCatch(apply(lee_day1_progeny_arg1pos_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
lee_day1_pathway_medians_arg1neg <- tryCatch(apply(lee_day1_progeny_arg1neg_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))

# Calculate log2 fold change (effect size)
lee_day1_pathway_shift <- pmax(0, -pmin(lee_day1_pathway_means_arg1pos, lee_day1_pathway_means_arg1neg, na.rm = TRUE)) + 1e-10
lee_day1_pathway_log2fc <- log2((lee_day1_pathway_means_arg1pos + lee_day1_pathway_shift) / (lee_day1_pathway_means_arg1neg + lee_day1_pathway_shift))

# Create comparison dataframe with effect sizes
lee_day1_pathway_comparison <- tryCatch(data.frame(
  pathway = names(lee_day1_pathway_means_arg1pos),
  Arg1pos_mean = lee_day1_pathway_means_arg1pos,
  Arg1neg_mean = lee_day1_pathway_means_arg1neg,
  Arg1pos_median = lee_day1_pathway_medians_arg1pos,
  Arg1neg_median = lee_day1_pathway_medians_arg1neg,
  mean_difference = lee_day1_pathway_means_arg1pos - lee_day1_pathway_means_arg1neg,
  median_difference = lee_day1_pathway_medians_arg1pos - lee_day1_pathway_medians_arg1neg,
  log2FC = lee_day1_pathway_log2fc,
  stringsAsFactors = FALSE
), error = function(e) data.frame())

# Exploratory pathway analysis: Use all PROGENy pathways (no pre-classification)
# Get all unique pathways from PROGENy scores matrix (use full matrix, not subset matrices)
# If scores matrix is empty, fall back to PROGENy network pathways
lee_day1_progeny_scores_mat_rownames <- rownames(lee_day1_progeny_scores_mat)
lee_day1_progeny_scores_mat_has_rownames <- !is.null(lee_day1_progeny_scores_mat_rownames) && length(lee_day1_progeny_scores_mat_rownames) > 0
lee_day1_all_pathways_from_scores <- lee_day1_progeny_scores_mat_rownames
lee_day1_all_pathways_from_network <- unique(progeny_network_lee_d1$source)
lee_day1_all_pathways_options <- list(lee_day1_all_pathways_from_network, lee_day1_all_pathways_from_scores)
lee_day1_all_pathways_index <- as.numeric(lee_day1_progeny_scores_mat_has_rownames) + 1
lee_day1_all_pathways <- lee_day1_all_pathways_options[[lee_day1_all_pathways_index]]
lee_day1_all_pathways <- unique(lee_day1_all_pathways)
lee_day1_pathway_results <- lee_day1_pathway_comparison

# Report sample sizes
lee_day1_arg1pos_count <- length(lee_day1_pos)
lee_day1_arg1neg_count <- length(lee_day1_neg)
print(paste("Lee Day 1: Arg1pos cells =", lee_day1_arg1pos_count, ", Arg1neg cells =", lee_day1_arg1neg_count))

# Statistical test: Wilcoxon test for pathway differences (all pathways)
# With sample size adequacy checks, effect sizes, and confidence intervals
lee_day1_pathway_pvals <- numeric(length(lee_day1_all_pathways))
names(lee_day1_pathway_pvals) <- lee_day1_all_pathways
lee_day1_pathway_effect_sizes <- numeric(length(lee_day1_all_pathways))
names(lee_day1_pathway_effect_sizes) <- lee_day1_all_pathways
lee_day1_pathway_ci_lower <- numeric(length(lee_day1_all_pathways))
names(lee_day1_pathway_ci_lower) <- lee_day1_all_pathways
lee_day1_pathway_ci_upper <- numeric(length(lee_day1_all_pathways))
names(lee_day1_pathway_ci_upper) <- lee_day1_all_pathways
lee_day1_pathway_adequate_n <- logical(length(lee_day1_all_pathways))
names(lee_day1_pathway_adequate_n) <- lee_day1_all_pathways

# Wilcoxon per pathway (same pattern as Lee Day 3 / Wang Day 3): all pathways in lee_day1_all_pathways
for (pw_idx in seq_along(lee_day1_all_pathways)) {
  pw <- lee_day1_all_pathways[pw_idx]
  arg1pos_scores_pw <- tryCatch(as.numeric(lee_day1_progeny_arg1pos_scores[pw, ]), error = function(e) numeric(0))
  arg1neg_scores_pw <- tryCatch(as.numeric(lee_day1_progeny_arg1neg_scores[pw, ]), error = function(e) numeric(0))
  arg1pos_scores_pw <- arg1pos_scores_pw[!is.na(arg1pos_scores_pw)]
  arg1neg_scores_pw <- arg1neg_scores_pw[!is.na(arg1neg_scores_pw)]
  n_arg1pos <- length(arg1pos_scores_pw)
  n_arg1neg <- length(arg1neg_scores_pw)
  lee_day1_pathway_adequate_n[pw] <- (n_arg1pos >= 5) & (n_arg1neg >= 5)
  median_arg1pos <- tryCatch(stats::median(arg1pos_scores_pw, na.rm = TRUE), error = function(e) 0)
  median_arg1neg <- tryCatch(stats::median(arg1neg_scores_pw, na.rm = TRUE), error = function(e) 0)
  lee_day1_pathway_effect_sizes[pw] <- median_arg1pos - median_arg1neg
  lee_day1_pathway_pvals[pw] <- 1.0
  lee_day1_pathway_ci_lower[pw] <- NA_real_
  lee_day1_pathway_ci_upper[pw] <- NA_real_
  wilcox_result <- tryCatch(stats::wilcox.test(arg1pos_scores_pw, arg1neg_scores_pw, conf.int = TRUE, conf.level = 0.95), error = function(e) NULL)
  wilcox_pvalue <- tryCatch(if (!is.null(wilcox_result)) wilcox_result$p.value else 1.0, error = function(e) 1.0)
  lee_day1_pathway_pvals[pw] <- if (length(wilcox_pvalue) > 0) wilcox_pvalue[1] else 1.0
  lee_day1_pathway_ci_lower[pw] <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[1] else NA_real_, error = function(e) NA_real_)
  lee_day1_pathway_ci_upper[pw] <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[2] else NA_real_, error = function(e) NA_real_)
  lee_day1_pathway_warning_msg <- paste("Warning: Pathway", pw, "has insufficient sample size (Arg1pos n =", n_arg1pos, ", Arg1neg n =", n_arg1neg, "). Skipping statistical test.")
  lee_day1_pathway_warning_vector <- c("", lee_day1_pathway_warning_msg)
  lee_day1_pathway_warning_index <- as.numeric(!lee_day1_pathway_adequate_n[pw]) + 1
  print(lee_day1_pathway_warning_vector[lee_day1_pathway_warning_index])
}

# Multiple testing correction across ALL pathways (FDR)
lee_day1_pathway_pvals_adj <- p.adjust(lee_day1_pathway_pvals, method = "BH")

# Add statistical results to comparison dataframe (p-values, effect sizes, CIs, sample size adequacy)
lee_day1_pathway_results$p_value <- tryCatch(lee_day1_pathway_pvals[lee_day1_pathway_results$pathway], error = function(e) rep(1.0, nrow(lee_day1_pathway_results)))
lee_day1_pathway_results$p_adj <- tryCatch(lee_day1_pathway_pvals_adj[lee_day1_pathway_results$pathway], error = function(e) rep(1.0, nrow(lee_day1_pathway_results)))
lee_day1_pathway_results$effect_size_median_diff <- tryCatch(lee_day1_pathway_effect_sizes[lee_day1_pathway_results$pathway], error = function(e) rep(0.0, nrow(lee_day1_pathway_results)))
lee_day1_pathway_results$ci_lower_95 <- tryCatch(lee_day1_pathway_ci_lower[lee_day1_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(lee_day1_pathway_results)))
lee_day1_pathway_results$ci_upper_95 <- tryCatch(lee_day1_pathway_ci_upper[lee_day1_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(lee_day1_pathway_results)))
lee_day1_pathway_results$adequate_sample_size <- tryCatch(lee_day1_pathway_adequate_n[lee_day1_pathway_results$pathway], error = function(e) rep(FALSE, nrow(lee_day1_pathway_results)))
lee_day1_pathway_results$significant <- tryCatch((lee_day1_pathway_results$p_adj < 0.05) & lee_day1_pathway_results$adequate_sample_size, error = function(e) rep(FALSE, nrow(lee_day1_pathway_results)))

# Save PROGENy results
write.csv(lee_day1_pathway_comparison, file.path(OUTPUT_DIR, "LeeDay1_PROGENy_PathwayComparison.csv"), row.names = FALSE)
write.csv(lee_day1_pathway_results, file.path(OUTPUT_DIR, "LeeDay1_PROGENy_PathwayResults.csv"), row.names = FALSE)
saveRDS(list(pathway_comparison = lee_day1_pathway_comparison, pathway_results = lee_day1_pathway_results, pathway_means_arg1pos = lee_day1_pathway_means_arg1pos, pathway_means_arg1neg = lee_day1_pathway_means_arg1neg, progeny_network = progeny_network_lee_d1), file.path(OUTPUT_DIR, "LeeDay1_PROGENy_Results.rds"))
# lee_day1_progeny_loaded <- readRDS(file.path(OUTPUT_DIR, "LeeDay1_PROGENy_Results.rds")); lee_day1_pathway_comparison <- lee_day1_progeny_loaded$pathway_comparison; lee_day1_pathway_results <- lee_day1_progeny_loaded$pathway_results; lee_day1_pathway_means_arg1pos <- lee_day1_progeny_loaded$pathway_means_arg1pos; lee_day1_pathway_means_arg1neg <- lee_day1_progeny_loaded$pathway_means_arg1neg; progeny_network_lee_d1 <- lee_day1_progeny_loaded$progeny_network

# Ligand-to-pathway mapping: Use top-ranked interactions for BOTH Arg1pos and Arg1neg
# Note: This uses top-ranked, not necessarily "specific" interactions
lee_day1_identified_ligands_arg1pos <- tryCatch(unique(lee_day1_liana_arg1pos_topranked$ligand.complex), error = function(e) character(0))
lee_day1_identified_ligands_arg1neg <- tryCatch(unique(lee_day1_liana_arg1neg_topranked$ligand.complex), error = function(e) character(0))
lee_day1_identified_ligands <- tryCatch(unique(c(lee_day1_identified_ligands_arg1pos, lee_day1_identified_ligands_arg1neg)), error = function(e) character(0))
lee_day1_ligand_source_population <- tryCatch(c(rep("Arg1pos", length(lee_day1_identified_ligands_arg1pos)), rep("Arg1neg", length(lee_day1_identified_ligands_arg1neg))), error = function(e) character(0))
# Exploratory ligand-to-pathway mapping: checks if ligand genes are PROGENy footprint targets (not receptor->pathway activation). Receptor->pathway links use integration section.
lee_day1_ligand_pathway_list <- list()
for (lig_idx in seq_along(lee_day1_identified_ligands)) {
  lig <- lee_day1_identified_ligands[lig_idx]
  lig_genes <- unlist(strsplit(lig, "_"))
  lig_genes <- unlist(strsplit(lig_genes, "[_+]"))
  lig_genes <- trimws(lig_genes[nzchar(lig_genes)])
  lig_pathways <- dplyr::filter(progeny_network_lee_d1, target %in% lig_genes)
  if (nrow(lig_pathways) > 0) lee_day1_ligand_pathway_list[[length(lee_day1_ligand_pathway_list) + 1L]] <- data.frame(ligand = lig, pathway = lig_pathways$source, weight = lig_pathways$weight, stringsAsFactors = FALSE)
}
lee_day1_ligand_pathway_map <- if (length(lee_day1_ligand_pathway_list) > 0) do.call(rbind, lee_day1_ligand_pathway_list) else data.frame(ligand = character(0), pathway = character(0), weight = numeric(0), stringsAsFactors = FALSE)
write.csv(lee_day1_ligand_pathway_map, file.path(OUTPUT_DIR, "LeeDay1_LigandToPathwayMapping.csv"), row.names = FALSE)

# Functional annotation: Exploratory annotation with population source (Arg1pos vs Arg1neg)
n_lig_d1 <- length(lee_day1_identified_ligands)
lee_day1_ligand_vec <- lee_day1_identified_ligands
lee_day1_ligand_genes_vec <- vapply(lee_day1_identified_ligands, function(lig) {
  lig_genes <- unlist(strsplit(lig, "[_+]"))
  paste(toupper(lig_genes), collapse = ";")
}, character(1))
ligand_in_arg1pos <- lee_day1_identified_ligands %in% lee_day1_identified_ligands_arg1pos
ligand_in_arg1neg <- lee_day1_identified_ligands %in% lee_day1_identified_ligands_arg1neg
lee_day1_target_pop_vec <- vapply(seq_len(n_lig_d1), function(i) {
  pops <- c("Arg1pos", "Arg1neg")[c(ligand_in_arg1pos[i], ligand_in_arg1neg[i])]
  paste(pops, collapse = ";")
}, character(1))
lee_day1_ligand_annotation <- data.frame(ligand = lee_day1_ligand_vec, ligand_genes = lee_day1_ligand_genes_vec, target_population = lee_day1_target_pop_vec, stringsAsFactors = FALSE)
write.csv(lee_day1_ligand_annotation, file.path(OUTPUT_DIR, "LeeDay1_LigandAnnotation.csv"), row.names = FALSE)

# PROGENy visualizations (Lee Day 1)
library(ggplot2)
lee_day1_pathway_long <- tryCatch(tidyr::pivot_longer(lee_day1_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
p_lee_d1_progeny1 <- tryCatch(ggplot(lee_day1_pathway_long, aes(x = pathway, y = Pathway_Score, fill = Group)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("Arg1pos_mean" = "red", "Arg1neg_mean" = "lightblue"), labels = c("Arg1pos", "Arg1neg")) + labs(title = "PROGENy Pathway Activity - Lee Day 1", x = "Pathway", y = "Pathway Activity Score") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
vals_pw <- numeric(0)
if (exists("lee_day1_pathway_long") && is.data.frame(lee_day1_pathway_long) && nrow(lee_day1_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day1_pathway_long)) vals_pw <- c(vals_pw, lee_day1_pathway_long$Pathway_Score)
if (exists("lee_day3_pathway_long") && is.data.frame(lee_day3_pathway_long) && nrow(lee_day3_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day3_pathway_long)) vals_pw <- c(vals_pw, lee_day3_pathway_long$Pathway_Score)
if (exists("wang_day3_pathway_long") && is.data.frame(wang_day3_pathway_long) && nrow(wang_day3_pathway_long) > 0 && "Pathway_Score" %in% names(wang_day3_pathway_long)) vals_pw <- c(vals_pw, wang_day3_pathway_long$Pathway_Score)
if (exists("qin_day1_pathway_long") && is.data.frame(qin_day1_pathway_long) && nrow(qin_day1_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day1_pathway_long)) vals_pw <- c(vals_pw, qin_day1_pathway_long$Pathway_Score)
if (exists("qin_day3_pathway_long") && is.data.frame(qin_day3_pathway_long) && nrow(qin_day3_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day3_pathway_long)) vals_pw <- c(vals_pw, qin_day3_pathway_long$Pathway_Score)
vals_pw <- vals_pw[is.finite(vals_pw)]
UNIFY_PW_Y <- if (length(vals_pw) > 0) range(vals_pw) else c(-1, 1)
if (UNIFY_PW_Y[1] == UNIFY_PW_Y[2]) UNIFY_PW_Y[2] <- UNIFY_PW_Y[1] + 1e-6
print(p_lee_d1_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
# Primary readout for sc: per-cell pathway score distributions (bar plot above = group means, supplementary)
top_pw_ld1 <- character(0)
if (nrow(lee_day1_pathway_comparison) > 0 && "mean_difference" %in% colnames(lee_day1_pathway_comparison)) {
  ord_pw_ld1 <- order(-abs(lee_day1_pathway_comparison$mean_difference))
  take_pw_ld1 <- min(6L, length(ord_pw_ld1))
  top_pw_ld1 <- as.character(lee_day1_pathway_comparison$pathway[ord_pw_ld1[seq_len(take_pw_ld1)]])
}
if (length(top_pw_ld1) > 0 && nrow(lee_day1_progeny_scores_mat) > 0) {
  n_cells_pw <- ncol(lee_day1_progeny_scores_mat)
  n_pw_sel <- length(top_pw_ld1)
  cell_vec_pw <- rep(colnames(lee_day1_progeny_scores_mat), each = n_pw_sel)
  pathway_vec_pw <- rep(top_pw_ld1, times = n_cells_pw)
  score_vec_pw <- as.vector(t(lee_day1_progeny_scores_mat[top_pw_ld1, , drop = FALSE]))
  pw_scores_long_ld1 <- data.frame(cell = cell_vec_pw, pathway = pathway_vec_pw, score = score_vec_pw, stringsAsFactors = FALSE)
  pw_scores_long_ld1$arg1_status <- ifelse(pw_scores_long_ld1$cell %in% lee_day1_arg1pos_cells, "Arg1pos", "Arg1neg")
  p_progeny_vln_ld1 <- ggplot(pw_scores_long_ld1, aes(x = arg1_status, y = score, fill = arg1_status)) + geom_violin(trim = FALSE) + geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) + facet_wrap(~pathway, scales = "free_y") + labs(title = "PROGENy pathway scores: Arg1+ vs Arg1- neutrophils (Lee Day 1)", x = "ARG1 status", y = "Activity score (norm_wmean)") + theme_minimal() + PLOT_TITLE_THEME
  print(p_progeny_vln_ld1)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "LeeDay1_PROGENy_TopPathways_Violin.png"), p_progeny_vln_ld1, width = 12, height = 8, dpi = 150)
}
progeny_heatmap_mat_ld1 <- tryCatch(rbind(Arg1pos = lee_day1_pathway_means_arg1pos, Arg1neg = lee_day1_pathway_means_arg1neg), error = function(e) matrix(0, nrow = 0, ncol = 0))
if (nrow(progeny_heatmap_mat_ld1) > 0 && ncol(progeny_heatmap_mat_ld1) > 0) {
  colors_progeny_ld1 <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors_use_progeny_ld1 <- grDevices::colorRampPalette(colors = colors_progeny_ld1)(100)
  p_progeny_heatmap_ld1 <- pheatmap::pheatmap(progeny_heatmap_mat_ld1, color = colors_use_progeny_ld1, border_color = "white", cellwidth = 20, cellheight = 20, main = "PROGENy Pathway Activity: Arg1+ vs Arg1- Neutrophils (Lee Day 1)")
  print(p_progeny_heatmap_ld1)
}
lee_day1_ligand_map_plot <- if (nrow(lee_day1_ligand_pathway_map) > 0 && "ligand" %in% names(lee_day1_ligand_pathway_map)) { top15_lig_d1 <- head(unique(lee_day1_ligand_pathway_map$ligand), 15); lee_day1_ligand_pathway_map[lee_day1_ligand_pathway_map$ligand %in% top15_lig_d1, ] } else data.frame()
p_lee_d1_progeny2 <- if (nrow(lee_day1_ligand_map_plot) > 0) ggplot(lee_day1_ligand_map_plot, aes(x = ligand, y = pathway, size = abs(weight), color = weight)) + geom_point(alpha = 0.7) + scale_color_gradient2(low = "blue", mid = "white", high = "red") + labs(title = "Top 15 Ligands Linked to PROGENy Pathways - Lee Day 1", x = "Ligand", y = "Pathway") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)) else ggplot() + theme_void() + labs(title = "No ligands from LIANA Day 1")
vals_w <- numeric(0)
map_tabnames <- c("lee_day1_ligand_pathway_map", "lee_day3_ligand_pathway_map", "wang_day3_ligand_pathway_map", "qin_day1_ligand_pathway_map", "qin_day3_ligand_pathway_map")
for (fn in map_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("weight" %in% names(d)) vals_w <- c(vals_w, d$weight)
}
vals_w <- vals_w[is.finite(vals_w)]
UNIFY_LIGW_ABS <- if (length(vals_w) > 0) range(abs(vals_w)) else c(0, 1)
if (UNIFY_LIGW_ABS[1] == UNIFY_LIGW_ABS[2]) UNIFY_LIGW_ABS[2] <- UNIFY_LIGW_ABS[1] + 1e-6
mxw <- if (length(vals_w) > 0) max(abs(vals_w)) else 1
if (!is.finite(mxw) || mxw <= 0) mxw <- 1
UNIFY_LIGW_COL <- c(-mxw, mxw)
print(p_lee_d1_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("✓ Lee Day 1 PROGENy pathway analysis complete")

# -------- Save environment and checkpoint RDS before BARTsc --------
# Restore full workspace: load(file.path(OUTPUT_DIR, "Workspace_BeforeBARTsc.RData"))
# Restore key objects only: checkpoint <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_BeforeBARTsc.rds")); list2env(checkpoint, envir = .GlobalEnv)
save.image(file.path(OUTPUT_DIR, "Workspace_BeforeBARTsc.RData"))
checkpoint_before_bartsc <- list(
  lee_day1 = lee_day1,
  lee_day1_neut_cells = lee_day1_neut_cells,
  lee_day1_liana_result_df = lee_day1_liana_result_df,
  lee_day1_liana_consensus = lee_day1_liana_consensus,
  lee_day1_arg1pos_received = lee_day1_arg1pos_received,
  lee_day1_pathway_means_arg1pos = lee_day1_pathway_means_arg1pos,
  lee_day1_pathway_means_arg1neg = lee_day1_pathway_means_arg1neg,
  lee_day1_pathway_results = lee_day1_pathway_results,
  progeny_network_lee_d1 = progeny_network_lee_d1,
  OUTPUT_DIR = OUTPUT_DIR
)
saveRDS(checkpoint_before_bartsc, file.path(OUTPUT_DIR, "Checkpoint_BeforeBARTsc.rds"))
# checkpoint_before_bartsc <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_BeforeBARTsc.rds")); list2env(checkpoint_before_bartsc, envir = .GlobalEnv)
print("✓ Environment and checkpoint saved (Workspace_BeforeBARTsc.RData, Checkpoint_BeforeBARTsc.rds)")

# ============================================================================
# FRESH RSTUDIO — RESUME LEE DAY 1 BARTsc + INTEGRATION (no full pipeline)
# ============================================================================
# You already ran through Lee Day 1 PROGENy once so these files exist:
#   MASTER_PIPELINE_RESULTS/Checkpoint_BeforeBARTsc.rds
#
# Step 1 — setwd() to the folder that CONTAINS MASTER_PIPELINE_RESULTS (not inside it).
#
# Step 2 — Run this block ONCE (uncomment or paste into console):
#
# OUTPUT_DIR <- "MASTER_PIPELINE_RESULTS"
# library(Seurat); library(ggplot2); library(BARTsc); library(viridis); library(scatterplot3d); library(reticulate)
# if (requireNamespace("OmnipathR", quietly = TRUE)) library(OmnipathR)
# if (requireNamespace("tidyr", quietly = TRUE)) library(tidyr)
# checkpoint_before_bartsc <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_BeforeBARTsc.rds"))
# list2env(checkpoint_before_bartsc, envir = .GlobalEnv)
# NEUTROPHIL_BASE_LABEL <- "Neutrophil"
# PLOT_TITLE_THEME <- theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12, margin = margin(t = 10)), axis.title.y = element_text(size = 12, margin = margin(r = 10)), legend.text = element_text(size = 12), legend.title = element_text(size = 12), plot.margin = margin(10, 10, 20, 10))
#
# Step 3 — In this script, select from the line:
#   bartsc_available <- requireNamespace("BARTsc", quietly = TRUE)
#   through the line that prints "Environment saved: Workspace_AfterLeeDay1_BARTsc"
#   (stop at the comment "# -------- Lee Day 3 LIANA --------"; do not run Lee Day 3 yet).
# ============================================================================

# -------- BARTsc TF Analysis: Lee Day 1 Neutrophils (Arg1pos vs Arg1neg; downstream TF support only) --------
# Official workflow (BARTsc vignettes/scRNA-seq.md): §1 load_bart2 — §3 bartsc + normalize_RNA — §4 find_signature_genes + find_pairwise_deg —
#   §5 run_signature_RNA then get_result("cell type signature") — §6 calc_crossCT_auc_RNA, crossCT_test, get_result("cross-cell-type") —
#   §7 find_key_regulators then get_result("Key regs ident"); saveRDS. Dot plots, heatmaps, and key_regulator_scatter run in the separate per-cohort BARTsc visualizations sections (reload *BARTsc_Object.rds).
# Pipeline gates crossCT_test/find_key_regulators when BARTSC_* + pairwise DEG counts fail.
# Runs after PROGENy so the cascade is: Ligand (LIANA) -> Pathway (PROGENy) -> TF (BARTsc) -> Gene expression.
# Target population (aligned with LIANA): NeutrophilArg1pos and NeutrophilArg1neg only.
# Official: https://github.com/hongpan-uva/BARTsc - install, then initialize() once, load_bart2() each session.
#
# DESIGN NOTE (tool choice): BARTsc is appropriate for SCI scRNA-seq here because it (1) infers cell-type-
# specific key TFs from expression, (2) compares two states (Arg1pos vs Arg1neg) as required, and (3) fits
# the ligand->pathway->TF cascade. Alternatives (SCENIC, DoRothEA/viper) could complement but do not
# replace BARTsc for this question. Key-regulator scatter uses key_regulator_scatter_unified (fixed axes + rank color scale; PNG + optional interactive redraw).
# crossCT_test + find_key_regulators: gated when BARTSC_* cell-count flags pass AND bart_deg_ok (enough pairwise DEG each direction); avoids BARTsc crossCT internal errors when Arg1neg signature is too weak. Order matches vignette §6–§7 (crossCT_test and cross-cell-type viz before find_key_regulators).
BARTSC_SKIP_CROSSCT_TEST <- FALSE
BARTSC_MIN_CELLS_CROSSCT <- 2L
BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N <- FALSE
# Shared BARTsc parameters for find_signature_genes + find_pairwise_deg (all datasets, all time points). Relaxed defaults improve sensitivity vs log2FC=1 / padj=0.05.
BART_MIN_PCT <- 0.1
BART_MIN_DIFF_PCT <- -Inf
BART_LOG2FC_THR <- 0.25
BART_PADJ_THR <- 0.1
BART_AUC_THR <- 0
# Minimum rows in each direction of @data$pairwise_DEG (Arg1pos::Arg1neg and Arg1neg::Arg1pos) required before crossCT_test / find_key_regulators. Set 0L to disable gate (not recommended for unstable cohorts).
BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT <- 5L
BARTSC_N_TF_DOTHEAT <- 5L
BARTSC_N_LABELED_TFS <- 5L
# Key-regulator 3D PNG size (key_regulator_scatter_unified; wider width for right-side TF labels).
BARTSC_KEYREG_SCATTER_W_IN <- 11
BARTSC_KEYREG_SCATTER_H_IN <- 7
# BARTsc default plot auto-scales axes and rank color; these fixed limits match across cohorts/datasets (adjust once for all figures).
BARTSC_KEYREG_XLIM <- c(-1, 1)
BARTSC_KEYREG_YLIM_SIG <- c(0, 7)
BARTSC_KEYREG_ZLIM_DMDR <- c(-2, 2)
BARTSC_KEYREG_RANK_MAX <- 60L
# Same logic as BARTsc::key_regulator_scatter (visualization.R) but: xlim/ylim_sig/zlim_dmdr fixed; rank color mapped to 1..rank_color_max; dashed leaders + numbered TF labels on the right.
key_regulator_scatter_unified <- function(object, mod, cell_type, tfs_labeled = NULL, pval.thr = 0.05,
  xlim = c(-1, 1), ylim_sig = c(0, 7), zlim_dmdr = c(-2, 2), rank_color_max = 60L,
  main = NULL, subtitle = "") {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) stop("Install package scatterplot3d (BARTsc dependency).")
  plot_df <- BARTsc::get_result(object, analysis = "Key regs ident", mod = mod)[[cell_type]]
  plot_df$role <- "other"
  if (!is.null(tfs_labeled) && length(tfs_labeled) > 0) {
    plot_df$role[which(plot_df$TF %in% tfs_labeled)] <- "labeled"
    plot_df$role <- factor(plot_df$role, levels = c("labeled", "other"))
  }
  n_keep <- length(which(plot_df$final_pvalue < pval.thr))
  print(paste0(n_keep, " significant key regulators were identified"))
  plot_df$label <- paste0(plot_df$final_rank, ":", plot_df$TF)
  plot_df$final_rank[which(plot_df$final_rank > n_keep)] <- NA
  plot_df$role2 <- rep("other", nrow(plot_df))
  plot_df$role2[which(plot_df$pval.thr < 0.05)] <- "key regulators"
  plot_df$role2 <- factor(plot_df$role2, levels = c("key regulators", "other"))
  label_df <- plot_df[which(plot_df$role == "labeled"), ]
  pal <- rev(viridis::plasma(50))
  rk <- plot_df$final_rank
  color_vec <- rep("grey", nrow(plot_df))
  fin <- is.finite(rk) & !is.na(rk)
  rk_c <- pmin(pmax(rk[fin], 1), rank_color_max)
  brks <- seq(1, rank_color_max, length.out = 51)
  bin <- findInterval(rk_c, brks, all.inside = TRUE)
  bin[bin < 1] <- 1L
  bin[bin > 50] <- 50L
  color_vec[fin] <- pal[bin]
  pm <- if (is.null(main)) cell_type else main
  p <- scatterplot3d::scatterplot3d(
    x = plot_df$MDR, y = plot_df$signature_score, z = plot_df$dMDR,
    xlab = "MDR", ylab = "signature_score", zlab = "dMDR",
    xlim = xlim, ylim = ylim_sig, zlim = zlim_dmdr,
    pch = 16, color = color_vec,
    cex.symbols = 1.5, cex.axis = 1.5, cex.lab = 1.5,
    main = pm
  )
  if (nzchar(subtitle)) title(sub = subtitle, cex.sub = 0.72)
  if (nrow(label_df) > 0) {
    label.coords <- p$xyz.convert(label_df$MDR, label_df$signature_score, label_df$dMDR)
    usr <- par("usr")
    seg_len <- (usr[2] - usr[1]) * 0.06
    i <- 1L
    while (i <= nrow(label_df)) {
      segments(label.coords$x[i], label.coords$y[i], label.coords$x[i] + seg_len, label.coords$y[i], lty = 2, col = "darkolivegreen3", lwd = 1.2)
      text(label.coords$x[i] + seg_len, label.coords$y[i], labels = paste0(i, ". ", label_df$TF[i]), pos = 4, cex = 1.35, col = "darkolivegreen4", xpd = NA)
      i <- i + 1L
    }
  }
  invisible(p)
}
# BARTsc get_result(..., "cross-cell-type") returns the deviation slot (named list: TF name -> deviation matrix). dot_plot/deviation_heatmap need object@resultsCrossCellType filled by crossCT_test (not only calc_crossCT_auc_RNA).
bartsc_available <- requireNamespace("BARTsc", quietly = TRUE)
if (!bartsc_available) {
  print("BARTsc not installed; entire BARTsc section skipped (no plots, no BARTsc_LeeDay1 outputs). Install: devtools::install_github('hongpan-uva/BARTsc')")
}
print("--- BARTsc TF Analysis: Lee Day 1 (Arg1pos vs Arg1neg neutrophils) ---")
# Secondary support only: BARTsc prioritizes downstream TF programs associated with the Arg1 state.
print(paste0("BARTsc shared thresholds (all datasets): min.pct=", BART_MIN_PCT, ", log2fc.thr=", BART_LOG2FC_THR, ", padj.thr=", BART_PADJ_THR, "; pairwise DEG asymmetry note threshold (per direction) = ", BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT))
bartsc_initialized <- FALSE
if (bartsc_available) {
  bartsc_initialized <- suppressWarnings(tryCatch({ BARTsc::load_bart2(); TRUE }, error = function(e) FALSE))
  if (!bartsc_initialized && exists("bart2", envir = .GlobalEnv)) bartsc_initialized <- TRUE
}
if (bartsc_available && !bartsc_initialized) print("BARTsc load_bart2() failed; skip BARTsc. Run BARTsc::initialize() once if first use.")
if (bartsc_initialized && exists("bart2", envir = .GlobalEnv)) {
  if (!exists("types", envir = .GlobalEnv)) types <<- reticulate::import("types")
  bart2_mods <- reticulate::py_to_r(reticulate::py_get_attr(bart2, "__all__"))
  for (m in bart2_mods) {
    if (!exists(m, envir = .GlobalEnv)) assign(m, reticulate::import(paste0("bart2.", m), delay_load = TRUE), envir = .GlobalEnv)
  }
}
if (bartsc_initialized) {
  options("mc.cores" = min(6L, parallel::detectCores()))
  lee_day1_neut_subset <- subset(lee_day1, cells = lee_day1_neut_cells)
  lee_day1_bart_rna_cnt <- Seurat::GetAssayData(lee_day1_neut_subset, layer = "counts")
  lee_day1_bart_label <- setNames(as.character(lee_day1_neut_subset$arg1_status), colnames(lee_day1_neut_subset))
  lee_day1_bart_label <- factor(lee_day1_bart_label, levels = c("Arg1pos", "Arg1neg"))
  lee_day1_bart_proj <- BARTsc::bartsc(name = "LeeDay1_Arg1", genome = "mm10", label = lee_day1_bart_label, cell_types_used = c("Arg1pos", "Arg1neg"), RNA_cnt_matrix = lee_day1_bart_rna_cnt)
  lee_day1_bart_proj <- BARTsc::normalize_RNA(lee_day1_bart_proj)
  lee_day1_bart_proj <- BARTsc::find_signature_genes(lee_day1_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  lee_day1_bart_proj <- BARTsc::find_pairwise_deg(lee_day1_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  lee_day1_pw <- lee_day1_bart_proj@data$pairwise_DEG
  n_ld1_posneg <- 0L
  n_ld1_negpos <- 0L
  if (!is.null(lee_day1_pw) && "Arg1pos::Arg1neg" %in% names(lee_day1_pw)) n_ld1_posneg <- { x <- lee_day1_pw[["Arg1pos::Arg1neg"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  if (!is.null(lee_day1_pw) && "Arg1neg::Arg1pos" %in% names(lee_day1_pw)) n_ld1_negpos <- { x <- lee_day1_pw[["Arg1neg::Arg1pos"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  print(paste0("Lee Day 1 BARTsc pairwise DEG count (@data$pairwise_DEG): Arg1pos::Arg1neg=", n_ld1_posneg, ", Arg1neg::Arg1pos=", n_ld1_negpos))
  bart_deg_ok_ld1 <- if (BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT > 0) n_ld1_posneg >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT && n_ld1_negpos >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT else TRUE
  bart_nt_ld1 <- table(lee_day1_bart_label)
  bart_unequal_ld1 <- length(bart_nt_ld1) == 2L && length(unique(as.vector(bart_nt_ld1))) > 1L
  bart_crossct_allowed_ld1 <- !isTRUE(BARTSC_SKIP_CROSSCT_TEST) && length(bart_nt_ld1) >= 2L && all(as.integer(bart_nt_ld1) >= BARTSC_MIN_CELLS_CROSSCT) && !(isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_ld1) && bart_deg_ok_ld1
  lee_day1_bart_proj <- BARTsc::run_signature_RNA(lee_day1_bart_proj)
  lee_day1_bart_sig <- BARTsc::get_result(lee_day1_bart_proj, analysis = "cell type signature", mod = "RNA")
  lee_day1_bart_proj <- BARTsc::calc_crossCT_auc_RNA(lee_day1_bart_proj)
  print(bart_nt_ld1)
  if (!bart_crossct_allowed_ld1) message(paste0(
    "BARTsc Lee Day 1: crossCT_test / find_key_regulators skipped (BARTSC_SKIP_CROSSCT_TEST=", isTRUE(BARTSC_SKIP_CROSSCT_TEST),
    "; min_cells_ok=", all(as.integer(bart_nt_ld1) >= BARTSC_MIN_CELLS_CROSSCT),
    "; skip_if_unequal_n=", isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_ld1,
    "; deg_ok=", bart_deg_ok_ld1,
    " [Arg1pos::Arg1neg=", n_ld1_posneg,
    ", Arg1neg::Arg1pos=", n_ld1_negpos,
    ", min_each_dir=", BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT,
    "]). Vignette scRNA-seq.md §6: dot_plot/deviation_heatmap follow crossCT_test."
  ))
  if (bart_crossct_allowed_ld1) lee_day1_bart_proj <- BARTsc::crossCT_test(lee_day1_bart_proj, mod = "RNA")
  lee_day1_bart_cross <- if (bart_crossct_allowed_ld1) BARTsc::get_result(lee_day1_bart_proj, analysis = "cross-cell-type", mod = "RNA") else list()
  lee_day1_bart_cross_dev <- if (is.null(lee_day1_bart_cross)) list() else if (is.list(lee_day1_bart_cross) && "deviation" %in% names(lee_day1_bart_cross) && is.list(lee_day1_bart_cross$deviation)) lee_day1_bart_cross$deviation else lee_day1_bart_cross
  bart_outdir <- file.path(OUTPUT_DIR, "BARTsc_LeeDay1")
  dir.create(bart_outdir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(lee_day1_bart_sig) && !is.null(lee_day1_bart_sig[["Arg1pos"]])) write.csv(lee_day1_bart_sig[["Arg1pos"]], file.path(bart_outdir, "LeeDay1_BARTsc_Signature_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(lee_day1_bart_sig) && !is.null(lee_day1_bart_sig[["Arg1neg"]])) write.csv(lee_day1_bart_sig[["Arg1neg"]], file.path(bart_outdir, "LeeDay1_BARTsc_Signature_Arg1neg.csv"), row.names = FALSE)
  lee_day1_bart_key <- NULL
  if (bart_crossct_allowed_ld1) {
    lee_day1_bart_proj <- BARTsc::find_key_regulators(lee_day1_bart_proj, mod = "RNA", min.N.profile = 3)
    lee_day1_bart_key <- BARTsc::get_result(lee_day1_bart_proj, analysis = "Key regs ident", mod = "RNA")
  }
  if (!is.null(lee_day1_bart_key) && !is.null(lee_day1_bart_key[["Arg1pos"]])) write.csv(lee_day1_bart_key[["Arg1pos"]], file.path(bart_outdir, "LeeDay1_BARTsc_KeyRegulators_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(lee_day1_bart_key) && !is.null(lee_day1_bart_key[["Arg1neg"]])) write.csv(lee_day1_bart_key[["Arg1neg"]], file.path(bart_outdir, "LeeDay1_BARTsc_KeyRegulators_Arg1neg.csv"), row.names = FALSE)
  saveRDS(lee_day1_bart_proj, file.path(bart_outdir, "LeeDay1_BARTsc_Object.rds"))
  print("✓ BARTsc TF analysis complete (Lee Day 1) — plots: run Lee Day 1 BARTsc visualizations section below")
}
if (!bartsc_initialized) print("BARTsc skipped (package not installed or load_bart2 failed). No BARTsc_LeeDay1/ analysis outputs; run BARTsc visualizations sections only after a successful analysis RDS exists.")

# -------- Lee Day 1: BARTsc visualizations (BARTsc::dot_plot, deviation_heatmap, key_regulator_scatter_unified; does not re-run bartsc/crossCT/find_key_regulators) --------
# Standalone: set OUTPUT_DIR, library(ggplot2), library(BARTsc), BARTsc::load_bart2(), then lee_day1_bart_proj <- readRDS(file.path(OUTPUT_DIR, "BARTsc_LeeDay1", "LeeDay1_BARTsc_Object.rds")) and run this block from bart_viz_ld1_rds onward.
# If you skip the main script's BARTsc parameter block, defaults below define PNG size and TF counts for this section only.
bart_viz_ld1_rds <- file.path(OUTPUT_DIR, "BARTsc_LeeDay1", "LeeDay1_BARTsc_Object.rds")
if (!requireNamespace("BARTsc", quietly = TRUE)) message("Lee Day 1 BARTsc visualizations skipped: package BARTsc not installed.")
if (requireNamespace("BARTsc", quietly = TRUE) && !file.exists(bart_viz_ld1_rds)) message(paste0("Lee Day 1 BARTsc visualizations skipped: missing ", bart_viz_ld1_rds))
if (requireNamespace("BARTsc", quietly = TRUE) && file.exists(bart_viz_ld1_rds)) {
  if (!exists("BARTSC_N_TF_DOTHEAT")) BARTSC_N_TF_DOTHEAT <- 5L
  if (!exists("BARTSC_N_LABELED_TFS")) BARTSC_N_LABELED_TFS <- 5L
  if (!exists("BARTSC_KEYREG_SCATTER_W_IN")) BARTSC_KEYREG_SCATTER_W_IN <- 11
  if (!exists("BARTSC_KEYREG_SCATTER_H_IN")) BARTSC_KEYREG_SCATTER_H_IN <- 7
  if (!exists("BARTSC_KEYREG_XLIM")) BARTSC_KEYREG_XLIM <- c(-1, 1)
  if (!exists("BARTSC_KEYREG_YLIM_SIG")) BARTSC_KEYREG_YLIM_SIG <- c(0, 7)
  if (!exists("BARTSC_KEYREG_ZLIM_DMDR")) BARTSC_KEYREG_ZLIM_DMDR <- c(-2, 2)
  if (!exists("BARTSC_KEYREG_RANK_MAX")) BARTSC_KEYREG_RANK_MAX <- 60L
  suppressWarnings(tryCatch({ BARTsc::load_bart2(); NULL }, error = function(e) NULL))
  if (exists("bart2", envir = .GlobalEnv)) {
    if (!exists("types", envir = .GlobalEnv)) types <<- reticulate::import("types")
    bart2_mods_viz_ld1 <- reticulate::py_to_r(reticulate::py_get_attr(bart2, "__all__"))
    for (m_viz_ld1 in bart2_mods_viz_ld1) {
      if (!exists(m_viz_ld1, envir = .GlobalEnv)) assign(m_viz_ld1, reticulate::import(paste0("bart2.", m_viz_ld1), delay_load = TRUE), envir = .GlobalEnv)
    }
  }
  if (!exists("lee_day1_bart_proj", envir = .GlobalEnv, inherits = FALSE)) lee_day1_bart_proj <- readRDS(bart_viz_ld1_rds)
  bart_outdir <- file.path(OUTPUT_DIR, "BARTsc_LeeDay1")
  dir.create(bart_outdir, showWarnings = FALSE, recursive = TRUE)
  if (interactive() && grDevices::dev.cur() == 1L) grDevices::dev.new()
  lee_day1_bart_cross_viz <- BARTsc::get_result(lee_day1_bart_proj, analysis = "cross-cell-type", mod = "RNA")
  lee_day1_bart_cross_dev <- if (is.null(lee_day1_bart_cross_viz)) list() else if (is.list(lee_day1_bart_cross_viz) && "deviation" %in% names(lee_day1_bart_cross_viz) && is.list(lee_day1_bart_cross_viz$deviation)) lee_day1_bart_cross_viz$deviation else lee_day1_bart_cross_viz
  lee_day1_bart_key <- BARTsc::get_result(lee_day1_bart_proj, analysis = "Key regs ident", mod = "RNA")
  bart_tf_names <- names(lee_day1_bart_cross_dev)
  tf_example <- if (length(bart_tf_names) > 0) bart_tf_names[1] else NULL
  if (length(lee_day1_bart_cross_dev) > 0 && !is.null(tf_example)) {
    p_bart_dot <- BARTsc::dot_plot(lee_day1_bart_proj, mod = "RNA", tf = tf_example, max_dot_size = 22)
    print(p_bart_dot)
    ggplot2::ggsave(file.path(bart_outdir, "LeeDay1_BARTsc_DotPlot.png"), p_bart_dot, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir, "LeeDay1_BARTsc_DotPlot.pdf"), p_bart_dot, width = 8, height = 6)
    p_bart_heat <- BARTsc::deviation_heatmap(lee_day1_bart_proj, mod = "RNA", tf = tf_example, tile_fontsize = 6)
    print(p_bart_heat)
    ggplot2::ggsave(file.path(bart_outdir, "LeeDay1_BARTsc_DeviationHeatmap.png"), p_bart_heat, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir, "LeeDay1_BARTsc_DeviationHeatmap.pdf"), p_bart_heat, width = 8, height = 6)
  }
  top_tf_pos <- if (length(lee_day1_bart_cross_dev) > 0 && !is.null(lee_day1_bart_key) && !is.null(lee_day1_bart_key[["Arg1pos"]]) && "TF" %in% colnames(lee_day1_bart_key[["Arg1pos"]])) head(lee_day1_bart_key[["Arg1pos"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_tf_neg <- if (length(lee_day1_bart_cross_dev) > 0 && !is.null(lee_day1_bart_key) && !is.null(lee_day1_bart_key[["Arg1neg"]]) && "TF" %in% colnames(lee_day1_bart_key[["Arg1neg"]])) head(lee_day1_bart_key[["Arg1neg"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  tfs_available <- names(lee_day1_bart_cross_dev)
  top_key_tfs <- intersect(unique(c(top_tf_pos, top_tf_neg)), tfs_available)
  if (length(lee_day1_bart_cross_dev) > 0 && length(top_key_tfs) > 0) {
    i_tf_ld1 <- 1L
    while (i_tf_ld1 <= length(top_key_tfs)) {
      tf_cur <- top_key_tfs[i_tf_ld1]
      p_dot_tf <- BARTsc::dot_plot(lee_day1_bart_proj, mod = "RNA", tf = tf_cur, max_dot_size = 22)
      print(p_dot_tf)
      ggplot2::ggsave(file.path(bart_outdir, paste0("LeeDay1_BARTsc_DotPlot_", tf_cur, ".png")), p_dot_tf, width = 8, height = 6)
      p_heat_tf <- BARTsc::deviation_heatmap(lee_day1_bart_proj, mod = "RNA", tf = tf_cur, tile_fontsize = 6)
      print(p_heat_tf)
      ggplot2::ggsave(file.path(bart_outdir, paste0("LeeDay1_BARTsc_DeviationHeatmap_", tf_cur, ".png")), p_heat_tf, width = 8, height = 6)
      i_tf_ld1 <- i_tf_ld1 + 1L
    }
  }
  tfs_labeled_arg1pos <- character(0)
  if (!is.null(lee_day1_bart_key) && !is.null(lee_day1_bart_key[["Arg1pos"]])) {
    df_kp <- lee_day1_bart_key[["Arg1pos"]]
    if (nrow(df_kp) > 0) {
      tf_cols_kp <- intersect(colnames(df_kp), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_kp <- if (length(tf_cols_kp) == 0) colnames(df_kp)[1] else tf_cols_kp[1]
      ord_kp <- if ("final_rank" %in% colnames(df_kp)) order(df_kp$final_rank, na.last = TRUE) else seq_len(nrow(df_kp))
      df_kp_ord <- df_kp[ord_kp, , drop = FALSE]
      tfs_labeled_arg1pos <- as.character(head(df_kp_ord[[tf_col_kp]], BARTSC_N_LABELED_TFS))
    }
  }
  tfs_labeled_arg1neg <- character(0)
  if (!is.null(lee_day1_bart_key) && !is.null(lee_day1_bart_key[["Arg1neg"]])) {
    df_kn <- lee_day1_bart_key[["Arg1neg"]]
    if (nrow(df_kn) > 0) {
      tf_cols_kn <- intersect(colnames(df_kn), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_kn <- if (length(tf_cols_kn) == 0) colnames(df_kn)[1] else tf_cols_kn[1]
      ord_kn <- if ("final_rank" %in% colnames(df_kn)) order(df_kn$final_rank, na.last = TRUE) else seq_len(nrow(df_kn))
      df_kn_ord <- df_kn[ord_kn, , drop = FALSE]
      tfs_labeled_arg1neg <- as.character(head(df_kn_ord[[tf_col_kn]], BARTSC_N_LABELED_TFS))
    }
  }
  if (length(tfs_labeled_arg1pos) > 0) {
    grDevices::png(file.path(bart_outdir, "LeeDay1_BARTsc_KeyRegScatter_Official_Arg1pos.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(lee_day1_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_labeled_arg1pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_labeled_arg1pos) > 0) message("Displaying official key_regulator_scatter for Arg1pos...")
  if (interactive() && length(tfs_labeled_arg1pos) > 0) key_regulator_scatter_unified(lee_day1_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_labeled_arg1pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
  if (length(tfs_labeled_arg1neg) > 0) {
    grDevices::png(file.path(bart_outdir, "LeeDay1_BARTsc_KeyRegScatter_Official_Arg1neg.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(lee_day1_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_labeled_arg1neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_labeled_arg1neg) > 0) message("Displaying official key_regulator_scatter for Arg1neg...")
  if (interactive() && length(tfs_labeled_arg1neg) > 0) key_regulator_scatter_unified(lee_day1_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_labeled_arg1neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
  print("✓ Lee Day 1 BARTsc visualizations complete")
}

# ============================================================================
# INTEGRATION: LIANA -> BARTsc + PROGENy (Lee Day 1) — DATA-DRIVEN
# ============================================================================
# Linear flow (read top to bottom): LIANA-expanded receptors → PROGENy footprint co-membership CSV → OmniPath receptor→TF (uses BARTsc key TFs when present) → full-chain CSV → sender×pathway pheatmap. BARTsc plots: section above (or reload RDS).
# OmniPath receptor->TF is the mechanistic overlay. PROGENy gives pathway activity Arg1pos vs Arg1neg only.
# Exploratory table below: LIANA receptor genes that are also PROGENy footprint *targets* for pathways with
# higher mean inferred activity in Arg1pos — co-membership, NOT "receptor activates pathway" (PROGENy is pathway->gene).
# ============================================================================

print("--- Integration: LIANA -> BARTsc + PROGENy (Lee Day 1, Data-Driven) ---")
# Exploratory overlay only: LIANA provides incoming-signal candidates; PROGENy/BARTsc/OmniPath add downstream support, not causal proof.

if (!exists("progeny_network_lee_d1") || !all(c("source", "target") %in% colnames(progeny_network_lee_d1))) {
  stop("Integration requires progeny_network_lee_d1 (source = pathway, target = gene). Run PROGENy section first.")
}
liana_receptors_ld1 <- unique(lee_day1_arg1pos_received$receptor.complex)
liana_receptors_ld1_split <- unique(trimws(unlist(strsplit(liana_receptors_ld1, "[_+]"))))
print(paste("Unique receptors received by Arg1+ neutrophils (LIANA):", length(liana_receptors_ld1_split)))

liana_expanded_ld1 <- lee_day1_arg1pos_received
liana_expanded_ld1$receptor_subunits <- lapply(strsplit(liana_expanded_ld1$receptor.complex, "[_+]"), function(x) trimws(x))
liana_expanded_ld1 <- tidyr::unnest_longer(liana_expanded_ld1, receptor_subunits)

pathways_active_ld1 <- names(lee_day1_pathway_means_arg1pos)[lee_day1_pathway_means_arg1pos > lee_day1_pathway_means_arg1neg]

# Exploratory footprint overlap (not causal receptor -> pathway activation)
receptor_pathway_links_ld1 <- progeny_network_lee_d1[
  progeny_network_lee_d1$target %in% liana_receptors_ld1_split &
  progeny_network_lee_d1$source %in% pathways_active_ld1,
]
receptor_pathway_links_ld1 <- receptor_pathway_links_ld1[, c("target", "source")]
colnames(receptor_pathway_links_ld1) <- c("receptor_gene", "pathway")

pathway_diff_ld1 <- lee_day1_pathway_means_arg1pos - lee_day1_pathway_means_arg1neg
receptor_pathway_links_ld1$pathway_activity_diff <- pathway_diff_ld1[receptor_pathway_links_ld1$pathway]

integration_receptor_pathway_ld1 <- merge(
  liana_expanded_ld1,
  receptor_pathway_links_ld1,
  by.x = "receptor_subunits",
  by.y = "receptor_gene",
  all.x = TRUE
)
integration_receptor_pathway_ld1 <- integration_receptor_pathway_ld1[!is.na(integration_receptor_pathway_ld1$pathway), ]

write.csv(integration_receptor_pathway_ld1, file.path(OUTPUT_DIR, "LeeDay1_LIANA_PROGENy_FootprintCoMembership_Exploratory.csv"), row.names = FALSE)
print(paste("LIANA x PROGENy footprint co-membership (exploratory):", nrow(integration_receptor_pathway_ld1), "rows (not receptor->pathway activation)"))

receptor_tf_links_ld1 <- data.frame(receptor_gene = character(0), tf = character(0), omnipath_via = character(0), omnipath_hops = integer(0), stringsAsFactors = FALSE)
active_tfs_ld1 <- character(0)

has_bart_key_ld1 <- exists("lee_day1_bart_key") &&
  !is.null(lee_day1_bart_key[["Arg1pos"]]) &&
  nrow(lee_day1_bart_key[["Arg1pos"]]) > 0

if (has_bart_key_ld1) {
  bart_tfs_ld1 <- lee_day1_bart_key[["Arg1pos"]]

  # robustly pick the TF column
  tf_col_ld1 <- intersect(colnames(bart_tfs_ld1),
                          c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
  if (length(tf_col_ld1) == 0) {
    tf_col_ld1 <- colnames(bart_tfs_ld1)[1]
  } else {
    tf_col_ld1 <- tf_col_ld1[1]
  }

  active_tfs_ld1 <- unique(as.character(bart_tfs_ld1[[tf_col_ld1]]))
}

omnipath_ok <- requireNamespace("OmnipathR", quietly = TRUE) && length(active_tfs_ld1) > 0

pathway_interactions_ld1 <- data.frame()
omni_ld1_first <- if (!omnipath_ok) list(dat = NULL, err1 = NA_character_) else tryCatch(list(dat = OmnipathR::import_omnipath_interactions(datasets = c("omnipath", "pathwayextra"), organism = 10090, genesymbols = TRUE), err1 = NA_character_), error = function(e1) list(dat = NULL, err1 = conditionMessage(e1)))
if (omnipath_ok) pathway_interactions_ld1 <- omni_ld1_first$dat
omnipath_ld1_err1 <- omni_ld1_first$err1
if (omnipath_ok && is.null(pathway_interactions_ld1)) pathway_interactions_ld1 <- tryCatch(OmnipathR::import_pathwayextra_interactions(organism = 10090, genesymbols = TRUE), error = function(e2) { message("OmniPath curated+pathwayextra import failed: ", omnipath_ld1_err1, "; fallback: ", conditionMessage(e2)); data.frame() })

omni_ld1_cols_ok <- omnipath_ok && nrow(pathway_interactions_ld1) > 0 && all(c("source_genesymbol", "target_genesymbol") %in% colnames(pathway_interactions_ld1))
omni_direct <- if (omni_ld1_cols_ok) subset(pathway_interactions_ld1, source_genesymbol %in% liana_receptors_ld1_split & target_genesymbol %in% active_tfs_ld1, select = c("source_genesymbol", "target_genesymbol")) else data.frame()
if (omni_ld1_cols_ok && nrow(omni_direct) > 0) receptor_tf_links_ld1 <- rbind(receptor_tf_links_ld1, data.frame(receptor_gene = omni_direct$source_genesymbol, tf = omni_direct$target_genesymbol, omnipath_via = NA_character_, omnipath_hops = 1L, stringsAsFactors = FALSE))

hop1 <- if (omni_ld1_cols_ok) pathway_interactions_ld1[pathway_interactions_ld1$source_genesymbol %in% liana_receptors_ld1_split, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_ld1_cols_ok) colnames(hop1) <- c("receptor_gene", "via")
hop2 <- if (omni_ld1_cols_ok) pathway_interactions_ld1[pathway_interactions_ld1$target_genesymbol %in% active_tfs_ld1, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_ld1_cols_ok) colnames(hop2) <- c("via", "tf")
omni_2hop <- if (omni_ld1_cols_ok) merge(hop1, hop2, by = "via") else data.frame()
if (omni_ld1_cols_ok && nrow(omni_2hop) > 0) omni_2hop <- unique(omni_2hop[, c("receptor_gene", "tf", "via")])
if (omni_ld1_cols_ok && nrow(omni_2hop) > 0) receptor_tf_links_ld1 <- rbind(receptor_tf_links_ld1, data.frame(receptor_gene = omni_2hop$receptor_gene, tf = omni_2hop$tf, omnipath_via = omni_2hop$via, omnipath_hops = 2L, stringsAsFactors = FALSE))
if (omni_ld1_cols_ok) receptor_tf_links_ld1 <- receptor_tf_links_ld1[!duplicated(paste(receptor_tf_links_ld1$receptor_gene, receptor_tf_links_ld1$tf)), ]

integration_receptor_tf_ld1 <- data.frame()
if (nrow(receptor_tf_links_ld1) > 0) {
  integration_receptor_tf_ld1 <- merge(
    liana_expanded_ld1,
    receptor_tf_links_ld1,
    by.x = "receptor_subunits",
    by.y = "receptor_gene",
    all.x = TRUE
  )
  integration_receptor_tf_ld1 <- integration_receptor_tf_ld1[!is.na(integration_receptor_tf_ld1$tf), ]
}

integration_tf_ld1_nonempty <- nrow(integration_receptor_tf_ld1) > 0
if (integration_tf_ld1_nonempty) write.csv(integration_receptor_tf_ld1, file.path(OUTPUT_DIR, "LeeDay1_LIANA_BARTsc_Integration_DataDriven.csv"), row.names = FALSE)
if (integration_tf_ld1_nonempty) n1 <- sum(integration_receptor_tf_ld1$omnipath_hops == 1L, na.rm = TRUE)
if (integration_tf_ld1_nonempty) n2 <- sum(integration_receptor_tf_ld1$omnipath_hops == 2L, na.rm = TRUE)
if (integration_tf_ld1_nonempty) print(paste0("LIANA -> BARTsc integration: ", nrow(integration_receptor_tf_ld1), " receptor-TF rows (OmniPath curated+pathwayextra; direct=", n1, ", two-hop=", n2, ")"))

receptor_tf_ld1_fallback <- nrow(receptor_tf_links_ld1) == 0 && length(active_tfs_ld1) > 0
if (receptor_tf_ld1_fallback) fallback_ld1 <- data.frame(receptor = liana_receptors_ld1_split, note = "Active TFs in Arg1pos (no OmniPath direct/two-hop link):", active_tfs = paste(active_tfs_ld1, collapse = "; "), stringsAsFactors = FALSE)
if (receptor_tf_ld1_fallback) write.csv(fallback_ld1, file.path(OUTPUT_DIR, "LeeDay1_LIANA_BARTsc_Receptors_and_TFs.csv"), row.names = FALSE)
if (receptor_tf_ld1_fallback) print("LIANA receptors and BARTsc TFs saved separately (no OmniPath link)")

if (nrow(receptor_tf_links_ld1) == 0) message("Note: No receptor-TF pairs found in OmniPath curated+pathwayextra (direct or two-hop via one intermediate). Signaling is often longer than two hops; the main evidence chain remains: Receptors (LIANA) -> Pathways (PROGENy) -> TFs (BARTsc)")

has_pathway_ld1 <- nrow(integration_receptor_pathway_ld1) > 0
has_tf_ld1 <- nrow(integration_receptor_tf_ld1) > 0
full_chain_ld1 <- data.frame()
if (has_pathway_ld1 && has_tf_ld1) full_chain_ld1 <- merge(integration_receptor_pathway_ld1, integration_receptor_tf_ld1[, c("source", "ligand.complex", "receptor_subunits", "aggregate_rank", "tf", "omnipath_via", "omnipath_hops")], by = c("source", "ligand.complex", "receptor_subunits", "aggregate_rank"), all = TRUE)
if (has_pathway_ld1 && !has_tf_ld1) { full_chain_ld1 <- integration_receptor_pathway_ld1; full_chain_ld1$tf <- NA_character_ }
if (!has_pathway_ld1 && has_tf_ld1) { full_chain_ld1 <- integration_receptor_tf_ld1; full_chain_ld1$pathway <- NA_character_; full_chain_ld1$pathway_activity_diff <- NA_real_ }
if (nrow(full_chain_ld1) > 0) full_chain_ld1$evidence_level <- ifelse(!is.na(full_chain_ld1$pathway) & !is.na(full_chain_ld1$tf), "STRONG (pathway + TF linked)", ifelse(!is.na(full_chain_ld1$pathway) | !is.na(full_chain_ld1$tf), "MODERATE (pathway or TF linked)", "WEAK (no link)"))
if (nrow(full_chain_ld1) > 0) full_chain_ld1 <- full_chain_ld1[order(full_chain_ld1$evidence_level, full_chain_ld1$aggregate_rank), ]
# Legacy filename retained for compatibility; contents are an exploratory overlay, not a causal signal chain.
if (nrow(full_chain_ld1) > 0) write.csv(full_chain_ld1, file.path(OUTPUT_DIR, "LeeDay1_FullSignalChain_DataDriven.csv"), row.names = FALSE)
if (nrow(full_chain_ld1) > 0) {
  cols_show_ld1 <- intersect(c("source", "receptor_subunits", "pathway", "tf", "omnipath_via", "omnipath_hops", "evidence_level"), colnames(full_chain_ld1))
  print(head(full_chain_ld1[, cols_show_ld1, drop = FALSE], 20))
}
heatmap_data_ld1 <- data.frame()
if (nrow(full_chain_ld1) > 0 && has_pathway_ld1) heatmap_data_ld1 <- as.data.frame.matrix(table(integration_receptor_pathway_ld1$source, integration_receptor_pathway_ld1$pathway))
if (nrow(heatmap_data_ld1) > 0 && ncol(heatmap_data_ld1) > 0) pheatmap::pheatmap(heatmap_data_ld1, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("white", "blue", "red"))(50), main = "Lee Day 1: Sender x pathway (PROGENy footprint co-membership, exploratory)")
if (nrow(full_chain_ld1) == 0) print("No integration rows (PROGENy/BARTsc may not overlap LIANA receptors)")
print("✓ Lee Day 1 LIANA <-> BARTsc <-> PROGENy integration complete (exploratory overlay)")
# Save environment after Lee Day 1 BARTsc + integration (resume from here: load(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay1_BARTsc.RData")))
save.image(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay1_BARTsc.RData"))
saveRDS(list(lee_day1 = lee_day1, lee_day1_liana_result_df = lee_day1_liana_result_df, lee_day1_pathway_results = lee_day1_pathway_results, progeny_network_lee_d1 = progeny_network_lee_d1, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay1_BARTsc.rds"))
# checkpoint_ld1_bartsc <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay1_BARTsc.rds")); list2env(checkpoint_ld1_bartsc, envir = .GlobalEnv)
print("✓ Environment saved: Workspace_AfterLeeDay1_BARTsc.RData, Checkpoint_AfterLeeDay1_BARTsc.rds")

# -------- Lee Day 3: LIANA → LIANA plots → CellChat/CellCall overlap → top-receptor Vln → PROGENy → BARTsc → integration (same cohort order as Lee Day 1) --------
# LIANA core: liana_wrap -> liana_aggregate; custom LR when LIANA_USE_CUSTOM_LR.
# STANDALONE: You can run only this section. It (1) loads Checkpoint_AfterLeeDay3_LIANA.rds if missing lee_day3 or if LEE_DAY3_LOAD_LIANA_CHECKPOINT is TRUE, else (2) builds lee_day3 from LEE_DAT_RDS_STANDALONE (default LeeDat.rds in getwd()). Setwd to your data folder; override paths below if needed.
if (!exists("OUTPUT_DIR")) OUTPUT_DIR <- "MASTER_PIPELINE_RESULTS"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!exists("LEE_DAY3_LOAD_LIANA_CHECKPOINT")) LEE_DAY3_LOAD_LIANA_CHECKPOINT <- FALSE
ck_ld3_liana_path <- file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_LIANA.rds")
load_ld3_ck <- file.exists(ck_ld3_liana_path) && (isTRUE(LEE_DAY3_LOAD_LIANA_CHECKPOINT) || !exists("lee_day3"))
if (load_ld3_ck) list2env(readRDS(ck_ld3_liana_path), envir = .GlobalEnv)
if (exists("lee_day3") && exists("lee_day3_pos") && exists("lee_day3_neg") && (!exists("lee_day3_arg1_status") || length(lee_day3_arg1_status) != ncol(lee_day3))) {
  lee_day3_arg1_status <- rep("Arg1neg", ncol(lee_day3))
  lee_day3_arg1_status[lee_day3_pos] <- "Arg1pos"
  lee_day3$arg1_status <- lee_day3_arg1_status
}
lee_d3_standalone_build <- !exists("lee_day3")
if (lee_d3_standalone_build && !exists("LEE_DAT_RDS_STANDALONE")) LEE_DAT_RDS_STANDALONE <- "LeeDat.rds"
if (lee_d3_standalone_build && !file.exists(LEE_DAT_RDS_STANDALONE)) stop(paste0("Lee Day 3 LIANA standalone: no lee_day3; checkpoint not found or incomplete at ", ck_ld3_liana_path, "; and Lee RDS not found: ", LEE_DAT_RDS_STANDALONE, ". Use setwd() to folder containing LeeDat.rds or set LEE_DAT_RDS_STANDALONE <- '/full/path/LeeDat.rds' before this section."))
if (lee_d3_standalone_build && !exists("NEUTROPHIL_BASE_LABEL")) NEUTROPHIL_BASE_LABEL <- "Neutrophil"
if (lee_d3_standalone_build && !exists("NEUTROPHIL_ALIASES")) NEUTROPHIL_ALIASES <- c("cNeutrophil", "Neutrophils")
if (lee_d3_standalone_build && !exists("GENE_CANDIDATES")) GENE_CANDIDATES <- c("Arg1", "ARG1")
if (lee_d3_standalone_build && !requireNamespace("Seurat", quietly = TRUE)) stop("Install Seurat for Lee Day 3 LIANA standalone")
if (lee_d3_standalone_build) library(Seurat)
if (lee_d3_standalone_build) lee_dat <- readRDS(LEE_DAT_RDS_STANDALONE)
if (lee_d3_standalone_build) Seurat::DefaultAssay(lee_dat) <- "RNA"
if (lee_d3_standalone_build) lee_dat$time <- as.character(lee_dat$time)
if (lee_d3_standalone_build) lee_dat$time <- gsub("Uninjured", "0", lee_dat$time)
if (lee_d3_standalone_build) lee_dat$time <- gsub("1dpi", "1", lee_dat$time)
if (lee_d3_standalone_build) lee_dat$time <- gsub("3dpi", "3", lee_dat$time)
if (lee_d3_standalone_build) lee_dat$time <- gsub("7dpi", "7", lee_dat$time)
if (lee_d3_standalone_build) lee_dat$time <- as.numeric(lee_dat$time)
if (lee_d3_standalone_build) lee_day3 <- subset(lee_dat, subset = time == 3)
if (lee_d3_standalone_build) Seurat::DefaultAssay(lee_day3) <- "RNA"
if (lee_d3_standalone_build) lee_day3_pruned_labels <- as.character(lee_day3$celltype)
if (lee_d3_standalone_build) lee_day3_pruned_labels <- gsub("-", "", lee_day3_pruned_labels)
if (lee_d3_standalone_build) lee_day3_pruned_labels[lee_day3_pruned_labels %in% NEUTROPHIL_ALIASES] <- NEUTROPHIL_BASE_LABEL
if (lee_d3_standalone_build) lee_day3_neut <- which(lee_day3_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL))
if (lee_d3_standalone_build) lee_day3_gene_arg1 <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(lee_day3)][1]
if (lee_d3_standalone_build) stopifnot(!is.na(lee_day3_gene_arg1) && nchar(lee_day3_gene_arg1) > 0)
if (lee_d3_standalone_build) lee_day3_expr_arg1 <- as.numeric(Seurat::GetAssayData(lee_day3, layer = "data")[lee_day3_gene_arg1, lee_day3_neut, drop = FALSE])
if (lee_d3_standalone_build) lee_day3_indicator <- lee_day3_expr_arg1 > 0
if (lee_d3_standalone_build) lee_day3_indicator[is.na(lee_day3_indicator)] <- FALSE
if (lee_d3_standalone_build) lee_day3_pos <- lee_day3_neut[lee_day3_indicator]
if (lee_d3_standalone_build) lee_day3_neg <- setdiff(lee_day3_neut, lee_day3_pos)
if (lee_d3_standalone_build) lee_day3_arg1_status <- rep("Arg1neg", ncol(lee_day3))
if (lee_d3_standalone_build) lee_day3_arg1_status[lee_day3_pos] <- "Arg1pos"
if (lee_d3_standalone_build) lee_day3$arg1_status <- lee_day3_arg1_status
if (lee_d3_standalone_build) Seurat::Idents(lee_day3) <- factor(lee_day3_pruned_labels)
if (lee_d3_standalone_build) lee_day3_neut_cells <- colnames(lee_day3)[c(lee_day3_pos, lee_day3_neg)]
if (lee_d3_standalone_build) print(paste0("Lee Day 3 LIANA standalone: built lee_day3 from ", LEE_DAT_RDS_STANDALONE))
if (!exists("NEUTROPHIL_BASE_LABEL")) NEUTROPHIL_BASE_LABEL <- "Neutrophil"
if (!exists("NEUTROPHIL_ALIASES")) NEUTROPHIL_ALIASES <- c("cNeutrophil", "Neutrophils")
if (!exists("LIANA_METHODS")) {
  LIANA_METHODS <- c("natmi", "connectome", "logfc", "sca", "cellphonedb")
  if (requireNamespace("CytoTalk", quietly = TRUE)) LIANA_METHODS <- c(LIANA_METHODS, "cytotalk")
}
if (!exists("NEUTROPHIL_STATES")) NEUTROPHIL_STATES <- c(paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos"), paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg"))
if (!exists("LIANA_MIN_CELLS")) LIANA_MIN_CELLS <- 3L
if (!exists("LIANA_TOP_RECEPTOR_VLN")) LIANA_TOP_RECEPTOR_VLN <- 5L
if (!exists("LIANA_USE_CUSTOM_LR")) LIANA_USE_CUSTOM_LR <- TRUE
if (!exists("CELLCHAT_LEE_DAY3_RDS")) CELLCHAT_LEE_DAY3_RDS <- "D:/Active Analysis/Mustafa/post-communication/cellchatThreeLee_ARG1_complete.rds"
if (!exists("CELLCALL_LEE_DAY3_RDS")) CELLCALL_LEE_DAY3_RDS <- "D:/Active Analysis/Mustafa/post-communication/CellCall_LEE_Day3_GLOBAL.rds"
if (!exists("LIANA_DOTPLOT_SIZE_RANGE")) LIANA_DOTPLOT_SIZE_RANGE <- c(2, 10)
if (!exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE <- c(2, 8)
if (!exists("PLOT_TITLE_THEME")) {
  PLOT_TITLE_THEME <- ggplot2::theme(
    plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 10)),
    axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 10)),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 12),
    plot.margin = ggplot2::margin(10, 10, 20, 10)
  )
}
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Install Seurat")
if (!requireNamespace("liana", quietly = TRUE)) stop("Install liana")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("Install SingleCellExperiment")
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap")
if (!requireNamespace("tidyr", quietly = TRUE)) stop("Install tidyr")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Install org.Mm.eg.db (Bioconductor)")
library(Seurat)
library(ggplot2)
library(liana)
library(SingleCellExperiment)
library(pheatmap)
library(grid)
library(org.Mm.eg.db)
library(tidyr)
lee_day3_liana_prereq <- c("lee_day3", "lee_day3_pruned_labels", "lee_day3_pos", "lee_day3_neg", "lee_day3_neut_cells")
lee_day3_liana_missing <- lee_day3_liana_prereq[!vapply(lee_day3_liana_prereq, exists, FUN.VALUE = logical(1))]
if (length(lee_day3_liana_missing) > 0) stop(paste0("Lee Day 3 LIANA internal error, still missing: ", paste(lee_day3_liana_missing, collapse = ", ")))
print("--- LIANA Analysis: Lee Day 3 ---")

## 1) Build LIANA labels (NeutrophilArg1pos / NeutrophilArg1neg / Other)
lee_day3_liana_labels <- as.character(lee_day3_pruned_labels)
lee_day3_liana_labels[lee_day3_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
lee_day3_liana_labels[lee_day3_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
lee_day3_liana_labels[is.na(lee_day3_liana_labels) | lee_day3_liana_labels == ""] <- "Other"

lee_day3$liana_label <- factor(lee_day3_liana_labels)
Seurat::Idents(lee_day3) <- lee_day3$liana_label

## 2) Clean Seurat for LIANA: RNA only, drop other assays/reductions/graphs
SeuratObject::DefaultAssay(lee_day3) <- "RNA"
for (assay_name in names(lee_day3@assays)) {
  if (assay_name != "RNA") lee_day3[[assay_name]] <- NULL
}
for (red_name in names(lee_day3@reductions)) lee_day3[[red_name]] <- NULL
for (graph_name in names(lee_day3@graphs)) lee_day3[[graph_name]] <- NULL
## 2b) Custom L-R from CellChat + CellCall (Lee Day 3) — same logic as Lee Day 1
lee_day3_liana_resource <- "MouseConsensus"
lee_day3_liana_external <- NULL
custom_lr_list_d3 <- list()
path_cc_d3_liana <- CELLCHAT_LEE_DAY3_RDS
path_ccall_d3_liana <- CELLCALL_LEE_DAY3_RDS
cc_d3_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_cc_d3_liana) && nzchar(trimws(path_cc_d3_liana)) && file.exists(path_cc_d3_liana) && requireNamespace("CellChat", quietly = TRUE)
cc_d3_liana <- NULL
if (cc_d3_ok) cc_d3_liana <- readRDS(path_cc_d3_liana)
comm_d3_liana <- NULL
if (!is.null(cc_d3_liana) && inherits(cc_d3_liana, "CellChat")) comm_d3_liana <- CellChat::subsetCommunication(cc_d3_liana)
has_comm_d3 <- !is.null(comm_d3_liana) && nrow(comm_d3_liana) > 0 && "ligand" %in% colnames(comm_d3_liana) && "receptor" %in% colnames(comm_d3_liana)
if (has_comm_d3) {
  cc_lr_d3 <- data.frame(source_genesymbol = comm_d3_liana$ligand, target_genesymbol = comm_d3_liana$receptor, stringsAsFactors = FALSE)
  cc_lr_d3 <- cc_lr_d3[!duplicated(cc_lr_d3), ]
  custom_lr_list_d3[["CellChat"]] <- cc_lr_d3
}
ccall_d3_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_ccall_d3_liana) && nzchar(trimws(path_ccall_d3_liana)) && file.exists(path_ccall_d3_liana)
ccall_d3_liana <- NULL
if (ccall_d3_ok) ccall_d3_liana <- readRDS(path_ccall_d3_liana)
has_ccall_lr_d3 <- !is.null(ccall_d3_liana) && !is.null(ccall_d3_liana@data$expr_l_r_log2_scale)
if (has_ccall_lr_d3) {
  lr_rownames_d3_lr <- rownames(ccall_d3_liana@data$expr_l_r_log2_scale)
  ccall_lig_d3_lr <- sub("-.*", "", lr_rownames_d3_lr)
  ccall_rec_d3_lr <- sub("^[^-]+-", "", lr_rownames_d3_lr)
  ccall_rec_d3_lr <- gsub("-", "_", ccall_rec_d3_lr)
  ccall_lr_d3 <- data.frame(source_genesymbol = ccall_lig_d3_lr, target_genesymbol = ccall_rec_d3_lr, stringsAsFactors = FALSE)
  ccall_lr_d3 <- ccall_lr_d3[nzchar(ccall_lr_d3$target_genesymbol), ]
  ccall_lr_d3 <- unique(ccall_lr_d3)
  custom_lr_list_d3[["CellCall"]] <- ccall_lr_d3
}
has_custom_lr_d3 <- length(custom_lr_list_d3) > 0
if (has_custom_lr_d3) {
  custom_lr_d3 <- do.call(rbind, custom_lr_list_d3)
  custom_lr_d3 <- custom_lr_d3[!duplicated(custom_lr_d3[, c("source_genesymbol", "target_genesymbol")]), ]
  lee_day3_liana_external <- data.frame(source_genesymbol = custom_lr_d3$source_genesymbol, target_genesymbol = custom_lr_d3$target_genesymbol, stringsAsFactors = FALSE)
  lee_day3_liana_resource <- "custom"
  print(paste("Loaded", nrow(lee_day3_liana_external), "custom L-R pairs from CellChat + CellCall (Lee Day 3)"))
}
## 3) SCE + liana_wrap (min_cells = LIANA_MIN_CELLS, base = exp(1); same as Lee Day 1)
lee_day3_counts <- tryCatch(Seurat::GetAssayData(lee_day3, slot = "counts"), error = function(e) Seurat::GetAssayData(lee_day3, layer = "counts"))
lee_day3_logcounts <- tryCatch(Seurat::GetAssayData(lee_day3, slot = "data"), error = function(e) Seurat::GetAssayData(lee_day3, layer = "data"))
lee_day3_sce <- SingleCellExperiment(assays = list(counts = lee_day3_counts, logcounts = lee_day3_logcounts))
lee_day3_sce$liana_label <- lee_day3$liana_label
SingleCellExperiment::colLabels(lee_day3_sce) <- lee_day3_sce$liana_label
liana_args_ld3 <- list(sce = lee_day3_sce, method = LIANA_METHODS, resource = lee_day3_liana_resource, idents_col = "liana_label", expr_prop = 0.05, verbose = TRUE, min_cells = LIANA_MIN_CELLS, base = exp(1))
if (!is.null(lee_day3_liana_external)) liana_args_ld3$external_resource <- lee_day3_liana_external
lee_day3_liana_result <- do.call(liana::liana_wrap, liana_args_ld3)
lee_day3_liana_result_df <- liana::liana_aggregate(lee_day3_liana_result)

lee_day3_liana_neutrophil_target <- dplyr::filter(lee_day3_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))
lee_day3_liana_neutrophil_source <- dplyr::filter(lee_day3_liana_result_df, source %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))
lee_day3_liana_neutrophil <- lee_day3_liana_neutrophil_target
lee_day3_arg1pos_received <- dplyr::filter(lee_day3_liana_result_df, target == "NeutrophilArg1pos")
lee_day3_arg1pos_received <- dplyr::arrange(lee_day3_arg1pos_received, aggregate_rank)
lee_day3_arg1pos_top_signals <- head(lee_day3_arg1pos_received, 50)
lee_day3_arg1neg_received <- dplyr::filter(lee_day3_liana_result_df, target == "NeutrophilArg1neg")
lee_day3_arg1neg_received <- dplyr::arrange(lee_day3_arg1neg_received, aggregate_rank)
topN <- 50
lee_day3_arg1pos_top <- head(lee_day3_arg1pos_received, topN)
lee_day3_arg1neg_top <- head(lee_day3_arg1neg_received, topN)
pos_pairs <- paste0(lee_day3_arg1pos_top$ligand.complex, "_", lee_day3_arg1pos_top$receptor.complex)
neg_pairs_all <- paste0(lee_day3_arg1neg_received$ligand.complex, "_", lee_day3_arg1neg_received$receptor.complex)
neg_rank_lookup <- setNames(lee_day3_arg1neg_received$aggregate_rank, neg_pairs_all)
neg_ranks_matched <- neg_rank_lookup[pos_pairs]
neg_ranks_matched[is.na(neg_ranks_matched)] <- 1.0
spec_index <- (neg_ranks_matched - lee_day3_arg1pos_top$aggregate_rank) / (neg_ranks_matched + lee_day3_arg1pos_top$aggregate_rank + 1e-10)
lee_day3_arg1pos_top$specificity_index <- spec_index
lee_day3_arg1pos_specific <- dplyr::filter(lee_day3_arg1pos_top, specificity_index > 0.3 | !(pos_pairs %in% neg_pairs_all))

lee_day3_liana_arg1_pos <- lee_day3_arg1pos_received
lee_day3_liana_arg1_neg <- lee_day3_arg1neg_received
lee_day3_liana_consensus <- head(lee_day3_arg1pos_received, 30)
lee_day3_liana_arg1pos_topranked <- head(lee_day3_arg1pos_received, 20)
lee_day3_liana_arg1neg_topranked <- head(lee_day3_arg1neg_received, 20)

## L-R reference with gene names (same as Lee Day 1)
lr_pairs_lee_d3 <- dplyr::distinct(lee_day3_liana_result_df, ligand.complex, receptor.complex, .keep_all = FALSE)
print("LIANA ligand.complex and receptor.complex columns - unique L-R pairs in full result (Lee Day 3):")
print(head(lr_pairs_lee_d3, 30))
all_sym_d3 <- unique(c(lr_pairs_lee_d3$ligand.complex, lr_pairs_lee_d3$receptor.complex))
all_sym_single_d3 <- unique(trimws(unlist(strsplit(all_sym_d3[grepl("[_+]", all_sym_d3)], "[_+]"))))
all_sym_single_d3 <- unique(c(all_sym_d3[!grepl("[_+]", all_sym_d3)], all_sym_single_d3))
tbl_lr_d3 <- data.frame(SYMBOL = character(0), GENENAME = character(0), stringsAsFactors = FALSE)
if (length(all_sym_single_d3) > 0) {
  tbl_lr_d3 <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = all_sym_single_d3, columns = "GENENAME", keytype = "SYMBOL"))
  tbl_lr_d3 <- dplyr::distinct(tbl_lr_d3, SYMBOL, .keep_all = TRUE)
}
sym_to_name_d3 <- if (nrow(tbl_lr_d3) > 0) stats::setNames(tbl_lr_d3$GENENAME, tbl_lr_d3$SYMBOL) else character(0)
lr_ref_d3 <- lr_pairs_lee_d3
lr_ref_d3$ligand_genename <- rep(NA_character_, nrow(lr_ref_d3))
lr_ref_d3$receptor_genename <- rep(NA_character_, nrow(lr_ref_d3))
has_sym_to_name_d3 <- length(sym_to_name_d3) > 0
if (!has_sym_to_name_d3) {
  lr_ref_d3$ligand_genename <- lr_ref_d3$ligand.complex
  lr_ref_d3$receptor_genename <- lr_ref_d3$receptor.complex
}
if (has_sym_to_name_d3) {
  ligand_is_complex_d3 <- grepl("[_+]", lr_ref_d3$ligand.complex)
  lr_ref_d3$ligand_genename[!ligand_is_complex_d3] <- ifelse(lr_ref_d3$ligand.complex[!ligand_is_complex_d3] %in% names(sym_to_name_d3), sym_to_name_d3[lr_ref_d3$ligand.complex[!ligand_is_complex_d3]], lr_ref_d3$ligand.complex[!ligand_is_complex_d3])
  complex_ligands_d3 <- lr_ref_d3$ligand.complex[ligand_is_complex_d3]
  mapped_complex_lig_d3 <- vapply(complex_ligands_d3, function(s) {
    parts <- trimws(strsplit(s, "[_+]")[[1]])
    mapped_parts <- ifelse(parts %in% names(sym_to_name_d3), sym_to_name_d3[parts], parts)
    paste(mapped_parts, collapse = "_")
  }, character(1))
  lr_ref_d3$ligand_genename[ligand_is_complex_d3] <- mapped_complex_lig_d3
  receptor_is_complex_d3 <- grepl("[_+]", lr_ref_d3$receptor.complex)
  lr_ref_d3$receptor_genename[!receptor_is_complex_d3] <- ifelse(lr_ref_d3$receptor.complex[!receptor_is_complex_d3] %in% names(sym_to_name_d3), sym_to_name_d3[lr_ref_d3$receptor.complex[!receptor_is_complex_d3]], lr_ref_d3$receptor.complex[!receptor_is_complex_d3])
  complex_receptors_d3 <- lr_ref_d3$receptor.complex[receptor_is_complex_d3]
  mapped_complex_rec_d3 <- vapply(complex_receptors_d3, function(s) {
    parts <- trimws(strsplit(s, "[_+]")[[1]])
    mapped_parts <- ifelse(parts %in% names(sym_to_name_d3), sym_to_name_d3[parts], parts)
    paste(mapped_parts, collapse = "_")
  }, character(1))
  lr_ref_d3$receptor_genename[receptor_is_complex_d3] <- mapped_complex_rec_d3
}
write.csv(lr_ref_d3, file.path(OUTPUT_DIR, "LeeDay3_LIANA_LigandReceptor_Reference.csv"), row.names = FALSE)
print(paste("Saved L-R reference with gene names to LeeDay3_LIANA_LigandReceptor_Reference.csv"))
print_gene_mapping_d3 <- has_sym_to_name_d3
if (print_gene_mapping_d3) {
  gene_map_df_d3 <- data.frame(symbol = names(sym_to_name_d3), full_name = unname(sym_to_name_d3), stringsAsFactors = FALSE)
  print("LIANA gene symbols in Lee Day 3 -> full names (org.Mm.eg.db):")
  print(gene_map_df_d3)
}

print(paste("Arg1pos received signals (all):", nrow(lee_day3_arg1pos_received)))
print(paste("Arg1pos top signals (top 50):", nrow(lee_day3_arg1pos_top_signals)))
print(paste("Arg1pos-specific signals:", nrow(lee_day3_arg1pos_specific)))

write.csv(lee_day3_liana_result_df, file.path(OUTPUT_DIR, "LeeDay3_LIANA_AllResults.csv"), row.names = FALSE)
write.csv(lee_day3_arg1pos_received, file.path(OUTPUT_DIR, "LeeDay3_Arg1pos_ReceivedSignals.csv"), row.names = FALSE)
write.csv(lee_day3_arg1pos_top_signals, file.path(OUTPUT_DIR, "LeeDay3_LIANA_Arg1pos_Top50Signals.csv"), row.names = FALSE)
write.csv(lee_day3_arg1pos_specific, file.path(OUTPUT_DIR, "LeeDay3_LIANA_Arg1pos_Specific.csv"), row.names = FALSE)
write.csv(lee_day3_liana_consensus, file.path(OUTPUT_DIR, "LeeDay3_LIANA_NeutrophilConsensus.csv"), row.names = FALSE)
write.csv(lee_day3_liana_arg1pos_topranked, file.path(OUTPUT_DIR, "LeeDay3_LIANA_Arg1pos_TopRanked.csv"), row.names = FALSE)
write.csv(lee_day3_liana_arg1neg_topranked, file.path(OUTPUT_DIR, "LeeDay3_LIANA_Arg1neg_TopRanked.csv"), row.names = FALSE)
saveRDS(list(liana_result = lee_day3_liana_result, liana_aggregated = lee_day3_liana_result_df, arg1pos_received = lee_day3_arg1pos_received, arg1pos_top_signals = lee_day3_arg1pos_top_signals, arg1pos_specific = lee_day3_arg1pos_specific, neutrophil_consensus = lee_day3_liana_consensus, arg1pos_topranked = lee_day3_liana_arg1pos_topranked, arg1neg_topranked = lee_day3_liana_arg1neg_topranked), file.path(OUTPUT_DIR, "LeeDay3_LIANA_Results.rds"))
# To skip re-run LIANA only: load LeeDay3_LIANA_Results.rds into list, assign slots below — you must still have lee_day3 in env (from full run or Checkpoint_AfterLeeDay3_LIANA.rds).
# Uncomment block below to load and skip re-running Lee Day 3 LIANA:
# lee_day3_liana_loaded <- readRDS(file.path(OUTPUT_DIR, "LeeDay3_LIANA_Results.rds"))
# lee_day3_liana_result <- lee_day3_liana_loaded$liana_result
# lee_day3_liana_result_df <- lee_day3_liana_loaded$liana_aggregated
# lee_day3_arg1pos_received <- lee_day3_liana_loaded$arg1pos_received
# lee_day3_arg1pos_top_signals <- lee_day3_liana_loaded$arg1pos_top_signals
# lee_day3_arg1pos_specific <- lee_day3_liana_loaded$arg1pos_specific
# lee_day3_liana_consensus <- lee_day3_liana_loaded$neutrophil_consensus
# lee_day3_liana_arg1pos_topranked <- lee_day3_liana_loaded$arg1pos_topranked
# lee_day3_liana_arg1neg_topranked <- lee_day3_liana_loaded$arg1neg_topranked

# -------- Lee Day 3: LIANA visualizations (official: dotplot, heat_freq, chord_freq, liana_heatmap, …) -> CellChat/CellCall overlap -> top LIANA receptors VlnPlot (same flow as Lee Day 1) --------
lee_d3_liana_viz_ok <- nrow(lee_day3_liana_result_df) > 0
liana_trunc_lee_d3 <- data.frame()
p_02f2 <- NULL
p_02g <- NULL
p_heat_freq_ld3 <- NULL
p_chord_freq_ld3 <- NULL
p_liana_heatmap_lee_d3 <- NULL
p_02o2 <- NULL
if (lee_d3_liana_viz_ok) liana_network_lee_d3_both <- lee_day3_liana_result_df |> dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |> dplyr::filter(!(source %in% NEUTROPHIL_STATES)) |> dplyr::group_by(target) |> dplyr::arrange(aggregate_rank) |> dplyr::slice_head(n = 20) |> dplyr::ungroup()
if (lee_d3_liana_viz_ok) liana_network_lee_d3_both$interaction_label <- paste0(liana_network_lee_d3_both$ligand.complex, "\u2013(", liana_network_lee_d3_both$receptor.complex, ")")
if (lee_d3_liana_viz_ok) p_02e <- tryCatch(ggplot(liana_network_lee_d3_both, aes(x = target, y = interaction_label, size = -log10(aggregate_rank + 1e-10), color = source)) + geom_point(alpha = 0.8) + scale_x_discrete(limits = NEUTROPHIL_STATES, drop = FALSE) + theme_minimal() + labs(title = "Top Signals Received: Arg1+ vs Arg1- (Lee Day 3)", subtitle = "Y = Ligand\u2013(Receptor); X = receiver; color = sender (neutrophil\u2192neutrophil excluded from data)", x = "Receiver (neutrophil state)", y = "Interaction", color = "Sender cell type", size = "Consensus support\n(-log10 aggregate rank)") + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), plot.margin = margin(10, 80, 10, 10)) + guides(color = guide_legend(override.aes = list(size = 3))), error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02e")) {
  vals_nlr <- numeric(0)
  if (exists("liana_network_lee_d1_both") && is.data.frame(liana_network_lee_d1_both) && nrow(liana_network_lee_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_lee_d3_both") && is.data.frame(liana_network_lee_d3_both) && nrow(liana_network_lee_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_wang_d3_both") && is.data.frame(liana_network_wang_d3_both) && nrow(liana_network_wang_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_wang_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d1_both") && is.data.frame(liana_network_qin_d1_both) && nrow(liana_network_qin_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d3_both") && is.data.frame(liana_network_qin_d3_both) && nrow(liana_network_qin_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d3_both$aggregate_rank + 1e-10))
  vals_nlr <- vals_nlr[is.finite(vals_nlr)]
  UNIFY_LIANA_NEGLOG10 <- if (length(vals_nlr) > 0) range(vals_nlr) else c(0, 10)
  if (length(UNIFY_LIANA_NEGLOG10) != 2 || !all(is.finite(UNIFY_LIANA_NEGLOG10))) UNIFY_LIANA_NEGLOG10 <- c(0, 10)
  if (UNIFY_LIANA_NEGLOG10[1] == UNIFY_LIANA_NEGLOG10[2]) UNIFY_LIANA_NEGLOG10[2] <- UNIFY_LIANA_NEGLOG10[1] + 1e-6
  print(p_02e + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = (if (exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE else c(2, 8))))
}
if (lee_d3_liana_viz_ok) lee_day3_external_to_neutrophils <- lee_day3_liana_result_df |> dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |> dplyr::filter(!(source %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))) |> dplyr::arrange(aggregate_rank)
if (lee_d3_liana_viz_ok) liana_top_lee_d3_external <- dplyr::slice_head(lee_day3_external_to_neutrophils, n = 20)
if (lee_d3_liana_viz_ok) liana_top_lee_d3_external$lr_label <- paste0(liana_top_lee_d3_external$source, " -> ", liana_top_lee_d3_external$target, "  ", liana_top_lee_d3_external$ligand.complex, "\u2013(", liana_top_lee_d3_external$receptor.complex, ")")
if (lee_d3_liana_viz_ok) p_02f0 <- tryCatch(ggplot(liana_top_lee_d3_external, aes(x = reorder(lr_label, aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 20 L-R Pairs - LIANA (Lee Day 3)\nExternal signals to Arg1+ / Arg1- neutrophils", x = "Source -> Target  Ligand\u2013(Receptor)", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02f0")) {
  vals_ext <- numeric(0)
  if (exists("lee_day1_external_to_neutrophils") && is.data.frame(lee_day1_external_to_neutrophils) && nrow(lee_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day1_external_to_neutrophils$aggregate_rank)
  if (exists("lee_day3_external_to_neutrophils") && is.data.frame(lee_day3_external_to_neutrophils) && nrow(lee_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day3_external_to_neutrophils$aggregate_rank)
  if (exists("wang_day3_external_to_neutrophils") && is.data.frame(wang_day3_external_to_neutrophils) && nrow(wang_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(wang_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, wang_day3_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day1_external_to_neutrophils") && is.data.frame(qin_day1_external_to_neutrophils) && nrow(qin_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day1_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day3_external_to_neutrophils") && is.data.frame(qin_day3_external_to_neutrophils) && nrow(qin_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day3_external_to_neutrophils$aggregate_rank)
  if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0 && "aggregate_rank" %in% names(liana_top_lee_d1)) vals_ext <- c(vals_ext, liana_top_lee_d1$aggregate_rank)
  if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_lee_d3_external)) vals_ext <- c(vals_ext, liana_top_lee_d3_external$aggregate_rank)
  if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_wang_d3_external)) vals_ext <- c(vals_ext, liana_top_wang_d3_external$aggregate_rank)
  if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d1_external)) vals_ext <- c(vals_ext, liana_top_qin_d1_external$aggregate_rank)
  if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d3_external)) vals_ext <- c(vals_ext, liana_top_qin_d3_external$aggregate_rank)
  vals_ext <- suppressWarnings(as.numeric(vals_ext))
  vals_ext <- vals_ext[is.finite(vals_ext)]
  plot_ext_ar <- suppressWarnings(as.numeric(liana_top_lee_d3_external[["aggregate_rank"]]))
  plot_ext_ar <- plot_ext_ar[is.finite(plot_ext_ar)]
  UNIFY_EXT_AR <- range(c(vals_ext, plot_ext_ar), na.rm = TRUE)
  if (length(plot_ext_ar) == 0 && length(vals_ext) == 0) UNIFY_EXT_AR <- c(0, 1)
  if (!all(is.finite(UNIFY_EXT_AR))) UNIFY_EXT_AR <- c(0, 1)
  if (UNIFY_EXT_AR[1] == UNIFY_EXT_AR[2]) UNIFY_EXT_AR[2] <- UNIFY_EXT_AR[1] + max(abs(UNIFY_EXT_AR[1]) * 1e-6, 1e-12)
  print(p_02f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
}
if (lee_d3_liana_viz_ok) liana_top_lee_d3 <- dplyr::slice_head(lee_day3_liana_consensus, n = 15)
if (lee_d3_liana_viz_ok) p_02f <- ggplot(liana_top_lee_d3, aes(x = reorder(paste0(source, " -> ", target), aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 15 L-R Pairs - LIANA (Lee Day 3)", x = "Source -> Target", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02f")) {
  vals_cons <- numeric(0)
  if (exists("liana_top_lee_d3") && is.data.frame(liana_top_lee_d3) && nrow(liana_top_lee_d3) > 0) vals_cons <- c(vals_cons, liana_top_lee_d3$aggregate_rank)
  if (exists("liana_top_wang_d3") && is.data.frame(liana_top_wang_d3) && nrow(liana_top_wang_d3) > 0) vals_cons <- c(vals_cons, liana_top_wang_d3$aggregate_rank)
  if (exists("liana_top_qin_d1") && is.data.frame(liana_top_qin_d1) && nrow(liana_top_qin_d1) > 0) vals_cons <- c(vals_cons, liana_top_qin_d1$aggregate_rank)
  if (exists("liana_top_qin_d3") && is.data.frame(liana_top_qin_d3) && nrow(liana_top_qin_d3) > 0) vals_cons <- c(vals_cons, liana_top_qin_d3$aggregate_rank)
  vals_cons <- vals_cons[is.finite(vals_cons)]
  UNIFY_CONS_AR <- if (length(vals_cons) > 0) range(vals_cons) else c(0, 1)
  if (UNIFY_CONS_AR[1] == UNIFY_CONS_AR[2]) UNIFY_CONS_AR[2] <- UNIFY_CONS_AR[1] + 1e-6
  print(p_02f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
}
if (lee_d3_liana_viz_ok) unique_sources_lee_d3 <- unique(lee_day3_liana_result_df$source)
if (lee_d3_liana_viz_ok) neutrophil_targets_ld3 <- NEUTROPHIL_STATES[NEUTROPHIL_STATES %in% unique(lee_day3_liana_result_df$target)]
if (lee_d3_liana_viz_ok) dot_sources_lee_d3 <- setdiff(unique_sources_lee_d3, neutrophil_targets_ld3)
if (lee_d3_liana_viz_ok && length(dot_sources_lee_d3) == 0) dot_sources_lee_d3 <- unique_sources_lee_d3
if (lee_d3_liana_viz_ok && length(neutrophil_targets_ld3) > 0 && length(dot_sources_lee_d3) > 0) p_02f2 <- tryCatch(liana::liana_dotplot(lee_day3_liana_result_df, source_groups = dot_sources_lee_d3, target_groups = neutrophil_targets_ld3, ntop = 20, size_range = LIANA_DOTPLOT_SIZE_RANGE), error = function(e) NULL)
if (lee_d3_liana_viz_ok && !is.null(p_02f2)) p_02f2 <- p_02f2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
if (lee_d3_liana_viz_ok && !is.null(p_02f2)) print(p_02f2)
if (lee_d3_liana_viz_ok) interaction_matrix_lee_d3 <- tryCatch(tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(lee_day3_liana_consensus, source, target), interaction_count = dplyr::n(), .groups = "drop"), names_from = target, values_from = interaction_count, values_fill = 0), error = function(e) data.frame())
if (lee_d3_liana_viz_ok) interaction_matrix_lee_d3_mat <- tryCatch({ cols_num <- setdiff(colnames(interaction_matrix_lee_d3), "source"); mat <- as.matrix(interaction_matrix_lee_d3[, cols_num, drop = FALSE]); rownames(mat) <- interaction_matrix_lee_d3$source; mat }, error = function(e) matrix(0, nrow = 0, ncol = 0))
if (lee_d3_liana_viz_ok) p_02g <- tryCatch(pheatmap::pheatmap(interaction_matrix_lee_d3_mat, color = colorRampPalette(c("white", "yellow", "orange", "red"))(100), main = "Cell Type Interaction Frequency - Lee Day 3", display_numbers = TRUE), error = function(e) NULL)
if (lee_d3_liana_viz_ok && !is.null(p_02g)) print(p_02g)
if (lee_d3_liana_viz_ok && !is.null(p_02g)) grid::grid.newpage()
if (lee_d3_liana_viz_ok && !is.null(p_02g)) grid::grid.draw(p_02g$gtable)
if (lee_d3_liana_viz_ok) receptor_freq_lee_d3 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(lee_day3_liana_consensus, receptor.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
if (lee_d3_liana_viz_ok) p_02h <- tryCatch(ggplot(receptor_freq_lee_d3, aes(x = reorder(receptor.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Receptors - Lee Day 3", x = "Receptor", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02h")) {
  vals_fc <- numeric(0)
  vals_fmr <- numeric(0)
  freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
  for (fn in freq_tabnames) {
    if (!exists(fn)) next
    d <- get(fn)
    if (!is.data.frame(d) || nrow(d) == 0) next
    if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
    if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
  }
  vals_fc <- vals_fc[is.finite(vals_fc)]
  vals_fmr <- vals_fmr[is.finite(vals_fmr)]
  UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
  UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
  if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
  if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
  print(p_02h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
}
if (lee_d3_liana_viz_ok) ligand_freq_lee_d3 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(lee_day3_liana_consensus, ligand.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
if (lee_d3_liana_viz_ok) p_02i <- tryCatch(ggplot(ligand_freq_lee_d3, aes(x = reorder(ligand.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Ligands - Lee Day 3", x = "Ligand", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02i")) {
  vals_fc <- numeric(0)
  vals_fmr <- numeric(0)
  freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
  for (fn in freq_tabnames) {
    if (!exists(fn)) next
    d <- get(fn)
    if (!is.data.frame(d) || nrow(d) == 0) next
    if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
    if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
  }
  vals_fc <- vals_fc[is.finite(vals_fc)]
  vals_fmr <- vals_fmr[is.finite(vals_fmr)]
  UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
  UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
  if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
  if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
  print(p_02i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
}
if (lee_d3_liana_viz_ok) liana_trunc_lee_d3 <- dplyr::filter(lee_day3_liana_result_df, aggregate_rank <= 0.01)
lee_d3_trunc_ok <- lee_d3_liana_viz_ok && nrow(liana_trunc_lee_d3) > 0
if (lee_d3_trunc_ok) liana_trunc_lee_d3$source <- as.character(liana_trunc_lee_d3$source)
if (lee_d3_trunc_ok) liana_trunc_lee_d3$target <- as.character(liana_trunc_lee_d3$target)
if (lee_d3_trunc_ok) p_heat_freq_ld3 <- tryCatch(liana::heat_freq(liana_trunc_lee_d3), error = function(e) { message("heat_freq Lee D3: ", conditionMessage(e)); NULL })
if (lee_d3_trunc_ok && !is.null(p_heat_freq_ld3)) print(p_heat_freq_ld3)
if (lee_d3_trunc_ok) unique_sources_chord_ld3 <- unique(liana_trunc_lee_d3$source)
if (lee_d3_trunc_ok) unique_targets_chord_ld3 <- unique(liana_trunc_lee_d3$target)
if (lee_d3_trunc_ok) grDevices::png(file.path(OUTPUT_DIR, "LeeDay3_LIANA_ChordFreq.png"), width = 1400, height = 1400, res = 150)
if (lee_d3_trunc_ok) tryCatch(liana::chord_freq(liana_trunc_lee_d3, source_groups = unique_sources_chord_ld3, target_groups = unique_targets_chord_ld3), error = function(e) message("chord_freq Lee D3 (PNG): ", conditionMessage(e)))
if (lee_d3_trunc_ok) grDevices::dev.off()
if (lee_d3_trunc_ok) p_chord_freq_ld3 <- tryCatch(liana::chord_freq(liana_trunc_lee_d3, source_groups = unique_sources_chord_ld3, target_groups = unique_targets_chord_ld3), error = function(e) { message("chord_freq Lee D3: ", conditionMessage(e)); NULL })
if (lee_d3_trunc_ok && !is.null(p_chord_freq_ld3)) print(p_chord_freq_ld3)
if (lee_d3_trunc_ok) liana_mat_lee_d3 <- as.matrix(table(liana_trunc_lee_d3$source, liana_trunc_lee_d3$target))
if (lee_d3_trunc_ok) p_liana_heatmap_lee_d3 <- NULL
if (lee_d3_trunc_ok && nrow(liana_mat_lee_d3) > 0 && ncol(liana_mat_lee_d3) > 0) p_liana_heatmap_lee_d3 <- tryCatch(liana::liana_heatmap(liana_mat_lee_d3), error = function(e) { message("liana_heatmap Lee D3: ", conditionMessage(e)); NULL })
if (lee_d3_trunc_ok && !is.null(p_liana_heatmap_lee_d3)) ComplexHeatmap::draw(p_liana_heatmap_lee_d3)
if (lee_d3_liana_viz_ok) source_importance_lee_d3 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(lee_day3_liana_consensus, source), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), importance_score = dplyr::n() * (1 - mean(aggregate_rank))), desc(importance_score)), error = function(e) data.frame())
if (lee_d3_liana_viz_ok) p_02k <- tryCatch(ggplot(source_importance_lee_d3, aes(x = reorder(source, importance_score), y = importance_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkred") + labs(title = "Source Cell Importance - Lee Day 3", x = "Cell Type", y = "Importance Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02k")) {
  vals_iy <- numeric(0)
  vals_imr <- numeric(0)
  imp_tabnames <- c("source_importance_lee_d1", "source_importance_lee_d3", "source_importance_wang_d3", "source_importance_qin_d1", "source_importance_qin_d3")
  for (fn in imp_tabnames) {
    if (!exists(fn)) next
    d <- get(fn)
    if (!is.data.frame(d) || nrow(d) == 0) next
    if ("importance_score" %in% names(d)) vals_iy <- c(vals_iy, d$importance_score)
    if ("mean_rank" %in% names(d)) vals_imr <- c(vals_imr, d$mean_rank)
  }
  vals_iy <- vals_iy[is.finite(vals_iy)]
  vals_imr <- vals_imr[is.finite(vals_imr)]
  UNIFY_SRC_IMP_Y <- if (length(vals_iy) > 0) range(vals_iy) else c(0, 1)
  UNIFY_SRC_IMP_MR <- if (length(vals_imr) > 0) range(vals_imr) else c(0, 1)
  if (UNIFY_SRC_IMP_Y[1] == UNIFY_SRC_IMP_Y[2]) UNIFY_SRC_IMP_Y[2] <- UNIFY_SRC_IMP_Y[1] + 1e-6
  if (UNIFY_SRC_IMP_MR[1] == UNIFY_SRC_IMP_MR[2]) UNIFY_SRC_IMP_MR[2] <- UNIFY_SRC_IMP_MR[1] + 1e-6
  print(p_02k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
}
if (lee_d3_liana_viz_ok) p_02n <- tryCatch(ggplot(head(lee_day3_liana_arg1pos_topranked, 15), aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), aggregate_rank), y = -aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top Arg1pos Interactions - Lee Day 3", x = "L-R Pair", y = "Rank Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02n")) {
  vals_a1 <- numeric(0)
  if (exists("lee_day1_liana_arg1pos_topranked") && is.data.frame(lee_day1_liana_arg1pos_topranked) && nrow(lee_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day1_liana_arg1pos_topranked$aggregate_rank, 15))
  if (exists("lee_day3_liana_arg1pos_topranked") && is.data.frame(lee_day3_liana_arg1pos_topranked) && nrow(lee_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day3_liana_arg1pos_topranked$aggregate_rank, 15))
  if (exists("wang_day3_liana_arg1pos_topranked") && is.data.frame(wang_day3_liana_arg1pos_topranked) && nrow(wang_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(wang_day3_liana_arg1pos_topranked$aggregate_rank, 15))
  if (exists("qin_day1_liana_arg1pos_topranked") && is.data.frame(qin_day1_liana_arg1pos_topranked) && nrow(qin_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day1_liana_arg1pos_topranked$aggregate_rank, 15))
  if (exists("qin_day3_liana_arg1pos_topranked") && is.data.frame(qin_day3_liana_arg1pos_topranked) && nrow(qin_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day3_liana_arg1pos_topranked$aggregate_rank, 15))
  vals_a1 <- vals_a1[is.finite(vals_a1)]
  UNIFY_ARG1_AR <- if (length(vals_a1) > 0) range(vals_a1) else c(0, 1)
  if (UNIFY_ARG1_AR[1] == UNIFY_ARG1_AR[2]) UNIFY_ARG1_AR[2] <- UNIFY_ARG1_AR[1] + 1e-6
  UNIFY_ARG1_NEGY <- c(-UNIFY_ARG1_AR[2], -UNIFY_ARG1_AR[1])
  print(p_02n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
}
if (lee_d3_liana_viz_ok) lee_day3_arg1pos_specific_plot <- head(lee_day3_arg1pos_specific, 15)
if (lee_d3_liana_viz_ok && nrow(lee_day3_arg1pos_specific_plot) > 0) p_02o2 <- tryCatch(ggplot(lee_day3_arg1pos_specific_plot, aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), specificity_index), y = specificity_index, fill = source)) + geom_bar(stat = "identity") + coord_flip() + labs(title = "Arg1+-Specific Signals (Lee Day 3)", x = "L-R Pair", y = "Specificity Index") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (lee_d3_liana_viz_ok && !is.null(p_02o2)) print(p_02o2)
if (lee_d3_liana_viz_ok) target_specificity_lee_d3 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(dplyr::filter(lee_day3_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")), target), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), specificity_score = dplyr::n() * (1 - mean(aggregate_rank)), .groups = "drop"), desc(specificity_score)), error = function(e) data.frame())
if (lee_d3_liana_viz_ok) p_02l <- tryCatch(ggplot(target_specificity_lee_d3, aes(x = reorder(target, specificity_score), y = specificity_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkgreen") + labs(title = "Target Cell Specificity - Lee Day 3", x = "Cell Type", y = "Specificity Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02l")) {
  vals_sy <- numeric(0)
  vals_smr <- numeric(0)
  spec_tabnames <- c("target_specificity_lee_d1", "target_specificity_lee_d3", "target_specificity_wang_d3", "target_specificity_qin_d1", "target_specificity_qin_d3")
  for (fn in spec_tabnames) {
    if (!exists(fn)) next
    d <- get(fn)
    if (!is.data.frame(d) || nrow(d) == 0) next
    if ("specificity_score" %in% names(d)) vals_sy <- c(vals_sy, d$specificity_score)
    if ("mean_rank" %in% names(d)) vals_smr <- c(vals_smr, d$mean_rank)
  }
  vals_sy <- vals_sy[is.finite(vals_sy)]
  vals_smr <- vals_smr[is.finite(vals_smr)]
  UNIFY_TGT_SPEC_Y <- if (length(vals_sy) > 0) range(vals_sy) else c(0, 1)
  UNIFY_TGT_SPEC_MR <- if (length(vals_smr) > 0) range(vals_smr) else c(0, 1)
  if (UNIFY_TGT_SPEC_Y[1] == UNIFY_TGT_SPEC_Y[2]) UNIFY_TGT_SPEC_Y[2] <- UNIFY_TGT_SPEC_Y[1] + 1e-6
  if (UNIFY_TGT_SPEC_MR[1] == UNIFY_TGT_SPEC_MR[2]) UNIFY_TGT_SPEC_MR[2] <- UNIFY_TGT_SPEC_MR[1] + 1e-6
  print(p_02l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
}
if (lee_d3_liana_viz_ok) p_02m <- tryCatch(ggplot(lee_day3_liana_consensus, aes(x = aggregate_rank)) + geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) + labs(title = "Distribution of L-R Pair Aggregate Ranks - Lee Day 3", x = "Aggregate Rank Score", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02m")) {
  vals_hist <- numeric(0)
  if (exists("lee_day1_liana_consensus") && is.data.frame(lee_day1_liana_consensus) && nrow(lee_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day1_liana_consensus)) vals_hist <- c(vals_hist, lee_day1_liana_consensus$aggregate_rank)
  if (exists("lee_day3_liana_consensus") && is.data.frame(lee_day3_liana_consensus) && nrow(lee_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day3_liana_consensus)) vals_hist <- c(vals_hist, lee_day3_liana_consensus$aggregate_rank)
  if (exists("wang_day3_liana_consensus") && is.data.frame(wang_day3_liana_consensus) && nrow(wang_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(wang_day3_liana_consensus)) vals_hist <- c(vals_hist, wang_day3_liana_consensus$aggregate_rank)
  if (exists("qin_day1_liana_consensus") && is.data.frame(qin_day1_liana_consensus) && nrow(qin_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day1_liana_consensus)) vals_hist <- c(vals_hist, qin_day1_liana_consensus$aggregate_rank)
  if (exists("qin_day3_liana_consensus") && is.data.frame(qin_day3_liana_consensus) && nrow(qin_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day3_liana_consensus)) vals_hist <- c(vals_hist, qin_day3_liana_consensus$aggregate_rank)
  vals_hist <- vals_hist[is.finite(vals_hist)]
  UNIFY_HIST_AR <- if (length(vals_hist) > 0) range(vals_hist) else c(0, 1)
  if (UNIFY_HIST_AR[1] == UNIFY_HIST_AR[2]) UNIFY_HIST_AR[2] <- UNIFY_HIST_AR[1] + 1e-6
  print(p_02m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
}
if (lee_d3_liana_viz_ok) method_cols_lee_d3 <- colnames(lee_day3_liana_result_df)[grep("pval_", colnames(lee_day3_liana_result_df))]
if (lee_d3_liana_viz_ok) consensus_lee_d3_step1 <- dplyr::slice_head(lee_day3_liana_consensus, n = 15)
if (lee_d3_liana_viz_ok) consensus_lee_d3 <- NULL
if (lee_d3_liana_viz_ok) consensus_lee_d3_has_cols <- length(method_cols_lee_d3) > 0 && nrow(consensus_lee_d3_step1) > 0
if (lee_d3_liana_viz_ok && consensus_lee_d3_has_cols) consensus_lee_d3 <- as.data.frame(dplyr::select(consensus_lee_d3_step1, dplyr::all_of(method_cols_lee_d3)))
if (lee_d3_liana_viz_ok && consensus_lee_d3_has_cols) consensus_lee_d3_rownames <- paste0(consensus_lee_d3_step1$source, " | ", consensus_lee_d3_step1$ligand.complex, " -> ", consensus_lee_d3_step1$receptor.complex)
if (lee_d3_liana_viz_ok && consensus_lee_d3_has_cols) rownames(consensus_lee_d3) <- make.unique(as.character(consensus_lee_d3_rownames))
if (lee_d3_liana_viz_ok) consensus_lee_d3_ready <- !is.null(consensus_lee_d3) && nrow(consensus_lee_d3) > 0
if (!lee_d3_liana_viz_ok) consensus_lee_d3_ready <- FALSE
if (consensus_lee_d3_ready) p_02o <- pheatmap::pheatmap(consensus_lee_d3, color = colorRampPalette(c("red", "white", "blue"))(100), main = "Method Consensus (p-values) - Lee Day 3")
if (consensus_lee_d3_ready) print(p_02o)
# -------- Lee Day 3: CellChat/CellCall downstream (integration = overlap with LIANA, not passed into LIANA) --------
# Subset CellChat communications to neutrophil source/target; parse CellCall LR pairs; write CSVs; overlap with LIANA downstream
neutrophil_labels_d3 <- c("NeutrophilArg1pos", "NeutrophilArg1neg")
path_cc_d3 <- CELLCHAT_LEE_DAY3_RDS
path_ccall_d3 <- CELLCALL_LEE_DAY3_RDS
cc_d3_down <- NULL
if (!is.na(path_cc_d3) && nzchar(path_cc_d3) && file.exists(path_cc_d3)) {
  cc_d3_down <- tryCatch(readRDS(path_cc_d3), error = function(e) { print(paste("Lee Day 3 CellChat (downstream) read failed:", conditionMessage(e))); NULL })
}
if (!is.null(cc_d3_down)) {
  comm_d3 <- CellChat::subsetCommunication(cc_d3_down)
  lee_day3_cellchat_neut <- dplyr::filter(comm_d3, source %in% neutrophil_labels_d3 | target %in% neutrophil_labels_d3)
  write.csv(lee_day3_cellchat_neut, file.path(OUTPUT_DIR, "LeeDay3_CellChat_NeutrophilCommunications.csv"), row.names = FALSE)
  print(paste("Lee Day 3: CellChat neutrophil communications:", nrow(lee_day3_cellchat_neut)))
}
if (!is.na(path_ccall_d3) && nzchar(path_ccall_d3) && file.exists(path_ccall_d3)) {
  ccall_d3 <- readRDS(path_ccall_d3)
  lr_rownames_d3 <- rownames(ccall_d3@data$expr_l_r_log2_scale)
  ccall_lig_d3 <- sub("-.*", "", lr_rownames_d3)
  ccall_rec_d3 <- sub("^[^-]+-", "", lr_rownames_d3)
  ccall_rec_d3 <- gsub("-", "_", ccall_rec_d3)
  lee_day3_cellcall_lr <- data.frame(ligand.complex = ccall_lig_d3, receptor.complex = ccall_rec_d3, stringsAsFactors = FALSE)
  lee_day3_cellcall_lr <- lee_day3_cellcall_lr[nzchar(lee_day3_cellcall_lr$receptor.complex), ]
  lee_day3_cellcall_lr <- dplyr::distinct(lee_day3_cellcall_lr, ligand.complex, receptor.complex)
  write.csv(lee_day3_cellcall_lr, file.path(OUTPUT_DIR, "LeeDay3_CellCall_LRPairs.csv"), row.names = FALSE)
  print(paste("Lee Day 3: CellCall LR pairs:", nrow(lee_day3_cellcall_lr)))
}
# Overlap with LIANA (LIANA uses "_" for complexes; CellChat/CellCall formatted to match)
lee_day3_liana_lr <- dplyr::distinct(lee_day3_liana_result_df, ligand.complex, receptor.complex)
if (nrow(lee_day3_liana_lr) > 0 && exists("lee_day3_cellcall_lr")) {
  lee_day3_lr_overlap_ccall <- dplyr::semi_join(lee_day3_liana_lr, lee_day3_cellcall_lr, by = c("ligand.complex", "receptor.complex"))
  write.csv(lee_day3_lr_overlap_ccall, file.path(OUTPUT_DIR, "LeeDay3_LIANA_CellCall_Overlap.csv"), row.names = FALSE)
  print(paste("Lee Day 3: LIANA-CellCall overlap:", nrow(lee_day3_lr_overlap_ccall), "LR pairs"))
} else if (nrow(lee_day3_liana_lr) == 0) {
  print("Lee Day 3: LIANA-CellCall overlap skipped (no LIANA results for Day 3)")
}
if (nrow(lee_day3_liana_lr) > 0 && exists("lee_day3_cellchat_neut") && nrow(lee_day3_cellchat_neut) > 0 && "ligand" %in% names(lee_day3_cellchat_neut) && "receptor" %in% names(lee_day3_cellchat_neut)) {
  lee_day3_cellchat_lr <- dplyr::distinct(lee_day3_cellchat_neut, ligand, receptor)
  lee_day3_cellchat_lr$ligand.complex <- lee_day3_cellchat_lr$ligand
  lee_day3_cellchat_lr$receptor.complex <- lee_day3_cellchat_lr$receptor
  lee_day3_lr_overlap_cc <- dplyr::semi_join(lee_day3_liana_lr, lee_day3_cellchat_lr, by = c("ligand.complex", "receptor.complex"))
  write.csv(lee_day3_lr_overlap_cc, file.path(OUTPUT_DIR, "LeeDay3_LIANA_CellChat_Overlap.csv"), row.names = FALSE)
  print(paste("Lee Day 3: LIANA-CellChat overlap:", nrow(lee_day3_lr_overlap_cc), "LR pairs"))
}

top_receptors_ld3 <- head(unique(lee_day3_arg1pos_received$receptor.complex), LIANA_TOP_RECEPTOR_VLN)
top_receptors_ld3_single <- top_receptors_ld3[!grepl("[_+]", top_receptors_ld3)]
top_receptors_ld3_in_data <- top_receptors_ld3_single[top_receptors_ld3_single %in% rownames(lee_day3)]
if (length(top_receptors_ld3_in_data) > 0) {
  lee_day3_neut_obj_vln <- subset(lee_day3, cells = lee_day3_neut_cells)
  p_receptor_vln_ld3 <- Seurat::VlnPlot(lee_day3_neut_obj_vln, features = top_receptors_ld3_in_data, group.by = "arg1_status", pt.size = 0.1, ncol = min(3L, length(top_receptors_ld3_in_data)))
  print(p_receptor_vln_ld3)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "LeeDay3_TopLIANA_Receptors_VlnPlot.png"), p_receptor_vln_ld3, width = 10, height = 6, dpi = 150)
}
print("✓ Lee Day 3 LIANA analysis complete")
## 4) Downstream: PROGENy -> BARTsc -> integration (same order as Lee Day 1)
# Checkpoint: resume from here if later section crashes. load(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay3_LIANA.RData"))
save.image(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay3_LIANA.RData"))
saveRDS(list(lee_day3 = lee_day3, lee_day3_pruned_labels = lee_day3_pruned_labels, lee_day3_pos = lee_day3_pos, lee_day3_neg = lee_day3_neg, lee_day3_neut_cells = lee_day3_neut_cells, lee_day3_arg1_status = lee_day3_arg1_status, lee_day3_liana_result_df = lee_day3_liana_result_df, lee_day3_arg1pos_received = lee_day3_arg1pos_received, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_LIANA.rds"))
# checkpoint_ld3 <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_LIANA.rds")); list2env(checkpoint_ld3, envir = .GlobalEnv)
print("✓ Checkpoint saved: Workspace_AfterLeeDay3_LIANA.RData, Checkpoint_AfterLeeDay3_LIANA.rds")

# -------- Lee Day 3 PROGENy Pathway Analysis (Exploratory; downstream pathway support only) --------
# Runs after LIANA: Ligand (LIANA) -> Pathway (PROGENy) -> TF (BARTsc). Same setup as Lee Day 1.
# Target population: NeutrophilArg1pos and NeutrophilArg1neg only; pathway activity compared Arg1pos vs Arg1neg.
print("--- PROGENy Pathway Analysis: Lee Day 3 (Exploratory) ---")
# Secondary support only: PROGENy summarizes downstream pathway state in Arg1pos vs Arg1neg neutrophils.

progeny_model_mouse_lee_d3 <- progeny::model_mouse_full
colnames(progeny_model_mouse_lee_d3) <- tolower(colnames(progeny_model_mouse_lee_d3))
colnames(progeny_model_mouse_lee_d3)[colnames(progeny_model_mouse_lee_d3) == "p.value"] <- "p_value"
stopifnot(all(c("gene", "pathway", "weight") %in% colnames(progeny_model_mouse_lee_d3)))
progeny_model_mouse_lee_d3_dim <- dim(progeny_model_mouse_lee_d3)
progeny_model_mouse_lee_d3_colnames <- colnames(progeny_model_mouse_lee_d3)
print(paste("PROGENy model dimensions (Lee Day 3):", paste(progeny_model_mouse_lee_d3_dim, collapse = " x ")))
print(paste("PROGENy model columns:", paste(progeny_model_mouse_lee_d3_colnames, collapse = ", ")))
progeny_network_lee_d3 <- data.frame(
  source = progeny_model_mouse_lee_d3$pathway,
  target = progeny_model_mouse_lee_d3$gene,
  weight = progeny_model_mouse_lee_d3$weight,
  stringsAsFactors = FALSE
)
progeny_network_lee_d3 <- progeny_network_lee_d3[progeny_network_lee_d3$weight != 0, ]
progeny_network_lee_d3_nrow <- nrow(progeny_network_lee_d3)
progeny_network_lee_d3_pathways_unique <- unique(progeny_network_lee_d3$source)
print(paste("PROGENy network (Lee Day 3): nrow =", progeny_network_lee_d3_nrow, ", unique pathways =", length(progeny_network_lee_d3_pathways_unique)))
print(paste("PROGENy pathways:", paste(progeny_network_lee_d3_pathways_unique, collapse = ", ")))
stopifnot(progeny_network_lee_d3_nrow > 0)
progeny_network_lee_d3_has_pvalue <- "p_value" %in% colnames(progeny_network_lee_d3)
progeny_network_lee_d3_cols <- colnames(progeny_network_lee_d3)
progeny_network_lee_d3_cols_no_pvalue <- progeny_network_lee_d3_cols[progeny_network_lee_d3_cols != "p_value"]
progeny_network_lee_d3 <- progeny_network_lee_d3[, progeny_network_lee_d3_cols_no_pvalue, drop = FALSE]
stopifnot("source" %in% colnames(progeny_network_lee_d3))
stopifnot("target" %in% colnames(progeny_network_lee_d3))

# Subset to neutrophils only (lee_day3_neut_cells defined above; same as Lee Day 1 pattern)
lee_day3_neut_expr_mat <- tryCatch(Seurat::GetAssayData(lee_day3, layer = "data")[, lee_day3_neut_cells, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
lee_day3_expr_mat <- lee_day3_neut_expr_mat

# Calculate pathway activity scores using weighted mean (WMEAN)
lee_day3_progeny_result <- tryCatch(decoupleR::run_wmean(mat = lee_day3_expr_mat, network = progeny_network_lee_d3, .source = "source", .target = "target", .mor = "weight", minsize = 5), error = function(e) data.frame())

# Convert result to wide format (pathways x cells). Use norm_wmean statistic only.
lee_day3_progeny_for_wide <- if (nrow(lee_day3_progeny_result) > 0 && "statistic" %in% colnames(lee_day3_progeny_result)) dplyr::filter(lee_day3_progeny_result, statistic == "norm_wmean") else lee_day3_progeny_result
if (nrow(lee_day3_progeny_for_wide) == 0 && nrow(lee_day3_progeny_result) > 0 && "statistic" %in% colnames(lee_day3_progeny_result)) lee_day3_progeny_for_wide <- dplyr::filter(lee_day3_progeny_result, statistic == "wmean")
lee_day3_pw_wide <- tryCatch(tidyr::pivot_wider(lee_day3_progeny_for_wide, names_from = "condition", values_from = "score", id_cols = "source"), error = function(e) data.frame())
lee_day3_pw_cols_num <- tryCatch(sapply(lee_day3_pw_wide[, -1, drop = FALSE], function(x) as.numeric(unlist(x))), error = function(e) matrix(0, nrow = 0, ncol = 0))
lee_day3_progeny_scores_mat <- tryCatch(as.matrix(lee_day3_pw_cols_num), error = function(e) matrix(0, nrow = 0, ncol = 0))
rownames(lee_day3_progeny_scores_mat) <- tryCatch(as.character(lee_day3_pw_wide$source), error = function(e) character(0))
colnames(lee_day3_progeny_scores_mat) <- tryCatch(colnames(lee_day3_pw_wide)[-1], error = function(e) character(0))

# Get cell names for Arg1pos and Arg1neg (use neutrophil indices directly, same as Day 1)
lee_day3_arg1pos_cells <- colnames(lee_day3)[lee_day3_pos]
lee_day3_arg1neg_cells <- colnames(lee_day3)[lee_day3_neg]

# Extract pathway scores for each group
lee_day3_arg1pos_cells_in_mat <- tryCatch(lee_day3_arg1pos_cells[lee_day3_arg1pos_cells %in% colnames(lee_day3_progeny_scores_mat)], error = function(e) character(0))
lee_day3_arg1neg_cells_in_mat <- tryCatch(lee_day3_arg1neg_cells[lee_day3_arg1neg_cells %in% colnames(lee_day3_progeny_scores_mat)], error = function(e) character(0))

lee_day3_progeny_arg1pos_scores <- tryCatch(lee_day3_progeny_scores_mat[, lee_day3_arg1pos_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
lee_day3_progeny_arg1neg_scores <- tryCatch(lee_day3_progeny_scores_mat[, lee_day3_arg1neg_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))

# Calculate mean and median pathway scores per group (effect sizes)
lee_day3_pathway_means_arg1pos <- tryCatch(rowMeans(lee_day3_progeny_arg1pos_scores, na.rm = TRUE), error = function(e) numeric(0))
lee_day3_pathway_means_arg1neg <- tryCatch(rowMeans(lee_day3_progeny_arg1neg_scores, na.rm = TRUE), error = function(e) numeric(0))
lee_day3_pathway_medians_arg1pos <- tryCatch(apply(lee_day3_progeny_arg1pos_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
lee_day3_pathway_medians_arg1neg <- tryCatch(apply(lee_day3_progeny_arg1neg_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))

lee_day3_pathway_shift <- pmax(0, -pmin(lee_day3_pathway_means_arg1pos, lee_day3_pathway_means_arg1neg, na.rm = TRUE)) + 1e-10
lee_day3_pathway_log2fc <- log2((lee_day3_pathway_means_arg1pos + lee_day3_pathway_shift) / (lee_day3_pathway_means_arg1neg + lee_day3_pathway_shift))

lee_day3_pathway_comparison <- tryCatch(data.frame(
  pathway = names(lee_day3_pathway_means_arg1pos),
  Arg1pos_mean = lee_day3_pathway_means_arg1pos,
  Arg1neg_mean = lee_day3_pathway_means_arg1neg,
  Arg1pos_median = lee_day3_pathway_medians_arg1pos,
  Arg1neg_median = lee_day3_pathway_medians_arg1neg,
  mean_difference = lee_day3_pathway_means_arg1pos - lee_day3_pathway_means_arg1neg,
  median_difference = lee_day3_pathway_medians_arg1pos - lee_day3_pathway_medians_arg1neg,
  log2FC = lee_day3_pathway_log2fc,
  stringsAsFactors = FALSE
), error = function(e) data.frame())

# Exploratory pathway analysis: Use all PROGENy pathways (no pre-classification)
# Get all unique pathways from PROGENy network
lee_day3_all_pathways <- tryCatch(unique(progeny_network_lee_d3$source), error = function(e) character(0))
lee_day3_pathway_results <- lee_day3_pathway_comparison

# Report sample sizes (neutrophils only, Arg1+ vs Arg1-)
lee_day3_arg1pos_count <- length(lee_day3_pos)
lee_day3_arg1neg_count <- length(lee_day3_neg)
print(paste("Lee Day 3: Arg1pos neutrophils =", lee_day3_arg1pos_count, ", Arg1neg neutrophils =", lee_day3_arg1neg_count))

# Statistical test: Wilcoxon test for pathway differences (all pathways)
# With sample size adequacy checks, effect sizes, and confidence intervals
lee_day3_pathway_pvals <- numeric(length(lee_day3_all_pathways))
names(lee_day3_pathway_pvals) <- lee_day3_all_pathways
lee_day3_pathway_effect_sizes <- numeric(length(lee_day3_all_pathways))
names(lee_day3_pathway_effect_sizes) <- lee_day3_all_pathways
lee_day3_pathway_ci_lower <- numeric(length(lee_day3_all_pathways))
names(lee_day3_pathway_ci_lower) <- lee_day3_all_pathways
lee_day3_pathway_ci_upper <- numeric(length(lee_day3_all_pathways))
names(lee_day3_pathway_ci_upper) <- lee_day3_all_pathways
lee_day3_pathway_adequate_n <- logical(length(lee_day3_all_pathways))
names(lee_day3_pathway_adequate_n) <- lee_day3_all_pathways

for (pw_idx in seq_along(lee_day3_all_pathways)) {
  pw <- lee_day3_all_pathways[pw_idx]
  arg1pos_scores_pw <- tryCatch(as.numeric(lee_day3_progeny_arg1pos_scores[pw, ]), error = function(e) numeric(0))
  arg1neg_scores_pw <- tryCatch(as.numeric(lee_day3_progeny_arg1neg_scores[pw, ]), error = function(e) numeric(0))
  arg1pos_scores_pw <- arg1pos_scores_pw[!is.na(arg1pos_scores_pw)]
  arg1neg_scores_pw <- arg1neg_scores_pw[!is.na(arg1neg_scores_pw)]
  
  # Sample size adequacy check (Wilcoxon requires n ≥ 5 per group)
  n_arg1pos <- length(arg1pos_scores_pw)
  n_arg1neg <- length(arg1neg_scores_pw)
  lee_day3_pathway_adequate_n[pw] <- (n_arg1pos >= 5) & (n_arg1neg >= 5)
  
  # Calculate effect size (median difference)
  median_arg1pos <- tryCatch(stats::median(arg1pos_scores_pw, na.rm = TRUE), error = function(e) 0)
  median_arg1neg <- tryCatch(stats::median(arg1neg_scores_pw, na.rm = TRUE), error = function(e) 0)
  lee_day3_pathway_effect_sizes[pw] <- median_arg1pos - median_arg1neg
  
  # Wilcoxon test with confidence interval (linear: always attempt, handle errors)
  lee_day3_pathway_pvals[pw] <- 1.0
  lee_day3_pathway_ci_lower[pw] <- NA_real_
  lee_day3_pathway_ci_upper[pw] <- NA_real_
  wilcox_result <- tryCatch(stats::wilcox.test(arg1pos_scores_pw, arg1neg_scores_pw, conf.int = TRUE, conf.level = 0.95), error = function(e) NULL)
  wilcox_pvalue <- tryCatch(if (!is.null(wilcox_result)) wilcox_result$p.value else 1.0, error = function(e) 1.0)
  wilcox_pvalue_length <- tryCatch(length(wilcox_pvalue), error = function(e) 0)
  wilcox_pvalue_final <- tryCatch(if (wilcox_pvalue_length > 0) wilcox_pvalue[1] else 1.0, error = function(e) 1.0)
  lee_day3_pathway_pvals[pw] <- wilcox_pvalue_final
  wilcox_ci_lower <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[1] else NA_real_, error = function(e) NA_real_)
  wilcox_ci_upper <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[2] else NA_real_, error = function(e) NA_real_)
  lee_day3_pathway_ci_lower[pw] <- wilcox_ci_lower
  lee_day3_pathway_ci_upper[pw] <- wilcox_ci_upper
  lee_day3_pathway_warning_msg <- tryCatch(paste("Warning: Pathway", pw, "has insufficient sample size (Arg1pos n =", n_arg1pos, ", Arg1neg n =", n_arg1neg, "). Skipping statistical test."), error = function(e) "")
  lee_day3_pathway_warning_vector <- c("", lee_day3_pathway_warning_msg)
  lee_day3_pathway_warning_index <- tryCatch(as.numeric(!lee_day3_pathway_adequate_n[pw]) + 1, error = function(e) 1)
  print(lee_day3_pathway_warning_vector[lee_day3_pathway_warning_index])
}

# Multiple testing correction across ALL pathways (FDR)
lee_day3_pathway_pvals_adj <- p.adjust(lee_day3_pathway_pvals, method = "BH")

# Add statistical results to comparison dataframe (p-values, effect sizes, CIs, sample size adequacy)
lee_day3_pathway_results$p_value <- tryCatch(lee_day3_pathway_pvals[lee_day3_pathway_results$pathway], error = function(e) rep(1.0, nrow(lee_day3_pathway_results)))
lee_day3_pathway_results$p_adj <- tryCatch(lee_day3_pathway_pvals_adj[lee_day3_pathway_results$pathway], error = function(e) rep(1.0, nrow(lee_day3_pathway_results)))
lee_day3_pathway_results$effect_size_median_diff <- tryCatch(lee_day3_pathway_effect_sizes[lee_day3_pathway_results$pathway], error = function(e) rep(0.0, nrow(lee_day3_pathway_results)))
lee_day3_pathway_results$ci_lower_95 <- tryCatch(lee_day3_pathway_ci_lower[lee_day3_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(lee_day3_pathway_results)))
lee_day3_pathway_results$ci_upper_95 <- tryCatch(lee_day3_pathway_ci_upper[lee_day3_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(lee_day3_pathway_results)))
lee_day3_pathway_results$adequate_sample_size <- tryCatch(lee_day3_pathway_adequate_n[lee_day3_pathway_results$pathway], error = function(e) rep(FALSE, nrow(lee_day3_pathway_results)))
lee_day3_pathway_results$significant <- tryCatch((lee_day3_pathway_results$p_adj < 0.05) & lee_day3_pathway_results$adequate_sample_size, error = function(e) rep(FALSE, nrow(lee_day3_pathway_results)))

# Save PROGENy results
write.csv(lee_day3_pathway_comparison, file.path(OUTPUT_DIR, "LeeDay3_PROGENy_PathwayComparison.csv"), row.names = FALSE)
write.csv(lee_day3_pathway_results, file.path(OUTPUT_DIR, "LeeDay3_PROGENy_PathwayResults.csv"), row.names = FALSE)
saveRDS(list(pathway_comparison = lee_day3_pathway_comparison, pathway_results = lee_day3_pathway_results, pathway_means_arg1pos = lee_day3_pathway_means_arg1pos, pathway_means_arg1neg = lee_day3_pathway_means_arg1neg, progeny_network = progeny_network_lee_d3), file.path(OUTPUT_DIR, "LeeDay3_PROGENy_Results.rds"))
# lee_day3_progeny_loaded <- readRDS(file.path(OUTPUT_DIR, "LeeDay3_PROGENy_Results.rds")); lee_day3_pathway_comparison <- lee_day3_progeny_loaded$pathway_comparison; lee_day3_pathway_results <- lee_day3_progeny_loaded$pathway_results; lee_day3_pathway_means_arg1pos <- lee_day3_progeny_loaded$pathway_means_arg1pos; lee_day3_pathway_means_arg1neg <- lee_day3_progeny_loaded$pathway_means_arg1neg; progeny_network_lee_d3 <- lee_day3_progeny_loaded$progeny_network

# Ligand-to-pathway mapping: top-ranked LIANA interactions (Arg1pos + Arg1neg); fallback CellCall ligands if empty
lee_day3_identified_ligands_arg1pos <- tryCatch(unique(lee_day3_liana_arg1pos_topranked$ligand.complex), error = function(e) character(0))
lee_day3_identified_ligands_arg1neg <- tryCatch(unique(lee_day3_liana_arg1neg_topranked$ligand.complex), error = function(e) character(0))
lee_day3_identified_ligands <- tryCatch(unique(c(lee_day3_identified_ligands_arg1pos, lee_day3_identified_ligands_arg1neg)), error = function(e) character(0))
if (length(lee_day3_identified_ligands) == 0 && exists("lee_day3_cellcall_lr") && !is.null(lee_day3_cellcall_lr) && nrow(lee_day3_cellcall_lr) > 0) {
  lee_day3_identified_ligands <- unique(lee_day3_cellcall_lr$ligand.complex)
  lee_day3_identified_ligands <- lee_day3_identified_ligands[nzchar(lee_day3_identified_ligands)]
}
lee_day3_ligand_source_population <- tryCatch(c(rep("Arg1pos", length(lee_day3_identified_ligands_arg1pos)), rep("Arg1neg", length(lee_day3_identified_ligands_arg1neg))), error = function(e) character(0))
lee_day3_ligand_pathway_map <- data.frame(ligand = character(0), pathway = character(0), weight = numeric(0), stringsAsFactors = FALSE)
for (lig_idx in seq_along(lee_day3_identified_ligands)) {
  lig <- lee_day3_identified_ligands[lig_idx]
  lig_genes <- tryCatch(unlist(strsplit(lig, "_")), error = function(e) character(0))
  lig_genes <- tryCatch(unlist(strsplit(lig_genes, "[_+]")), error = function(e) lig_genes)
  lig_pathways <- tryCatch(dplyr::filter(progeny_network_lee_d3, target %in% lig_genes), error = function(e) data.frame())
  lig_pathway_df <- tryCatch(data.frame(ligand = lig, pathway = unique(lig_pathways$source), weight = lig_pathways$weight, stringsAsFactors = FALSE), error = function(e) data.frame(ligand = character(0), pathway = character(0), weight = numeric(0)))
  if (nrow(lig_pathway_df) > 0) lee_day3_ligand_pathway_map <- rbind(lee_day3_ligand_pathway_map, lig_pathway_df)
}
write.csv(lee_day3_ligand_pathway_map, file.path(OUTPUT_DIR, "LeeDay3_LigandToPathwayMapping.csv"), row.names = FALSE)

# Functional annotation: Exploratory annotation with population source (Arg1pos vs Arg1neg)
n_lig_d3 <- length(lee_day3_identified_ligands)
lee_day3_ligand_vec <- character(n_lig_d3)
lee_day3_ligand_genes_vec <- character(n_lig_d3)
lee_day3_target_pop_vec <- character(n_lig_d3)
for (lig_idx in seq_len(n_lig_d3)) {
  lig <- lee_day3_identified_ligands[lig_idx]
  lig_genes <- unlist(strsplit(lig, "_"))
  lig_genes <- unlist(strsplit(lig_genes, "[_+]"))
  lig_upper <- toupper(lig_genes)
  lig_in_arg1pos <- lig %in% lee_day3_identified_ligands_arg1pos
  lig_in_arg1neg <- lig %in% lee_day3_identified_ligands_arg1neg
  lig_population_vector <- c("Arg1pos", "Arg1neg")[c(lig_in_arg1pos, lig_in_arg1neg)]
  lig_population_source <- paste(lig_population_vector, collapse = ";")
  lee_day3_ligand_vec[lig_idx] <- lig
  lee_day3_ligand_genes_vec[lig_idx] <- paste(lig_upper, collapse = ";")
  lee_day3_target_pop_vec[lig_idx] <- lig_population_source
}
lee_day3_ligand_annotation <- data.frame(ligand = lee_day3_ligand_vec, ligand_genes = lee_day3_ligand_genes_vec, target_population = lee_day3_target_pop_vec, stringsAsFactors = FALSE)
write.csv(lee_day3_ligand_annotation, file.path(OUTPUT_DIR, "LeeDay3_LigandAnnotation.csv"), row.names = FALSE)

# PROGENy visualizations (Lee Day 3)
lee_day3_pathway_long <- tryCatch(tidyr::pivot_longer(lee_day3_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
p_lee_d3_progeny1 <- tryCatch(ggplot(lee_day3_pathway_long, aes(x = pathway, y = Pathway_Score, fill = Group)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("Arg1pos_mean" = "red", "Arg1neg_mean" = "lightblue"), labels = c("Arg1pos", "Arg1neg")) + labs(title = "PROGENy Pathway Activity - Lee Day 3", x = "Pathway", y = "Pathway Activity Score") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
vals_pw <- numeric(0)
if (exists("lee_day1_pathway_long") && is.data.frame(lee_day1_pathway_long) && nrow(lee_day1_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day1_pathway_long)) vals_pw <- c(vals_pw, lee_day1_pathway_long$Pathway_Score)
if (exists("lee_day3_pathway_long") && is.data.frame(lee_day3_pathway_long) && nrow(lee_day3_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day3_pathway_long)) vals_pw <- c(vals_pw, lee_day3_pathway_long$Pathway_Score)
if (exists("wang_day3_pathway_long") && is.data.frame(wang_day3_pathway_long) && nrow(wang_day3_pathway_long) > 0 && "Pathway_Score" %in% names(wang_day3_pathway_long)) vals_pw <- c(vals_pw, wang_day3_pathway_long$Pathway_Score)
if (exists("qin_day1_pathway_long") && is.data.frame(qin_day1_pathway_long) && nrow(qin_day1_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day1_pathway_long)) vals_pw <- c(vals_pw, qin_day1_pathway_long$Pathway_Score)
if (exists("qin_day3_pathway_long") && is.data.frame(qin_day3_pathway_long) && nrow(qin_day3_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day3_pathway_long)) vals_pw <- c(vals_pw, qin_day3_pathway_long$Pathway_Score)
vals_pw <- vals_pw[is.finite(vals_pw)]
UNIFY_PW_Y <- if (length(vals_pw) > 0) range(vals_pw) else c(-1, 1)
if (UNIFY_PW_Y[1] == UNIFY_PW_Y[2]) UNIFY_PW_Y[2] <- UNIFY_PW_Y[1] + 1e-6
print(p_lee_d3_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
progeny_heatmap_mat_ld3 <- tryCatch(rbind(Arg1pos = lee_day3_pathway_means_arg1pos, Arg1neg = lee_day3_pathway_means_arg1neg), error = function(e) matrix(0, nrow = 0, ncol = 0))
if (nrow(progeny_heatmap_mat_ld3) > 0 && ncol(progeny_heatmap_mat_ld3) > 0) {
  colors_progeny_ld3 <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors_use_progeny_ld3 <- grDevices::colorRampPalette(colors = colors_progeny_ld3)(100)
  p_progeny_heatmap_ld3 <- pheatmap::pheatmap(progeny_heatmap_mat_ld3, color = colors_use_progeny_ld3, border_color = "white", cellwidth = 20, cellheight = 20, main = "PROGENy Pathway Activity: Arg1+ vs Arg1- Neutrophils (Lee Day 3)")
  print(p_progeny_heatmap_ld3)
}
top15_lig_d3 <- if (nrow(lee_day3_ligand_pathway_map) > 0 && "ligand" %in% names(lee_day3_ligand_pathway_map)) head(unique(lee_day3_ligand_pathway_map$ligand), 15) else character(0)
lee_day3_ligand_map_plot <- if (length(top15_lig_d3) > 0) lee_day3_ligand_pathway_map[lee_day3_ligand_pathway_map$ligand %in% top15_lig_d3, ] else data.frame()
p_lee_d3_progeny2 <- if (nrow(lee_day3_ligand_map_plot) > 0) ggplot(lee_day3_ligand_map_plot, aes(x = ligand, y = pathway, size = abs(weight), color = weight)) + geom_point(alpha = 0.7) + scale_color_gradient2(low = "blue", mid = "white", high = "red") + labs(title = "Top 15 Ligands Linked to PROGENy Pathways - Lee Day 3", x = "Ligand", y = "Pathway") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)) else ggplot() + theme_void() + labs(title = "No ligands for pathway mapping (Lee Day 3)")
vals_w <- numeric(0)
map_tabnames <- c("lee_day1_ligand_pathway_map", "lee_day3_ligand_pathway_map", "wang_day3_ligand_pathway_map", "qin_day1_ligand_pathway_map", "qin_day3_ligand_pathway_map")
for (fn in map_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("weight" %in% names(d)) vals_w <- c(vals_w, d$weight)
}
vals_w <- vals_w[is.finite(vals_w)]
UNIFY_LIGW_ABS <- if (length(vals_w) > 0) range(abs(vals_w)) else c(0, 1)
if (UNIFY_LIGW_ABS[1] == UNIFY_LIGW_ABS[2]) UNIFY_LIGW_ABS[2] <- UNIFY_LIGW_ABS[1] + 1e-6
mxw <- if (length(vals_w) > 0) max(abs(vals_w)) else 1
if (!is.finite(mxw) || mxw <= 0) mxw <- 1
UNIFY_LIGW_COL <- c(-mxw, mxw)
print(p_lee_d3_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("✓ Lee Day 3 PROGENy pathway analysis complete")
# Save environment after Lee Day 3 PROGENy (resume: load(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay3_PROGENy.RData")))
stopifnot(exists("lee_day3"), exists("lee_day3_neut_cells"), exists("lee_day3_pathway_results"), exists("progeny_network_lee_d3"))
save.image(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay3_PROGENy.RData"))
saveRDS(list(lee_day3 = lee_day3, lee_day3_neut_cells = lee_day3_neut_cells, lee_day3_pos = lee_day3_pos, lee_day3_neg = lee_day3_neg, lee_day3_pathway_results = lee_day3_pathway_results, progeny_network_lee_d3 = progeny_network_lee_d3, lee_day3_liana_result_df = lee_day3_liana_result_df, lee_day3_arg1pos_received = lee_day3_arg1pos_received, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_PROGENy.rds"))
# checkpoint_ld3_progeny <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_PROGENy.rds")); list2env(checkpoint_ld3_progeny, envir = .GlobalEnv)
print("✓ Environment saved: Workspace_AfterLeeDay3_PROGENy.RData, Checkpoint_AfterLeeDay3_PROGENy.rds")

# FRESH RSTUDIO — RESUME LEE DAY 3 BARTsc: libraries as Lee Day 1 Step 2; then:
#   ck <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_PROGENy.rds")); list2env(ck, envir = .GlobalEnv)
#   BARTSC_N_LABELED_TFS <- 5L
# Fresh R session: load checkpoint after Lee Day 3 PROGENy + libraries; BARTsc block uses inline TF column detection (same pattern as Lee Day 1).
# Run from "Lee Day 3 BARTsc TF Analysis" through end of that BARTsc block (before "PROGENy Day 1 vs Day 3 diagnostic").

# -------- Lee Day 3 BARTsc TF Analysis (Arg1pos vs Arg1neg; downstream TF support only) --------
# Same BART_* thresholds as Lee Day 1. Single linear block: official vignette scRNA-seq.md §3–§7 (signature result before calc_crossCT_auc_RNA; dot_plot/deviation_heatmap after crossCT_test; find_key_regulators last).
print("--- BARTsc TF Analysis: Lee Day 3 (Arg1pos vs Arg1neg neutrophils) ---")
# Secondary support only: BARTsc prioritizes downstream TF programs associated with the Arg1 state.
if (bartsc_initialized) {
  options("mc.cores" = min(6L, parallel::detectCores()))
  lee_day3_neut_subset <- subset(lee_day3, cells = lee_day3_neut_cells)
  lee_day3_bart_rna_cnt <- Seurat::GetAssayData(lee_day3_neut_subset, layer = "counts")
  lee_day3_bart_label <- setNames(as.character(lee_day3_neut_subset$arg1_status), colnames(lee_day3_neut_subset))
  lee_day3_bart_label <- factor(lee_day3_bart_label, levels = c("Arg1pos", "Arg1neg"))
  print(table(lee_day3_bart_label))
  lee_day3_bart_proj <- BARTsc::bartsc(name = "LeeDay3_Arg1", genome = "mm10", label = lee_day3_bart_label, cell_types_used = c("Arg1pos", "Arg1neg"), RNA_cnt_matrix = lee_day3_bart_rna_cnt)
  lee_day3_bart_proj <- BARTsc::normalize_RNA(lee_day3_bart_proj)
  lee_day3_bart_proj <- BARTsc::find_signature_genes(lee_day3_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  lee_day3_bart_proj <- BARTsc::find_pairwise_deg(lee_day3_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  lee_day3_pw <- lee_day3_bart_proj@data$pairwise_DEG
  n_ld3_posneg <- 0L
  n_ld3_negpos <- 0L
  if (!is.null(lee_day3_pw) && "Arg1pos::Arg1neg" %in% names(lee_day3_pw)) n_ld3_posneg <- { x <- lee_day3_pw[["Arg1pos::Arg1neg"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  if (!is.null(lee_day3_pw) && "Arg1neg::Arg1pos" %in% names(lee_day3_pw)) n_ld3_negpos <- { x <- lee_day3_pw[["Arg1neg::Arg1pos"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  print(paste0("Lee Day 3 BARTsc pairwise DEG count (@data$pairwise_DEG): Arg1pos::Arg1neg=", n_ld3_posneg, ", Arg1neg::Arg1pos=", n_ld3_negpos))
  bart_deg_ok_ld3 <- if (BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT > 0) n_ld3_posneg >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT && n_ld3_negpos >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT else TRUE
  bart_nt_ld3 <- table(lee_day3_bart_label)
  bart_unequal_ld3 <- length(bart_nt_ld3) == 2L && length(unique(as.vector(bart_nt_ld3))) > 1L
  bart_crossct_allowed_ld3 <- !isTRUE(BARTSC_SKIP_CROSSCT_TEST) && length(bart_nt_ld3) >= 2L && all(as.integer(bart_nt_ld3) >= BARTSC_MIN_CELLS_CROSSCT) && !(isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_ld3) && bart_deg_ok_ld3
  lee_day3_bart_proj <- BARTsc::run_signature_RNA(lee_day3_bart_proj)
  lee_day3_bart_sig <- BARTsc::get_result(lee_day3_bart_proj, analysis = "cell type signature", mod = "RNA")
  lee_day3_bart_proj <- BARTsc::calc_crossCT_auc_RNA(lee_day3_bart_proj)
  print(bart_nt_ld3)
  if (!bart_crossct_allowed_ld3) message(paste0(
    "BARTsc Lee Day 3: crossCT_test / find_key_regulators skipped (BARTSC_SKIP_CROSSCT_TEST=", isTRUE(BARTSC_SKIP_CROSSCT_TEST),
    "; min_cells_ok=", all(as.integer(bart_nt_ld3) >= BARTSC_MIN_CELLS_CROSSCT),
    "; skip_if_unequal_n=", isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_ld3,
    "; deg_ok=", bart_deg_ok_ld3,
    " [Arg1pos::Arg1neg=", n_ld3_posneg,
    ", Arg1neg::Arg1pos=", n_ld3_negpos,
    ", min_each_dir=", BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT,
    "]). Vignette scRNA-seq.md §6: dot_plot/deviation_heatmap follow crossCT_test."
  ))
  if (bart_crossct_allowed_ld3) lee_day3_bart_proj <- BARTsc::crossCT_test(lee_day3_bart_proj, mod = "RNA")
  lee_day3_bart_cross <- if (bart_crossct_allowed_ld3) BARTsc::get_result(lee_day3_bart_proj, analysis = "cross-cell-type", mod = "RNA") else list()
  lee_day3_bart_cross_dev <- if (is.null(lee_day3_bart_cross)) list() else if (is.list(lee_day3_bart_cross) && "deviation" %in% names(lee_day3_bart_cross) && is.list(lee_day3_bart_cross$deviation)) lee_day3_bart_cross$deviation else lee_day3_bart_cross
  bart_outdir_lee_d3 <- file.path(OUTPUT_DIR, "BARTsc_LeeDay3")
  dir.create(bart_outdir_lee_d3, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(lee_day3_bart_sig) && !is.null(lee_day3_bart_sig[["Arg1pos"]])) write.csv(lee_day3_bart_sig[["Arg1pos"]], file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_Signature_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(lee_day3_bart_sig) && !is.null(lee_day3_bart_sig[["Arg1neg"]])) write.csv(lee_day3_bart_sig[["Arg1neg"]], file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_Signature_Arg1neg.csv"), row.names = FALSE)
  lee_day3_bart_key <- NULL
  if (bart_crossct_allowed_ld3) {
    lee_day3_bart_proj <- BARTsc::find_key_regulators(lee_day3_bart_proj, mod = "RNA", min.N.profile = 3)
    lee_day3_bart_key <- BARTsc::get_result(lee_day3_bart_proj, analysis = "Key regs ident", mod = "RNA")
  }
  if (!is.null(lee_day3_bart_key) && !is.null(lee_day3_bart_key[["Arg1pos"]])) write.csv(lee_day3_bart_key[["Arg1pos"]], file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_KeyRegulators_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(lee_day3_bart_key) && !is.null(lee_day3_bart_key[["Arg1neg"]])) write.csv(lee_day3_bart_key[["Arg1neg"]], file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_KeyRegulators_Arg1neg.csv"), row.names = FALSE)
  saveRDS(lee_day3_bart_proj, file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_Object.rds"))
  print("✓ BARTsc TF analysis complete (Lee Day 3) — plots: Lee Day 3 BARTsc visualizations section below")
}

# -------- Lee Day 3: BARTsc visualizations (reload LeeDay3_BARTsc_Object.rds; no re-run of bartsc pipeline) --------
bart_viz_ld3_rds <- file.path(OUTPUT_DIR, "BARTsc_LeeDay3", "LeeDay3_BARTsc_Object.rds")
if (!requireNamespace("BARTsc", quietly = TRUE)) message("Lee Day 3 BARTsc visualizations skipped: package BARTsc not installed.")
if (requireNamespace("BARTsc", quietly = TRUE) && !file.exists(bart_viz_ld3_rds)) message(paste0("Lee Day 3 BARTsc visualizations skipped: missing ", bart_viz_ld3_rds))
if (requireNamespace("BARTsc", quietly = TRUE) && file.exists(bart_viz_ld3_rds)) {
  if (!exists("BARTSC_N_TF_DOTHEAT")) BARTSC_N_TF_DOTHEAT <- 5L
  if (!exists("BARTSC_N_LABELED_TFS")) BARTSC_N_LABELED_TFS <- 5L
  if (!exists("BARTSC_KEYREG_SCATTER_W_IN")) BARTSC_KEYREG_SCATTER_W_IN <- 11
  if (!exists("BARTSC_KEYREG_SCATTER_H_IN")) BARTSC_KEYREG_SCATTER_H_IN <- 7
  if (!exists("BARTSC_KEYREG_XLIM")) BARTSC_KEYREG_XLIM <- c(-1, 1)
  if (!exists("BARTSC_KEYREG_YLIM_SIG")) BARTSC_KEYREG_YLIM_SIG <- c(0, 7)
  if (!exists("BARTSC_KEYREG_ZLIM_DMDR")) BARTSC_KEYREG_ZLIM_DMDR <- c(-2, 2)
  if (!exists("BARTSC_KEYREG_RANK_MAX")) BARTSC_KEYREG_RANK_MAX <- 60L
  suppressWarnings(tryCatch({ BARTsc::load_bart2(); NULL }, error = function(e) NULL))
  if (exists("bart2", envir = .GlobalEnv)) {
    if (!exists("types", envir = .GlobalEnv)) types <<- reticulate::import("types")
    bart2_mods_viz_ld3 <- reticulate::py_to_r(reticulate::py_get_attr(bart2, "__all__"))
    for (m_viz_ld3 in bart2_mods_viz_ld3) {
      if (!exists(m_viz_ld3, envir = .GlobalEnv)) assign(m_viz_ld3, reticulate::import(paste0("bart2.", m_viz_ld3), delay_load = TRUE), envir = .GlobalEnv)
    }
  }
  if (!exists("lee_day3_bart_proj", envir = .GlobalEnv, inherits = FALSE)) lee_day3_bart_proj <- readRDS(bart_viz_ld3_rds)
  bart_outdir_lee_d3 <- file.path(OUTPUT_DIR, "BARTsc_LeeDay3")
  dir.create(bart_outdir_lee_d3, showWarnings = FALSE, recursive = TRUE)
  if (interactive() && grDevices::dev.cur() == 1L) grDevices::dev.new()
  lee_day3_bart_cross_viz <- BARTsc::get_result(lee_day3_bart_proj, analysis = "cross-cell-type", mod = "RNA")
  lee_day3_bart_cross_dev <- if (is.null(lee_day3_bart_cross_viz)) list() else if (is.list(lee_day3_bart_cross_viz) && "deviation" %in% names(lee_day3_bart_cross_viz) && is.list(lee_day3_bart_cross_viz$deviation)) lee_day3_bart_cross_viz$deviation else lee_day3_bart_cross_viz
  lee_day3_bart_key <- BARTsc::get_result(lee_day3_bart_proj, analysis = "Key regs ident", mod = "RNA")
  bart_tf_names_lee_d3 <- names(lee_day3_bart_cross_dev)
  tf_example_lee_d3 <- if (length(bart_tf_names_lee_d3) > 0) bart_tf_names_lee_d3[1] else NULL
  if (length(lee_day3_bart_cross_dev) > 0 && !is.null(tf_example_lee_d3)) {
    p_lee_d3_bart_dot <- BARTsc::dot_plot(lee_day3_bart_proj, mod = "RNA", tf = tf_example_lee_d3, max_dot_size = 22)
    print(p_lee_d3_bart_dot)
    ggplot2::ggsave(file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_DotPlot.png"), p_lee_d3_bart_dot, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_DotPlot.pdf"), p_lee_d3_bart_dot, width = 8, height = 6)
    p_lee_d3_bart_heat <- BARTsc::deviation_heatmap(lee_day3_bart_proj, mod = "RNA", tf = tf_example_lee_d3, tile_fontsize = 6)
    print(p_lee_d3_bart_heat)
    ggplot2::ggsave(file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_DeviationHeatmap.png"), p_lee_d3_bart_heat, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_DeviationHeatmap.pdf"), p_lee_d3_bart_heat, width = 8, height = 6)
  }
  top_tf_pos_lee_d3 <- if (length(lee_day3_bart_cross_dev) > 0 && !is.null(lee_day3_bart_key) && !is.null(lee_day3_bart_key[["Arg1pos"]]) && "TF" %in% colnames(lee_day3_bart_key[["Arg1pos"]])) head(lee_day3_bart_key[["Arg1pos"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_tf_neg_lee_d3 <- if (length(lee_day3_bart_cross_dev) > 0 && !is.null(lee_day3_bart_key) && !is.null(lee_day3_bart_key[["Arg1neg"]]) && "TF" %in% colnames(lee_day3_bart_key[["Arg1neg"]])) head(lee_day3_bart_key[["Arg1neg"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  tfs_available_lee_d3 <- names(lee_day3_bart_cross_dev)
  top_key_tfs_lee_d3 <- intersect(unique(c(top_tf_pos_lee_d3, top_tf_neg_lee_d3)), tfs_available_lee_d3)
  if (length(lee_day3_bart_cross_dev) > 0 && length(top_key_tfs_lee_d3) > 0) {
    i_tf_ld3 <- 1L
    while (i_tf_ld3 <= length(top_key_tfs_lee_d3)) {
      tf_cur_lee_d3 <- top_key_tfs_lee_d3[i_tf_ld3]
      p_dot_tf_lee_d3 <- BARTsc::dot_plot(lee_day3_bart_proj, mod = "RNA", tf = tf_cur_lee_d3, max_dot_size = 22)
      print(p_dot_tf_lee_d3)
      ggplot2::ggsave(file.path(bart_outdir_lee_d3, paste0("LeeDay3_BARTsc_DotPlot_", tf_cur_lee_d3, ".png")), p_dot_tf_lee_d3, width = 8, height = 6)
      p_heat_tf_lee_d3 <- BARTsc::deviation_heatmap(lee_day3_bart_proj, mod = "RNA", tf = tf_cur_lee_d3, tile_fontsize = 6)
      print(p_heat_tf_lee_d3)
      ggplot2::ggsave(file.path(bart_outdir_lee_d3, paste0("LeeDay3_BARTsc_DeviationHeatmap_", tf_cur_lee_d3, ".png")), p_heat_tf_lee_d3, width = 8, height = 6)
      i_tf_ld3 <- i_tf_ld3 + 1L
    }
  }
  tfs_labeled_lee_d3_pos <- character(0)
  if (!is.null(lee_day3_bart_key) && !is.null(lee_day3_bart_key[["Arg1pos"]])) {
    df_ld3p <- lee_day3_bart_key[["Arg1pos"]]
    if (nrow(df_ld3p) > 0) {
      tf_c_ld3p <- intersect(colnames(df_ld3p), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_ld3p <- if (length(tf_c_ld3p) == 0) colnames(df_ld3p)[1] else tf_c_ld3p[1]
      ord_ld3p <- if ("final_rank" %in% colnames(df_ld3p)) order(df_ld3p$final_rank, na.last = TRUE) else seq_len(nrow(df_ld3p))
      tfs_labeled_lee_d3_pos <- as.character(head(df_ld3p[ord_ld3p, , drop = FALSE][[tf_col_ld3p]], BARTSC_N_LABELED_TFS))
    }
  }
  tfs_labeled_lee_d3_neg <- character(0)
  if (!is.null(lee_day3_bart_key) && !is.null(lee_day3_bart_key[["Arg1neg"]])) {
    df_ld3n <- lee_day3_bart_key[["Arg1neg"]]
    if (nrow(df_ld3n) > 0) {
      tf_c_ld3n <- intersect(colnames(df_ld3n), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_ld3n <- if (length(tf_c_ld3n) == 0) colnames(df_ld3n)[1] else tf_c_ld3n[1]
      ord_ld3n <- if ("final_rank" %in% colnames(df_ld3n)) order(df_ld3n$final_rank, na.last = TRUE) else seq_len(nrow(df_ld3n))
      tfs_labeled_lee_d3_neg <- as.character(head(df_ld3n[ord_ld3n, , drop = FALSE][[tf_col_ld3n]], BARTSC_N_LABELED_TFS))
    }
  }
  if (length(tfs_labeled_lee_d3_pos) > 0) {
    grDevices::png(file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_KeyRegScatter_Official_Arg1pos.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(lee_day3_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_labeled_lee_d3_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_labeled_lee_d3_pos) > 0) message("Displaying official key_regulator_scatter for Arg1pos (Lee Day 3)...")
  if (interactive() && length(tfs_labeled_lee_d3_pos) > 0) key_regulator_scatter_unified(lee_day3_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_labeled_lee_d3_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
  if (length(tfs_labeled_lee_d3_neg) > 0) {
    grDevices::png(file.path(bart_outdir_lee_d3, "LeeDay3_BARTsc_KeyRegScatter_Official_Arg1neg.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(lee_day3_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_labeled_lee_d3_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_labeled_lee_d3_neg) > 0) message("Displaying official key_regulator_scatter for Arg1neg (Lee Day 3)...")
  if (interactive() && length(tfs_labeled_lee_d3_neg) > 0) key_regulator_scatter_unified(lee_day3_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_labeled_lee_d3_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
  print("✓ Lee Day 3 BARTsc visualizations complete")
}

# ============================================================================
# INTEGRATION: LIANA -> BARTsc + PROGENy (Lee Day 3) — same logic as Lee Day 1
# ============================================================================
# Linear: footprint co-membership CSV → OmniPath receptor→TF → full-chain CSV → sender×pathway heatmap (BARTsc plots: section above).
print("--- Integration: LIANA -> BARTsc + PROGENy (Lee Day 3, Data-Driven) ---")
# Exploratory overlay only: LIANA provides incoming-signal candidates; PROGENy/BARTsc/OmniPath add downstream support, not causal proof.

if (!exists("progeny_network_lee_d3") || !all(c("source", "target") %in% colnames(progeny_network_lee_d3))) {
  stop("Lee Day 3 integration requires progeny_network_lee_d3. Run Lee Day 3 PROGENy section first.")
}
liana_receptors_ld3 <- unique(lee_day3_arg1pos_received$receptor.complex)
liana_receptors_ld3_split <- unique(trimws(unlist(strsplit(liana_receptors_ld3, "[_+]"))))
print(paste("Unique receptors received by Arg1+ neutrophils (LIANA, Lee Day 3):", length(liana_receptors_ld3_split)))

liana_expanded_ld3 <- lee_day3_arg1pos_received
liana_expanded_ld3$receptor_subunits <- lapply(strsplit(liana_expanded_ld3$receptor.complex, "[_+]"), function(x) trimws(x))
liana_expanded_ld3 <- tidyr::unnest_longer(liana_expanded_ld3, receptor_subunits)

pathways_active_ld3 <- names(lee_day3_pathway_means_arg1pos)[lee_day3_pathway_means_arg1pos > lee_day3_pathway_means_arg1neg]

# Exploratory footprint overlap (not causal receptor -> pathway activation)
receptor_pathway_links_ld3 <- progeny_network_lee_d3[
  progeny_network_lee_d3$target %in% liana_receptors_ld3_split &
    progeny_network_lee_d3$source %in% pathways_active_ld3,
]
receptor_pathway_links_ld3 <- receptor_pathway_links_ld3[, c("target", "source")]
colnames(receptor_pathway_links_ld3) <- c("receptor_gene", "pathway")

pathway_diff_ld3 <- lee_day3_pathway_means_arg1pos - lee_day3_pathway_means_arg1neg
receptor_pathway_links_ld3$pathway_activity_diff <- pathway_diff_ld3[receptor_pathway_links_ld3$pathway]

integration_receptor_pathway_ld3 <- merge(
  liana_expanded_ld3,
  receptor_pathway_links_ld3,
  by.x = "receptor_subunits",
  by.y = "receptor_gene",
  all.x = TRUE
)
integration_receptor_pathway_ld3 <- integration_receptor_pathway_ld3[!is.na(integration_receptor_pathway_ld3$pathway), ]

write.csv(integration_receptor_pathway_ld3, file.path(OUTPUT_DIR, "LeeDay3_LIANA_PROGENy_FootprintCoMembership_Exploratory.csv"), row.names = FALSE)
print(paste("LIANA x PROGENy footprint co-membership (exploratory, Lee Day 3):", nrow(integration_receptor_pathway_ld3), "rows"))

receptor_tf_links_ld3 <- data.frame(receptor_gene = character(0), tf = character(0), omnipath_via = character(0), omnipath_hops = integer(0), stringsAsFactors = FALSE)
active_tfs_ld3 <- character(0)

has_bart_key_ld3 <- exists("lee_day3_bart_key") &&
  !is.null(lee_day3_bart_key) &&
  !is.null(lee_day3_bart_key[["Arg1pos"]]) &&
  nrow(lee_day3_bart_key[["Arg1pos"]]) > 0

if (has_bart_key_ld3) {
  bart_tfs_ld3 <- lee_day3_bart_key[["Arg1pos"]]
  tf_col_ld3 <- intersect(colnames(bart_tfs_ld3), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
  if (length(tf_col_ld3) == 0) {
    tf_col_ld3 <- colnames(bart_tfs_ld3)[1]
  } else {
    tf_col_ld3 <- tf_col_ld3[1]
  }
  active_tfs_ld3 <- unique(as.character(bart_tfs_ld3[[tf_col_ld3]]))
}

omnipath_ok_ld3 <- requireNamespace("OmnipathR", quietly = TRUE) && length(active_tfs_ld3) > 0

pathway_interactions_ld3 <- data.frame()
omni_ld3_first <- if (!omnipath_ok_ld3) list(dat = NULL, err1 = NA_character_) else tryCatch(list(dat = OmnipathR::import_omnipath_interactions(datasets = c("omnipath", "pathwayextra"), organism = 10090, genesymbols = TRUE), err1 = NA_character_), error = function(e1) list(dat = NULL, err1 = conditionMessage(e1)))
if (omnipath_ok_ld3) pathway_interactions_ld3 <- omni_ld3_first$dat
omnipath_ld3_err1 <- omni_ld3_first$err1
if (omnipath_ok_ld3 && is.null(pathway_interactions_ld3)) pathway_interactions_ld3 <- tryCatch(OmnipathR::import_pathwayextra_interactions(organism = 10090, genesymbols = TRUE), error = function(e2) { message("OmniPath curated+pathwayextra import failed (Lee Day 3): ", omnipath_ld3_err1, "; fallback: ", conditionMessage(e2)); data.frame() })

omni_ld3_cols_ok <- omnipath_ok_ld3 && nrow(pathway_interactions_ld3) > 0 && all(c("source_genesymbol", "target_genesymbol") %in% colnames(pathway_interactions_ld3))
omni_direct_ld3 <- if (omni_ld3_cols_ok) subset(pathway_interactions_ld3, source_genesymbol %in% liana_receptors_ld3_split & target_genesymbol %in% active_tfs_ld3, select = c("source_genesymbol", "target_genesymbol")) else data.frame()
if (omni_ld3_cols_ok && nrow(omni_direct_ld3) > 0) receptor_tf_links_ld3 <- rbind(receptor_tf_links_ld3, data.frame(receptor_gene = omni_direct_ld3$source_genesymbol, tf = omni_direct_ld3$target_genesymbol, omnipath_via = NA_character_, omnipath_hops = 1L, stringsAsFactors = FALSE))

hop1_ld3 <- if (omni_ld3_cols_ok) pathway_interactions_ld3[pathway_interactions_ld3$source_genesymbol %in% liana_receptors_ld3_split, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_ld3_cols_ok) colnames(hop1_ld3) <- c("receptor_gene", "via")
hop2_ld3 <- if (omni_ld3_cols_ok) pathway_interactions_ld3[pathway_interactions_ld3$target_genesymbol %in% active_tfs_ld3, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_ld3_cols_ok) colnames(hop2_ld3) <- c("via", "tf")
omni_2hop_ld3 <- if (omni_ld3_cols_ok) merge(hop1_ld3, hop2_ld3, by = "via") else data.frame()
if (omni_ld3_cols_ok && nrow(omni_2hop_ld3) > 0) omni_2hop_ld3 <- unique(omni_2hop_ld3[, c("receptor_gene", "tf", "via")])
if (omni_ld3_cols_ok && nrow(omni_2hop_ld3) > 0) receptor_tf_links_ld3 <- rbind(receptor_tf_links_ld3, data.frame(receptor_gene = omni_2hop_ld3$receptor_gene, tf = omni_2hop_ld3$tf, omnipath_via = omni_2hop_ld3$via, omnipath_hops = 2L, stringsAsFactors = FALSE))
if (omni_ld3_cols_ok) receptor_tf_links_ld3 <- receptor_tf_links_ld3[!duplicated(paste(receptor_tf_links_ld3$receptor_gene, receptor_tf_links_ld3$tf)), ]

integration_receptor_tf_ld3 <- data.frame()
if (nrow(receptor_tf_links_ld3) > 0) {
  integration_receptor_tf_ld3 <- merge(
    liana_expanded_ld3,
    receptor_tf_links_ld3,
    by.x = "receptor_subunits",
    by.y = "receptor_gene",
    all.x = TRUE
  )
  integration_receptor_tf_ld3 <- integration_receptor_tf_ld3[!is.na(integration_receptor_tf_ld3$tf), ]
}

if (nrow(integration_receptor_tf_ld3) > 0) {
  write.csv(integration_receptor_tf_ld3, file.path(OUTPUT_DIR, "LeeDay3_LIANA_BARTsc_Integration_DataDriven.csv"), row.names = FALSE)
  n1_ld3 <- sum(integration_receptor_tf_ld3$omnipath_hops == 1L, na.rm = TRUE)
  n2_ld3 <- sum(integration_receptor_tf_ld3$omnipath_hops == 2L, na.rm = TRUE)
  print(paste0(
    "LIANA -> BARTsc integration (Lee Day 3): ",
    nrow(integration_receptor_tf_ld3),
    " receptor-TF rows (OmniPath curated+pathwayextra; direct=",
    n1_ld3,
    ", two-hop=",
    n2_ld3,
    ")"
  ))
}

if (nrow(receptor_tf_links_ld3) == 0 && length(active_tfs_ld3) > 0) {
  fallback_ld3 <- data.frame(
    receptor   = liana_receptors_ld3_split,
    note       = "Active TFs in Arg1pos (no OmniPath direct/two-hop link):",
    active_tfs = paste(active_tfs_ld3, collapse = "; "),
    stringsAsFactors = FALSE
  )
  write.csv(fallback_ld3, file.path(OUTPUT_DIR, "LeeDay3_LIANA_BARTsc_Receptors_and_TFs.csv"), row.names = FALSE)
  print("LIANA receptors and BARTsc TFs saved separately (Lee Day 3, no OmniPath link)")
}

if (nrow(receptor_tf_links_ld3) == 0) {
  message("
  Note (Lee Day 3): No receptor-TF pairs in OmniPath curated+pathwayextra (direct or two-hop).
  Main chain: Receptors (LIANA) -> Pathways (PROGENy) -> TFs (BARTsc)
  ")
}

has_pathway_ld3 <- nrow(integration_receptor_pathway_ld3) > 0
has_tf_ld3 <- nrow(integration_receptor_tf_ld3) > 0
full_chain_ld3 <- data.frame()
if (has_pathway_ld3 && has_tf_ld3) full_chain_ld3 <- merge(integration_receptor_pathway_ld3, integration_receptor_tf_ld3[, c("source", "ligand.complex", "receptor_subunits", "aggregate_rank", "tf", "omnipath_via", "omnipath_hops")], by = c("source", "ligand.complex", "receptor_subunits", "aggregate_rank"), all = TRUE)
if (has_pathway_ld3 && !has_tf_ld3) { full_chain_ld3 <- integration_receptor_pathway_ld3; full_chain_ld3$tf <- NA_character_ }
if (!has_pathway_ld3 && has_tf_ld3) { full_chain_ld3 <- integration_receptor_tf_ld3; full_chain_ld3$pathway <- NA_character_; full_chain_ld3$pathway_activity_diff <- NA_real_ }
if (nrow(full_chain_ld3) > 0) full_chain_ld3$evidence_level <- ifelse(!is.na(full_chain_ld3$pathway) & !is.na(full_chain_ld3$tf), "STRONG (pathway + TF linked)", ifelse(!is.na(full_chain_ld3$pathway) | !is.na(full_chain_ld3$tf), "MODERATE (pathway or TF linked)", "WEAK (no link)"))
if (nrow(full_chain_ld3) > 0) full_chain_ld3 <- full_chain_ld3[order(full_chain_ld3$evidence_level, full_chain_ld3$aggregate_rank), ]
# Legacy filename retained for compatibility; contents are an exploratory overlay, not a causal signal chain.
if (nrow(full_chain_ld3) > 0) write.csv(full_chain_ld3, file.path(OUTPUT_DIR, "LeeDay3_FullSignalChain_DataDriven.csv"), row.names = FALSE)
if (nrow(full_chain_ld3) > 0) {
  cols_show_ld3 <- intersect(c("source", "receptor_subunits", "pathway", "tf", "omnipath_via", "omnipath_hops", "evidence_level"), colnames(full_chain_ld3))
  print(head(full_chain_ld3[, cols_show_ld3, drop = FALSE], 20))
}
heatmap_data_ld3 <- data.frame()
if (nrow(full_chain_ld3) > 0 && has_pathway_ld3) heatmap_data_ld3 <- as.data.frame.matrix(table(integration_receptor_pathway_ld3$source, integration_receptor_pathway_ld3$pathway))
if (nrow(heatmap_data_ld3) > 0 && ncol(heatmap_data_ld3) > 0) pheatmap::pheatmap(heatmap_data_ld3, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("white", "blue", "red"))(50), main = "Lee Day 3: Sender x pathway (PROGENy footprint co-membership, exploratory)")
if (nrow(full_chain_ld3) == 0) print("No integration rows Lee Day 3 (PROGENy/BARTsc may not overlap LIANA receptors)")
print("✓ Lee Day 3 LIANA <-> BARTsc <-> PROGENy integration complete (exploratory overlay)")
save.image(file.path(OUTPUT_DIR, "Workspace_AfterLeeDay3_BARTscIntegration.RData"))
checkpoint_ld3_integration <- list(lee_day3 = lee_day3, lee_day3_liana_result_df = lee_day3_liana_result_df, lee_day3_pathway_results = lee_day3_pathway_results, progeny_network_lee_d3 = progeny_network_lee_d3, OUTPUT_DIR = OUTPUT_DIR)
if (exists("lee_day3_bart_key")) checkpoint_ld3_integration$lee_day3_bart_key <- lee_day3_bart_key
saveRDS(checkpoint_ld3_integration, file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_BARTscIntegration.rds"))
# checkpoint_ld3_int <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterLeeDay3_BARTscIntegration.rds")); list2env(checkpoint_ld3_int, envir = .GlobalEnv)
# If resuming BARTsc figures only: run the cohort's "BARTsc visualizations" block; ensure OUTPUT_DIR, BARTSC_N_TF_DOTHEAT, BARTSC_N_LABELED_TFS, BARTSC_KEYREG_SCATTER_* from Lee Day 1 parameter block (or set in session).
print("✓ Saved: Workspace_AfterLeeDay3_BARTscIntegration.RData, Checkpoint_AfterLeeDay3_BARTscIntegration.rds")

# Diagnostic: verify Day 1 vs Day 3 use different data (print sample sizes and sample values)
print("--- PROGENy Day 1 vs Day 3 diagnostic ---")
if (!exists("lee_day1_arg1pos_count") && exists("lee_day1_pos")) lee_day1_arg1pos_count <- length(lee_day1_pos)
if (!exists("lee_day1_arg1neg_count") && exists("lee_day1_neg")) lee_day1_arg1neg_count <- length(lee_day1_neg)
if (!exists("lee_day3_arg1pos_count") && exists("lee_day3_pos")) lee_day3_arg1pos_count <- length(lee_day3_pos)
if (!exists("lee_day3_arg1neg_count") && exists("lee_day3_neg")) lee_day3_arg1neg_count <- length(lee_day3_neg)
print(paste("Day 1: Arg1pos n =", if (exists("lee_day1_arg1pos_count")) lee_day1_arg1pos_count else "?", ", Arg1neg n =", if (exists("lee_day1_arg1neg_count")) lee_day1_arg1neg_count else "?"))
print(paste("Day 3: Arg1pos n =", if (exists("lee_day3_arg1pos_count")) lee_day3_arg1pos_count else "?", ", Arg1neg n =", if (exists("lee_day3_arg1neg_count")) lee_day3_arg1neg_count else "?"))
nfkb_d1_pos <- lee_day1_pathway_results$Arg1pos_mean[lee_day1_pathway_results$pathway == "NFkB"]
nfkb_d1_neg <- lee_day1_pathway_results$Arg1neg_mean[lee_day1_pathway_results$pathway == "NFkB"]
nfkb_d3_pos <- lee_day3_pathway_results$Arg1pos_mean[lee_day3_pathway_results$pathway == "NFkB"]
nfkb_d3_neg <- lee_day3_pathway_results$Arg1neg_mean[lee_day3_pathway_results$pathway == "NFkB"]
print(paste("NFkB: Day1 Arg1pos =", round(nfkb_d1_pos, 4), ", Day1 Arg1neg =", round(nfkb_d1_neg, 4)))
print(paste("NFkB: Day3 Arg1pos =", round(nfkb_d3_pos, 4), ", Day3 Arg1neg =", round(nfkb_d3_neg, 4)))
print(paste("Day1 vs Day3 identical?", identical(lee_day1_pathway_results, lee_day3_pathway_results)))

# Combined Day 1 vs Day 3 PROGENy comparison (shows temporal differences)
lee_d1_long <- tryCatch(tidyr::pivot_longer(lee_day1_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
lee_d1_long$Timepoint <- "Day 1"
lee_d1_long$Group <- gsub("_mean", "", lee_d1_long$Group)
lee_d3_long <- tryCatch(tidyr::pivot_longer(lee_day3_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
lee_d3_long$Timepoint <- "Day 3"
lee_d3_long$Group <- gsub("_mean", "", lee_d3_long$Group)
lee_d1d3_combined <- tryCatch(rbind(lee_d1_long, lee_d3_long), error = function(e) data.frame())
lee_d1d3_combined$Timepoint_Group <- paste0(lee_d1d3_combined$Timepoint, " ", lee_d1d3_combined$Group)
if (nrow(lee_d1d3_combined) > 0) {
  p_lee_d1d3_progeny <- ggplot(lee_d1d3_combined, aes(x = pathway, y = Pathway_Score, fill = Timepoint_Group)) + geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) + scale_fill_manual(values = c("Day 1 Arg1pos" = "#87CEEB", "Day 1 Arg1neg" = "#CD5C5C", "Day 3 Arg1pos" = "#4682B4", "Day 3 Arg1neg" = "#8B0000")) + labs(title = "PROGENy Pathway Activity - Lee Day 1 vs Day 3", x = "Pathway", y = "Pathway Activity Score") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
  print(p_lee_d1d3_progeny)
}

# -------- Wang Day 3: LIANA → LIANA plots → top-receptor Vln → PROGENy → BARTsc → integration (same cohort order as Lee Day 1; no separate CellChat CSV overlap block) --------
# Run LIANA on Seurat. Use custom LR from CellChat+CellCall when available; else Consensus.
print("--- LIANA Analysis: Wang Day 3 ---")

# Preamble: when resuming from Checkpoint_AfterLeeDay3_BARTscIntegration.rds, Wang objects are missing. Rebuild from wang_dat.
if (!exists("wang_day3_pruned_labels")) {
  if (!exists("NEUTROPHIL_BASE_LABEL")) NEUTROPHIL_BASE_LABEL <- "Neutrophil"
  if (!exists("NEUTROPHIL_ALIASES")) NEUTROPHIL_ALIASES <- c("cNeutrophil", "Neutrophils")
  if (!exists("GENE_CANDIDATES")) GENE_CANDIDATES <- c("Arg1", "ARG1")
  if (!exists("wang_dat")) {
    wang_dat <- readRDS("WangDat.rds")
    DefaultAssay(wang_dat) <- "RNA"
    wang_dat$time <- as.character(wang_dat$time)
    wang_dat$time <- gsub("Uninjured", "0", wang_dat$time)
    wang_dat$time <- gsub("1dpi", "1", wang_dat$time)
    wang_dat$time <- gsub("3dpi", "3", wang_dat$time)
    wang_dat$time <- gsub("7dpi", "7", wang_dat$time)
    wang_dat$time <- as.numeric(wang_dat$time)
  }
  wang_day3 <- subset(wang_dat, time == 3)
  wang_day3_label_col <- if ("celltype" %in% colnames(wang_day3@meta.data)) wang_day3$celltype else wang_day3$pruned_labels
  wang_day3_pruned_labels <- as.character(wang_day3_label_col)
  wang_day3_pruned_labels <- gsub("-", "", wang_day3_pruned_labels)
  wang_day3_pruned_labels[wang_day3_pruned_labels %in% c(NEUTROPHIL_ALIASES, NEUTROPHIL_BASE_LABEL)] <- NEUTROPHIL_BASE_LABEL
  wang_day3_pruned_labels[is.na(wang_day3_pruned_labels) | wang_day3_pruned_labels == ""] <- "Other"
  wang_day3_neut <- which(wang_day3_pruned_labels == NEUTROPHIL_BASE_LABEL)
  wang_day3_gene_candidates <- GENE_CANDIDATES[GENE_CANDIDATES %in% rownames(wang_day3)]
  wang_day3_gene <- head(wang_day3_gene_candidates, 1)
  wang_day3_expr <- tryCatch(Seurat::GetAssayData(wang_day3, layer = "data")[wang_day3_gene, wang_day3_neut, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
  wang_day3_expr_vec <- as.numeric(wang_day3_expr)
  wang_day3_indicator <- wang_day3_expr_vec > 0
  wang_day3_indicator[is.na(wang_day3_indicator)] <- FALSE
  wang_day3_pos <- wang_day3_neut[wang_day3_indicator]
  wang_day3_neg <- setdiff(wang_day3_neut, wang_day3_pos)
  wang_day3_arg1_status <- rep("Arg1neg", ncol(wang_day3))
  wang_day3_arg1_status[wang_day3_pos] <- "Arg1pos"
  wang_day3$arg1_status <- wang_day3_arg1_status
  Idents(wang_day3) <- factor(wang_day3_pruned_labels)
  wang_day3_neut_cells <- colnames(wang_day3)[c(wang_day3_pos, wang_day3_neg)]
  print("Wang Day 3 objects rebuilt from wang_dat (resumed from checkpoint).")
}

wang_day3_liana_labels <- wang_day3_pruned_labels
wang_day3_liana_labels[wang_day3_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
wang_day3_liana_labels[wang_day3_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
wang_day3_liana_labels[is.na(wang_day3_liana_labels) | wang_day3_liana_labels == ""] <- "Other"
Seurat::Idents(wang_day3) <- factor(wang_day3_liana_labels)

# Wang Day 3: set labels, clean Seurat, build SCE (bypass validObject bug), run LIANA
wang_day3$liana_label <- factor(wang_day3_liana_labels)
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
for (assay_name in names(wang_day3@assays)) {
  if (assay_name != "RNA") wang_day3[[assay_name]] <- NULL
}
for (red_name in names(wang_day3@reductions)) wang_day3[[red_name]] <- NULL
for (graph_name in names(wang_day3@graphs)) wang_day3[[graph_name]] <- NULL
# Custom L-R from CellChat + CellCall (Wang Day 3) — same logic as Lee Day 1 / Lee Day 3; resource = custom + external_resource when pairs load.
wang_day3_liana_resource <- "MouseConsensus"
wang_day3_liana_external <- NULL
custom_lr_list_wang_d3 <- list()
path_cc_wang_d3 <- CELLCHAT_WANG_DAY3_RDS
path_ccall_wang_d3 <- CELLCALL_WANG_DAY3_RDS
cc_wang_d3_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_cc_wang_d3) && nzchar(trimws(path_cc_wang_d3)) && file.exists(path_cc_wang_d3) && requireNamespace("CellChat", quietly = TRUE)
cc_wang_d3 <- NULL
if (cc_wang_d3_ok) cc_wang_d3 <- readRDS(path_cc_wang_d3)
comm_wang_d3 <- NULL
if (!is.null(cc_wang_d3) && inherits(cc_wang_d3, "CellChat")) comm_wang_d3 <- CellChat::subsetCommunication(cc_wang_d3)
has_comm_wang_d3 <- !is.null(comm_wang_d3) && nrow(comm_wang_d3) > 0 && "ligand" %in% colnames(comm_wang_d3) && "receptor" %in% colnames(comm_wang_d3)
if (has_comm_wang_d3) {
  cc_lr_wang_d3 <- data.frame(source_genesymbol = comm_wang_d3$ligand, target_genesymbol = comm_wang_d3$receptor, stringsAsFactors = FALSE)
  cc_lr_wang_d3 <- cc_lr_wang_d3[!duplicated(cc_lr_wang_d3), ]
  custom_lr_list_wang_d3[["CellChat"]] <- cc_lr_wang_d3
}
ccall_wang_d3_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_ccall_wang_d3) && nzchar(trimws(path_ccall_wang_d3)) && file.exists(path_ccall_wang_d3)
ccall_wang_d3 <- NULL
if (ccall_wang_d3_ok) ccall_wang_d3 <- readRDS(path_ccall_wang_d3)
has_ccall_lr_wang_d3 <- !is.null(ccall_wang_d3) && !is.null(ccall_wang_d3@data$expr_l_r_log2_scale)
if (has_ccall_lr_wang_d3) {
  lr_rownames_wang_d3 <- rownames(ccall_wang_d3@data$expr_l_r_log2_scale)
  ccall_lig_wang_d3 <- sub("-.*", "", lr_rownames_wang_d3)
  ccall_rec_wang_d3 <- sub("^[^-]+-", "", lr_rownames_wang_d3)
  ccall_rec_wang_d3 <- gsub("-", "_", ccall_rec_wang_d3)
  ccall_lr_wang_d3 <- data.frame(source_genesymbol = ccall_lig_wang_d3, target_genesymbol = ccall_rec_wang_d3, stringsAsFactors = FALSE)
  ccall_lr_wang_d3 <- ccall_lr_wang_d3[nzchar(ccall_lr_wang_d3$target_genesymbol), ]
  ccall_lr_wang_d3 <- unique(ccall_lr_wang_d3)
  custom_lr_list_wang_d3[["CellCall"]] <- ccall_lr_wang_d3
}
has_custom_lr_wang_d3 <- length(custom_lr_list_wang_d3) > 0
if (has_custom_lr_wang_d3) {
  custom_lr_wang_d3 <- do.call(rbind, custom_lr_list_wang_d3)
  custom_lr_wang_d3 <- custom_lr_wang_d3[!duplicated(custom_lr_wang_d3[, c("source_genesymbol", "target_genesymbol")]), ]
  wang_day3_liana_external <- data.frame(source_genesymbol = custom_lr_wang_d3$source_genesymbol, target_genesymbol = custom_lr_wang_d3$target_genesymbol, stringsAsFactors = FALSE)
  wang_day3_liana_resource <- "custom"
  print(paste("Loaded", nrow(wang_day3_liana_external), "custom L-R pairs from CellChat + CellCall (Wang Day 3)"))
}
if (isTRUE(LIANA_USE_CUSTOM_LR) && !has_custom_lr_wang_d3) print("Wang Day 3 LIANA: LIANA_USE_CUSTOM_LR is TRUE but no CellChat/CellCall pairs loaded — using MouseConsensus only (check RDS paths and CellChat package).")
wang_day3_counts <- tryCatch(Seurat::GetAssayData(wang_day3, slot = "counts"), error = function(e) Seurat::GetAssayData(wang_day3, layer = "counts"))
wang_day3_logcounts <- tryCatch(Seurat::GetAssayData(wang_day3, slot = "data"), error = function(e) Seurat::GetAssayData(wang_day3, layer = "data"))
wang_day3_sce <- SingleCellExperiment(assays = list(counts = wang_day3_counts, logcounts = wang_day3_logcounts))
wang_day3_sce$liana_label <- wang_day3$liana_label
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_sce$liana_label
liana_args_wang_d3 <- list(sce = wang_day3_sce, method = LIANA_METHODS, resource = wang_day3_liana_resource, idents_col = "liana_label", expr_prop = 0.05, verbose = TRUE, min_cells = LIANA_MIN_CELLS, base = exp(1))
if (!is.null(wang_day3_liana_external)) liana_args_wang_d3$external_resource <- wang_day3_liana_external
wang_day3_liana_result <- do.call(liana::liana_wrap, liana_args_wang_d3)
wang_day3_liana_result_df <- liana::liana_aggregate(wang_day3_liana_result)

# Neutrophil targets only: signals RECEIVED by NeutrophilArg1pos and NeutrophilArg1neg (downstream plots/CSVs still emphasize Arg1pos unless noted)
wang_day3_liana_neutrophil <- dplyr::filter(wang_day3_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))
wang_day3_arg1pos_received <- dplyr::filter(wang_day3_liana_neutrophil, target == "NeutrophilArg1pos")
wang_day3_arg1pos_received <- dplyr::arrange(wang_day3_arg1pos_received, aggregate_rank)
wang_day3_arg1pos_top_signals <- head(wang_day3_arg1pos_received, 50)
wang_day3_arg1neg_received <- dplyr::filter(wang_day3_liana_neutrophil, target == "NeutrophilArg1neg")
wang_day3_arg1neg_received <- dplyr::arrange(wang_day3_arg1neg_received, aggregate_rank)
wang_day3_arg1neg_top_signals <- head(wang_day3_arg1neg_received, 50)
topN <- 50
wang_day3_arg1pos_top <- head(wang_day3_arg1pos_received, topN)
wang_day3_arg1neg_top <- head(wang_day3_arg1neg_received, topN)
pos_pairs <- paste0(wang_day3_arg1pos_top$ligand.complex, "_", wang_day3_arg1pos_top$receptor.complex)
neg_pairs_all <- paste0(wang_day3_arg1neg_received$ligand.complex, "_", wang_day3_arg1neg_received$receptor.complex)
neg_rank_lookup <- setNames(wang_day3_arg1neg_received$aggregate_rank, neg_pairs_all)
pos_ranks <- wang_day3_arg1pos_top$aggregate_rank
neg_ranks_matched <- neg_rank_lookup[pos_pairs]
neg_ranks_matched[is.na(neg_ranks_matched)] <- 1.0
spec_index <- (neg_ranks_matched - pos_ranks) / (neg_ranks_matched + pos_ranks + 1e-10)
wang_day3_arg1pos_top$specificity_index <- spec_index
wang_day3_arg1pos_specific <- dplyr::filter(wang_day3_arg1pos_top, specificity_index > 0.3 | !(pos_pairs %in% neg_pairs_all))
pos_pairs_all <- paste0(wang_day3_arg1pos_received$ligand.complex, "_", wang_day3_arg1pos_received$receptor.complex)
pos_rank_lookup_for_neg <- setNames(wang_day3_arg1pos_received$aggregate_rank, pos_pairs_all)
neg_lr_keys <- paste0(wang_day3_arg1neg_top$ligand.complex, "_", wang_day3_arg1neg_top$receptor.complex)
pos_ranks_for_neg <- pos_rank_lookup_for_neg[neg_lr_keys]
pos_ranks_for_neg[is.na(pos_ranks_for_neg)] <- 1.0
spec_index_neg <- (pos_ranks_for_neg - wang_day3_arg1neg_top$aggregate_rank) / (pos_ranks_for_neg + wang_day3_arg1neg_top$aggregate_rank + 1e-10)
wang_day3_arg1neg_top$specificity_index <- spec_index_neg
wang_day3_arg1neg_specific <- wang_day3_arg1neg_top[wang_day3_arg1neg_top$specificity_index > 0.3 | !(neg_lr_keys %in% pos_pairs_all), , drop = FALSE]
wang_day3_liana_arg1_pos <- wang_day3_arg1pos_received
wang_day3_liana_arg1_neg <- wang_day3_arg1neg_received
wang_day3_liana_consensus <- head(wang_day3_arg1pos_received, 30)
wang_day3_liana_consensus_arg1neg <- head(wang_day3_arg1neg_received, 30)
wang_day3_liana_arg1pos_topranked <- head(wang_day3_arg1pos_received, 20)
wang_day3_liana_arg1neg_topranked <- head(wang_day3_arg1neg_received, 20)

print(paste("Arg1pos received signals (all):", nrow(wang_day3_arg1pos_received)))
print(paste("Arg1neg received signals (all):", nrow(wang_day3_arg1neg_received)))
print(paste("Arg1pos top signals (top 50):", nrow(wang_day3_arg1pos_top_signals)))
print(paste("Arg1neg top signals (top 50):", nrow(wang_day3_arg1neg_top_signals)))
print(paste("Arg1pos-specific signals (specificity > 0.3 or unique to Arg1pos):", nrow(wang_day3_arg1pos_specific)))
print(paste("Arg1neg-specific signals (specificity > 0.3 or unique to Arg1neg):", nrow(wang_day3_arg1neg_specific)))

write.csv(wang_day3_liana_result_df, file.path(OUTPUT_DIR, "WangDay3_LIANA_AllResults.csv"), row.names = FALSE)
write.csv(wang_day3_arg1pos_received, file.path(OUTPUT_DIR, "WangDay3_Arg1pos_ReceivedSignals.csv"), row.names = FALSE)
write.csv(wang_day3_arg1neg_received, file.path(OUTPUT_DIR, "WangDay3_Arg1neg_ReceivedSignals.csv"), row.names = FALSE)
write.csv(wang_day3_arg1pos_top_signals, file.path(OUTPUT_DIR, "WangDay3_LIANA_Arg1pos_Top50Signals.csv"), row.names = FALSE)
write.csv(wang_day3_arg1neg_top_signals, file.path(OUTPUT_DIR, "WangDay3_LIANA_Arg1neg_Top50Signals.csv"), row.names = FALSE)
write.csv(wang_day3_arg1pos_specific, file.path(OUTPUT_DIR, "WangDay3_LIANA_Arg1pos_Specific.csv"), row.names = FALSE)
write.csv(wang_day3_arg1neg_specific, file.path(OUTPUT_DIR, "WangDay3_LIANA_Arg1neg_Specific.csv"), row.names = FALSE)
# NeutrophilConsensus.csv = top 30 by aggregate_rank for Arg1pos only (legacy filename; Arg1neg counterpart = NeutrophilConsensus_Arg1neg.csv)
write.csv(wang_day3_liana_consensus, file.path(OUTPUT_DIR, "WangDay3_LIANA_NeutrophilConsensus.csv"), row.names = FALSE)
write.csv(wang_day3_liana_consensus_arg1neg, file.path(OUTPUT_DIR, "WangDay3_LIANA_NeutrophilConsensus_Arg1neg.csv"), row.names = FALSE)
write.csv(wang_day3_liana_arg1pos_topranked, file.path(OUTPUT_DIR, "WangDay3_LIANA_Arg1pos_TopRanked.csv"), row.names = FALSE)
write.csv(wang_day3_liana_arg1neg_topranked, file.path(OUTPUT_DIR, "WangDay3_LIANA_Arg1neg_TopRanked.csv"), row.names = FALSE)
saveRDS(list(liana_result = wang_day3_liana_result, liana_aggregated = wang_day3_liana_result_df, arg1pos_received = wang_day3_arg1pos_received, arg1neg_received = wang_day3_arg1neg_received, arg1pos_top_signals = wang_day3_arg1pos_top_signals, arg1neg_top_signals = wang_day3_arg1neg_top_signals, arg1pos_specific = wang_day3_arg1pos_specific, arg1neg_specific = wang_day3_arg1neg_specific, neutrophil_consensus = wang_day3_liana_consensus, neutrophil_consensus_arg1neg = wang_day3_liana_consensus_arg1neg, arg1pos_topranked = wang_day3_liana_arg1pos_topranked, arg1neg_topranked = wang_day3_liana_arg1neg_topranked), file.path(OUTPUT_DIR, "WangDay3_LIANA_Results.rds"))
# Uncomment block below to load and skip re-running Wang Day 3 LIANA:
# wang_day3_liana_loaded <- readRDS(file.path(OUTPUT_DIR, "WangDay3_LIANA_Results.rds"))
# wang_day3_liana_result <- wang_day3_liana_loaded$liana_result
# wang_day3_liana_result_df <- wang_day3_liana_loaded$liana_aggregated
# wang_day3_arg1pos_received <- wang_day3_liana_loaded$arg1pos_received
# wang_day3_arg1neg_received <- wang_day3_liana_loaded$arg1neg_received
# wang_day3_arg1pos_top_signals <- wang_day3_liana_loaded$arg1pos_top_signals
# wang_day3_arg1neg_top_signals <- wang_day3_liana_loaded$arg1neg_top_signals
# wang_day3_arg1pos_specific <- wang_day3_liana_loaded$arg1pos_specific
# wang_day3_arg1neg_specific <- wang_day3_liana_loaded$arg1neg_specific
# wang_day3_liana_consensus <- wang_day3_liana_loaded$neutrophil_consensus
# wang_day3_liana_consensus_arg1neg <- wang_day3_liana_loaded$neutrophil_consensus_arg1neg
# wang_day3_liana_arg1pos_topranked <- wang_day3_liana_loaded$arg1pos_topranked
# wang_day3_liana_arg1neg_topranked <- wang_day3_liana_loaded$arg1neg_topranked

# -------- Wang Day 3: LIANA visualizations -> top LIANA receptors VlnPlot (same flow as Lee Day 1) --------
if (nrow(wang_day3_liana_result_df) > 0) {
  liana_network_wang_d3_both <- wang_day3_liana_result_df |>
    dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |>
    dplyr::filter(!(source %in% NEUTROPHIL_STATES)) |>
    dplyr::group_by(target) |>
    dplyr::arrange(aggregate_rank) |>
    dplyr::slice_head(n = 20) |>
    dplyr::ungroup()
  liana_network_wang_d3_both$interaction_label <- paste0(liana_network_wang_d3_both$ligand.complex, "\u2013(", liana_network_wang_d3_both$receptor.complex, ")")
  p_03e <- ggplot(liana_network_wang_d3_both, aes(x = target, y = interaction_label, size = -log10(aggregate_rank + 1e-10), color = source)) + geom_point(alpha = 0.8) + scale_x_discrete(limits = NEUTROPHIL_STATES, drop = FALSE) + theme_minimal() + labs(title = "Top Signals Received: Arg1+ vs Arg1- (Wang Day 3)", subtitle = "Y = Ligand\u2013(Receptor); X = receiver; color = sender (neutrophil\u2192neutrophil excluded from data)", x = "Receiver (neutrophil state)", y = "Interaction", color = "Sender cell type", size = "Consensus support\n(-log10 aggregate rank)") + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), plot.margin = margin(10, 80, 10, 10)) + guides(color = guide_legend(override.aes = list(size = 3)))
  vals_nlr <- numeric(0)
  if (exists("liana_network_lee_d1_both") && is.data.frame(liana_network_lee_d1_both) && nrow(liana_network_lee_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_lee_d3_both") && is.data.frame(liana_network_lee_d3_both) && nrow(liana_network_lee_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_wang_d3_both") && is.data.frame(liana_network_wang_d3_both) && nrow(liana_network_wang_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_wang_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d1_both") && is.data.frame(liana_network_qin_d1_both) && nrow(liana_network_qin_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d3_both") && is.data.frame(liana_network_qin_d3_both) && nrow(liana_network_qin_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d3_both$aggregate_rank + 1e-10))
  vals_nlr <- vals_nlr[is.finite(vals_nlr)]
  UNIFY_LIANA_NEGLOG10 <- if (length(vals_nlr) > 0) range(vals_nlr) else c(0, 10)
  if (length(UNIFY_LIANA_NEGLOG10) != 2 || !all(is.finite(UNIFY_LIANA_NEGLOG10))) UNIFY_LIANA_NEGLOG10 <- c(0, 10)
  if (UNIFY_LIANA_NEGLOG10[1] == UNIFY_LIANA_NEGLOG10[2]) UNIFY_LIANA_NEGLOG10[2] <- UNIFY_LIANA_NEGLOG10[1] + 1e-6
  print(p_03e + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = (if (exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE else c(2, 8))))
  wang_day3_external_to_neutrophils <- wang_day3_liana_result_df |> dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |> dplyr::filter(!(source %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))) |> dplyr::arrange(aggregate_rank)
  liana_top_wang_d3_external <- dplyr::slice_head(wang_day3_external_to_neutrophils, n = 20)
  liana_top_wang_d3_external$lr_label <- paste0(liana_top_wang_d3_external$source, " -> ", liana_top_wang_d3_external$target, "  ", liana_top_wang_d3_external$ligand.complex, "\u2013(", liana_top_wang_d3_external$receptor.complex, ")")
  p_03f0 <- ggplot(liana_top_wang_d3_external, aes(x = reorder(lr_label, aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 20 L-R Pairs - LIANA (Wang Day 3)\nExternal signals to Arg1+ / Arg1- neutrophils", x = "Source -> Target  Ligand\u2013(Receptor)", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME
  vals_ext <- numeric(0)
  if (exists("lee_day1_external_to_neutrophils") && is.data.frame(lee_day1_external_to_neutrophils) && nrow(lee_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day1_external_to_neutrophils$aggregate_rank)
  if (exists("lee_day3_external_to_neutrophils") && is.data.frame(lee_day3_external_to_neutrophils) && nrow(lee_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day3_external_to_neutrophils$aggregate_rank)
  if (exists("wang_day3_external_to_neutrophils") && is.data.frame(wang_day3_external_to_neutrophils) && nrow(wang_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(wang_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, wang_day3_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day1_external_to_neutrophils") && is.data.frame(qin_day1_external_to_neutrophils) && nrow(qin_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day1_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day3_external_to_neutrophils") && is.data.frame(qin_day3_external_to_neutrophils) && nrow(qin_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day3_external_to_neutrophils$aggregate_rank)
  if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0 && "aggregate_rank" %in% names(liana_top_lee_d1)) vals_ext <- c(vals_ext, liana_top_lee_d1$aggregate_rank)
  if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_lee_d3_external)) vals_ext <- c(vals_ext, liana_top_lee_d3_external$aggregate_rank)
  if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_wang_d3_external)) vals_ext <- c(vals_ext, liana_top_wang_d3_external$aggregate_rank)
  if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d1_external)) vals_ext <- c(vals_ext, liana_top_qin_d1_external$aggregate_rank)
  if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d3_external)) vals_ext <- c(vals_ext, liana_top_qin_d3_external$aggregate_rank)
  vals_ext <- suppressWarnings(as.numeric(vals_ext))
  vals_ext <- vals_ext[is.finite(vals_ext)]
  plot_ext_ar <- suppressWarnings(as.numeric(liana_top_wang_d3_external[["aggregate_rank"]]))
  plot_ext_ar <- plot_ext_ar[is.finite(plot_ext_ar)]
  UNIFY_EXT_AR <- range(c(vals_ext, plot_ext_ar), na.rm = TRUE)
  if (length(plot_ext_ar) == 0 && length(vals_ext) == 0) UNIFY_EXT_AR <- c(0, 1)
  if (!all(is.finite(UNIFY_EXT_AR))) UNIFY_EXT_AR <- c(0, 1)
  if (UNIFY_EXT_AR[1] == UNIFY_EXT_AR[2]) UNIFY_EXT_AR[2] <- UNIFY_EXT_AR[1] + max(abs(UNIFY_EXT_AR[1]) * 1e-6, 1e-12)
  print(p_03f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
}
liana_top_wang_d3 <- dplyr::slice_head(wang_day3_liana_consensus, n = 15)
p_03f <- ggplot(liana_top_wang_d3, aes(x = reorder(paste0(source, " -> ", target), aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 15 L-R Pairs - LIANA (Wang Day 3)", x = "Source -> Target", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME
vals_cons <- numeric(0)
if (exists("liana_top_lee_d3") && is.data.frame(liana_top_lee_d3) && nrow(liana_top_lee_d3) > 0) vals_cons <- c(vals_cons, liana_top_lee_d3$aggregate_rank)
if (exists("liana_top_wang_d3") && is.data.frame(liana_top_wang_d3) && nrow(liana_top_wang_d3) > 0) vals_cons <- c(vals_cons, liana_top_wang_d3$aggregate_rank)
if (exists("liana_top_qin_d1") && is.data.frame(liana_top_qin_d1) && nrow(liana_top_qin_d1) > 0) vals_cons <- c(vals_cons, liana_top_qin_d1$aggregate_rank)
if (exists("liana_top_qin_d3") && is.data.frame(liana_top_qin_d3) && nrow(liana_top_qin_d3) > 0) vals_cons <- c(vals_cons, liana_top_qin_d3$aggregate_rank)
vals_cons <- vals_cons[is.finite(vals_cons)]
UNIFY_CONS_AR <- if (length(vals_cons) > 0) range(vals_cons) else c(0, 1)
if (UNIFY_CONS_AR[1] == UNIFY_CONS_AR[2]) UNIFY_CONS_AR[2] <- UNIFY_CONS_AR[1] + 1e-6
print(p_03f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
unique_sources_wang_d3 <- unique(wang_day3_liana_result_df$source)
neutrophil_targets_wang <- NEUTROPHIL_STATES[NEUTROPHIL_STATES %in% unique(wang_day3_liana_result_df$target)]
dot_sources_wang_d3 <- setdiff(unique_sources_wang_d3, neutrophil_targets_wang)
if (length(dot_sources_wang_d3) == 0) dot_sources_wang_d3 <- unique_sources_wang_d3
p_03f2 <- NULL
if (length(neutrophil_targets_wang) > 0 && length(dot_sources_wang_d3) > 0) p_03f2 <- liana::liana_dotplot(wang_day3_liana_result_df, source_groups = dot_sources_wang_d3, target_groups = neutrophil_targets_wang, ntop = 20, size_range = LIANA_DOTPLOT_SIZE_RANGE)
if (!is.null(p_03f2)) { p_03f2 <- p_03f2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)); print(p_03f2) }
interaction_matrix_wang_d3 <- tryCatch(tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(wang_day3_liana_consensus, source, target), interaction_count = dplyr::n(), .groups = "drop"), names_from = target, values_from = interaction_count, values_fill = 0), error = function(e) data.frame())
interaction_matrix_wang_d3_mat <- tryCatch({ cols_num <- setdiff(colnames(interaction_matrix_wang_d3), "source"); mat <- as.matrix(interaction_matrix_wang_d3[, cols_num, drop = FALSE]); rownames(mat) <- interaction_matrix_wang_d3$source; mat }, error = function(e) matrix(0, nrow = 0, ncol = 0))
p_03g <- tryCatch(pheatmap::pheatmap(interaction_matrix_wang_d3_mat, color = colorRampPalette(c("white", "yellow", "orange", "red"))(100), main = "Cell Type Interaction Frequency - Wang Day 3", display_numbers = TRUE), error = function(e) NULL)
if (!is.null(p_03g)) print(p_03g)
if (!is.null(p_03g)) { grid::grid.newpage(); grid::grid.draw(p_03g$gtable) }
receptor_freq_wang_d3 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(wang_day3_liana_consensus, receptor.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
p_03h <- tryCatch(ggplot(receptor_freq_wang_d3, aes(x = reorder(receptor.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Receptors - Wang Day 3", x = "Receptor", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_03h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
ligand_freq_wang_d3 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(wang_day3_liana_consensus, ligand.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
p_03i <- tryCatch(ggplot(ligand_freq_wang_d3, aes(x = reorder(ligand.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Ligands - Wang Day 3", x = "Ligand", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_03i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
liana_trunc_wang_d3 <- dplyr::filter(wang_day3_liana_result_df, aggregate_rank <= 0.01)
if (nrow(liana_trunc_wang_d3) > 0) {
  liana_trunc_wang_d3$source <- as.character(liana_trunc_wang_d3$source)
  liana_trunc_wang_d3$target <- as.character(liana_trunc_wang_d3$target)
  p_heat_freq_wang <- tryCatch(liana::heat_freq(liana_trunc_wang_d3), error = function(e) { message("heat_freq Wang D3: ", conditionMessage(e)); NULL })
  if (!is.null(p_heat_freq_wang)) print(p_heat_freq_wang)
  unique_sources_chord_wang <- unique(liana_trunc_wang_d3$source)
  unique_targets_chord_wang <- unique(liana_trunc_wang_d3$target)
  grDevices::png(file.path(OUTPUT_DIR, "WangDay3_LIANA_ChordFreq.png"), width = 1400, height = 1400, res = 150)
  tryCatch(liana::chord_freq(liana_trunc_wang_d3, source_groups = unique_sources_chord_wang, target_groups = unique_targets_chord_wang), error = function(e) message("chord_freq Wang D3 (PNG): ", conditionMessage(e)))
  grDevices::dev.off()
  p_chord_freq_wang <- tryCatch(liana::chord_freq(liana_trunc_wang_d3, source_groups = unique_sources_chord_wang, target_groups = unique_targets_chord_wang), error = function(e) { message("chord_freq Wang D3: ", conditionMessage(e)); NULL })
  if (!is.null(p_chord_freq_wang)) print(p_chord_freq_wang)
  liana_mat_wang_d3 <- as.matrix(table(liana_trunc_wang_d3$source, liana_trunc_wang_d3$target))
  p_liana_heatmap_wang_d3 <- NULL
  if (nrow(liana_mat_wang_d3) > 0 && ncol(liana_mat_wang_d3) > 0) p_liana_heatmap_wang_d3 <- tryCatch(liana::liana_heatmap(liana_mat_wang_d3), error = function(e) { message("liana_heatmap Wang D3: ", conditionMessage(e)); NULL })
  if (!is.null(p_liana_heatmap_wang_d3)) ComplexHeatmap::draw(p_liana_heatmap_wang_d3)
}
source_importance_wang_d3 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(wang_day3_liana_consensus, source), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), importance_score = dplyr::n() * (1 - mean(aggregate_rank))), desc(importance_score)), error = function(e) data.frame())
p_03k <- tryCatch(ggplot(source_importance_wang_d3, aes(x = reorder(source, importance_score), y = importance_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkred") + labs(title = "Source Cell Importance - Wang Day 3", x = "Cell Type", y = "Importance Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_iy <- numeric(0)
vals_imr <- numeric(0)
imp_tabnames <- c("source_importance_lee_d1", "source_importance_lee_d3", "source_importance_wang_d3", "source_importance_qin_d1", "source_importance_qin_d3")
for (fn in imp_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("importance_score" %in% names(d)) vals_iy <- c(vals_iy, d$importance_score)
  if ("mean_rank" %in% names(d)) vals_imr <- c(vals_imr, d$mean_rank)
}
vals_iy <- vals_iy[is.finite(vals_iy)]
vals_imr <- vals_imr[is.finite(vals_imr)]
UNIFY_SRC_IMP_Y <- if (length(vals_iy) > 0) range(vals_iy) else c(0, 1)
UNIFY_SRC_IMP_MR <- if (length(vals_imr) > 0) range(vals_imr) else c(0, 1)
if (UNIFY_SRC_IMP_Y[1] == UNIFY_SRC_IMP_Y[2]) UNIFY_SRC_IMP_Y[2] <- UNIFY_SRC_IMP_Y[1] + 1e-6
if (UNIFY_SRC_IMP_MR[1] == UNIFY_SRC_IMP_MR[2]) UNIFY_SRC_IMP_MR[2] <- UNIFY_SRC_IMP_MR[1] + 1e-6
print(p_03k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
p_03n <- tryCatch(ggplot(head(wang_day3_liana_arg1pos_topranked, 15), aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), aggregate_rank), y = -aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top Arg1pos Interactions - Wang Day 3", x = "L-R Pair", y = "Rank Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_a1 <- numeric(0)
if (exists("lee_day1_liana_arg1pos_topranked") && is.data.frame(lee_day1_liana_arg1pos_topranked) && nrow(lee_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("lee_day3_liana_arg1pos_topranked") && is.data.frame(lee_day3_liana_arg1pos_topranked) && nrow(lee_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("wang_day3_liana_arg1pos_topranked") && is.data.frame(wang_day3_liana_arg1pos_topranked) && nrow(wang_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(wang_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day1_liana_arg1pos_topranked") && is.data.frame(qin_day1_liana_arg1pos_topranked) && nrow(qin_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day3_liana_arg1pos_topranked") && is.data.frame(qin_day3_liana_arg1pos_topranked) && nrow(qin_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day3_liana_arg1pos_topranked$aggregate_rank, 15))
vals_a1 <- vals_a1[is.finite(vals_a1)]
UNIFY_ARG1_AR <- if (length(vals_a1) > 0) range(vals_a1) else c(0, 1)
if (UNIFY_ARG1_AR[1] == UNIFY_ARG1_AR[2]) UNIFY_ARG1_AR[2] <- UNIFY_ARG1_AR[1] + 1e-6
UNIFY_ARG1_NEGY <- c(-UNIFY_ARG1_AR[2], -UNIFY_ARG1_AR[1])
print(p_03n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
wang_day3_arg1pos_specific_plot <- head(wang_day3_arg1pos_specific, 15)
p_03o2 <- tryCatch(ggplot(wang_day3_arg1pos_specific_plot, aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), specificity_index), y = specificity_index, fill = source)) + geom_bar(stat = "identity") + coord_flip() + labs(title = "Arg1+-Specific Signals (Wang Day 3)", x = "L-R Pair", y = "Specificity Index") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
print(p_03o2)
target_specificity_wang_d3 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(dplyr::filter(wang_day3_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")), target), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), specificity_score = dplyr::n() * (1 - mean(aggregate_rank)), .groups = "drop"), desc(specificity_score)), error = function(e) data.frame())
p_03l <- tryCatch(ggplot(target_specificity_wang_d3, aes(x = reorder(target, specificity_score), y = specificity_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkgreen") + labs(title = "Target Cell Specificity - Wang Day 3", x = "Cell Type", y = "Specificity Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_sy <- numeric(0)
vals_smr <- numeric(0)
spec_tabnames <- c("target_specificity_lee_d1", "target_specificity_lee_d3", "target_specificity_wang_d3", "target_specificity_qin_d1", "target_specificity_qin_d3")
for (fn in spec_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("specificity_score" %in% names(d)) vals_sy <- c(vals_sy, d$specificity_score)
  if ("mean_rank" %in% names(d)) vals_smr <- c(vals_smr, d$mean_rank)
}
vals_sy <- vals_sy[is.finite(vals_sy)]
vals_smr <- vals_smr[is.finite(vals_smr)]
UNIFY_TGT_SPEC_Y <- if (length(vals_sy) > 0) range(vals_sy) else c(0, 1)
UNIFY_TGT_SPEC_MR <- if (length(vals_smr) > 0) range(vals_smr) else c(0, 1)
if (UNIFY_TGT_SPEC_Y[1] == UNIFY_TGT_SPEC_Y[2]) UNIFY_TGT_SPEC_Y[2] <- UNIFY_TGT_SPEC_Y[1] + 1e-6
if (UNIFY_TGT_SPEC_MR[1] == UNIFY_TGT_SPEC_MR[2]) UNIFY_TGT_SPEC_MR[2] <- UNIFY_TGT_SPEC_MR[1] + 1e-6
print(p_03l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
p_03m_rank <- tryCatch(ggplot(wang_day3_liana_consensus, aes(x = aggregate_rank)) + geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) + labs(title = "Distribution of L-R Pair Aggregate Ranks - Wang Day 3", x = "Aggregate Rank Score", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_hist <- numeric(0)
if (exists("lee_day1_liana_consensus") && is.data.frame(lee_day1_liana_consensus) && nrow(lee_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day1_liana_consensus)) vals_hist <- c(vals_hist, lee_day1_liana_consensus$aggregate_rank)
if (exists("lee_day3_liana_consensus") && is.data.frame(lee_day3_liana_consensus) && nrow(lee_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day3_liana_consensus)) vals_hist <- c(vals_hist, lee_day3_liana_consensus$aggregate_rank)
if (exists("wang_day3_liana_consensus") && is.data.frame(wang_day3_liana_consensus) && nrow(wang_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(wang_day3_liana_consensus)) vals_hist <- c(vals_hist, wang_day3_liana_consensus$aggregate_rank)
if (exists("qin_day1_liana_consensus") && is.data.frame(qin_day1_liana_consensus) && nrow(qin_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day1_liana_consensus)) vals_hist <- c(vals_hist, qin_day1_liana_consensus$aggregate_rank)
if (exists("qin_day3_liana_consensus") && is.data.frame(qin_day3_liana_consensus) && nrow(qin_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day3_liana_consensus)) vals_hist <- c(vals_hist, qin_day3_liana_consensus$aggregate_rank)
vals_hist <- vals_hist[is.finite(vals_hist)]
UNIFY_HIST_AR <- if (length(vals_hist) > 0) range(vals_hist) else c(0, 1)
if (UNIFY_HIST_AR[1] == UNIFY_HIST_AR[2]) UNIFY_HIST_AR[2] <- UNIFY_HIST_AR[1] + 1e-6
print(p_03m_rank + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
method_cols_wang_d3 <- colnames(wang_day3_liana_result_df)[grep("pval_", colnames(wang_day3_liana_result_df))]
consensus_wang_d3_step1 <- dplyr::slice_head(wang_day3_liana_consensus, n = 15)
consensus_wang_d3 <- NULL
consensus_wang_d3_has_cols <- length(method_cols_wang_d3) > 0 && nrow(consensus_wang_d3_step1) > 0
if (consensus_wang_d3_has_cols) consensus_wang_d3 <- as.data.frame(dplyr::select(consensus_wang_d3_step1, dplyr::all_of(method_cols_wang_d3)))
if (consensus_wang_d3_has_cols) consensus_wang_d3_rownames <- paste0(consensus_wang_d3_step1$source, " | ", consensus_wang_d3_step1$ligand.complex, " -> ", consensus_wang_d3_step1$receptor.complex)
if (consensus_wang_d3_has_cols) rownames(consensus_wang_d3) <- make.unique(as.character(consensus_wang_d3_rownames))
consensus_wang_d3_ready <- !is.null(consensus_wang_d3) && nrow(consensus_wang_d3) > 0
if (consensus_wang_d3_ready) { p_03o <- pheatmap::pheatmap(consensus_wang_d3, color = colorRampPalette(c("red", "white", "blue"))(100), main = "Method Consensus (p-values) - Wang Day 3"); print(p_03o) }


wang_day3_neut_cells_vln <- colnames(wang_day3)[c(wang_day3_pos, wang_day3_neg)]
top_receptors_wang_d3 <- head(unique(wang_day3_arg1pos_received$receptor.complex), LIANA_TOP_RECEPTOR_VLN)
top_receptors_wang_d3_single <- top_receptors_wang_d3[!grepl("[_+]", top_receptors_wang_d3)]
top_receptors_wang_d3_in_data <- top_receptors_wang_d3_single[top_receptors_wang_d3_single %in% rownames(wang_day3)]
if (length(top_receptors_wang_d3_in_data) > 0) {
  wang_day3_neut_obj_vln <- subset(wang_day3, cells = wang_day3_neut_cells_vln)
  p_receptor_vln_wang_d3 <- Seurat::VlnPlot(wang_day3_neut_obj_vln, features = top_receptors_wang_d3_in_data, group.by = "arg1_status", pt.size = 0.1, ncol = min(3L, length(top_receptors_wang_d3_in_data)))
  print(p_receptor_vln_wang_d3)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "WangDay3_TopLIANA_Receptors_VlnPlot.png"), p_receptor_vln_wang_d3, width = 10, height = 6, dpi = 150)
}
print("✓ Wang Day 3 LIANA analysis complete")
# Checkpoint: Wang Day 3 LIANA only (RDS). For Qin Day 1 LIANA next, use Workspace_AfterWangDay3_LIANA.RData (full env) or run Section 3 so qin_day1* exist. save.image may warn depth(NULL) from a NULL slot in the workspace (often harmless).
save.image(file.path(OUTPUT_DIR, "Workspace_AfterWangDay3_LIANA.RData"))
saveRDS(list(OUTPUT_DIR = OUTPUT_DIR, NEUTROPHIL_BASE_LABEL = NEUTROPHIL_BASE_LABEL, NEUTROPHIL_STATES = NEUTROPHIL_STATES, wang_day3 = wang_day3, wang_day3_pruned_labels = wang_day3_pruned_labels, wang_day3_pos = wang_day3_pos, wang_day3_neg = wang_day3_neg, wang_day3_arg1_status = wang_day3_arg1_status, wang_day3_liana_result_df = wang_day3_liana_result_df, wang_day3_arg1pos_received = wang_day3_arg1pos_received), file.path(OUTPUT_DIR, "Checkpoint_AfterWangDay3_LIANA.rds"))
# checkpoint_wd3 <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterWangDay3_LIANA.rds")); list2env(checkpoint_wd3, envir = .GlobalEnv)
print("✓ Checkpoint saved: Workspace_AfterWangDay3_LIANA.RData, Checkpoint_AfterWangDay3_LIANA.rds")
# -------- Wang Day 3 PROGENy Pathway Analysis (Exploratory; downstream pathway support only) --------
# Target population (aligned with LIANA): NeutrophilArg1pos and NeutrophilArg1neg only. Expression subset to neutrophils; pathway activity compared Arg1pos vs Arg1neg.
# Note: wang_day3_arg1_status is already defined above and will be reused here
print("--- PROGENy Pathway Analysis: Wang Day 3 (Exploratory) ---")
# Secondary support only: PROGENy summarizes downstream pathway state in Arg1pos vs Arg1neg neutrophils.

# Get PROGENy pathway gene sets for mouse
# Use progeny package directly (avoids OmnipathR logging issues)
# Note: progeny::model_mouse_full is a dataframe with columns: gene, pathway, weight, p.value
progeny_model_mouse_wang_d3 <- progeny::model_mouse_full
colnames(progeny_model_mouse_wang_d3) <- tolower(colnames(progeny_model_mouse_wang_d3))
colnames(progeny_model_mouse_wang_d3)[colnames(progeny_model_mouse_wang_d3) == "p.value"] <- "p_value"
stopifnot(all(c("gene", "pathway", "weight") %in% colnames(progeny_model_mouse_wang_d3)))
progeny_network_wang_d3 <- data.frame(
  source = progeny_model_mouse_wang_d3$pathway,
  target = progeny_model_mouse_wang_d3$gene,
  weight = progeny_model_mouse_wang_d3$weight,
  stringsAsFactors = FALSE
)
progeny_network_wang_d3 <- progeny_network_wang_d3[progeny_network_wang_d3$weight != 0, ]
progeny_network_wang_d3_pathways_unique <- unique(progeny_network_wang_d3$source)
print(paste("PROGENy loaded successfully: nrow =", nrow(progeny_network_wang_d3), ", pathways =", length(progeny_network_wang_d3_pathways_unique)))

# Subset to neutrophils only for pathway analysis (Arg1+ and Arg1- neutrophils)
wang_day3_neut_cells <- colnames(wang_day3)[c(wang_day3_pos, wang_day3_neg)]
wang_day3_neut_expr_mat <- tryCatch(Seurat::GetAssayData(wang_day3, layer = "data")[, wang_day3_neut_cells, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
wang_day3_expr_mat <- wang_day3_neut_expr_mat

# Calculate pathway activity scores using weighted mean (WMEAN)
wang_day3_progeny_result <- tryCatch(decoupleR::run_wmean(mat = wang_day3_expr_mat, network = progeny_network_wang_d3, .source = "source", .target = "target", .mor = "weight", minsize = 5), error = function(e) data.frame())

# Convert result to wide format (pathways x cells). Use norm_wmean statistic only.
wang_day3_progeny_for_wide <- if (nrow(wang_day3_progeny_result) > 0 && "statistic" %in% colnames(wang_day3_progeny_result)) dplyr::filter(wang_day3_progeny_result, statistic == "norm_wmean") else wang_day3_progeny_result
if (nrow(wang_day3_progeny_for_wide) == 0 && nrow(wang_day3_progeny_result) > 0 && "statistic" %in% colnames(wang_day3_progeny_result)) wang_day3_progeny_for_wide <- dplyr::filter(wang_day3_progeny_result, statistic == "wmean")
wang_day3_pw_wide <- tryCatch(tidyr::pivot_wider(wang_day3_progeny_for_wide, names_from = "condition", values_from = "score", id_cols = "source"), error = function(e) data.frame())
wang_day3_pw_cols_num <- tryCatch(sapply(wang_day3_pw_wide[, -1, drop = FALSE], function(x) as.numeric(unlist(x))), error = function(e) matrix(0, nrow = 0, ncol = 0))
wang_day3_progeny_scores_mat <- tryCatch(as.matrix(wang_day3_pw_cols_num), error = function(e) matrix(0, nrow = 0, ncol = 0))
rownames(wang_day3_progeny_scores_mat) <- tryCatch(as.character(wang_day3_pw_wide$source), error = function(e) character(0))
colnames(wang_day3_progeny_scores_mat) <- tryCatch(colnames(wang_day3_pw_wide)[-1], error = function(e) character(0))

# Get cell names for Arg1pos and Arg1neg (neutrophils only)
wang_day3_arg1pos_cells <- colnames(wang_day3)[wang_day3_pos]
wang_day3_arg1neg_cells <- colnames(wang_day3)[wang_day3_neg]

# Extract pathway scores for each group
wang_day3_arg1pos_cells_in_mat <- tryCatch(wang_day3_arg1pos_cells[wang_day3_arg1pos_cells %in% colnames(wang_day3_progeny_scores_mat)], error = function(e) character(0))
wang_day3_arg1neg_cells_in_mat <- tryCatch(wang_day3_arg1neg_cells[wang_day3_arg1neg_cells %in% colnames(wang_day3_progeny_scores_mat)], error = function(e) character(0))

wang_day3_progeny_arg1pos_scores <- tryCatch(wang_day3_progeny_scores_mat[, wang_day3_arg1pos_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
wang_day3_progeny_arg1neg_scores <- tryCatch(wang_day3_progeny_scores_mat[, wang_day3_arg1neg_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))

# Calculate mean and median pathway scores per group (effect sizes)
wang_day3_pathway_means_arg1pos <- tryCatch(rowMeans(wang_day3_progeny_arg1pos_scores, na.rm = TRUE), error = function(e) numeric(0))
wang_day3_pathway_means_arg1neg <- tryCatch(rowMeans(wang_day3_progeny_arg1neg_scores, na.rm = TRUE), error = function(e) numeric(0))
wang_day3_pathway_medians_arg1pos <- tryCatch(apply(wang_day3_progeny_arg1pos_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
wang_day3_pathway_medians_arg1neg <- tryCatch(apply(wang_day3_progeny_arg1neg_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))

wang_day3_pathway_shift <- pmax(0, -pmin(wang_day3_pathway_means_arg1pos, wang_day3_pathway_means_arg1neg, na.rm = TRUE)) + 1e-10
wang_day3_pathway_log2fc <- log2((wang_day3_pathway_means_arg1pos + wang_day3_pathway_shift) / (wang_day3_pathway_means_arg1neg + wang_day3_pathway_shift))

wang_day3_pathway_comparison <- tryCatch(data.frame(
  pathway = names(wang_day3_pathway_means_arg1pos),
  Arg1pos_mean = wang_day3_pathway_means_arg1pos,
  Arg1neg_mean = wang_day3_pathway_means_arg1neg,
  Arg1pos_median = wang_day3_pathway_medians_arg1pos,
  Arg1neg_median = wang_day3_pathway_medians_arg1neg,
  mean_difference = wang_day3_pathway_means_arg1pos - wang_day3_pathway_means_arg1neg,
  median_difference = wang_day3_pathway_medians_arg1pos - wang_day3_pathway_medians_arg1neg,
  log2FC = wang_day3_pathway_log2fc,
  stringsAsFactors = FALSE
), error = function(e) data.frame())

# Exploratory pathway analysis: Use all PROGENy pathways (no pre-classification)
# Get all unique pathways from PROGENy network
wang_day3_all_pathways <- tryCatch(unique(progeny_network_wang_d3$source), error = function(e) character(0))
wang_day3_pathway_results <- wang_day3_pathway_comparison

# Report sample sizes
wang_day3_arg1pos_count <- tryCatch(sum(wang_day3_arg1_status == "Arg1pos"), error = function(e) 0)
wang_day3_arg1neg_count <- tryCatch(sum(wang_day3_arg1_status == "Arg1neg"), error = function(e) 0)
print(paste("Wang Day 3: Arg1pos cells =", wang_day3_arg1pos_count, ", Arg1neg cells =", wang_day3_arg1neg_count))

# Statistical test: Wilcoxon test for pathway differences (all pathways)
# With sample size adequacy checks, effect sizes, and confidence intervals
wang_day3_pathway_pvals <- numeric(length(wang_day3_all_pathways))
names(wang_day3_pathway_pvals) <- wang_day3_all_pathways
wang_day3_pathway_effect_sizes <- numeric(length(wang_day3_all_pathways))
names(wang_day3_pathway_effect_sizes) <- wang_day3_all_pathways
wang_day3_pathway_ci_lower <- numeric(length(wang_day3_all_pathways))
names(wang_day3_pathway_ci_lower) <- wang_day3_all_pathways
wang_day3_pathway_ci_upper <- numeric(length(wang_day3_all_pathways))
names(wang_day3_pathway_ci_upper) <- wang_day3_all_pathways
wang_day3_pathway_adequate_n <- logical(length(wang_day3_all_pathways))
names(wang_day3_pathway_adequate_n) <- wang_day3_all_pathways

for (pw_idx in seq_along(wang_day3_all_pathways)) {
  pw <- wang_day3_all_pathways[pw_idx]
  arg1pos_scores_pw <- tryCatch(as.numeric(wang_day3_progeny_arg1pos_scores[pw, ]), error = function(e) numeric(0))
  arg1neg_scores_pw <- tryCatch(as.numeric(wang_day3_progeny_arg1neg_scores[pw, ]), error = function(e) numeric(0))
  arg1pos_scores_pw <- arg1pos_scores_pw[!is.na(arg1pos_scores_pw)]
  arg1neg_scores_pw <- arg1neg_scores_pw[!is.na(arg1neg_scores_pw)]
  
  # Sample size adequacy check (Wilcoxon requires n ≥ 5 per group)
  n_arg1pos <- length(arg1pos_scores_pw)
  n_arg1neg <- length(arg1neg_scores_pw)
  wang_day3_pathway_adequate_n[pw] <- (n_arg1pos >= 5) & (n_arg1neg >= 5)
  
  # Calculate effect size (median difference)
  median_arg1pos <- tryCatch(stats::median(arg1pos_scores_pw, na.rm = TRUE), error = function(e) 0)
  median_arg1neg <- tryCatch(stats::median(arg1neg_scores_pw, na.rm = TRUE), error = function(e) 0)
  wang_day3_pathway_effect_sizes[pw] <- median_arg1pos - median_arg1neg
  
  # Wilcoxon test with confidence interval (linear: always attempt, handle errors)
  wang_day3_pathway_pvals[pw] <- 1.0
  wang_day3_pathway_ci_lower[pw] <- NA_real_
  wang_day3_pathway_ci_upper[pw] <- NA_real_
  wilcox_result <- tryCatch(stats::wilcox.test(arg1pos_scores_pw, arg1neg_scores_pw, conf.int = TRUE, conf.level = 0.95), error = function(e) NULL)
  wilcox_pvalue <- tryCatch(if (!is.null(wilcox_result)) wilcox_result$p.value else 1.0, error = function(e) 1.0)
  wilcox_pvalue_length <- tryCatch(length(wilcox_pvalue), error = function(e) 0)
  wilcox_pvalue_final <- tryCatch(if (wilcox_pvalue_length > 0) wilcox_pvalue[1] else 1.0, error = function(e) 1.0)
  wang_day3_pathway_pvals[pw] <- wilcox_pvalue_final
  wilcox_ci_lower <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[1] else NA_real_, error = function(e) NA_real_)
  wilcox_ci_upper <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[2] else NA_real_, error = function(e) NA_real_)
  wang_day3_pathway_ci_lower[pw] <- wilcox_ci_lower
  wang_day3_pathway_ci_upper[pw] <- wilcox_ci_upper
  wang_day3_pathway_warning_msg <- tryCatch(paste("Warning: Pathway", pw, "has insufficient sample size (Arg1pos n =", n_arg1pos, ", Arg1neg n =", n_arg1neg, "). Skipping statistical test."), error = function(e) "")
  wang_day3_pathway_warning_vector <- c("", wang_day3_pathway_warning_msg)
  wang_day3_pathway_warning_index <- tryCatch(as.numeric(!wang_day3_pathway_adequate_n[pw]) + 1, error = function(e) 1)
  print(wang_day3_pathway_warning_vector[wang_day3_pathway_warning_index])
}

# Multiple testing correction across ALL pathways (FDR)
wang_day3_pathway_pvals_adj <- p.adjust(wang_day3_pathway_pvals, method = "BH")

# Add statistical results to comparison dataframe (p-values, effect sizes, CIs, sample size adequacy)
wang_day3_pathway_results$p_value <- tryCatch(wang_day3_pathway_pvals[wang_day3_pathway_results$pathway], error = function(e) rep(1.0, nrow(wang_day3_pathway_results)))
wang_day3_pathway_results$p_adj <- tryCatch(wang_day3_pathway_pvals_adj[wang_day3_pathway_results$pathway], error = function(e) rep(1.0, nrow(wang_day3_pathway_results)))
wang_day3_pathway_results$effect_size_median_diff <- tryCatch(wang_day3_pathway_effect_sizes[wang_day3_pathway_results$pathway], error = function(e) rep(0.0, nrow(wang_day3_pathway_results)))
wang_day3_pathway_results$ci_lower_95 <- tryCatch(wang_day3_pathway_ci_lower[wang_day3_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(wang_day3_pathway_results)))
wang_day3_pathway_results$ci_upper_95 <- tryCatch(wang_day3_pathway_ci_upper[wang_day3_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(wang_day3_pathway_results)))
wang_day3_pathway_results$adequate_sample_size <- tryCatch(wang_day3_pathway_adequate_n[wang_day3_pathway_results$pathway], error = function(e) rep(FALSE, nrow(wang_day3_pathway_results)))
wang_day3_pathway_results$significant <- tryCatch((wang_day3_pathway_results$p_adj < 0.05) & wang_day3_pathway_results$adequate_sample_size, error = function(e) rep(FALSE, nrow(wang_day3_pathway_results)))

# Save PROGENy results
write.csv(wang_day3_pathway_comparison, file.path(OUTPUT_DIR, "WangDay3_PROGENy_PathwayComparison.csv"), row.names = FALSE)
write.csv(wang_day3_pathway_results, file.path(OUTPUT_DIR, "WangDay3_PROGENy_PathwayResults.csv"), row.names = FALSE)
saveRDS(list(pathway_comparison = wang_day3_pathway_comparison, pathway_results = wang_day3_pathway_results, pathway_means_arg1pos = wang_day3_pathway_means_arg1pos, pathway_means_arg1neg = wang_day3_pathway_means_arg1neg, progeny_network = progeny_network_wang_d3), file.path(OUTPUT_DIR, "WangDay3_PROGENy_Results.rds"))
# wang_day3_progeny_loaded <- readRDS(file.path(OUTPUT_DIR, "WangDay3_PROGENy_Results.rds")); wang_day3_pathway_comparison <- wang_day3_progeny_loaded$pathway_comparison; wang_day3_pathway_results <- wang_day3_progeny_loaded$pathway_results; wang_day3_pathway_means_arg1pos <- wang_day3_progeny_loaded$pathway_means_arg1pos; wang_day3_pathway_means_arg1neg <- wang_day3_progeny_loaded$pathway_means_arg1neg; progeny_network_wang_d3 <- wang_day3_progeny_loaded$progeny_network

# Ligand-to-pathway mapping: Use top-ranked interactions for BOTH Arg1pos and Arg1neg
# Note: This uses top-ranked, not necessarily "specific" interactions
wang_day3_identified_ligands_arg1pos <- tryCatch(unique(wang_day3_liana_arg1pos_topranked$ligand.complex), error = function(e) character(0))
wang_day3_identified_ligands_arg1neg <- tryCatch(unique(wang_day3_liana_arg1neg_topranked$ligand.complex), error = function(e) character(0))
wang_day3_identified_ligands <- tryCatch(unique(c(wang_day3_identified_ligands_arg1pos, wang_day3_identified_ligands_arg1neg)), error = function(e) character(0))
wang_day3_ligand_source_population <- tryCatch(c(rep("Arg1pos", length(wang_day3_identified_ligands_arg1pos)), rep("Arg1neg", length(wang_day3_identified_ligands_arg1neg))), error = function(e) character(0))
wang_day3_ligand_pathway_map <- data.frame(ligand = character(0), pathway = character(0), weight = numeric(0), stringsAsFactors = FALSE)
for (lig_idx in seq_along(wang_day3_identified_ligands)) {
  lig <- wang_day3_identified_ligands[lig_idx]
  lig_genes <- tryCatch(unlist(strsplit(lig, "_")), error = function(e) character(0))
  lig_genes <- tryCatch(unlist(strsplit(lig_genes, "[_+]")), error = function(e) lig_genes)
  lig_pathways <- tryCatch(dplyr::filter(progeny_network_wang_d3, target %in% lig_genes), error = function(e) data.frame())
  lig_pathway_df <- tryCatch(data.frame(ligand = lig, pathway = unique(lig_pathways$source), weight = lig_pathways$weight, stringsAsFactors = FALSE), error = function(e) data.frame())
  wang_day3_ligand_pathway_map <- tryCatch(rbind(wang_day3_ligand_pathway_map, lig_pathway_df), error = function(e) wang_day3_ligand_pathway_map)
}
write.csv(wang_day3_ligand_pathway_map, file.path(OUTPUT_DIR, "WangDay3_LigandToPathwayMapping.csv"), row.names = FALSE)

# Functional annotation: Exploratory annotation with population source (Arg1pos vs Arg1neg)
n_lig_wang <- length(wang_day3_identified_ligands)
wang_day3_ligand_vec <- character(n_lig_wang)
wang_day3_ligand_genes_vec <- character(n_lig_wang)
wang_day3_target_pop_vec <- character(n_lig_wang)
for (lig_idx in seq_len(n_lig_wang)) {
  lig <- wang_day3_identified_ligands[lig_idx]
  lig_genes <- unlist(strsplit(lig, "_"))
  lig_genes <- unlist(strsplit(lig_genes, "[_+]"))
  lig_upper <- toupper(lig_genes)
  lig_in_arg1pos <- lig %in% wang_day3_identified_ligands_arg1pos
  lig_in_arg1neg <- lig %in% wang_day3_identified_ligands_arg1neg
  lig_population_vector <- c("Arg1pos", "Arg1neg")[c(lig_in_arg1pos, lig_in_arg1neg)]
  lig_population_source <- paste(lig_population_vector, collapse = ";")
  wang_day3_ligand_vec[lig_idx] <- lig
  wang_day3_ligand_genes_vec[lig_idx] <- paste(lig_upper, collapse = ";")
  wang_day3_target_pop_vec[lig_idx] <- lig_population_source
}
wang_day3_ligand_annotation <- data.frame(ligand = wang_day3_ligand_vec, ligand_genes = wang_day3_ligand_genes_vec, target_population = wang_day3_target_pop_vec, stringsAsFactors = FALSE)
write.csv(wang_day3_ligand_annotation, file.path(OUTPUT_DIR, "WangDay3_LigandAnnotation.csv"), row.names = FALSE)

# PROGENy visualizations (Wang Day 3)
wang_day3_pathway_long <- tryCatch(tidyr::pivot_longer(wang_day3_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
p_wang_d3_progeny1 <- tryCatch(ggplot(wang_day3_pathway_long, aes(x = pathway, y = Pathway_Score, fill = Group)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("Arg1pos_mean" = "red", "Arg1neg_mean" = "lightblue"), labels = c("Arg1pos", "Arg1neg")) + labs(title = "PROGENy Pathway Activity - Wang Day 3", x = "Pathway", y = "Pathway Activity Score") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
vals_pw <- numeric(0)
if (exists("lee_day1_pathway_long") && is.data.frame(lee_day1_pathway_long) && nrow(lee_day1_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day1_pathway_long)) vals_pw <- c(vals_pw, lee_day1_pathway_long$Pathway_Score)
if (exists("lee_day3_pathway_long") && is.data.frame(lee_day3_pathway_long) && nrow(lee_day3_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day3_pathway_long)) vals_pw <- c(vals_pw, lee_day3_pathway_long$Pathway_Score)
if (exists("wang_day3_pathway_long") && is.data.frame(wang_day3_pathway_long) && nrow(wang_day3_pathway_long) > 0 && "Pathway_Score" %in% names(wang_day3_pathway_long)) vals_pw <- c(vals_pw, wang_day3_pathway_long$Pathway_Score)
if (exists("qin_day1_pathway_long") && is.data.frame(qin_day1_pathway_long) && nrow(qin_day1_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day1_pathway_long)) vals_pw <- c(vals_pw, qin_day1_pathway_long$Pathway_Score)
if (exists("qin_day3_pathway_long") && is.data.frame(qin_day3_pathway_long) && nrow(qin_day3_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day3_pathway_long)) vals_pw <- c(vals_pw, qin_day3_pathway_long$Pathway_Score)
vals_pw <- vals_pw[is.finite(vals_pw)]
UNIFY_PW_Y <- if (length(vals_pw) > 0) range(vals_pw) else c(-1, 1)
if (UNIFY_PW_Y[1] == UNIFY_PW_Y[2]) UNIFY_PW_Y[2] <- UNIFY_PW_Y[1] + 1e-6
print(p_wang_d3_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
progeny_heatmap_mat_wang_d3 <- tryCatch(rbind(Arg1pos = wang_day3_pathway_means_arg1pos, Arg1neg = wang_day3_pathway_means_arg1neg), error = function(e) matrix(0, nrow = 0, ncol = 0))
if (nrow(progeny_heatmap_mat_wang_d3) > 0 && ncol(progeny_heatmap_mat_wang_d3) > 0) {
  colors_progeny_wang <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors_use_progeny_wang <- grDevices::colorRampPalette(colors = colors_progeny_wang)(100)
  p_progeny_heatmap_wang_d3 <- pheatmap::pheatmap(progeny_heatmap_mat_wang_d3, color = colors_use_progeny_wang, border_color = "white", cellwidth = 20, cellheight = 20, main = "PROGENy Pathway Activity: Arg1+ vs Arg1- Neutrophils (Wang Day 3)")
  print(p_progeny_heatmap_wang_d3)
}
p_wang_d3_progeny2 <- tryCatch(ggplot(wang_day3_ligand_pathway_map, aes(x = ligand, y = pathway, size = abs(weight), color = weight)) + geom_point(alpha = 0.7) + scale_color_gradient2(low = "blue", mid = "white", high = "red") + labs(title = "Identified Ligands Linked to PROGENy Pathways - Wang Day 3", x = "Ligand", y = "Pathway") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
vals_w <- numeric(0)
map_tabnames <- c("lee_day1_ligand_pathway_map", "lee_day3_ligand_pathway_map", "wang_day3_ligand_pathway_map", "qin_day1_ligand_pathway_map", "qin_day3_ligand_pathway_map")
for (fn in map_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("weight" %in% names(d)) vals_w <- c(vals_w, d$weight)
}
vals_w <- vals_w[is.finite(vals_w)]
UNIFY_LIGW_ABS <- if (length(vals_w) > 0) range(abs(vals_w)) else c(0, 1)
if (UNIFY_LIGW_ABS[1] == UNIFY_LIGW_ABS[2]) UNIFY_LIGW_ABS[2] <- UNIFY_LIGW_ABS[1] + 1e-6
mxw <- if (length(vals_w) > 0) max(abs(vals_w)) else 1
if (!is.finite(mxw) || mxw <= 0) mxw <- 1
UNIFY_LIGW_COL <- c(-mxw, mxw)
print(p_wang_d3_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("✓ Wang Day 3 PROGENy pathway analysis complete")
# Save environment after Wang Day 3 PROGENy (resume: load(file.path(OUTPUT_DIR, "Workspace_AfterWangDay3_PROGENy.RData")))
save.image(file.path(OUTPUT_DIR, "Workspace_AfterWangDay3_PROGENy.RData"))
saveRDS(list(wang_day3 = wang_day3, wang_day3_neut_cells = wang_day3_neut_cells, wang_day3_pos = wang_day3_pos, wang_day3_neg = wang_day3_neg, wang_day3_pathway_results = wang_day3_pathway_results, progeny_network_wang_d3 = progeny_network_wang_d3, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterWangDay3_PROGENy.rds"))
# checkpoint_wang_d3_progeny <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterWangDay3_PROGENy.rds")); list2env(checkpoint_wang_d3_progeny, envir = .GlobalEnv)
print("✓ Environment saved: Workspace_AfterWangDay3_PROGENy.RData, Checkpoint_AfterWangDay3_PROGENy.rds")

# FRESH RSTUDIO — RESUME WANG DAY 3 BARTsc: readRDS Checkpoint_AfterWangDay3_PROGENy.rds + libraries; run Wang Day 3 BARTsc block only.

# -------- Wang Day 3 BARTsc TF Analysis (Arg1pos vs Arg1neg; downstream TF support only) --------
if (bartsc_initialized) {
  bartsc_initialized <- suppressWarnings(tryCatch({ BARTsc::load_bart2(); TRUE }, error = function(e) FALSE))
  if (!bartsc_initialized && exists("bart2", envir = .GlobalEnv)) bartsc_initialized <- TRUE
}
print("--- BARTsc TF Analysis: Wang Day 3 (Arg1pos vs Arg1neg neutrophils) ---")
if (bartsc_initialized) {
  options("mc.cores" = min(6L, parallel::detectCores()))
  wang_day3_neut_subset <- subset(wang_day3, cells = wang_day3_neut_cells)
  wang_day3_bart_label <- setNames(ifelse(colnames(wang_day3_neut_subset) %in% colnames(wang_day3)[wang_day3_pos], "Arg1pos", "Arg1neg"), colnames(wang_day3_neut_subset))
  wang_day3_bart_label <- factor(wang_day3_bart_label, levels = c("Arg1pos", "Arg1neg"))
  wang_day3_bart_proj <- BARTsc::bartsc(name = "WangDay3_Arg1", genome = "mm10", label = wang_day3_bart_label, cell_types_used = c("Arg1pos", "Arg1neg"), RNA_cnt_matrix = Seurat::GetAssayData(wang_day3_neut_subset, layer = "counts"))
  wang_day3_bart_proj <- BARTsc::normalize_RNA(wang_day3_bart_proj)
  wang_day3_bart_proj <- BARTsc::find_signature_genes(wang_day3_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  wang_day3_bart_proj <- BARTsc::find_pairwise_deg(wang_day3_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  n_wang_posneg <- 0L
  n_wang_negpos <- 0L
  wang_pw <- wang_day3_bart_proj@data$pairwise_DEG
  if (!is.null(wang_pw) && "Arg1pos::Arg1neg" %in% names(wang_pw)) n_wang_posneg <- { x <- wang_pw[["Arg1pos::Arg1neg"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  if (!is.null(wang_pw) && "Arg1neg::Arg1pos" %in% names(wang_pw)) n_wang_negpos <- { x <- wang_pw[["Arg1neg::Arg1pos"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  print(paste0("Wang Day 3 BARTsc pairwise DEG count (@data$pairwise_DEG): Arg1pos::Arg1neg=", n_wang_posneg, ", Arg1neg::Arg1pos=", n_wang_negpos))
  bart_deg_ok_wang_d3 <- if (BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT > 0) n_wang_posneg >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT && n_wang_negpos >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT else TRUE
  bart_nt_wang_d3 <- table(wang_day3_bart_label)
  bart_unequal_wang_d3 <- length(bart_nt_wang_d3) == 2L && length(unique(as.vector(bart_nt_wang_d3))) > 1L
  bart_crossct_allowed_wang_d3 <- !isTRUE(BARTSC_SKIP_CROSSCT_TEST) && length(bart_nt_wang_d3) >= 2L && all(as.integer(bart_nt_wang_d3) >= BARTSC_MIN_CELLS_CROSSCT) && !(isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_wang_d3) && bart_deg_ok_wang_d3
  wang_day3_bart_proj <- BARTsc::run_signature_RNA(wang_day3_bart_proj)
  wang_day3_bart_sig <- BARTsc::get_result(wang_day3_bart_proj, analysis = "cell type signature", mod = "RNA")
  wang_day3_bart_proj <- BARTsc::calc_crossCT_auc_RNA(wang_day3_bart_proj)
  print(bart_nt_wang_d3)
  if (!bart_crossct_allowed_wang_d3) message(paste0(
    "BARTsc Wang Day 3: crossCT_test / find_key_regulators skipped (BARTSC_SKIP_CROSSCT_TEST=", isTRUE(BARTSC_SKIP_CROSSCT_TEST),
    "; min_cells_ok=", all(as.integer(bart_nt_wang_d3) >= BARTSC_MIN_CELLS_CROSSCT),
    "; skip_if_unequal_n=", isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_wang_d3,
    "; deg_ok=", bart_deg_ok_wang_d3,
    " [Arg1pos::Arg1neg=", n_wang_posneg,
    ", Arg1neg::Arg1pos=", n_wang_negpos,
    ", min_each_dir=", BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT,
    "]). Vignette scRNA-seq.md §6: dot_plot/deviation_heatmap follow crossCT_test."
  ))
  if (bart_crossct_allowed_wang_d3) wang_day3_bart_proj <- BARTsc::crossCT_test(wang_day3_bart_proj, mod = "RNA")
  wang_day3_bart_cross <- if (bart_crossct_allowed_wang_d3) BARTsc::get_result(wang_day3_bart_proj, analysis = "cross-cell-type", mod = "RNA") else list()
  wang_day3_bart_cross_dev <- if (is.null(wang_day3_bart_cross)) list() else if (is.list(wang_day3_bart_cross) && "deviation" %in% names(wang_day3_bart_cross) && is.list(wang_day3_bart_cross$deviation)) wang_day3_bart_cross$deviation else wang_day3_bart_cross
  bart_outdir_wang_d3 <- file.path(OUTPUT_DIR, "BARTsc_WangDay3")
  dir.create(bart_outdir_wang_d3, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(wang_day3_bart_sig) && !is.null(wang_day3_bart_sig[["Arg1pos"]])) write.csv(wang_day3_bart_sig[["Arg1pos"]], file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_Signature_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(wang_day3_bart_sig) && !is.null(wang_day3_bart_sig[["Arg1neg"]])) write.csv(wang_day3_bart_sig[["Arg1neg"]], file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_Signature_Arg1neg.csv"), row.names = FALSE)
  wang_day3_bart_key <- NULL
  if (bart_crossct_allowed_wang_d3) {
    wang_day3_bart_proj <- BARTsc::find_key_regulators(wang_day3_bart_proj, mod = "RNA", min.N.profile = 3)
    wang_day3_bart_key <- BARTsc::get_result(wang_day3_bart_proj, analysis = "Key regs ident", mod = "RNA")
  }
  if (!is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1pos"]])) write.csv(wang_day3_bart_key[["Arg1pos"]], file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_KeyRegulators_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1neg"]])) write.csv(wang_day3_bart_key[["Arg1neg"]], file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_KeyRegulators_Arg1neg.csv"), row.names = FALSE)
  saveRDS(wang_day3_bart_proj, file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_Object.rds"))
  print("✓ Wang Day 3 BARTsc TF analysis complete — plots: Wang Day 3 BARTsc visualizations section below")
}

# -------- Wang Day 3: BARTsc visualizations (reload WangDay3_BARTsc_Object.rds; no re-run of bartsc pipeline) --------
bart_viz_wd3_rds <- file.path(OUTPUT_DIR, "BARTsc_WangDay3", "WangDay3_BARTsc_Object.rds")
if (!requireNamespace("BARTsc", quietly = TRUE)) message("Wang Day 3 BARTsc visualizations skipped: package BARTsc not installed.")
if (requireNamespace("BARTsc", quietly = TRUE) && !file.exists(bart_viz_wd3_rds)) message(paste0("Wang Day 3 BARTsc visualizations skipped: missing ", bart_viz_wd3_rds))
if (requireNamespace("BARTsc", quietly = TRUE) && file.exists(bart_viz_wd3_rds)) {
  if (!exists("BARTSC_N_TF_DOTHEAT")) BARTSC_N_TF_DOTHEAT <- 5L
  if (!exists("BARTSC_N_LABELED_TFS")) BARTSC_N_LABELED_TFS <- 5L
  if (!exists("BARTSC_KEYREG_SCATTER_W_IN")) BARTSC_KEYREG_SCATTER_W_IN <- 11
  if (!exists("BARTSC_KEYREG_SCATTER_H_IN")) BARTSC_KEYREG_SCATTER_H_IN <- 7
  if (!exists("BARTSC_KEYREG_XLIM")) BARTSC_KEYREG_XLIM <- c(-1, 1)
  if (!exists("BARTSC_KEYREG_YLIM_SIG")) BARTSC_KEYREG_YLIM_SIG <- c(0, 7)
  if (!exists("BARTSC_KEYREG_ZLIM_DMDR")) BARTSC_KEYREG_ZLIM_DMDR <- c(-2, 2)
  if (!exists("BARTSC_KEYREG_RANK_MAX")) BARTSC_KEYREG_RANK_MAX <- 60L
  suppressWarnings(tryCatch({ BARTsc::load_bart2(); NULL }, error = function(e) NULL))
  if (exists("bart2", envir = .GlobalEnv)) {
    if (!exists("types", envir = .GlobalEnv)) types <<- reticulate::import("types")
    bart2_mods_viz_wd3 <- reticulate::py_to_r(reticulate::py_get_attr(bart2, "__all__"))
    for (m_viz_wd3 in bart2_mods_viz_wd3) {
      if (!exists(m_viz_wd3, envir = .GlobalEnv)) assign(m_viz_wd3, reticulate::import(paste0("bart2.", m_viz_wd3), delay_load = TRUE), envir = .GlobalEnv)
    }
  }
  if (!exists("wang_day3_bart_proj", envir = .GlobalEnv, inherits = FALSE)) wang_day3_bart_proj <- readRDS(bart_viz_wd3_rds)
  bart_outdir_wang_d3 <- file.path(OUTPUT_DIR, "BARTsc_WangDay3")
  dir.create(bart_outdir_wang_d3, showWarnings = FALSE, recursive = TRUE)
  if (interactive() && grDevices::dev.cur() == 1L) grDevices::dev.new()
  wang_day3_bart_cross_viz <- BARTsc::get_result(wang_day3_bart_proj, analysis = "cross-cell-type", mod = "RNA")
  wang_day3_bart_cross_dev <- if (is.null(wang_day3_bart_cross_viz)) list() else if (is.list(wang_day3_bart_cross_viz) && "deviation" %in% names(wang_day3_bart_cross_viz) && is.list(wang_day3_bart_cross_viz$deviation)) wang_day3_bart_cross_viz$deviation else wang_day3_bart_cross_viz
  wang_day3_bart_key <- BARTsc::get_result(wang_day3_bart_proj, analysis = "Key regs ident", mod = "RNA")
  bart_tf_names_wang_d3 <- names(wang_day3_bart_cross_dev)
  tf_example_wang_d3 <- if (length(bart_tf_names_wang_d3) > 0) bart_tf_names_wang_d3[1] else NULL
  if (length(wang_day3_bart_cross_dev) > 0 && !is.null(tf_example_wang_d3)) {
    p_wang_d3_bart_dot <- BARTsc::dot_plot(wang_day3_bart_proj, mod = "RNA", tf = tf_example_wang_d3, max_dot_size = 22)
    print(p_wang_d3_bart_dot)
    ggplot2::ggsave(file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_DotPlot.png"), p_wang_d3_bart_dot, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_DotPlot.pdf"), p_wang_d3_bart_dot, width = 8, height = 6)
    p_wang_d3_bart_heat <- BARTsc::deviation_heatmap(wang_day3_bart_proj, mod = "RNA", tf = tf_example_wang_d3, tile_fontsize = 6)
    print(p_wang_d3_bart_heat)
    ggplot2::ggsave(file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_DeviationHeatmap.png"), p_wang_d3_bart_heat, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_DeviationHeatmap.pdf"), p_wang_d3_bart_heat, width = 8, height = 6)
  }
  top_tf_wang_pos <- if (length(wang_day3_bart_cross_dev) > 0 && !is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1pos"]]) && "TF" %in% colnames(wang_day3_bart_key[["Arg1pos"]])) head(wang_day3_bart_key[["Arg1pos"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_tf_wang_neg <- if (length(wang_day3_bart_cross_dev) > 0 && !is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1neg"]]) && "TF" %in% colnames(wang_day3_bart_key[["Arg1neg"]])) head(wang_day3_bart_key[["Arg1neg"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_key_tfs_wang <- intersect(unique(c(top_tf_wang_pos, top_tf_wang_neg)), names(wang_day3_bart_cross_dev))
  if (length(wang_day3_bart_cross_dev) > 0 && length(top_key_tfs_wang) > 0) {
    i_tf_wd3 <- 1L
    while (i_tf_wd3 <= length(top_key_tfs_wang)) {
      tf_cur_wang <- top_key_tfs_wang[i_tf_wd3]
      p_dot_wang <- BARTsc::dot_plot(wang_day3_bart_proj, mod = "RNA", tf = tf_cur_wang, max_dot_size = 22)
      print(p_dot_wang)
      ggplot2::ggsave(file.path(bart_outdir_wang_d3, paste0("WangDay3_BARTsc_DotPlot_", tf_cur_wang, ".png")), p_dot_wang, width = 8, height = 6)
      p_heat_wang <- BARTsc::deviation_heatmap(wang_day3_bart_proj, mod = "RNA", tf = tf_cur_wang, tile_fontsize = 6)
      print(p_heat_wang)
      ggplot2::ggsave(file.path(bart_outdir_wang_d3, paste0("WangDay3_BARTsc_DeviationHeatmap_", tf_cur_wang, ".png")), p_heat_wang, width = 8, height = 6)
      i_tf_wd3 <- i_tf_wd3 + 1L
    }
  }
  tfs_wang_d3_pos <- character(0)
  if (!is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1pos"]])) {
    df_wp <- wang_day3_bart_key[["Arg1pos"]]
    if (nrow(df_wp) > 0) {
      tf_c_wp <- intersect(colnames(df_wp), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_wp <- if (length(tf_c_wp) == 0) colnames(df_wp)[1] else tf_c_wp[1]
      ord_wp <- if ("final_rank" %in% colnames(df_wp)) order(df_wp$final_rank, na.last = TRUE) else seq_len(nrow(df_wp))
      tfs_wang_d3_pos <- as.character(head(df_wp[ord_wp, , drop = FALSE][[tf_col_wp]], BARTSC_N_LABELED_TFS))
    }
  }
  tfs_wang_d3_neg <- character(0)
  if (!is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1neg"]])) {
    df_wn <- wang_day3_bart_key[["Arg1neg"]]
    if (nrow(df_wn) > 0) {
      tf_c_wn <- intersect(colnames(df_wn), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_wn <- if (length(tf_c_wn) == 0) colnames(df_wn)[1] else tf_c_wn[1]
      ord_wn <- if ("final_rank" %in% colnames(df_wn)) order(df_wn$final_rank, na.last = TRUE) else seq_len(nrow(df_wn))
      tfs_wang_d3_neg <- as.character(head(df_wn[ord_wn, , drop = FALSE][[tf_col_wn]], BARTSC_N_LABELED_TFS))
    }
  }
  if (length(tfs_wang_d3_pos) > 0) {
    grDevices::png(file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_KeyRegScatter_Official_Arg1pos.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(wang_day3_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_wang_d3_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_wang_d3_pos) > 0) message("Displaying official key_regulator_scatter for Arg1pos (Wang Day 3)...")
  if (interactive() && length(tfs_wang_d3_pos) > 0) key_regulator_scatter_unified(wang_day3_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_wang_d3_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
  if (length(tfs_wang_d3_neg) > 0) {
    grDevices::png(file.path(bart_outdir_wang_d3, "WangDay3_BARTsc_KeyRegScatter_Official_Arg1neg.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(wang_day3_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_wang_d3_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_wang_d3_neg) > 0) message("Displaying official key_regulator_scatter for Arg1neg (Wang Day 3)...")
  if (interactive() && length(tfs_wang_d3_neg) > 0) key_regulator_scatter_unified(wang_day3_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_wang_d3_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
  print("✓ Wang Day 3 BARTsc visualizations complete")
}

# ============================================================================
# INTEGRATION: LIANA -> BARTsc + PROGENy (Wang Day 3) — same logic as Lee Day 3
# ============================================================================
# Linear: footprint co-membership CSV → OmniPath receptor→TF → full-chain CSV → sender×pathway heatmap (BARTsc plots: section above).
print("--- Integration: LIANA -> BARTsc + PROGENy (Wang Day 3, Data-Driven) ---")
# Exploratory overlay only: LIANA provides incoming-signal candidates; PROGENy/BARTsc/OmniPath add downstream support, not causal proof.
if (!exists("progeny_network_wang_d3") || !all(c("source", "target") %in% colnames(progeny_network_wang_d3)) || !exists("wang_day3_arg1pos_received") || nrow(wang_day3_arg1pos_received) == 0) {
  message("Wang Day 3 integration skipped (need PROGENy network + LIANA Arg1pos received with rows).")
} else {
liana_receptors_wang_d3 <- unique(wang_day3_arg1pos_received$receptor.complex)
liana_receptors_wang_d3_split <- unique(trimws(unlist(strsplit(liana_receptors_wang_d3, "[_+]"))))
print(paste("Unique receptors received by Arg1+ neutrophils (LIANA, Wang Day 3):", length(liana_receptors_wang_d3_split)))
liana_expanded_wang_d3 <- wang_day3_arg1pos_received
liana_expanded_wang_d3$receptor_subunits <- lapply(strsplit(liana_expanded_wang_d3$receptor.complex, "[_+]"), function(x) trimws(x))
liana_expanded_wang_d3 <- tidyr::unnest_longer(liana_expanded_wang_d3, receptor_subunits)
pathways_active_wang_d3 <- names(wang_day3_pathway_means_arg1pos)[wang_day3_pathway_means_arg1pos > wang_day3_pathway_means_arg1neg]
# Exploratory footprint overlap (not causal receptor -> pathway activation)
receptor_pathway_links_wang_d3 <- progeny_network_wang_d3[progeny_network_wang_d3$target %in% liana_receptors_wang_d3_split & progeny_network_wang_d3$source %in% pathways_active_wang_d3, ]
receptor_pathway_links_wang_d3 <- receptor_pathway_links_wang_d3[, c("target", "source")]
colnames(receptor_pathway_links_wang_d3) <- c("receptor_gene", "pathway")
pathway_diff_wang_d3 <- wang_day3_pathway_means_arg1pos - wang_day3_pathway_means_arg1neg
receptor_pathway_links_wang_d3$pathway_activity_diff <- pathway_diff_wang_d3[receptor_pathway_links_wang_d3$pathway]
integration_receptor_pathway_wang_d3 <- merge(liana_expanded_wang_d3, receptor_pathway_links_wang_d3, by.x = "receptor_subunits", by.y = "receptor_gene", all.x = TRUE)
integration_receptor_pathway_wang_d3 <- integration_receptor_pathway_wang_d3[!is.na(integration_receptor_pathway_wang_d3$pathway), ]
write.csv(integration_receptor_pathway_wang_d3, file.path(OUTPUT_DIR, "WangDay3_LIANA_PROGENy_FootprintCoMembership_Exploratory.csv"), row.names = FALSE)
print(paste("LIANA x PROGENy footprint co-membership (exploratory, Wang Day 3):", nrow(integration_receptor_pathway_wang_d3), "rows"))
receptor_tf_links_wang_d3 <- data.frame(receptor_gene = character(0), tf = character(0), omnipath_via = character(0), omnipath_hops = integer(0), stringsAsFactors = FALSE)
active_tfs_wang_d3 <- character(0)
has_bart_key_wang_d3 <- exists("wang_day3_bart_key") && !is.null(wang_day3_bart_key) && !is.null(wang_day3_bart_key[["Arg1pos"]]) && nrow(wang_day3_bart_key[["Arg1pos"]]) > 0
if (has_bart_key_wang_d3) {
  bart_tfs_wang_d3 <- wang_day3_bart_key[["Arg1pos"]]
  tf_col_wang_d3 <- intersect(colnames(bart_tfs_wang_d3), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
  if (length(tf_col_wang_d3) == 0) tf_col_wang_d3 <- colnames(bart_tfs_wang_d3)[1] else tf_col_wang_d3 <- tf_col_wang_d3[1]
  active_tfs_wang_d3 <- unique(as.character(bart_tfs_wang_d3[[tf_col_wang_d3]]))
}
omnipath_ok_wang_d3 <- requireNamespace("OmnipathR", quietly = TRUE) && length(active_tfs_wang_d3) > 0

pathway_interactions_wang_d3 <- data.frame()
omni_wang_d3_first <- if (!omnipath_ok_wang_d3) list(dat = NULL, err1 = NA_character_) else tryCatch(list(dat = OmnipathR::import_omnipath_interactions(datasets = c("omnipath", "pathwayextra"), organism = 10090, genesymbols = TRUE), err1 = NA_character_), error = function(e1) list(dat = NULL, err1 = conditionMessage(e1)))
if (omnipath_ok_wang_d3) pathway_interactions_wang_d3 <- omni_wang_d3_first$dat
omnipath_wang_d3_err1 <- omni_wang_d3_first$err1
if (omnipath_ok_wang_d3 && is.null(pathway_interactions_wang_d3)) pathway_interactions_wang_d3 <- tryCatch(OmnipathR::import_pathwayextra_interactions(organism = 10090, genesymbols = TRUE), error = function(e2) { message("OmniPath curated+pathwayextra import failed (Wang Day 3): ", omnipath_wang_d3_err1, "; fallback: ", conditionMessage(e2)); data.frame() })

omni_wang_d3_cols_ok <- omnipath_ok_wang_d3 && nrow(pathway_interactions_wang_d3) > 0 && all(c("source_genesymbol", "target_genesymbol") %in% colnames(pathway_interactions_wang_d3))
omni_direct_wang_d3 <- if (omni_wang_d3_cols_ok) subset(pathway_interactions_wang_d3, source_genesymbol %in% liana_receptors_wang_d3_split & target_genesymbol %in% active_tfs_wang_d3, select = c("source_genesymbol", "target_genesymbol")) else data.frame()
if (omni_wang_d3_cols_ok && nrow(omni_direct_wang_d3) > 0) receptor_tf_links_wang_d3 <- rbind(receptor_tf_links_wang_d3, data.frame(receptor_gene = omni_direct_wang_d3$source_genesymbol, tf = omni_direct_wang_d3$target_genesymbol, omnipath_via = NA_character_, omnipath_hops = 1L, stringsAsFactors = FALSE))

hop1_wang_d3 <- if (omni_wang_d3_cols_ok) pathway_interactions_wang_d3[pathway_interactions_wang_d3$source_genesymbol %in% liana_receptors_wang_d3_split, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_wang_d3_cols_ok) colnames(hop1_wang_d3) <- c("receptor_gene", "via")
hop2_wang_d3 <- if (omni_wang_d3_cols_ok) pathway_interactions_wang_d3[pathway_interactions_wang_d3$target_genesymbol %in% active_tfs_wang_d3, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_wang_d3_cols_ok) colnames(hop2_wang_d3) <- c("via", "tf")
omni_2hop_wang_d3 <- if (omni_wang_d3_cols_ok) merge(hop1_wang_d3, hop2_wang_d3, by = "via") else data.frame()
if (omni_wang_d3_cols_ok && nrow(omni_2hop_wang_d3) > 0) omni_2hop_wang_d3 <- unique(omni_2hop_wang_d3[, c("receptor_gene", "tf", "via")])
if (omni_wang_d3_cols_ok && nrow(omni_2hop_wang_d3) > 0) receptor_tf_links_wang_d3 <- rbind(receptor_tf_links_wang_d3, data.frame(receptor_gene = omni_2hop_wang_d3$receptor_gene, tf = omni_2hop_wang_d3$tf, omnipath_via = omni_2hop_wang_d3$via, omnipath_hops = 2L, stringsAsFactors = FALSE))
if (omni_wang_d3_cols_ok) receptor_tf_links_wang_d3 <- receptor_tf_links_wang_d3[!duplicated(paste(receptor_tf_links_wang_d3$receptor_gene, receptor_tf_links_wang_d3$tf)), ]

integration_receptor_tf_wang_d3 <- data.frame()
if (nrow(receptor_tf_links_wang_d3) > 0) {
  integration_receptor_tf_wang_d3 <- merge(liana_expanded_wang_d3, receptor_tf_links_wang_d3, by.x = "receptor_subunits", by.y = "receptor_gene", all.x = TRUE)
  integration_receptor_tf_wang_d3 <- integration_receptor_tf_wang_d3[!is.na(integration_receptor_tf_wang_d3$tf), ]
}
if (nrow(integration_receptor_tf_wang_d3) > 0) {
  write.csv(integration_receptor_tf_wang_d3, file.path(OUTPUT_DIR, "WangDay3_LIANA_BARTsc_Integration_DataDriven.csv"), row.names = FALSE)
  n1_wang <- sum(integration_receptor_tf_wang_d3$omnipath_hops == 1L, na.rm = TRUE)
  n2_wang <- sum(integration_receptor_tf_wang_d3$omnipath_hops == 2L, na.rm = TRUE)
  print(paste0("LIANA -> BARTsc integration (Wang Day 3): ", nrow(integration_receptor_tf_wang_d3), " receptor-TF rows (OmniPath curated+pathwayextra; direct=", n1_wang, ", two-hop=", n2_wang, ")"))
}
if (nrow(receptor_tf_links_wang_d3) == 0 && length(active_tfs_wang_d3) > 0) {
  write.csv(data.frame(receptor = liana_receptors_wang_d3_split, note = "Active TFs in Arg1pos (no OmniPath direct/two-hop link):", active_tfs = paste(active_tfs_wang_d3, collapse = "; "), stringsAsFactors = FALSE), file.path(OUTPUT_DIR, "WangDay3_LIANA_BARTsc_Receptors_and_TFs.csv"), row.names = FALSE)
  print("LIANA receptors and BARTsc TFs saved separately (Wang Day 3, no OmniPath link)")
}
has_pathway_wang_d3 <- nrow(integration_receptor_pathway_wang_d3) > 0
has_tf_wang_d3 <- nrow(integration_receptor_tf_wang_d3) > 0
full_chain_wang_d3 <- data.frame()
if (has_pathway_wang_d3 && has_tf_wang_d3) full_chain_wang_d3 <- merge(integration_receptor_pathway_wang_d3, integration_receptor_tf_wang_d3[, c("source", "ligand.complex", "receptor_subunits", "aggregate_rank", "tf", "omnipath_via", "omnipath_hops")], by = c("source", "ligand.complex", "receptor_subunits", "aggregate_rank"), all = TRUE)
if (has_pathway_wang_d3 && !has_tf_wang_d3) { full_chain_wang_d3 <- integration_receptor_pathway_wang_d3; full_chain_wang_d3$tf <- NA_character_ }
if (!has_pathway_wang_d3 && has_tf_wang_d3) { full_chain_wang_d3 <- integration_receptor_tf_wang_d3; full_chain_wang_d3$pathway <- NA_character_; full_chain_wang_d3$pathway_activity_diff <- NA_real_ }
if (nrow(full_chain_wang_d3) > 0) full_chain_wang_d3$evidence_level <- ifelse(!is.na(full_chain_wang_d3$pathway) & !is.na(full_chain_wang_d3$tf), "STRONG (pathway + TF linked)", ifelse(!is.na(full_chain_wang_d3$pathway) | !is.na(full_chain_wang_d3$tf), "MODERATE (pathway or TF linked)", "WEAK (no link)"))
if (nrow(full_chain_wang_d3) > 0) full_chain_wang_d3 <- full_chain_wang_d3[order(full_chain_wang_d3$evidence_level, full_chain_wang_d3$aggregate_rank), ]
# Legacy filename retained for compatibility; contents are an exploratory overlay, not a causal signal chain.
if (nrow(full_chain_wang_d3) > 0) write.csv(full_chain_wang_d3, file.path(OUTPUT_DIR, "WangDay3_FullSignalChain_DataDriven.csv"), row.names = FALSE)
if (nrow(full_chain_wang_d3) > 0) {
  cols_show_wang_d3 <- intersect(c("source", "receptor_subunits", "pathway", "tf", "omnipath_via", "omnipath_hops", "evidence_level"), colnames(full_chain_wang_d3))
  print(head(full_chain_wang_d3[, cols_show_wang_d3, drop = FALSE], 20))
}
heatmap_data_wang_d3 <- data.frame()
if (nrow(full_chain_wang_d3) > 0 && has_pathway_wang_d3) heatmap_data_wang_d3 <- as.data.frame.matrix(table(integration_receptor_pathway_wang_d3$source, integration_receptor_pathway_wang_d3$pathway))
if (nrow(heatmap_data_wang_d3) > 0 && ncol(heatmap_data_wang_d3) > 0) pheatmap::pheatmap(heatmap_data_wang_d3, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("white", "blue", "red"))(50), main = "Wang Day 3: Sender x pathway (PROGENy footprint co-membership, exploratory)")
if (nrow(full_chain_wang_d3) == 0) print("No integration rows Wang Day 3 (PROGENy/BARTsc may not overlap LIANA receptors)")
print("✓ Wang Day 3 LIANA <-> BARTsc <-> PROGENy integration complete (exploratory overlay)")
}

# -------- Qin Day 1: LIANA → LIANA plots → top-receptor Vln → PROGENy → BARTsc → integration (same cohort order as Lee Day 1) --------
print("--- LIANA Analysis: Qin Day 1 ---")

qin_day1_liana_labels <- qin_day1_pruned_labels
qin_day1_liana_labels[qin_day1_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
qin_day1_liana_labels[qin_day1_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
qin_day1_liana_labels[is.na(qin_day1_liana_labels) | qin_day1_liana_labels == ""] <- "Other"
Seurat::Idents(qin_day1) <- factor(qin_day1_liana_labels)

# Qin Day 1: set labels, clean Seurat, build SCE (bypass validObject bug), run LIANA
qin_day1$liana_label <- factor(qin_day1_liana_labels)
SeuratObject::DefaultAssay(qin_day1) <- "RNA"
for (assay_name in names(qin_day1@assays)) {
  if (assay_name != "RNA") qin_day1[[assay_name]] <- NULL
}
for (red_name in names(qin_day1@reductions)) qin_day1[[red_name]] <- NULL
for (graph_name in names(qin_day1@graphs)) qin_day1[[graph_name]] <- NULL
# Custom L-R from CellChat + CellCall (Qin Day 1) — same logic as Lee Day 1.
qin_day1_liana_resource <- "MouseConsensus"
qin_day1_liana_external <- NULL
custom_lr_list_qin_d1 <- list()
path_cc_qin_d1 <- CELLCHAT_QIN_DAY1_RDS
path_ccall_qin_d1 <- CELLCALL_QIN_DAY1_RDS
cc_qin_d1_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_cc_qin_d1) && nzchar(trimws(path_cc_qin_d1)) && file.exists(path_cc_qin_d1) && requireNamespace("CellChat", quietly = TRUE)
cc_qin_d1 <- NULL
if (cc_qin_d1_ok) cc_qin_d1 <- readRDS(path_cc_qin_d1)
comm_qin_d1 <- NULL
if (!is.null(cc_qin_d1) && inherits(cc_qin_d1, "CellChat")) comm_qin_d1 <- CellChat::subsetCommunication(cc_qin_d1)
has_comm_qin_d1 <- !is.null(comm_qin_d1) && nrow(comm_qin_d1) > 0 && "ligand" %in% colnames(comm_qin_d1) && "receptor" %in% colnames(comm_qin_d1)
if (has_comm_qin_d1) {
  cc_lr_qin_d1 <- data.frame(source_genesymbol = comm_qin_d1$ligand, target_genesymbol = comm_qin_d1$receptor, stringsAsFactors = FALSE)
  cc_lr_qin_d1 <- cc_lr_qin_d1[!duplicated(cc_lr_qin_d1), ]
  custom_lr_list_qin_d1[["CellChat"]] <- cc_lr_qin_d1
}
ccall_qin_d1_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_ccall_qin_d1) && nzchar(trimws(path_ccall_qin_d1)) && file.exists(path_ccall_qin_d1)
ccall_qin_d1 <- NULL
if (ccall_qin_d1_ok) ccall_qin_d1 <- readRDS(path_ccall_qin_d1)
has_ccall_lr_qin_d1 <- !is.null(ccall_qin_d1) && !is.null(ccall_qin_d1@data$expr_l_r_log2_scale)
if (has_ccall_lr_qin_d1) {
  lr_rownames_qin_d1 <- rownames(ccall_qin_d1@data$expr_l_r_log2_scale)
  ccall_lig_qin_d1 <- sub("-.*", "", lr_rownames_qin_d1)
  ccall_rec_qin_d1 <- sub("^[^-]+-", "", lr_rownames_qin_d1)
  ccall_rec_qin_d1 <- gsub("-", "_", ccall_rec_qin_d1)
  ccall_lr_qin_d1 <- data.frame(source_genesymbol = ccall_lig_qin_d1, target_genesymbol = ccall_rec_qin_d1, stringsAsFactors = FALSE)
  ccall_lr_qin_d1 <- ccall_lr_qin_d1[nzchar(ccall_lr_qin_d1$target_genesymbol), ]
  ccall_lr_qin_d1 <- unique(ccall_lr_qin_d1)
  custom_lr_list_qin_d1[["CellCall"]] <- ccall_lr_qin_d1
}
has_custom_lr_qin_d1 <- length(custom_lr_list_qin_d1) > 0
if (has_custom_lr_qin_d1) {
  custom_lr_qin_d1 <- do.call(rbind, custom_lr_list_qin_d1)
  custom_lr_qin_d1 <- custom_lr_qin_d1[!duplicated(custom_lr_qin_d1[, c("source_genesymbol", "target_genesymbol")]), ]
  qin_day1_liana_external <- data.frame(source_genesymbol = custom_lr_qin_d1$source_genesymbol, target_genesymbol = custom_lr_qin_d1$target_genesymbol, stringsAsFactors = FALSE)
  qin_day1_liana_resource <- "custom"
  print(paste("Loaded", nrow(qin_day1_liana_external), "custom L-R pairs from CellChat + CellCall (Qin Day 1)"))
}
if (isTRUE(LIANA_USE_CUSTOM_LR) && !has_custom_lr_qin_d1) print("Qin Day 1 LIANA: LIANA_USE_CUSTOM_LR is TRUE but no CellChat/CellCall pairs loaded — using MouseConsensus only (check RDS paths and CellChat package).")
qin_day1_counts <- tryCatch(Seurat::GetAssayData(qin_day1, slot = "counts"), error = function(e) Seurat::GetAssayData(qin_day1, layer = "counts"))
qin_day1_logcounts <- tryCatch(Seurat::GetAssayData(qin_day1, slot = "data"), error = function(e) Seurat::GetAssayData(qin_day1, layer = "data"))
qin_day1_sce <- SingleCellExperiment(assays = list(counts = qin_day1_counts, logcounts = qin_day1_logcounts))
qin_day1_sce$liana_label <- qin_day1$liana_label
SingleCellExperiment::colLabels(qin_day1_sce) <- qin_day1_sce$liana_label
liana_args_qin_d1 <- list(sce = qin_day1_sce, method = LIANA_METHODS, resource = qin_day1_liana_resource, idents_col = "liana_label", expr_prop = 0.05, verbose = TRUE, min_cells = LIANA_MIN_CELLS, base = exp(1))
if (!is.null(qin_day1_liana_external)) liana_args_qin_d1$external_resource <- qin_day1_liana_external
qin_day1_liana_result <- do.call(liana::liana_wrap, liana_args_qin_d1)
qin_day1_liana_result_df <- liana::liana_aggregate(qin_day1_liana_result)

# Keep only interactions where target is Arg1+ neutrophils (signals RECEIVED by Arg1+ neutrophils)
qin_day1_arg1pos_received <- dplyr::filter(qin_day1_liana_result_df, target == "NeutrophilArg1pos")
qin_day1_arg1pos_received <- dplyr::arrange(qin_day1_arg1pos_received, aggregate_rank)
qin_day1_arg1pos_top_signals <- head(qin_day1_arg1pos_received, 50)
# Arg1- neutrophils received (for specificity contrast)
qin_day1_arg1neg_received <- dplyr::filter(qin_day1_liana_result_df, target == "NeutrophilArg1neg")
qin_day1_arg1neg_received <- dplyr::arrange(qin_day1_arg1neg_received, aggregate_rank)

# Specificity: signals preferentially received by Arg1+ vs Arg1- neutrophils
topN <- 50
qin_day1_arg1pos_top <- head(qin_day1_arg1pos_received, topN)
qin_day1_arg1neg_top <- head(qin_day1_arg1neg_received, topN)
pos_pairs <- paste0(qin_day1_arg1pos_top$ligand.complex, "_", qin_day1_arg1pos_top$receptor.complex)
neg_pairs_all <- paste0(qin_day1_arg1neg_received$ligand.complex, "_", qin_day1_arg1neg_received$receptor.complex)
neg_rank_lookup <- setNames(qin_day1_arg1neg_received$aggregate_rank, neg_pairs_all)
pos_ranks <- qin_day1_arg1pos_top$aggregate_rank
neg_ranks_matched <- neg_rank_lookup[pos_pairs]
neg_ranks_matched[is.na(neg_ranks_matched)] <- 1.0
spec_index <- (neg_ranks_matched - pos_ranks) / (neg_ranks_matched + pos_ranks + 1e-10)
qin_day1_arg1pos_top$specificity_index <- spec_index
qin_day1_arg1pos_specific <- dplyr::filter(qin_day1_arg1pos_top, specificity_index > 0.3 | !(pos_pairs %in% neg_pairs_all))

# Backward-compatible names for downstream code
qin_day1_liana_arg1_pos <- qin_day1_arg1pos_received
qin_day1_liana_arg1_neg <- qin_day1_arg1neg_received
qin_day1_liana_consensus <- head(qin_day1_arg1pos_received, 30)
qin_day1_liana_arg1pos_topranked <- head(qin_day1_arg1pos_received, 20)
qin_day1_liana_arg1neg_topranked <- head(qin_day1_arg1neg_received, 20)

print(paste("Arg1pos received signals (all):", nrow(qin_day1_arg1pos_received)))
print(paste("Arg1pos top signals (top 50):", nrow(qin_day1_arg1pos_top_signals)))
print(paste("Arg1pos-specific signals (specificity > 0.3 or unique to Arg1pos):", nrow(qin_day1_arg1pos_specific)))

write.csv(qin_day1_liana_result_df, file.path(OUTPUT_DIR, "QinDay1_LIANA_AllResults.csv"), row.names = FALSE)
write.csv(qin_day1_arg1pos_received, file.path(OUTPUT_DIR, "QinDay1_Arg1pos_ReceivedSignals.csv"), row.names = FALSE)
write.csv(qin_day1_arg1pos_top_signals, file.path(OUTPUT_DIR, "QinDay1_LIANA_Arg1pos_Top50Signals.csv"), row.names = FALSE)
write.csv(qin_day1_arg1pos_specific, file.path(OUTPUT_DIR, "QinDay1_LIANA_Arg1pos_Specific.csv"), row.names = FALSE)
write.csv(qin_day1_liana_consensus, file.path(OUTPUT_DIR, "QinDay1_LIANA_NeutrophilConsensus.csv"), row.names = FALSE)
write.csv(qin_day1_liana_arg1pos_topranked, file.path(OUTPUT_DIR, "QinDay1_LIANA_Arg1pos_TopRanked.csv"), row.names = FALSE)
write.csv(qin_day1_liana_arg1neg_topranked, file.path(OUTPUT_DIR, "QinDay1_LIANA_Arg1neg_TopRanked.csv"), row.names = FALSE)
saveRDS(list(liana_result = qin_day1_liana_result, liana_aggregated = qin_day1_liana_result_df, arg1pos_received = qin_day1_arg1pos_received, arg1pos_top_signals = qin_day1_arg1pos_top_signals, arg1pos_specific = qin_day1_arg1pos_specific, neutrophil_consensus = qin_day1_liana_consensus, arg1pos_topranked = qin_day1_liana_arg1pos_topranked, arg1neg_topranked = qin_day1_liana_arg1neg_topranked), file.path(OUTPUT_DIR, "QinDay1_LIANA_Results.rds"))
# Uncomment block below to load and skip re-running Qin Day 1 LIANA:
# qin_day1_liana_loaded <- readRDS(file.path(OUTPUT_DIR, "QinDay1_LIANA_Results.rds"))
# qin_day1_liana_result <- qin_day1_liana_loaded$liana_result
# qin_day1_liana_result_df <- qin_day1_liana_loaded$liana_aggregated
# qin_day1_arg1pos_received <- qin_day1_liana_loaded$arg1pos_received
# qin_day1_arg1pos_top_signals <- qin_day1_liana_loaded$arg1pos_top_signals
# qin_day1_arg1pos_specific <- qin_day1_liana_loaded$arg1pos_specific
# qin_day1_liana_consensus <- qin_day1_liana_loaded$neutrophil_consensus
# qin_day1_liana_arg1pos_topranked <- qin_day1_liana_loaded$arg1pos_topranked
# qin_day1_liana_arg1neg_topranked <- qin_day1_liana_loaded$arg1neg_topranked

# -------- Qin Day 1: LIANA visualizations -> top LIANA receptors VlnPlot (same flow as Lee Day 1) --------
if (nrow(qin_day1_liana_result_df) > 0) {
  liana_network_qin_d1_both <- qin_day1_liana_result_df |>
    dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |>
    dplyr::filter(!(source %in% NEUTROPHIL_STATES)) |>
    dplyr::group_by(target) |>
    dplyr::arrange(aggregate_rank) |>
    dplyr::slice_head(n = 20) |>
    dplyr::ungroup()
  liana_network_qin_d1_both$interaction_label <- paste0(liana_network_qin_d1_both$ligand.complex, "\u2013(", liana_network_qin_d1_both$receptor.complex, ")")
  p_qin_d1_network <- tryCatch(ggplot(liana_network_qin_d1_both, aes(x = target, y = interaction_label, size = -log10(aggregate_rank + 1e-10), color = source)) + geom_point(alpha = 0.8) + scale_x_discrete(limits = NEUTROPHIL_STATES, drop = FALSE) + theme_minimal() + labs(title = "Top Signals Received: Arg1+ vs Arg1- (Qin Day 1)", subtitle = "Y = Ligand\u2013(Receptor); X = receiver; color = sender (neutrophil\u2192neutrophil excluded from data)", x = "Receiver (neutrophil state)", y = "Interaction", color = "Sender cell type", size = "Consensus support\n(-log10 aggregate rank)") + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), plot.margin = margin(10, 80, 10, 10)) + guides(color = guide_legend(override.aes = list(size = 3))), error = function(e) ggplot() + theme_void())
  vals_nlr <- numeric(0)
  if (exists("liana_network_lee_d1_both") && is.data.frame(liana_network_lee_d1_both) && nrow(liana_network_lee_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_lee_d3_both") && is.data.frame(liana_network_lee_d3_both) && nrow(liana_network_lee_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_wang_d3_both") && is.data.frame(liana_network_wang_d3_both) && nrow(liana_network_wang_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_wang_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d1_both") && is.data.frame(liana_network_qin_d1_both) && nrow(liana_network_qin_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d3_both") && is.data.frame(liana_network_qin_d3_both) && nrow(liana_network_qin_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d3_both$aggregate_rank + 1e-10))
  vals_nlr <- vals_nlr[is.finite(vals_nlr)]
  UNIFY_LIANA_NEGLOG10 <- if (length(vals_nlr) > 0) range(vals_nlr) else c(0, 10)
  if (length(UNIFY_LIANA_NEGLOG10) != 2 || !all(is.finite(UNIFY_LIANA_NEGLOG10))) UNIFY_LIANA_NEGLOG10 <- c(0, 10)
  if (UNIFY_LIANA_NEGLOG10[1] == UNIFY_LIANA_NEGLOG10[2]) UNIFY_LIANA_NEGLOG10[2] <- UNIFY_LIANA_NEGLOG10[1] + 1e-6
  print(p_qin_d1_network + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = (if (exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE else c(2, 8))))
  qin_day1_external_to_neutrophils <- qin_day1_liana_result_df |> dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |> dplyr::filter(!(source %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))) |> dplyr::arrange(aggregate_rank)
  liana_top_qin_d1_external <- dplyr::slice_head(qin_day1_external_to_neutrophils, n = 20)
  liana_top_qin_d1_external$lr_label <- paste0(liana_top_qin_d1_external$source, " -> ", liana_top_qin_d1_external$target, "  ", liana_top_qin_d1_external$ligand.complex, "\u2013(", liana_top_qin_d1_external$receptor.complex, ")")
  p_qin_d1_f0 <- tryCatch(ggplot(liana_top_qin_d1_external, aes(x = reorder(lr_label, aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 20 L-R Pairs - LIANA (Qin Day 1)\nExternal signals to Arg1+ / Arg1- neutrophils", x = "Source -> Target  Ligand\u2013(Receptor)", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
  vals_ext <- numeric(0)
  if (exists("lee_day1_external_to_neutrophils") && is.data.frame(lee_day1_external_to_neutrophils) && nrow(lee_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day1_external_to_neutrophils$aggregate_rank)
  if (exists("lee_day3_external_to_neutrophils") && is.data.frame(lee_day3_external_to_neutrophils) && nrow(lee_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day3_external_to_neutrophils$aggregate_rank)
  if (exists("wang_day3_external_to_neutrophils") && is.data.frame(wang_day3_external_to_neutrophils) && nrow(wang_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(wang_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, wang_day3_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day1_external_to_neutrophils") && is.data.frame(qin_day1_external_to_neutrophils) && nrow(qin_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day1_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day3_external_to_neutrophils") && is.data.frame(qin_day3_external_to_neutrophils) && nrow(qin_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day3_external_to_neutrophils$aggregate_rank)
  if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0 && "aggregate_rank" %in% names(liana_top_lee_d1)) vals_ext <- c(vals_ext, liana_top_lee_d1$aggregate_rank)
  if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_lee_d3_external)) vals_ext <- c(vals_ext, liana_top_lee_d3_external$aggregate_rank)
  if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_wang_d3_external)) vals_ext <- c(vals_ext, liana_top_wang_d3_external$aggregate_rank)
  if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d1_external)) vals_ext <- c(vals_ext, liana_top_qin_d1_external$aggregate_rank)
  if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d3_external)) vals_ext <- c(vals_ext, liana_top_qin_d3_external$aggregate_rank)
  vals_ext <- suppressWarnings(as.numeric(vals_ext))
  vals_ext <- vals_ext[is.finite(vals_ext)]
  plot_ext_ar <- suppressWarnings(as.numeric(liana_top_qin_d1_external[["aggregate_rank"]]))
  plot_ext_ar <- plot_ext_ar[is.finite(plot_ext_ar)]
  UNIFY_EXT_AR <- range(c(vals_ext, plot_ext_ar), na.rm = TRUE)
  if (length(plot_ext_ar) == 0 && length(vals_ext) == 0) UNIFY_EXT_AR <- c(0, 1)
  if (!all(is.finite(UNIFY_EXT_AR))) UNIFY_EXT_AR <- c(0, 1)
  if (UNIFY_EXT_AR[1] == UNIFY_EXT_AR[2]) UNIFY_EXT_AR[2] <- UNIFY_EXT_AR[1] + max(abs(UNIFY_EXT_AR[1]) * 1e-6, 1e-12)
  print(p_qin_d1_f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
}
liana_top_qin_d1 <- dplyr::slice_head(qin_day1_liana_consensus, n = 15)
if (nrow(liana_top_qin_d1) > 0) {
  p_qin_d1_f <- ggplot(liana_top_qin_d1, aes(x = reorder(paste0(source, " -> ", target), aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 15 L-R Pairs - LIANA (Qin Day 1)", x = "Source -> Target", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME
  vals_cons <- numeric(0)
  if (exists("liana_top_lee_d3") && is.data.frame(liana_top_lee_d3) && nrow(liana_top_lee_d3) > 0) vals_cons <- c(vals_cons, liana_top_lee_d3$aggregate_rank)
  if (exists("liana_top_wang_d3") && is.data.frame(liana_top_wang_d3) && nrow(liana_top_wang_d3) > 0) vals_cons <- c(vals_cons, liana_top_wang_d3$aggregate_rank)
  if (exists("liana_top_qin_d1") && is.data.frame(liana_top_qin_d1) && nrow(liana_top_qin_d1) > 0) vals_cons <- c(vals_cons, liana_top_qin_d1$aggregate_rank)
  if (exists("liana_top_qin_d3") && is.data.frame(liana_top_qin_d3) && nrow(liana_top_qin_d3) > 0) vals_cons <- c(vals_cons, liana_top_qin_d3$aggregate_rank)
  vals_cons <- vals_cons[is.finite(vals_cons)]
  UNIFY_CONS_AR <- if (length(vals_cons) > 0) range(vals_cons) else c(0, 1)
  if (UNIFY_CONS_AR[1] == UNIFY_CONS_AR[2]) UNIFY_CONS_AR[2] <- UNIFY_CONS_AR[1] + 1e-6
  print(p_qin_d1_f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
}
unique_sources_qin_d1 <- unique(qin_day1_liana_result_df$source)
neutrophil_targets_qin_d1 <- NEUTROPHIL_STATES[NEUTROPHIL_STATES %in% unique(qin_day1_liana_result_df$target)]
dot_sources_qin_d1 <- setdiff(unique_sources_qin_d1, neutrophil_targets_qin_d1)
if (length(dot_sources_qin_d1) == 0) dot_sources_qin_d1 <- unique_sources_qin_d1
p_qin_d1_dot <- NULL
if (length(neutrophil_targets_qin_d1) > 0 && length(dot_sources_qin_d1) > 0) p_qin_d1_dot <- tryCatch(liana::liana_dotplot(qin_day1_liana_result_df, source_groups = dot_sources_qin_d1, target_groups = neutrophil_targets_qin_d1, ntop = 20, size_range = LIANA_DOTPLOT_SIZE_RANGE), error = function(e) NULL)
if (!is.null(p_qin_d1_dot)) { p_qin_d1_dot <- p_qin_d1_dot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)); print(p_qin_d1_dot) }
receptor_freq_qin_d1 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(qin_day1_liana_consensus, receptor.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
p_qin_d1_h <- tryCatch(ggplot(receptor_freq_qin_d1, aes(x = reorder(receptor.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Receptors - Qin Day 1", x = "Receptor", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_qin_d1_h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
ligand_freq_qin_d1 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(qin_day1_liana_consensus, ligand.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
p_qin_d1_i <- tryCatch(ggplot(ligand_freq_qin_d1, aes(x = reorder(ligand.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Ligands - Qin Day 1", x = "Ligand", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_qin_d1_i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
interaction_matrix_qin_d1 <- tryCatch(tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(qin_day1_liana_consensus, source, target), interaction_count = dplyr::n(), .groups = "drop"), names_from = target, values_from = interaction_count, values_fill = 0), error = function(e) data.frame())
interaction_matrix_qin_d1_mat <- tryCatch({ cols_num <- setdiff(colnames(interaction_matrix_qin_d1), "source"); mat <- as.matrix(interaction_matrix_qin_d1[, cols_num, drop = FALSE]); rownames(mat) <- interaction_matrix_qin_d1$source; mat }, error = function(e) matrix(0, nrow = 0, ncol = 0))
p_qin_d1_j <- tryCatch(pheatmap::pheatmap(interaction_matrix_qin_d1_mat, color = colorRampPalette(c("white", "yellow", "orange", "red"))(100), main = "Cell Type Interaction Frequency - Qin Day 1", display_numbers = TRUE), error = function(e) NULL)
if (!is.null(p_qin_d1_j)) { print(p_qin_d1_j); grid::grid.newpage(); grid::grid.draw(p_qin_d1_j$gtable) }
liana_trunc_qin_d1 <- dplyr::filter(qin_day1_liana_result_df, aggregate_rank <= 0.01)
if (nrow(liana_trunc_qin_d1) > 0) {
  liana_trunc_qin_d1$source <- as.character(liana_trunc_qin_d1$source)
  liana_trunc_qin_d1$target <- as.character(liana_trunc_qin_d1$target)
  p_heat_freq_qin_d1 <- tryCatch(liana::heat_freq(liana_trunc_qin_d1), error = function(e) { message("heat_freq Qin D1: ", conditionMessage(e)); NULL })
  if (!is.null(p_heat_freq_qin_d1)) print(p_heat_freq_qin_d1)
  unique_sources_chord_qin_d1 <- unique(liana_trunc_qin_d1$source)
  unique_targets_chord_qin_d1 <- unique(liana_trunc_qin_d1$target)
  grDevices::png(file.path(OUTPUT_DIR, "QinDay1_LIANA_ChordFreq.png"), width = 1400, height = 1400, res = 150)
  tryCatch(liana::chord_freq(liana_trunc_qin_d1, source_groups = unique_sources_chord_qin_d1, target_groups = unique_targets_chord_qin_d1), error = function(e) message("chord_freq Qin D1 (PNG): ", conditionMessage(e)))
  grDevices::dev.off()
  p_chord_freq_qin_d1 <- tryCatch(liana::chord_freq(liana_trunc_qin_d1, source_groups = unique_sources_chord_qin_d1, target_groups = unique_targets_chord_qin_d1), error = function(e) { message("chord_freq Qin D1: ", conditionMessage(e)); NULL })
  if (!is.null(p_chord_freq_qin_d1)) print(p_chord_freq_qin_d1)
  liana_mat_qin_d1 <- as.matrix(table(liana_trunc_qin_d1$source, liana_trunc_qin_d1$target))
  p_liana_heatmap_qin_d1 <- NULL
  if (nrow(liana_mat_qin_d1) > 0 && ncol(liana_mat_qin_d1) > 0) p_liana_heatmap_qin_d1 <- tryCatch(liana::liana_heatmap(liana_mat_qin_d1), error = function(e) { message("liana_heatmap Qin D1: ", conditionMessage(e)); NULL })
  if (!is.null(p_liana_heatmap_qin_d1)) ComplexHeatmap::draw(p_liana_heatmap_qin_d1)
}
source_importance_qin_d1 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(qin_day1_liana_consensus, source), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), importance_score = dplyr::n() * (1 - mean(aggregate_rank))), desc(importance_score)), error = function(e) data.frame())
p_qin_d1_k <- tryCatch(ggplot(source_importance_qin_d1, aes(x = reorder(source, importance_score), y = importance_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkred") + labs(title = "Source Cell Importance - Qin Day 1", x = "Cell Type", y = "Importance Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_iy <- numeric(0)
vals_imr <- numeric(0)
imp_tabnames <- c("source_importance_lee_d1", "source_importance_lee_d3", "source_importance_wang_d3", "source_importance_qin_d1", "source_importance_qin_d3")
for (fn in imp_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("importance_score" %in% names(d)) vals_iy <- c(vals_iy, d$importance_score)
  if ("mean_rank" %in% names(d)) vals_imr <- c(vals_imr, d$mean_rank)
}
vals_iy <- vals_iy[is.finite(vals_iy)]
vals_imr <- vals_imr[is.finite(vals_imr)]
UNIFY_SRC_IMP_Y <- if (length(vals_iy) > 0) range(vals_iy) else c(0, 1)
UNIFY_SRC_IMP_MR <- if (length(vals_imr) > 0) range(vals_imr) else c(0, 1)
if (UNIFY_SRC_IMP_Y[1] == UNIFY_SRC_IMP_Y[2]) UNIFY_SRC_IMP_Y[2] <- UNIFY_SRC_IMP_Y[1] + 1e-6
if (UNIFY_SRC_IMP_MR[1] == UNIFY_SRC_IMP_MR[2]) UNIFY_SRC_IMP_MR[2] <- UNIFY_SRC_IMP_MR[1] + 1e-6
print(p_qin_d1_k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
p_qin_d1_n <- tryCatch(ggplot(head(qin_day1_liana_arg1pos_topranked, 15), aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), aggregate_rank), y = -aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top Arg1pos Interactions - Qin Day 1", x = "L-R Pair", y = "Rank Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_a1 <- numeric(0)
if (exists("lee_day1_liana_arg1pos_topranked") && is.data.frame(lee_day1_liana_arg1pos_topranked) && nrow(lee_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("lee_day3_liana_arg1pos_topranked") && is.data.frame(lee_day3_liana_arg1pos_topranked) && nrow(lee_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("wang_day3_liana_arg1pos_topranked") && is.data.frame(wang_day3_liana_arg1pos_topranked) && nrow(wang_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(wang_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day1_liana_arg1pos_topranked") && is.data.frame(qin_day1_liana_arg1pos_topranked) && nrow(qin_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day3_liana_arg1pos_topranked") && is.data.frame(qin_day3_liana_arg1pos_topranked) && nrow(qin_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day3_liana_arg1pos_topranked$aggregate_rank, 15))
vals_a1 <- vals_a1[is.finite(vals_a1)]
UNIFY_ARG1_AR <- if (length(vals_a1) > 0) range(vals_a1) else c(0, 1)
if (UNIFY_ARG1_AR[1] == UNIFY_ARG1_AR[2]) UNIFY_ARG1_AR[2] <- UNIFY_ARG1_AR[1] + 1e-6
UNIFY_ARG1_NEGY <- c(-UNIFY_ARG1_AR[2], -UNIFY_ARG1_AR[1])
print(p_qin_d1_n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
qin_day1_arg1pos_specific_plot <- head(qin_day1_arg1pos_specific, 15)
p_qin_d1_o2 <- tryCatch(ggplot(qin_day1_arg1pos_specific_plot, aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), specificity_index), y = specificity_index, fill = source)) + geom_bar(stat = "identity") + coord_flip() + labs(title = "Arg1+-Specific Signals (Qin Day 1)", x = "L-R Pair", y = "Specificity Index") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
print(p_qin_d1_o2)
target_specificity_qin_d1 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(dplyr::filter(qin_day1_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")), target), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), specificity_score = dplyr::n() * (1 - mean(aggregate_rank)), .groups = "drop"), desc(specificity_score)), error = function(e) data.frame())
p_qin_d1_l <- tryCatch(ggplot(target_specificity_qin_d1, aes(x = reorder(target, specificity_score), y = specificity_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkgreen") + labs(title = "Target Cell Specificity - Qin Day 1", x = "Cell Type", y = "Specificity Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_sy <- numeric(0)
vals_smr <- numeric(0)
spec_tabnames <- c("target_specificity_lee_d1", "target_specificity_lee_d3", "target_specificity_wang_d3", "target_specificity_qin_d1", "target_specificity_qin_d3")
for (fn in spec_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("specificity_score" %in% names(d)) vals_sy <- c(vals_sy, d$specificity_score)
  if ("mean_rank" %in% names(d)) vals_smr <- c(vals_smr, d$mean_rank)
}
vals_sy <- vals_sy[is.finite(vals_sy)]
vals_smr <- vals_smr[is.finite(vals_smr)]
UNIFY_TGT_SPEC_Y <- if (length(vals_sy) > 0) range(vals_sy) else c(0, 1)
UNIFY_TGT_SPEC_MR <- if (length(vals_smr) > 0) range(vals_smr) else c(0, 1)
if (UNIFY_TGT_SPEC_Y[1] == UNIFY_TGT_SPEC_Y[2]) UNIFY_TGT_SPEC_Y[2] <- UNIFY_TGT_SPEC_Y[1] + 1e-6
if (UNIFY_TGT_SPEC_MR[1] == UNIFY_TGT_SPEC_MR[2]) UNIFY_TGT_SPEC_MR[2] <- UNIFY_TGT_SPEC_MR[1] + 1e-6
print(p_qin_d1_l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
p_qin_d1_m <- tryCatch(ggplot(qin_day1_liana_consensus, aes(x = aggregate_rank)) + geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) + labs(title = "Distribution of L-R Pair Aggregate Ranks - Qin Day 1", x = "Aggregate Rank Score", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_hist <- numeric(0)
if (exists("lee_day1_liana_consensus") && is.data.frame(lee_day1_liana_consensus) && nrow(lee_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day1_liana_consensus)) vals_hist <- c(vals_hist, lee_day1_liana_consensus$aggregate_rank)
if (exists("lee_day3_liana_consensus") && is.data.frame(lee_day3_liana_consensus) && nrow(lee_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day3_liana_consensus)) vals_hist <- c(vals_hist, lee_day3_liana_consensus$aggregate_rank)
if (exists("wang_day3_liana_consensus") && is.data.frame(wang_day3_liana_consensus) && nrow(wang_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(wang_day3_liana_consensus)) vals_hist <- c(vals_hist, wang_day3_liana_consensus$aggregate_rank)
if (exists("qin_day1_liana_consensus") && is.data.frame(qin_day1_liana_consensus) && nrow(qin_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day1_liana_consensus)) vals_hist <- c(vals_hist, qin_day1_liana_consensus$aggregate_rank)
if (exists("qin_day3_liana_consensus") && is.data.frame(qin_day3_liana_consensus) && nrow(qin_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day3_liana_consensus)) vals_hist <- c(vals_hist, qin_day3_liana_consensus$aggregate_rank)
vals_hist <- vals_hist[is.finite(vals_hist)]
UNIFY_HIST_AR <- if (length(vals_hist) > 0) range(vals_hist) else c(0, 1)
if (UNIFY_HIST_AR[1] == UNIFY_HIST_AR[2]) UNIFY_HIST_AR[2] <- UNIFY_HIST_AR[1] + 1e-6
print(p_qin_d1_m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
method_cols_qin_d1 <- colnames(qin_day1_liana_result_df)[grep("pval_", colnames(qin_day1_liana_result_df))]
consensus_qin_d1_step1 <- dplyr::slice_head(qin_day1_liana_consensus, n = 15)
consensus_qin_d1 <- NULL
consensus_qin_d1_has_cols <- length(method_cols_qin_d1) > 0 && nrow(consensus_qin_d1_step1) > 0
if (consensus_qin_d1_has_cols) consensus_qin_d1 <- as.data.frame(dplyr::select(consensus_qin_d1_step1, dplyr::all_of(method_cols_qin_d1)))
if (consensus_qin_d1_has_cols) consensus_qin_d1_rownames <- paste0(consensus_qin_d1_step1$source, " | ", consensus_qin_d1_step1$ligand.complex, " -> ", consensus_qin_d1_step1$receptor.complex)
if (consensus_qin_d1_has_cols) rownames(consensus_qin_d1) <- make.unique(as.character(consensus_qin_d1_rownames))
consensus_qin_d1_ready <- !is.null(consensus_qin_d1) && nrow(consensus_qin_d1) > 0
if (consensus_qin_d1_ready) { p_qin_d1_o <- pheatmap::pheatmap(consensus_qin_d1, color = colorRampPalette(c("red", "white", "blue"))(100), main = "Method Consensus (p-values) - Qin Day 1"); print(p_qin_d1_o) }


qin_day1_neut_cells_vln <- colnames(qin_day1)[c(qin_day1_pos, qin_day1_neg)]
top_receptors_qin_d1 <- head(unique(qin_day1_arg1pos_received$receptor.complex), LIANA_TOP_RECEPTOR_VLN)
top_receptors_qin_d1_single <- top_receptors_qin_d1[!grepl("[_+]", top_receptors_qin_d1)]
top_receptors_qin_d1_in_data <- top_receptors_qin_d1_single[top_receptors_qin_d1_single %in% rownames(qin_day1)]
if (length(top_receptors_qin_d1_in_data) > 0) {
  qin_day1_neut_obj_vln <- subset(qin_day1, cells = qin_day1_neut_cells_vln)
  p_receptor_vln_qin_d1 <- Seurat::VlnPlot(qin_day1_neut_obj_vln, features = top_receptors_qin_d1_in_data, group.by = "arg1_status", pt.size = 0.1, ncol = min(3L, length(top_receptors_qin_d1_in_data)))
  print(p_receptor_vln_qin_d1)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "QinDay1_TopLIANA_Receptors_VlnPlot.png"), p_receptor_vln_qin_d1, width = 10, height = 6, dpi = 150)
}
print("✓ Qin Day 1 LIANA analysis complete")
# Checkpoint: resume from here to run Qin Day 3 LIANA. load(file.path(OUTPUT_DIR, "Workspace_AfterQinDay1_LIANA.RData")) then run from Qin Day 3 section.
save.image(file.path(OUTPUT_DIR, "Workspace_AfterQinDay1_LIANA.RData"))
saveRDS(list(qin_day1_liana_result_df = qin_day1_liana_result_df, qin_day3 = qin_day3, qin_day3_pruned_labels = qin_day3_pruned_labels, qin_day3_pos = qin_day3_pos, qin_day3_neg = qin_day3_neg, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay1_LIANA.rds"))
# checkpoint_qd1 <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay1_LIANA.rds")); list2env(checkpoint_qd1, envir = .GlobalEnv)
print("✓ Checkpoint saved: Workspace_AfterQinDay1_LIANA.RData, Checkpoint_AfterQinDay1_LIANA.rds")

# -------- Qin Day 1 PROGENy Pathway Analysis (Exploratory; downstream pathway support only) --------
print("--- PROGENy Pathway Analysis: Qin Day 1 (Exploratory) ---")
# Secondary support only: PROGENy summarizes downstream pathway state in Arg1pos vs Arg1neg neutrophils.
progeny_model_mouse_qin_d1 <- progeny::model_mouse_full
colnames(progeny_model_mouse_qin_d1) <- tolower(colnames(progeny_model_mouse_qin_d1))
colnames(progeny_model_mouse_qin_d1)[colnames(progeny_model_mouse_qin_d1) == "p.value"] <- "p_value"
stopifnot(all(c("gene", "pathway", "weight") %in% colnames(progeny_model_mouse_qin_d1)))
progeny_network_qin_d1 <- data.frame(source = progeny_model_mouse_qin_d1$pathway, target = progeny_model_mouse_qin_d1$gene, weight = progeny_model_mouse_qin_d1$weight, stringsAsFactors = FALSE)
progeny_network_qin_d1 <- progeny_network_qin_d1[progeny_network_qin_d1$weight != 0, ]
qin_day1_neut_cells <- colnames(qin_day1)[c(qin_day1_pos, qin_day1_neg)]
qin_day1_neut_expr_mat <- tryCatch(Seurat::GetAssayData(qin_day1, layer = "data")[, qin_day1_neut_cells, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day1_progeny_result <- tryCatch(decoupleR::run_wmean(mat = qin_day1_neut_expr_mat, network = progeny_network_qin_d1, .source = "source", .target = "target", .mor = "weight", minsize = 5), error = function(e) data.frame())
qin_day1_progeny_for_wide <- if (nrow(qin_day1_progeny_result) > 0 && "statistic" %in% colnames(qin_day1_progeny_result)) dplyr::filter(qin_day1_progeny_result, statistic == "norm_wmean") else qin_day1_progeny_result
if (nrow(qin_day1_progeny_for_wide) == 0 && nrow(qin_day1_progeny_result) > 0 && "statistic" %in% colnames(qin_day1_progeny_result)) qin_day1_progeny_for_wide <- dplyr::filter(qin_day1_progeny_result, statistic == "wmean")
qin_day1_pw_wide <- tryCatch(tidyr::pivot_wider(qin_day1_progeny_for_wide, names_from = "condition", values_from = "score", id_cols = "source"), error = function(e) data.frame())
qin_day1_pw_cols_num <- tryCatch(sapply(qin_day1_pw_wide[, -1, drop = FALSE], function(x) as.numeric(unlist(x))), error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day1_progeny_scores_mat <- tryCatch(as.matrix(qin_day1_pw_cols_num), error = function(e) matrix(0, nrow = 0, ncol = 0))
rownames(qin_day1_progeny_scores_mat) <- tryCatch(as.character(qin_day1_pw_wide$source), error = function(e) character(0))
colnames(qin_day1_progeny_scores_mat) <- tryCatch(colnames(qin_day1_pw_wide)[-1], error = function(e) character(0))
qin_day1_arg1pos_cells <- colnames(qin_day1)[qin_day1_pos]
qin_day1_arg1neg_cells <- colnames(qin_day1)[qin_day1_neg]
qin_day1_arg1pos_cells_in_mat <- qin_day1_arg1pos_cells[qin_day1_arg1pos_cells %in% colnames(qin_day1_progeny_scores_mat)]
qin_day1_arg1neg_cells_in_mat <- qin_day1_arg1neg_cells[qin_day1_arg1neg_cells %in% colnames(qin_day1_progeny_scores_mat)]
qin_day1_progeny_arg1pos_scores <- tryCatch(qin_day1_progeny_scores_mat[, qin_day1_arg1pos_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day1_progeny_arg1neg_scores <- tryCatch(qin_day1_progeny_scores_mat[, qin_day1_arg1neg_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day1_pathway_means_arg1pos <- tryCatch(rowMeans(qin_day1_progeny_arg1pos_scores, na.rm = TRUE), error = function(e) numeric(0))
qin_day1_pathway_means_arg1neg <- tryCatch(rowMeans(qin_day1_progeny_arg1neg_scores, na.rm = TRUE), error = function(e) numeric(0))
qin_day1_pathway_medians_arg1pos <- tryCatch(apply(qin_day1_progeny_arg1pos_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
qin_day1_pathway_medians_arg1neg <- tryCatch(apply(qin_day1_progeny_arg1neg_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
qin_day1_pathway_shift <- pmax(0, -pmin(qin_day1_pathway_means_arg1pos, qin_day1_pathway_means_arg1neg, na.rm = TRUE)) + 1e-10
qin_day1_pathway_log2fc <- log2((qin_day1_pathway_means_arg1pos + qin_day1_pathway_shift) / (qin_day1_pathway_means_arg1neg + qin_day1_pathway_shift))
qin_day1_pathway_comparison <- tryCatch(data.frame(
  pathway = names(qin_day1_pathway_means_arg1pos),
  Arg1pos_mean = qin_day1_pathway_means_arg1pos,
  Arg1neg_mean = qin_day1_pathway_means_arg1neg,
  Arg1pos_median = qin_day1_pathway_medians_arg1pos,
  Arg1neg_median = qin_day1_pathway_medians_arg1neg,
  mean_difference = qin_day1_pathway_means_arg1pos - qin_day1_pathway_means_arg1neg,
  median_difference = qin_day1_pathway_medians_arg1pos - qin_day1_pathway_medians_arg1neg,
  log2FC = qin_day1_pathway_log2fc,
  stringsAsFactors = FALSE
), error = function(e) data.frame())
qin_day1_all_pathways <- tryCatch(unique(progeny_network_qin_d1$source), error = function(e) character(0))
qin_day1_pathway_results <- qin_day1_pathway_comparison
qin_day1_arg1pos_count <- length(qin_day1_pos)
qin_day1_arg1neg_count <- length(qin_day1_neg)
print(paste("Qin Day 1: Arg1pos neutrophils =", qin_day1_arg1pos_count, ", Arg1neg neutrophils =", qin_day1_arg1neg_count))
qin_day1_pathway_pvals <- numeric(length(qin_day1_all_pathways))
names(qin_day1_pathway_pvals) <- qin_day1_all_pathways
qin_day1_pathway_effect_sizes <- numeric(length(qin_day1_all_pathways))
names(qin_day1_pathway_effect_sizes) <- qin_day1_all_pathways
qin_day1_pathway_ci_lower <- numeric(length(qin_day1_all_pathways))
names(qin_day1_pathway_ci_lower) <- qin_day1_all_pathways
qin_day1_pathway_ci_upper <- numeric(length(qin_day1_all_pathways))
names(qin_day1_pathway_ci_upper) <- qin_day1_all_pathways
qin_day1_pathway_adequate_n <- logical(length(qin_day1_all_pathways))
names(qin_day1_pathway_adequate_n) <- qin_day1_all_pathways
for (pw_idx in seq_along(qin_day1_all_pathways)) {
  pw <- qin_day1_all_pathways[pw_idx]
  arg1pos_scores_pw <- tryCatch(as.numeric(qin_day1_progeny_arg1pos_scores[pw, ]), error = function(e) numeric(0))
  arg1neg_scores_pw <- tryCatch(as.numeric(qin_day1_progeny_arg1neg_scores[pw, ]), error = function(e) numeric(0))
  arg1pos_scores_pw <- arg1pos_scores_pw[!is.na(arg1pos_scores_pw)]
  arg1neg_scores_pw <- arg1neg_scores_pw[!is.na(arg1neg_scores_pw)]
  n_arg1pos <- length(arg1pos_scores_pw)
  n_arg1neg <- length(arg1neg_scores_pw)
  qin_day1_pathway_adequate_n[pw] <- (n_arg1pos >= 5) & (n_arg1neg >= 5)
  median_arg1pos <- tryCatch(stats::median(arg1pos_scores_pw, na.rm = TRUE), error = function(e) 0)
  median_arg1neg <- tryCatch(stats::median(arg1neg_scores_pw, na.rm = TRUE), error = function(e) 0)
  qin_day1_pathway_effect_sizes[pw] <- median_arg1pos - median_arg1neg
  qin_day1_pathway_pvals[pw] <- 1.0
  qin_day1_pathway_ci_lower[pw] <- NA_real_
  qin_day1_pathway_ci_upper[pw] <- NA_real_
  wilcox_result <- tryCatch(stats::wilcox.test(arg1pos_scores_pw, arg1neg_scores_pw, conf.int = TRUE, conf.level = 0.95), error = function(e) NULL)
  wilcox_pvalue <- tryCatch(if (!is.null(wilcox_result)) wilcox_result$p.value else 1.0, error = function(e) 1.0)
  wilcox_pvalue_length <- tryCatch(length(wilcox_pvalue), error = function(e) 0)
  wilcox_pvalue_final <- tryCatch(if (wilcox_pvalue_length > 0) wilcox_pvalue[1] else 1.0, error = function(e) 1.0)
  qin_day1_pathway_pvals[pw] <- wilcox_pvalue_final
  wilcox_ci_lower <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[1] else NA_real_, error = function(e) NA_real_)
  wilcox_ci_upper <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[2] else NA_real_, error = function(e) NA_real_)
  qin_day1_pathway_ci_lower[pw] <- wilcox_ci_lower
  qin_day1_pathway_ci_upper[pw] <- wilcox_ci_upper
  qin_day1_pathway_warning_msg <- tryCatch(paste("Warning: Pathway", pw, "has insufficient sample size (Arg1pos n =", n_arg1pos, ", Arg1neg n =", n_arg1neg, "). Skipping statistical test."), error = function(e) "")
  qin_day1_pathway_warning_vector <- c("", qin_day1_pathway_warning_msg)
  qin_day1_pathway_warning_index <- tryCatch(as.numeric(!qin_day1_pathway_adequate_n[pw]) + 1, error = function(e) 1)
  print(qin_day1_pathway_warning_vector[qin_day1_pathway_warning_index])
}
qin_day1_pathway_pvals_adj <- p.adjust(qin_day1_pathway_pvals, method = "BH")
qin_day1_pathway_results$p_value <- tryCatch(qin_day1_pathway_pvals[qin_day1_pathway_results$pathway], error = function(e) rep(1.0, nrow(qin_day1_pathway_results)))
qin_day1_pathway_results$p_adj <- tryCatch(qin_day1_pathway_pvals_adj[qin_day1_pathway_results$pathway], error = function(e) rep(1.0, nrow(qin_day1_pathway_results)))
qin_day1_pathway_results$effect_size_median_diff <- tryCatch(qin_day1_pathway_effect_sizes[qin_day1_pathway_results$pathway], error = function(e) rep(0.0, nrow(qin_day1_pathway_results)))
qin_day1_pathway_results$ci_lower_95 <- tryCatch(qin_day1_pathway_ci_lower[qin_day1_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(qin_day1_pathway_results)))
qin_day1_pathway_results$ci_upper_95 <- tryCatch(qin_day1_pathway_ci_upper[qin_day1_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(qin_day1_pathway_results)))
qin_day1_pathway_results$adequate_sample_size <- tryCatch(qin_day1_pathway_adequate_n[qin_day1_pathway_results$pathway], error = function(e) rep(FALSE, nrow(qin_day1_pathway_results)))
qin_day1_pathway_results$significant <- tryCatch((qin_day1_pathway_results$p_adj < 0.05) & qin_day1_pathway_results$adequate_sample_size, error = function(e) rep(FALSE, nrow(qin_day1_pathway_results)))
write.csv(qin_day1_pathway_comparison, file.path(OUTPUT_DIR, "QinDay1_PROGENy_PathwayComparison.csv"), row.names = FALSE)
write.csv(qin_day1_pathway_results, file.path(OUTPUT_DIR, "QinDay1_PROGENy_PathwayResults.csv"), row.names = FALSE)
saveRDS(list(pathway_comparison = qin_day1_pathway_comparison, pathway_results = qin_day1_pathway_results, pathway_means_arg1pos = qin_day1_pathway_means_arg1pos, pathway_means_arg1neg = qin_day1_pathway_means_arg1neg, progeny_network = progeny_network_qin_d1), file.path(OUTPUT_DIR, "QinDay1_PROGENy_Results.rds"))
qin_day1_identified_ligands <- tryCatch(unique(c(qin_day1_liana_arg1pos_topranked$ligand.complex, qin_day1_liana_arg1neg_topranked$ligand.complex)), error = function(e) character(0))
qin_day1_ligand_pathway_map <- data.frame(ligand = character(0), pathway = character(0), weight = numeric(0), stringsAsFactors = FALSE)
for (lig in qin_day1_identified_ligands) { lig_genes <- trimws(unlist(strsplit(lig, "[_+]"))); lig_pathways <- dplyr::filter(progeny_network_qin_d1, target %in% lig_genes); if (nrow(lig_pathways) > 0) qin_day1_ligand_pathway_map <- rbind(qin_day1_ligand_pathway_map, data.frame(ligand = lig, pathway = lig_pathways$source, weight = lig_pathways$weight, stringsAsFactors = FALSE)) }
qin_day1_pathway_long <- tryCatch(tidyr::pivot_longer(qin_day1_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
p_qin_d1_progeny1 <- tryCatch(ggplot(qin_day1_pathway_long, aes(x = pathway, y = Pathway_Score, fill = Group)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("Arg1pos_mean" = "red", "Arg1neg_mean" = "lightblue"), labels = c("Arg1pos", "Arg1neg")) + labs(title = "PROGENy Pathway Activity - Qin Day 1", x = "Pathway", y = "Pathway Activity Score") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
vals_pw <- numeric(0)
if (exists("lee_day1_pathway_long") && is.data.frame(lee_day1_pathway_long) && nrow(lee_day1_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day1_pathway_long)) vals_pw <- c(vals_pw, lee_day1_pathway_long$Pathway_Score)
if (exists("lee_day3_pathway_long") && is.data.frame(lee_day3_pathway_long) && nrow(lee_day3_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day3_pathway_long)) vals_pw <- c(vals_pw, lee_day3_pathway_long$Pathway_Score)
if (exists("wang_day3_pathway_long") && is.data.frame(wang_day3_pathway_long) && nrow(wang_day3_pathway_long) > 0 && "Pathway_Score" %in% names(wang_day3_pathway_long)) vals_pw <- c(vals_pw, wang_day3_pathway_long$Pathway_Score)
if (exists("qin_day1_pathway_long") && is.data.frame(qin_day1_pathway_long) && nrow(qin_day1_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day1_pathway_long)) vals_pw <- c(vals_pw, qin_day1_pathway_long$Pathway_Score)
if (exists("qin_day3_pathway_long") && is.data.frame(qin_day3_pathway_long) && nrow(qin_day3_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day3_pathway_long)) vals_pw <- c(vals_pw, qin_day3_pathway_long$Pathway_Score)
vals_pw <- vals_pw[is.finite(vals_pw)]
UNIFY_PW_Y <- if (length(vals_pw) > 0) range(vals_pw) else c(-1, 1)
if (UNIFY_PW_Y[1] == UNIFY_PW_Y[2]) UNIFY_PW_Y[2] <- UNIFY_PW_Y[1] + 1e-6
print(p_qin_d1_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
top_pw_qin_d1 <- character(0)
if (nrow(qin_day1_pathway_comparison) > 0 && "mean_difference" %in% colnames(qin_day1_pathway_comparison)) {
  ord_pw_qin_d1 <- order(-abs(qin_day1_pathway_comparison$mean_difference))
  take_pw_qin_d1 <- min(6L, length(ord_pw_qin_d1))
  top_pw_qin_d1 <- as.character(qin_day1_pathway_comparison$pathway[ord_pw_qin_d1[seq_len(take_pw_qin_d1)]])
}
if (length(top_pw_qin_d1) > 0 && nrow(qin_day1_progeny_scores_mat) > 0) {
  n_cells_pw_q1 <- ncol(qin_day1_progeny_scores_mat)
  n_pw_sel_q1 <- length(top_pw_qin_d1)
  cell_vec_pw_q1 <- rep(colnames(qin_day1_progeny_scores_mat), each = n_pw_sel_q1)
  pathway_vec_pw_q1 <- rep(top_pw_qin_d1, times = n_cells_pw_q1)
  score_vec_pw_q1 <- as.vector(t(qin_day1_progeny_scores_mat[top_pw_qin_d1, , drop = FALSE]))
  pw_scores_long_qin_d1 <- data.frame(cell = cell_vec_pw_q1, pathway = pathway_vec_pw_q1, score = score_vec_pw_q1, stringsAsFactors = FALSE)
  pw_scores_long_qin_d1$arg1_status <- ifelse(pw_scores_long_qin_d1$cell %in% qin_day1_arg1pos_cells, "Arg1pos", "Arg1neg")
  p_progeny_vln_qin_d1 <- ggplot(pw_scores_long_qin_d1, aes(x = arg1_status, y = score, fill = arg1_status)) + geom_violin(trim = FALSE) + geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) + facet_wrap(~pathway, scales = "free_y") + labs(title = "PROGENy pathway scores: Arg1+ vs Arg1- neutrophils (Qin Day 1)", x = "ARG1 status", y = "Activity score (norm_wmean)") + theme_minimal() + PLOT_TITLE_THEME
  print(p_progeny_vln_qin_d1)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "QinDay1_PROGENy_TopPathways_Violin.png"), p_progeny_vln_qin_d1, width = 12, height = 8, dpi = 150)
}
progeny_heatmap_mat_qin_d1 <- tryCatch(rbind(Arg1pos = qin_day1_pathway_means_arg1pos, Arg1neg = qin_day1_pathway_means_arg1neg), error = function(e) matrix(0, nrow = 0, ncol = 0))
if (nrow(progeny_heatmap_mat_qin_d1) > 0 && ncol(progeny_heatmap_mat_qin_d1) > 0) {
  colors_progeny_qin_d1 <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors_use_progeny_qin_d1 <- grDevices::colorRampPalette(colors = colors_progeny_qin_d1)(100)
  p_qin_d1_progeny_heat <- pheatmap::pheatmap(progeny_heatmap_mat_qin_d1, color = colors_use_progeny_qin_d1, border_color = "white", cellwidth = 20, cellheight = 20, main = "PROGENy Pathway Activity: Arg1+ vs Arg1- Neutrophils (Qin Day 1)")
  print(p_qin_d1_progeny_heat)
}
qin_day1_ligand_map_plot <- if (nrow(qin_day1_ligand_pathway_map) > 0 && "ligand" %in% names(qin_day1_ligand_pathway_map)) { top15_lig_qin_d1 <- head(unique(qin_day1_ligand_pathway_map$ligand), 15); qin_day1_ligand_pathway_map[qin_day1_ligand_pathway_map$ligand %in% top15_lig_qin_d1, ] } else data.frame()
p_qin_d1_progeny2 <- if (nrow(qin_day1_ligand_map_plot) > 0) tryCatch(ggplot(qin_day1_ligand_map_plot, aes(x = ligand, y = pathway, size = abs(weight), color = weight)) + geom_point(alpha = 0.7) + scale_color_gradient2(low = "blue", mid = "white", high = "red") + labs(title = "Top 15 Ligands Linked to PROGENy Pathways - Qin Day 1", x = "Ligand", y = "Pathway") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available")) else ggplot() + theme_void() + labs(title = "No ligands from LIANA (Qin Day 1)")
vals_w <- numeric(0)
map_tabnames <- c("lee_day1_ligand_pathway_map", "lee_day3_ligand_pathway_map", "wang_day3_ligand_pathway_map", "qin_day1_ligand_pathway_map", "qin_day3_ligand_pathway_map")
for (fn in map_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("weight" %in% names(d)) vals_w <- c(vals_w, d$weight)
}
vals_w <- vals_w[is.finite(vals_w)]
UNIFY_LIGW_ABS <- if (length(vals_w) > 0) range(abs(vals_w)) else c(0, 1)
if (UNIFY_LIGW_ABS[1] == UNIFY_LIGW_ABS[2]) UNIFY_LIGW_ABS[2] <- UNIFY_LIGW_ABS[1] + 1e-6
mxw <- if (length(vals_w) > 0) max(abs(vals_w)) else 1
if (!is.finite(mxw) || mxw <= 0) mxw <- 1
UNIFY_LIGW_COL <- c(-mxw, mxw)
print(p_qin_d1_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))
print("✓ Qin Day 1 PROGENy pathway analysis complete")
saveRDS(list(qin_day1 = qin_day1, qin_day1_neut_cells = qin_day1_neut_cells, qin_day1_pos = qin_day1_pos, qin_day1_neg = qin_day1_neg, qin_day1_pathway_results = qin_day1_pathway_results, progeny_network_qin_d1 = progeny_network_qin_d1, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay1_PROGENy.rds"))
# ck_qin_d1_p <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay1_PROGENy.rds")); list2env(ck_qin_d1_p, envir = .GlobalEnv)
# FRESH RSTUDIO — Qin Day 1 BARTsc: load Checkpoint_AfterQinDay1_PROGENy.rds + libraries; BARTSC_* from Lee Day 1 block if missing.

# -------- Qin Day 1 BARTsc TF Analysis (Arg1pos vs Arg1neg; downstream TF support only) --------
if (bartsc_initialized) {
  bartsc_initialized <- suppressWarnings(tryCatch({ BARTsc::load_bart2(); TRUE }, error = function(e) FALSE))
  if (!bartsc_initialized && exists("bart2", envir = .GlobalEnv)) bartsc_initialized <- TRUE
}
print("--- BARTsc TF Analysis: Qin Day 1 (Arg1pos vs Arg1neg neutrophils) ---")
if (bartsc_initialized) {
  options("mc.cores" = min(6L, parallel::detectCores()))
  qin_day1_neut_subset <- subset(qin_day1, cells = qin_day1_neut_cells)
  qin_day1_bart_label <- setNames(ifelse(colnames(qin_day1_neut_subset) %in% colnames(qin_day1)[qin_day1_pos], "Arg1pos", "Arg1neg"), colnames(qin_day1_neut_subset))
  qin_day1_bart_label <- factor(qin_day1_bart_label, levels = c("Arg1pos", "Arg1neg"))
  qin_day1_bart_proj <- BARTsc::bartsc(name = "QinDay1_Arg1", genome = "mm10", label = qin_day1_bart_label, cell_types_used = c("Arg1pos", "Arg1neg"), RNA_cnt_matrix = Seurat::GetAssayData(qin_day1_neut_subset, layer = "counts"))
  qin_day1_bart_proj <- BARTsc::normalize_RNA(qin_day1_bart_proj)
  qin_day1_bart_proj <- BARTsc::find_signature_genes(qin_day1_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  qin_day1_bart_proj <- BARTsc::find_pairwise_deg(qin_day1_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  n_qin_d1_posneg <- 0L
  n_qin_d1_negpos <- 0L
  qin_d1_pw <- qin_day1_bart_proj@data$pairwise_DEG
  if (!is.null(qin_d1_pw) && "Arg1pos::Arg1neg" %in% names(qin_d1_pw)) n_qin_d1_posneg <- { x <- qin_d1_pw[["Arg1pos::Arg1neg"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  if (!is.null(qin_d1_pw) && "Arg1neg::Arg1pos" %in% names(qin_d1_pw)) n_qin_d1_negpos <- { x <- qin_d1_pw[["Arg1neg::Arg1pos"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  print(paste0("Qin Day 1 BARTsc pairwise DEG count (@data$pairwise_DEG): Arg1pos::Arg1neg=", n_qin_d1_posneg, ", Arg1neg::Arg1pos=", n_qin_d1_negpos))
  bart_deg_ok_qin_d1 <- if (BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT > 0) n_qin_d1_posneg >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT && n_qin_d1_negpos >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT else TRUE
  bart_nt_qin_d1 <- table(qin_day1_bart_label)
  bart_unequal_qin_d1 <- length(bart_nt_qin_d1) == 2L && length(unique(as.vector(bart_nt_qin_d1))) > 1L
  bart_crossct_allowed_qin_d1 <- !isTRUE(BARTSC_SKIP_CROSSCT_TEST) && length(bart_nt_qin_d1) >= 2L && all(as.integer(bart_nt_qin_d1) >= BARTSC_MIN_CELLS_CROSSCT) && !(isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_qin_d1) && bart_deg_ok_qin_d1
  qin_day1_bart_proj <- BARTsc::run_signature_RNA(qin_day1_bart_proj)
  qin_day1_bart_sig <- BARTsc::get_result(qin_day1_bart_proj, analysis = "cell type signature", mod = "RNA")
  qin_day1_bart_proj <- BARTsc::calc_crossCT_auc_RNA(qin_day1_bart_proj)
  print(bart_nt_qin_d1)
  if (!bart_crossct_allowed_qin_d1) message(paste0(
    "BARTsc Qin Day 1: crossCT_test / find_key_regulators skipped (BARTSC_SKIP_CROSSCT_TEST=", isTRUE(BARTSC_SKIP_CROSSCT_TEST),
    "; min_cells_ok=", all(as.integer(bart_nt_qin_d1) >= BARTSC_MIN_CELLS_CROSSCT),
    "; skip_if_unequal_n=", isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_qin_d1,
    "; deg_ok=", bart_deg_ok_qin_d1,
    " [Arg1pos::Arg1neg=", n_qin_d1_posneg,
    ", Arg1neg::Arg1pos=", n_qin_d1_negpos,
    ", min_each_dir=", BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT,
    "]). Vignette scRNA-seq.md §6: dot_plot/deviation_heatmap follow crossCT_test."
  ))
  if (bart_crossct_allowed_qin_d1) qin_day1_bart_proj <- BARTsc::crossCT_test(qin_day1_bart_proj, mod = "RNA")
  qin_day1_bart_cross <- if (bart_crossct_allowed_qin_d1) BARTsc::get_result(qin_day1_bart_proj, analysis = "cross-cell-type", mod = "RNA") else list()
  qin_day1_bart_cross_dev <- if (is.null(qin_day1_bart_cross)) list() else if (is.list(qin_day1_bart_cross) && "deviation" %in% names(qin_day1_bart_cross) && is.list(qin_day1_bart_cross$deviation)) qin_day1_bart_cross$deviation else qin_day1_bart_cross
  bart_outdir_qin_d1 <- file.path(OUTPUT_DIR, "BARTsc_QinDay1")
  dir.create(bart_outdir_qin_d1, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(qin_day1_bart_sig) && !is.null(qin_day1_bart_sig[["Arg1pos"]])) write.csv(qin_day1_bart_sig[["Arg1pos"]], file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_Signature_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(qin_day1_bart_sig) && !is.null(qin_day1_bart_sig[["Arg1neg"]])) write.csv(qin_day1_bart_sig[["Arg1neg"]], file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_Signature_Arg1neg.csv"), row.names = FALSE)
  qin_day1_bart_key <- NULL
  if (bart_crossct_allowed_qin_d1) {
    qin_day1_bart_proj <- BARTsc::find_key_regulators(qin_day1_bart_proj, mod = "RNA", min.N.profile = 3)
    qin_day1_bart_key <- BARTsc::get_result(qin_day1_bart_proj, analysis = "Key regs ident", mod = "RNA")
  }
  if (!is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1pos"]])) write.csv(qin_day1_bart_key[["Arg1pos"]], file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_KeyRegulators_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1neg"]])) write.csv(qin_day1_bart_key[["Arg1neg"]], file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_KeyRegulators_Arg1neg.csv"), row.names = FALSE)
  saveRDS(qin_day1_bart_proj, file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_Object.rds"))
  print("✓ Qin Day 1 BARTsc TF analysis complete — plots: Qin Day 1 BARTsc visualizations section below")
}

# -------- Qin Day 1: BARTsc visualizations (reload QinDay1_BARTsc_Object.rds; no re-run of bartsc pipeline) --------
bart_viz_qd1_rds <- file.path(OUTPUT_DIR, "BARTsc_QinDay1", "QinDay1_BARTsc_Object.rds")
if (!requireNamespace("BARTsc", quietly = TRUE)) message("Qin Day 1 BARTsc visualizations skipped: package BARTsc not installed.")
if (requireNamespace("BARTsc", quietly = TRUE) && !file.exists(bart_viz_qd1_rds)) message(paste0("Qin Day 1 BARTsc visualizations skipped: missing ", bart_viz_qd1_rds))
if (requireNamespace("BARTsc", quietly = TRUE) && file.exists(bart_viz_qd1_rds)) {
  if (!exists("BARTSC_N_TF_DOTHEAT")) BARTSC_N_TF_DOTHEAT <- 5L
  if (!exists("BARTSC_N_LABELED_TFS")) BARTSC_N_LABELED_TFS <- 5L
  if (!exists("BARTSC_KEYREG_SCATTER_W_IN")) BARTSC_KEYREG_SCATTER_W_IN <- 11
  if (!exists("BARTSC_KEYREG_SCATTER_H_IN")) BARTSC_KEYREG_SCATTER_H_IN <- 7
  if (!exists("BARTSC_KEYREG_XLIM")) BARTSC_KEYREG_XLIM <- c(-1, 1)
  if (!exists("BARTSC_KEYREG_YLIM_SIG")) BARTSC_KEYREG_YLIM_SIG <- c(0, 7)
  if (!exists("BARTSC_KEYREG_ZLIM_DMDR")) BARTSC_KEYREG_ZLIM_DMDR <- c(-2, 2)
  if (!exists("BARTSC_KEYREG_RANK_MAX")) BARTSC_KEYREG_RANK_MAX <- 60L
  suppressWarnings(tryCatch({ BARTsc::load_bart2(); NULL }, error = function(e) NULL))
  if (exists("bart2", envir = .GlobalEnv)) {
    if (!exists("types", envir = .GlobalEnv)) types <<- reticulate::import("types")
    bart2_mods_viz_qd1 <- reticulate::py_to_r(reticulate::py_get_attr(bart2, "__all__"))
    for (m_viz_qd1 in bart2_mods_viz_qd1) {
      if (!exists(m_viz_qd1, envir = .GlobalEnv)) assign(m_viz_qd1, reticulate::import(paste0("bart2.", m_viz_qd1), delay_load = TRUE), envir = .GlobalEnv)
    }
  }
  if (!exists("qin_day1_bart_proj", envir = .GlobalEnv, inherits = FALSE)) qin_day1_bart_proj <- readRDS(bart_viz_qd1_rds)
  bart_outdir_qin_d1 <- file.path(OUTPUT_DIR, "BARTsc_QinDay1")
  dir.create(bart_outdir_qin_d1, showWarnings = FALSE, recursive = TRUE)
  if (interactive() && grDevices::dev.cur() == 1L) grDevices::dev.new()
  qin_day1_bart_cross_viz <- BARTsc::get_result(qin_day1_bart_proj, analysis = "cross-cell-type", mod = "RNA")
  qin_day1_bart_cross_dev <- if (is.null(qin_day1_bart_cross_viz)) list() else if (is.list(qin_day1_bart_cross_viz) && "deviation" %in% names(qin_day1_bart_cross_viz) && is.list(qin_day1_bart_cross_viz$deviation)) qin_day1_bart_cross_viz$deviation else qin_day1_bart_cross_viz
  qin_day1_bart_key <- BARTsc::get_result(qin_day1_bart_proj, analysis = "Key regs ident", mod = "RNA")
  bart_tf_names_qin_d1 <- names(qin_day1_bart_cross_dev)
  tf_example_qin_d1 <- if (length(bart_tf_names_qin_d1) > 0) bart_tf_names_qin_d1[1] else NULL
  if (length(qin_day1_bart_cross_dev) > 0 && !is.null(tf_example_qin_d1)) {
    p_qin_d1_bart_dot <- BARTsc::dot_plot(qin_day1_bart_proj, mod = "RNA", tf = tf_example_qin_d1, max_dot_size = 22)
    print(p_qin_d1_bart_dot)
    ggplot2::ggsave(file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_DotPlot.png"), p_qin_d1_bart_dot, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_DotPlot.pdf"), p_qin_d1_bart_dot, width = 8, height = 6)
    p_qin_d1_bart_heat <- BARTsc::deviation_heatmap(qin_day1_bart_proj, mod = "RNA", tf = tf_example_qin_d1, tile_fontsize = 6)
    print(p_qin_d1_bart_heat)
    ggplot2::ggsave(file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_DeviationHeatmap.png"), p_qin_d1_bart_heat, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_DeviationHeatmap.pdf"), p_qin_d1_bart_heat, width = 8, height = 6)
  }
  top_tf_qin_d1_pos <- if (length(qin_day1_bart_cross_dev) > 0 && !is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1pos"]]) && "TF" %in% colnames(qin_day1_bart_key[["Arg1pos"]])) head(qin_day1_bart_key[["Arg1pos"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_tf_qin_d1_neg <- if (length(qin_day1_bart_cross_dev) > 0 && !is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1neg"]]) && "TF" %in% colnames(qin_day1_bart_key[["Arg1neg"]])) head(qin_day1_bart_key[["Arg1neg"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_key_tfs_qin_d1 <- intersect(unique(c(top_tf_qin_d1_pos, top_tf_qin_d1_neg)), names(qin_day1_bart_cross_dev))
  if (length(qin_day1_bart_cross_dev) > 0 && length(top_key_tfs_qin_d1) > 0) {
    i_tf_qd1 <- 1L
    while (i_tf_qd1 <= length(top_key_tfs_qin_d1)) {
      tf_cur_qin_d1 <- top_key_tfs_qin_d1[i_tf_qd1]
      p_dot_qd1 <- BARTsc::dot_plot(qin_day1_bart_proj, mod = "RNA", tf = tf_cur_qin_d1, max_dot_size = 22)
      print(p_dot_qd1)
      ggplot2::ggsave(file.path(bart_outdir_qin_d1, paste0("QinDay1_BARTsc_DotPlot_", tf_cur_qin_d1, ".png")), p_dot_qd1, width = 8, height = 6)
      p_heat_qd1 <- BARTsc::deviation_heatmap(qin_day1_bart_proj, mod = "RNA", tf = tf_cur_qin_d1, tile_fontsize = 6)
      print(p_heat_qd1)
      ggplot2::ggsave(file.path(bart_outdir_qin_d1, paste0("QinDay1_BARTsc_DeviationHeatmap_", tf_cur_qin_d1, ".png")), p_heat_qd1, width = 8, height = 6)
      i_tf_qd1 <- i_tf_qd1 + 1L
    }
  }
  tfs_qin_d1_pos <- character(0)
  if (!is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1pos"]])) {
    df_q1p <- qin_day1_bart_key[["Arg1pos"]]
    if (nrow(df_q1p) > 0) {
      tf_c_q1p <- intersect(colnames(df_q1p), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_q1p <- if (length(tf_c_q1p) == 0) colnames(df_q1p)[1] else tf_c_q1p[1]
      ord_q1p <- if ("final_rank" %in% colnames(df_q1p)) order(df_q1p$final_rank, na.last = TRUE) else seq_len(nrow(df_q1p))
      tfs_qin_d1_pos <- as.character(head(df_q1p[ord_q1p, , drop = FALSE][[tf_col_q1p]], BARTSC_N_LABELED_TFS))
    }
  }
  tfs_qin_d1_neg <- character(0)
  if (!is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1neg"]])) {
    df_q1n <- qin_day1_bart_key[["Arg1neg"]]
    if (nrow(df_q1n) > 0) {
      tf_c_q1n <- intersect(colnames(df_q1n), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_q1n <- if (length(tf_c_q1n) == 0) colnames(df_q1n)[1] else tf_c_q1n[1]
      ord_q1n <- if ("final_rank" %in% colnames(df_q1n)) order(df_q1n$final_rank, na.last = TRUE) else seq_len(nrow(df_q1n))
      tfs_qin_d1_neg <- as.character(head(df_q1n[ord_q1n, , drop = FALSE][[tf_col_q1n]], BARTSC_N_LABELED_TFS))
    }
  }
  if (length(tfs_qin_d1_pos) > 0) {
    grDevices::png(file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_KeyRegScatter_Official_Arg1pos.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(qin_day1_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_qin_d1_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_qin_d1_pos) > 0) message("Displaying official key_regulator_scatter for Arg1pos (Qin Day 1)...")
  if (interactive() && length(tfs_qin_d1_pos) > 0) key_regulator_scatter_unified(qin_day1_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_qin_d1_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
  if (length(tfs_qin_d1_neg) > 0) {
    grDevices::png(file.path(bart_outdir_qin_d1, "QinDay1_BARTsc_KeyRegScatter_Official_Arg1neg.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(qin_day1_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_qin_d1_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_qin_d1_neg) > 0) message("Displaying official key_regulator_scatter for Arg1neg (Qin Day 1)...")
  if (interactive() && length(tfs_qin_d1_neg) > 0) key_regulator_scatter_unified(qin_day1_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_qin_d1_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
  print("✓ Qin Day 1 BARTsc visualizations complete")
}

# ============================================================================
# INTEGRATION: LIANA + PROGENy footprint co-membership (exploratory) + OmniPath receptor→TF (Qin Day 1)
# ============================================================================
# Linear: footprint co-membership CSV → OmniPath receptor→TF → full-chain CSV → sender×pathway heatmap (BARTsc plots: section above).
print("--- Integration: LIANA + PROGENy footprint overlap (exploratory) + BARTsc/OmniPath (Qin Day 1) ---")
# Exploratory overlay only: LIANA provides incoming-signal candidates; PROGENy/BARTsc/OmniPath add downstream support, not causal proof.
if (!exists("progeny_network_qin_d1") || !all(c("source", "target") %in% colnames(progeny_network_qin_d1)) || !exists("qin_day1_arg1pos_received") || nrow(qin_day1_arg1pos_received) == 0) {
  message("Qin Day 1 integration skipped (need PROGENy network + LIANA Arg1pos received with rows).")
} else {
liana_receptors_qin_d1 <- unique(qin_day1_arg1pos_received$receptor.complex)
liana_receptors_qin_d1_split <- unique(trimws(unlist(strsplit(liana_receptors_qin_d1, "[_+]"))))
print(paste("Unique receptors received by Arg1+ neutrophils (LIANA, Qin Day 1):", length(liana_receptors_qin_d1_split)))
liana_expanded_qin_d1 <- qin_day1_arg1pos_received
liana_expanded_qin_d1$receptor_subunits <- lapply(strsplit(liana_expanded_qin_d1$receptor.complex, "[_+]"), function(x) trimws(x))
liana_expanded_qin_d1 <- tidyr::unnest_longer(liana_expanded_qin_d1, receptor_subunits)
pathways_active_qin_d1 <- names(qin_day1_pathway_means_arg1pos)[qin_day1_pathway_means_arg1pos > qin_day1_pathway_means_arg1neg]
# Exploratory only: PROGENy links pathway -> gene targets (footprints), not receptor -> pathway activation.
receptor_pathway_links_qin_d1 <- progeny_network_qin_d1[progeny_network_qin_d1$target %in% liana_receptors_qin_d1_split & progeny_network_qin_d1$source %in% pathways_active_qin_d1, ]
receptor_pathway_links_qin_d1 <- receptor_pathway_links_qin_d1[, c("target", "source")]
colnames(receptor_pathway_links_qin_d1) <- c("receptor_gene", "pathway")
pathway_diff_qin_d1 <- qin_day1_pathway_means_arg1pos - qin_day1_pathway_means_arg1neg
receptor_pathway_links_qin_d1$pathway_activity_diff <- pathway_diff_qin_d1[receptor_pathway_links_qin_d1$pathway]
integration_receptor_pathway_qin_d1 <- merge(liana_expanded_qin_d1, receptor_pathway_links_qin_d1, by.x = "receptor_subunits", by.y = "receptor_gene", all.x = TRUE)
integration_receptor_pathway_qin_d1 <- integration_receptor_pathway_qin_d1[!is.na(integration_receptor_pathway_qin_d1$pathway), ]
write.csv(integration_receptor_pathway_qin_d1, file.path(OUTPUT_DIR, "QinDay1_LIANA_PROGENy_FootprintCoMembership_Exploratory.csv"), row.names = FALSE)
print(paste("LIANA x PROGENy footprint co-membership (exploratory, Qin Day 1):", nrow(integration_receptor_pathway_qin_d1), "rows (not causal receptor->pathway)"))
receptor_tf_links_qin_d1 <- data.frame(receptor_gene = character(0), tf = character(0), omnipath_via = character(0), omnipath_hops = integer(0), stringsAsFactors = FALSE)
active_tfs_qin_d1 <- character(0)
has_bart_key_qin_d1 <- exists("qin_day1_bart_key") && !is.null(qin_day1_bart_key) && !is.null(qin_day1_bart_key[["Arg1pos"]]) && nrow(qin_day1_bart_key[["Arg1pos"]]) > 0
if (has_bart_key_qin_d1) {
  bart_tfs_qin_d1 <- qin_day1_bart_key[["Arg1pos"]]
  tf_col_qin_d1 <- intersect(colnames(bart_tfs_qin_d1), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
  if (length(tf_col_qin_d1) == 0) tf_col_qin_d1 <- colnames(bart_tfs_qin_d1)[1] else tf_col_qin_d1 <- tf_col_qin_d1[1]
  active_tfs_qin_d1 <- unique(as.character(bart_tfs_qin_d1[[tf_col_qin_d1]]))
}
omnipath_ok_qin_d1 <- requireNamespace("OmnipathR", quietly = TRUE) && length(active_tfs_qin_d1) > 0

pathway_interactions_qin_d1 <- data.frame()
omni_qin_d1_first <- if (!omnipath_ok_qin_d1) list(dat = NULL, err1 = NA_character_) else tryCatch(list(dat = OmnipathR::import_omnipath_interactions(datasets = c("omnipath", "pathwayextra"), organism = 10090, genesymbols = TRUE), err1 = NA_character_), error = function(e1) list(dat = NULL, err1 = conditionMessage(e1)))
if (omnipath_ok_qin_d1) pathway_interactions_qin_d1 <- omni_qin_d1_first$dat
omnipath_qin_d1_err1 <- omni_qin_d1_first$err1
if (omnipath_ok_qin_d1 && is.null(pathway_interactions_qin_d1)) pathway_interactions_qin_d1 <- tryCatch(OmnipathR::import_pathwayextra_interactions(organism = 10090, genesymbols = TRUE), error = function(e2) { message("OmniPath curated+pathwayextra import failed (Qin Day 1): ", omnipath_qin_d1_err1, "; fallback: ", conditionMessage(e2)); data.frame() })

omni_qin_d1_cols_ok <- omnipath_ok_qin_d1 && nrow(pathway_interactions_qin_d1) > 0 && all(c("source_genesymbol", "target_genesymbol") %in% colnames(pathway_interactions_qin_d1))
omni_direct_qin_d1 <- if (omni_qin_d1_cols_ok) subset(pathway_interactions_qin_d1, source_genesymbol %in% liana_receptors_qin_d1_split & target_genesymbol %in% active_tfs_qin_d1, select = c("source_genesymbol", "target_genesymbol")) else data.frame()
if (omni_qin_d1_cols_ok && nrow(omni_direct_qin_d1) > 0) receptor_tf_links_qin_d1 <- rbind(receptor_tf_links_qin_d1, data.frame(receptor_gene = omni_direct_qin_d1$source_genesymbol, tf = omni_direct_qin_d1$target_genesymbol, omnipath_via = NA_character_, omnipath_hops = 1L, stringsAsFactors = FALSE))

hop1_qin_d1 <- if (omni_qin_d1_cols_ok) pathway_interactions_qin_d1[pathway_interactions_qin_d1$source_genesymbol %in% liana_receptors_qin_d1_split, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_qin_d1_cols_ok) colnames(hop1_qin_d1) <- c("receptor_gene", "via")
hop2_qin_d1 <- if (omni_qin_d1_cols_ok) pathway_interactions_qin_d1[pathway_interactions_qin_d1$target_genesymbol %in% active_tfs_qin_d1, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_qin_d1_cols_ok) colnames(hop2_qin_d1) <- c("via", "tf")
omni_2hop_qin_d1 <- if (omni_qin_d1_cols_ok) merge(hop1_qin_d1, hop2_qin_d1, by = "via") else data.frame()
if (omni_qin_d1_cols_ok && nrow(omni_2hop_qin_d1) > 0) omni_2hop_qin_d1 <- unique(omni_2hop_qin_d1[, c("receptor_gene", "tf", "via")])
if (omni_qin_d1_cols_ok && nrow(omni_2hop_qin_d1) > 0) receptor_tf_links_qin_d1 <- rbind(receptor_tf_links_qin_d1, data.frame(receptor_gene = omni_2hop_qin_d1$receptor_gene, tf = omni_2hop_qin_d1$tf, omnipath_via = omni_2hop_qin_d1$via, omnipath_hops = 2L, stringsAsFactors = FALSE))
if (omni_qin_d1_cols_ok) receptor_tf_links_qin_d1 <- receptor_tf_links_qin_d1[!duplicated(paste(receptor_tf_links_qin_d1$receptor_gene, receptor_tf_links_qin_d1$tf)), ]

integration_receptor_tf_qin_d1 <- data.frame()
if (nrow(receptor_tf_links_qin_d1) > 0) {
  integration_receptor_tf_qin_d1 <- merge(liana_expanded_qin_d1, receptor_tf_links_qin_d1, by.x = "receptor_subunits", by.y = "receptor_gene", all.x = TRUE)
  integration_receptor_tf_qin_d1 <- integration_receptor_tf_qin_d1[!is.na(integration_receptor_tf_qin_d1$tf), ]
}
if (nrow(integration_receptor_tf_qin_d1) > 0) {
  write.csv(integration_receptor_tf_qin_d1, file.path(OUTPUT_DIR, "QinDay1_LIANA_BARTsc_Integration_DataDriven.csv"), row.names = FALSE)
  n1_qd1 <- sum(integration_receptor_tf_qin_d1$omnipath_hops == 1L, na.rm = TRUE)
  n2_qd1 <- sum(integration_receptor_tf_qin_d1$omnipath_hops == 2L, na.rm = TRUE)
  print(paste0("LIANA -> BARTsc integration (Qin Day 1): ", nrow(integration_receptor_tf_qin_d1), " receptor-TF rows (OmniPath curated+pathwayextra; direct=", n1_qd1, ", two-hop=", n2_qd1, ")"))
}
if (nrow(receptor_tf_links_qin_d1) == 0 && length(active_tfs_qin_d1) > 0) {
  write.csv(data.frame(receptor = liana_receptors_qin_d1_split, note = "Active TFs in Arg1pos (no OmniPath direct/two-hop link):", active_tfs = paste(active_tfs_qin_d1, collapse = "; "), stringsAsFactors = FALSE), file.path(OUTPUT_DIR, "QinDay1_LIANA_BARTsc_Receptors_and_TFs.csv"), row.names = FALSE)
  print("LIANA receptors and BARTsc TFs saved separately (Qin Day 1, no OmniPath link)")
}
has_pathway_qin_d1 <- nrow(integration_receptor_pathway_qin_d1) > 0
has_tf_qin_d1 <- nrow(integration_receptor_tf_qin_d1) > 0
full_chain_qin_d1 <- data.frame()
if (has_pathway_qin_d1 && has_tf_qin_d1) full_chain_qin_d1 <- merge(integration_receptor_pathway_qin_d1, integration_receptor_tf_qin_d1[, c("source", "ligand.complex", "receptor_subunits", "aggregate_rank", "tf", "omnipath_via", "omnipath_hops")], by = c("source", "ligand.complex", "receptor_subunits", "aggregate_rank"), all = TRUE)
if (has_pathway_qin_d1 && !has_tf_qin_d1) { full_chain_qin_d1 <- integration_receptor_pathway_qin_d1; full_chain_qin_d1$tf <- NA_character_ }
if (!has_pathway_qin_d1 && has_tf_qin_d1) { full_chain_qin_d1 <- integration_receptor_tf_qin_d1; full_chain_qin_d1$pathway <- NA_character_; full_chain_qin_d1$pathway_activity_diff <- NA_real_ }
if (nrow(full_chain_qin_d1) > 0) full_chain_qin_d1$evidence_level <- ifelse(!is.na(full_chain_qin_d1$pathway) & !is.na(full_chain_qin_d1$tf), "STRONG (pathway + TF linked)", ifelse(!is.na(full_chain_qin_d1$pathway) | !is.na(full_chain_qin_d1$tf), "MODERATE (pathway or TF linked)", "WEAK (no link)"))
if (nrow(full_chain_qin_d1) > 0) full_chain_qin_d1 <- full_chain_qin_d1[order(full_chain_qin_d1$evidence_level, full_chain_qin_d1$aggregate_rank), ]
# Legacy filename retained for compatibility; contents are an exploratory overlay, not a causal signal chain.
if (nrow(full_chain_qin_d1) > 0) write.csv(full_chain_qin_d1, file.path(OUTPUT_DIR, "QinDay1_FullSignalChain_DataDriven.csv"), row.names = FALSE)
if (nrow(full_chain_qin_d1) > 0) {
  cols_show_qin_d1 <- intersect(c("source", "receptor_subunits", "pathway", "tf", "omnipath_via", "omnipath_hops", "evidence_level"), colnames(full_chain_qin_d1))
  print(head(full_chain_qin_d1[, cols_show_qin_d1, drop = FALSE], 20))
}
heatmap_data_qin_d1 <- data.frame()
if (nrow(full_chain_qin_d1) > 0 && has_pathway_qin_d1) heatmap_data_qin_d1 <- as.data.frame.matrix(table(integration_receptor_pathway_qin_d1$source, integration_receptor_pathway_qin_d1$pathway))
if (nrow(heatmap_data_qin_d1) > 0 && ncol(heatmap_data_qin_d1) > 0) pheatmap::pheatmap(heatmap_data_qin_d1, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("white", "blue", "red"))(50), main = "Qin Day 1: Sender x Pathway (LIANA receptor x PROGENy footprint co-membership; exploratory)")
if (nrow(full_chain_qin_d1) == 0) print("No integration rows Qin Day 1 (PROGENy/BARTsc may not overlap LIANA receptors)")
print("✓ Qin Day 1 LIANA <-> BARTsc <-> PROGENy integration complete (exploratory overlay)")
}

# -------- Qin Day 3: LIANA → LIANA plots → top-receptor Vln → PROGENy → BARTsc → integration (same cohort order as Lee Day 1) --------
print("--- LIANA Analysis: Qin Day 3 ---")

qin_day3_liana_labels <- qin_day3_pruned_labels
qin_day3_liana_labels[qin_day3_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
qin_day3_liana_labels[qin_day3_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
qin_day3_liana_labels[is.na(qin_day3_liana_labels) | qin_day3_liana_labels == ""] <- "Other"
Seurat::Idents(qin_day3) <- factor(qin_day3_liana_labels)

# Qin Day 3: set labels, clean Seurat, build SCE (bypass validObject bug), run LIANA
qin_day3$liana_label <- factor(qin_day3_liana_labels)
SeuratObject::DefaultAssay(qin_day3) <- "RNA"
for (assay_name in names(qin_day3@assays)) {
  if (assay_name != "RNA") qin_day3[[assay_name]] <- NULL
}
for (red_name in names(qin_day3@reductions)) qin_day3[[red_name]] <- NULL
for (graph_name in names(qin_day3@graphs)) qin_day3[[graph_name]] <- NULL
# Custom L-R from CellChat + CellCall (Qin Day 3) — same logic as Lee Day 1.
qin_day3_liana_resource <- "MouseConsensus"
qin_day3_liana_external <- NULL
custom_lr_list_qin_d3 <- list()
path_cc_qin_d3 <- CELLCHAT_QIN_DAY3_RDS
path_ccall_qin_d3 <- CELLCALL_QIN_DAY3_RDS
cc_qin_d3_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_cc_qin_d3) && nzchar(trimws(path_cc_qin_d3)) && file.exists(path_cc_qin_d3) && requireNamespace("CellChat", quietly = TRUE)
cc_qin_d3 <- NULL
if (cc_qin_d3_ok) cc_qin_d3 <- readRDS(path_cc_qin_d3)
comm_qin_d3 <- NULL
if (!is.null(cc_qin_d3) && inherits(cc_qin_d3, "CellChat")) comm_qin_d3 <- CellChat::subsetCommunication(cc_qin_d3)
has_comm_qin_d3 <- !is.null(comm_qin_d3) && nrow(comm_qin_d3) > 0 && "ligand" %in% colnames(comm_qin_d3) && "receptor" %in% colnames(comm_qin_d3)
if (has_comm_qin_d3) {
  cc_lr_qin_d3 <- data.frame(source_genesymbol = comm_qin_d3$ligand, target_genesymbol = comm_qin_d3$receptor, stringsAsFactors = FALSE)
  cc_lr_qin_d3 <- cc_lr_qin_d3[!duplicated(cc_lr_qin_d3), ]
  custom_lr_list_qin_d3[["CellChat"]] <- cc_lr_qin_d3
}
ccall_qin_d3_ok <- isTRUE(LIANA_USE_CUSTOM_LR) && !is.na(path_ccall_qin_d3) && nzchar(trimws(path_ccall_qin_d3)) && file.exists(path_ccall_qin_d3)
ccall_qin_d3 <- NULL
if (ccall_qin_d3_ok) ccall_qin_d3 <- readRDS(path_ccall_qin_d3)
has_ccall_lr_qin_d3 <- !is.null(ccall_qin_d3) && !is.null(ccall_qin_d3@data$expr_l_r_log2_scale)
if (has_ccall_lr_qin_d3) {
  lr_rownames_qin_d3 <- rownames(ccall_qin_d3@data$expr_l_r_log2_scale)
  ccall_lig_qin_d3 <- sub("-.*", "", lr_rownames_qin_d3)
  ccall_rec_qin_d3 <- sub("^[^-]+-", "", lr_rownames_qin_d3)
  ccall_rec_qin_d3 <- gsub("-", "_", ccall_rec_qin_d3)
  ccall_lr_qin_d3 <- data.frame(source_genesymbol = ccall_lig_qin_d3, target_genesymbol = ccall_rec_qin_d3, stringsAsFactors = FALSE)
  ccall_lr_qin_d3 <- ccall_lr_qin_d3[nzchar(ccall_lr_qin_d3$target_genesymbol), ]
  ccall_lr_qin_d3 <- unique(ccall_lr_qin_d3)
  custom_lr_list_qin_d3[["CellCall"]] <- ccall_lr_qin_d3
}
has_custom_lr_qin_d3 <- length(custom_lr_list_qin_d3) > 0
if (has_custom_lr_qin_d3) {
  custom_lr_qin_d3 <- do.call(rbind, custom_lr_list_qin_d3)
  custom_lr_qin_d3 <- custom_lr_qin_d3[!duplicated(custom_lr_qin_d3[, c("source_genesymbol", "target_genesymbol")]), ]
  qin_day3_liana_external <- data.frame(source_genesymbol = custom_lr_qin_d3$source_genesymbol, target_genesymbol = custom_lr_qin_d3$target_genesymbol, stringsAsFactors = FALSE)
  qin_day3_liana_resource <- "custom"
  print(paste("Loaded", nrow(qin_day3_liana_external), "custom L-R pairs from CellChat + CellCall (Qin Day 3)"))
}
if (isTRUE(LIANA_USE_CUSTOM_LR) && !has_custom_lr_qin_d3) print("Qin Day 3 LIANA: LIANA_USE_CUSTOM_LR is TRUE but no CellChat/CellCall pairs loaded — using MouseConsensus only (check RDS paths and CellChat package).")
qin_day3_counts <- tryCatch(Seurat::GetAssayData(qin_day3, slot = "counts"), error = function(e) Seurat::GetAssayData(qin_day3, layer = "counts"))
qin_day3_logcounts <- tryCatch(Seurat::GetAssayData(qin_day3, slot = "data"), error = function(e) Seurat::GetAssayData(qin_day3, layer = "data"))
qin_day3_sce <- SingleCellExperiment(assays = list(counts = qin_day3_counts, logcounts = qin_day3_logcounts))
qin_day3_sce$liana_label <- qin_day3$liana_label
SingleCellExperiment::colLabels(qin_day3_sce) <- qin_day3_sce$liana_label
liana_args_qin_d3 <- list(sce = qin_day3_sce, method = LIANA_METHODS, resource = qin_day3_liana_resource, idents_col = "liana_label", expr_prop = 0.05, verbose = TRUE, min_cells = LIANA_MIN_CELLS, base = exp(1))
if (!is.null(qin_day3_liana_external)) liana_args_qin_d3$external_resource <- qin_day3_liana_external
qin_day3_liana_result <- do.call(liana::liana_wrap, liana_args_qin_d3)
qin_day3_liana_result_df <- liana::liana_aggregate(qin_day3_liana_result)

# Keep only interactions where target is Arg1+ neutrophils (signals RECEIVED by Arg1+ neutrophils)
qin_day3_arg1pos_received <- dplyr::filter(qin_day3_liana_result_df, target == "NeutrophilArg1pos")
qin_day3_arg1pos_received <- dplyr::arrange(qin_day3_arg1pos_received, aggregate_rank)
qin_day3_arg1pos_top_signals <- head(qin_day3_arg1pos_received, 50)
# Arg1- neutrophils received (for specificity contrast)
qin_day3_arg1neg_received <- dplyr::filter(qin_day3_liana_result_df, target == "NeutrophilArg1neg")
qin_day3_arg1neg_received <- dplyr::arrange(qin_day3_arg1neg_received, aggregate_rank)

# Specificity: signals preferentially received by Arg1+ vs Arg1- neutrophils
topN <- 50
qin_day3_arg1pos_top <- head(qin_day3_arg1pos_received, topN)
qin_day3_arg1neg_top <- head(qin_day3_arg1neg_received, topN)
pos_pairs <- paste0(qin_day3_arg1pos_top$ligand.complex, "_", qin_day3_arg1pos_top$receptor.complex)
neg_pairs_all <- paste0(qin_day3_arg1neg_received$ligand.complex, "_", qin_day3_arg1neg_received$receptor.complex)
neg_rank_lookup <- setNames(qin_day3_arg1neg_received$aggregate_rank, neg_pairs_all)
pos_ranks <- qin_day3_arg1pos_top$aggregate_rank
neg_ranks_matched <- neg_rank_lookup[pos_pairs]
neg_ranks_matched[is.na(neg_ranks_matched)] <- 1.0
spec_index <- (neg_ranks_matched - pos_ranks) / (neg_ranks_matched + pos_ranks + 1e-10)
qin_day3_arg1pos_top$specificity_index <- spec_index
qin_day3_arg1pos_specific <- dplyr::filter(qin_day3_arg1pos_top, specificity_index > 0.3 | !(pos_pairs %in% neg_pairs_all))

# Backward-compatible names for downstream code
qin_day3_liana_arg1_pos <- qin_day3_arg1pos_received
qin_day3_liana_arg1_neg <- qin_day3_arg1neg_received
qin_day3_liana_consensus <- head(qin_day3_arg1pos_received, 30)
qin_day3_liana_arg1pos_topranked <- head(qin_day3_arg1pos_received, 20)
qin_day3_liana_arg1neg_topranked <- head(qin_day3_arg1neg_received, 20)

print(paste("Arg1pos received signals (all):", nrow(qin_day3_arg1pos_received)))
print(paste("Arg1pos top signals (top 50):", nrow(qin_day3_arg1pos_top_signals)))
print(paste("Arg1pos-specific signals (specificity > 0.3 or unique to Arg1pos):", nrow(qin_day3_arg1pos_specific)))

# L-R symbols -> full gene names (org.Mm.eg.db), same logic as Lee Day 1 / Lee Day 3
lr_pairs_qin_d3 <- dplyr::distinct(qin_day3_liana_result_df, ligand.complex, receptor.complex, .keep_all = FALSE)
print("Qin Day 3 LIANA: unique L-R pairs in full result (ligand.complex / receptor.complex):")
print(head(lr_pairs_qin_d3, 30))
all_sym_qin_d3 <- unique(c(lr_pairs_qin_d3$ligand.complex, lr_pairs_qin_d3$receptor.complex))
all_sym_single_qin_d3 <- unique(trimws(unlist(strsplit(all_sym_qin_d3[grepl("[_+]", all_sym_qin_d3)], "[_+]"))))
all_sym_single_qin_d3 <- unique(c(all_sym_qin_d3[!grepl("[_+]", all_sym_qin_d3)], all_sym_single_qin_d3))
tbl_lr_qin_d3 <- data.frame(SYMBOL = character(0), GENENAME = character(0), stringsAsFactors = FALSE)
if (length(all_sym_single_qin_d3) > 0) {
  tbl_lr_qin_d3 <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = all_sym_single_qin_d3, columns = "GENENAME", keytype = "SYMBOL"))
  tbl_lr_qin_d3 <- dplyr::distinct(tbl_lr_qin_d3, SYMBOL, .keep_all = TRUE)
}
sym_to_name_qin_d3 <- if (nrow(tbl_lr_qin_d3) > 0) stats::setNames(tbl_lr_qin_d3$GENENAME, tbl_lr_qin_d3$SYMBOL) else character(0)
lr_ref_qin_d3 <- lr_pairs_qin_d3
lr_ref_qin_d3$ligand_genename <- rep(NA_character_, nrow(lr_ref_qin_d3))
lr_ref_qin_d3$receptor_genename <- rep(NA_character_, nrow(lr_ref_qin_d3))
has_sym_to_name_qin_d3 <- length(sym_to_name_qin_d3) > 0
if (!has_sym_to_name_qin_d3) {
  lr_ref_qin_d3$ligand_genename <- lr_ref_qin_d3$ligand.complex
  lr_ref_qin_d3$receptor_genename <- lr_ref_qin_d3$receptor.complex
}
if (has_sym_to_name_qin_d3) {
  ligand_is_complex_q3 <- grepl("[_+]", lr_ref_qin_d3$ligand.complex)
  lr_ref_qin_d3$ligand_genename[!ligand_is_complex_q3] <- ifelse(lr_ref_qin_d3$ligand.complex[!ligand_is_complex_q3] %in% names(sym_to_name_qin_d3), sym_to_name_qin_d3[lr_ref_qin_d3$ligand.complex[!ligand_is_complex_q3]], lr_ref_qin_d3$ligand.complex[!ligand_is_complex_q3])
  complex_ligands_q3 <- lr_ref_qin_d3$ligand.complex[ligand_is_complex_q3]
  mapped_complex_lig_q3 <- vapply(complex_ligands_q3, function(s) {
    parts <- trimws(strsplit(s, "[_+]")[[1]])
    mapped_parts <- ifelse(parts %in% names(sym_to_name_qin_d3), sym_to_name_qin_d3[parts], parts)
    paste(mapped_parts, collapse = "_")
  }, character(1))
  lr_ref_qin_d3$ligand_genename[ligand_is_complex_q3] <- mapped_complex_lig_q3
  receptor_is_complex_q3 <- grepl("[_+]", lr_ref_qin_d3$receptor.complex)
  lr_ref_qin_d3$receptor_genename[!receptor_is_complex_q3] <- ifelse(lr_ref_qin_d3$receptor.complex[!receptor_is_complex_q3] %in% names(sym_to_name_qin_d3), sym_to_name_qin_d3[lr_ref_qin_d3$receptor.complex[!receptor_is_complex_q3]], lr_ref_qin_d3$receptor.complex[!receptor_is_complex_q3])
  complex_receptors_q3 <- lr_ref_qin_d3$receptor.complex[receptor_is_complex_q3]
  mapped_complex_rec_q3 <- vapply(complex_receptors_q3, function(s) {
    parts <- trimws(strsplit(s, "[_+]")[[1]])
    mapped_parts <- ifelse(parts %in% names(sym_to_name_qin_d3), sym_to_name_qin_d3[parts], parts)
    paste(mapped_parts, collapse = "_")
  }, character(1))
  lr_ref_qin_d3$receptor_genename[receptor_is_complex_q3] <- mapped_complex_rec_q3
}
write.csv(lr_ref_qin_d3, file.path(OUTPUT_DIR, "QinDay3_LIANA_LigandReceptor_Reference.csv"), row.names = FALSE)
print(paste("Saved Qin Day 3 L-R reference with gene names to QinDay3_LIANA_LigandReceptor_Reference.csv"))
if (has_sym_to_name_qin_d3) {
  gene_map_df_qin_d3 <- data.frame(symbol = names(sym_to_name_qin_d3), full_name = unname(sym_to_name_qin_d3), stringsAsFactors = FALSE)
  print("LIANA gene symbols in Qin Day 3 -> full names (org.Mm.eg.db):")
  print(gene_map_df_qin_d3)
}

write.csv(qin_day3_liana_result_df, file.path(OUTPUT_DIR, "QinDay3_LIANA_AllResults.csv"), row.names = FALSE)
write.csv(qin_day3_arg1pos_received, file.path(OUTPUT_DIR, "QinDay3_Arg1pos_ReceivedSignals.csv"), row.names = FALSE)
write.csv(qin_day3_arg1pos_top_signals, file.path(OUTPUT_DIR, "QinDay3_LIANA_Arg1pos_Top50Signals.csv"), row.names = FALSE)
write.csv(qin_day3_arg1pos_specific, file.path(OUTPUT_DIR, "QinDay3_LIANA_Arg1pos_Specific.csv"), row.names = FALSE)
write.csv(qin_day3_liana_consensus, file.path(OUTPUT_DIR, "QinDay3_LIANA_NeutrophilConsensus.csv"), row.names = FALSE)
write.csv(qin_day3_liana_arg1pos_topranked, file.path(OUTPUT_DIR, "QinDay3_LIANA_Arg1pos_TopRanked.csv"), row.names = FALSE)
write.csv(qin_day3_liana_arg1neg_topranked, file.path(OUTPUT_DIR, "QinDay3_LIANA_Arg1neg_TopRanked.csv"), row.names = FALSE)
saveRDS(list(liana_result = qin_day3_liana_result, liana_aggregated = qin_day3_liana_result_df, arg1pos_received = qin_day3_arg1pos_received, arg1pos_top_signals = qin_day3_arg1pos_top_signals, arg1pos_specific = qin_day3_arg1pos_specific, neutrophil_consensus = qin_day3_liana_consensus, arg1pos_topranked = qin_day3_liana_arg1pos_topranked, arg1neg_topranked = qin_day3_liana_arg1neg_topranked), file.path(OUTPUT_DIR, "QinDay3_LIANA_Results.rds"))
# Uncomment block below to load and skip re-running Qin Day 3 LIANA:
# qin_day3_liana_loaded <- readRDS(file.path(OUTPUT_DIR, "QinDay3_LIANA_Results.rds"))
# qin_day3_liana_result <- qin_day3_liana_loaded$liana_result
# qin_day3_liana_result_df <- qin_day3_liana_loaded$liana_aggregated
# qin_day3_arg1pos_received <- qin_day3_liana_loaded$arg1pos_received
# qin_day3_arg1pos_top_signals <- qin_day3_liana_loaded$arg1pos_top_signals
# qin_day3_arg1pos_specific <- qin_day3_liana_loaded$arg1pos_specific
# qin_day3_liana_consensus <- qin_day3_liana_loaded$neutrophil_consensus
# qin_day3_liana_arg1pos_topranked <- qin_day3_liana_loaded$arg1pos_topranked
# qin_day3_liana_arg1neg_topranked <- qin_day3_liana_loaded$arg1neg_topranked

# -------- Qin Day 3: LIANA visualizations -> top LIANA receptors VlnPlot (same flow as Lee Day 1) --------
if (nrow(qin_day3_liana_result_df) > 0) {
  liana_network_qin_d3_both <- qin_day3_liana_result_df |>
    dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |>
    dplyr::filter(!(source %in% NEUTROPHIL_STATES)) |>
    dplyr::group_by(target) |>
    dplyr::arrange(aggregate_rank) |>
    dplyr::slice_head(n = 20) |>
    dplyr::ungroup()
  liana_network_qin_d3_both$interaction_label <- paste0(liana_network_qin_d3_both$ligand.complex, "\u2013(", liana_network_qin_d3_both$receptor.complex, ")")
  p_qin_d3_network <- tryCatch(ggplot(liana_network_qin_d3_both, aes(x = target, y = interaction_label, size = -log10(aggregate_rank + 1e-10), color = source)) + geom_point(alpha = 0.8) + scale_x_discrete(limits = NEUTROPHIL_STATES, drop = FALSE) + theme_minimal() + labs(title = "Top Signals Received: Arg1+ vs Arg1- (Qin Day 3)", subtitle = "Y = Ligand\u2013(Receptor); X = receiver; color = sender (neutrophil\u2192neutrophil excluded from data)", x = "Receiver (neutrophil state)", y = "Interaction", color = "Sender cell type", size = "Consensus support\n(-log10 aggregate rank)") + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), plot.margin = margin(10, 80, 10, 10)) + guides(color = guide_legend(override.aes = list(size = 3))), error = function(e) ggplot() + theme_void())
  vals_nlr <- numeric(0)
  if (exists("liana_network_lee_d1_both") && is.data.frame(liana_network_lee_d1_both) && nrow(liana_network_lee_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_lee_d3_both") && is.data.frame(liana_network_lee_d3_both) && nrow(liana_network_lee_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_wang_d3_both") && is.data.frame(liana_network_wang_d3_both) && nrow(liana_network_wang_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_wang_d3_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d1_both") && is.data.frame(liana_network_qin_d1_both) && nrow(liana_network_qin_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d1_both$aggregate_rank + 1e-10))
  if (exists("liana_network_qin_d3_both") && is.data.frame(liana_network_qin_d3_both) && nrow(liana_network_qin_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d3_both$aggregate_rank + 1e-10))
  vals_nlr <- vals_nlr[is.finite(vals_nlr)]
  UNIFY_LIANA_NEGLOG10 <- if (length(vals_nlr) > 0) range(vals_nlr) else c(0, 10)
  if (length(UNIFY_LIANA_NEGLOG10) != 2 || !all(is.finite(UNIFY_LIANA_NEGLOG10))) UNIFY_LIANA_NEGLOG10 <- c(0, 10)
  if (UNIFY_LIANA_NEGLOG10[1] == UNIFY_LIANA_NEGLOG10[2]) UNIFY_LIANA_NEGLOG10[2] <- UNIFY_LIANA_NEGLOG10[1] + 1e-6
  print(p_qin_d3_network + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = (if (exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE else c(2, 8))))
  qin_day3_external_to_neutrophils <- qin_day3_liana_result_df |> dplyr::filter(target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")) |> dplyr::filter(!(source %in% c("NeutrophilArg1pos", "NeutrophilArg1neg"))) |> dplyr::arrange(aggregate_rank)
  liana_top_qin_d3_external <- dplyr::slice_head(qin_day3_external_to_neutrophils, n = 20)
  liana_top_qin_d3_external$lr_label <- paste0(liana_top_qin_d3_external$source, " -> ", liana_top_qin_d3_external$target, "  ", liana_top_qin_d3_external$ligand.complex, "\u2013(", liana_top_qin_d3_external$receptor.complex, ")")
  p_qin_d3_f0 <- tryCatch(ggplot(liana_top_qin_d3_external, aes(x = reorder(lr_label, aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 20 L-R Pairs - LIANA (Qin Day 3)\nExternal signals to Arg1+ / Arg1- neutrophils", x = "Source -> Target  Ligand\u2013(Receptor)", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
  vals_ext <- numeric(0)
  if (exists("lee_day1_external_to_neutrophils") && is.data.frame(lee_day1_external_to_neutrophils) && nrow(lee_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day1_external_to_neutrophils$aggregate_rank)
  if (exists("lee_day3_external_to_neutrophils") && is.data.frame(lee_day3_external_to_neutrophils) && nrow(lee_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day3_external_to_neutrophils$aggregate_rank)
  if (exists("wang_day3_external_to_neutrophils") && is.data.frame(wang_day3_external_to_neutrophils) && nrow(wang_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(wang_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, wang_day3_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day1_external_to_neutrophils") && is.data.frame(qin_day1_external_to_neutrophils) && nrow(qin_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day1_external_to_neutrophils$aggregate_rank)
  if (exists("qin_day3_external_to_neutrophils") && is.data.frame(qin_day3_external_to_neutrophils) && nrow(qin_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day3_external_to_neutrophils$aggregate_rank)
  if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0 && "aggregate_rank" %in% names(liana_top_lee_d1)) vals_ext <- c(vals_ext, liana_top_lee_d1$aggregate_rank)
  if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_lee_d3_external)) vals_ext <- c(vals_ext, liana_top_lee_d3_external$aggregate_rank)
  if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_wang_d3_external)) vals_ext <- c(vals_ext, liana_top_wang_d3_external$aggregate_rank)
  if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d1_external)) vals_ext <- c(vals_ext, liana_top_qin_d1_external$aggregate_rank)
  if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d3_external)) vals_ext <- c(vals_ext, liana_top_qin_d3_external$aggregate_rank)
  vals_ext <- suppressWarnings(as.numeric(vals_ext))
  vals_ext <- vals_ext[is.finite(vals_ext)]
  plot_ext_ar <- suppressWarnings(as.numeric(liana_top_qin_d3_external[["aggregate_rank"]]))
  plot_ext_ar <- plot_ext_ar[is.finite(plot_ext_ar)]
  UNIFY_EXT_AR <- range(c(vals_ext, plot_ext_ar), na.rm = TRUE)
  if (length(plot_ext_ar) == 0 && length(vals_ext) == 0) UNIFY_EXT_AR <- c(0, 1)
  if (!all(is.finite(UNIFY_EXT_AR))) UNIFY_EXT_AR <- c(0, 1)
  if (UNIFY_EXT_AR[1] == UNIFY_EXT_AR[2]) UNIFY_EXT_AR[2] <- UNIFY_EXT_AR[1] + max(abs(UNIFY_EXT_AR[1]) * 1e-6, 1e-12)
  print(p_qin_d3_f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
}
liana_top_qin_d3 <- dplyr::slice_head(qin_day3_liana_consensus, n = 15)
if (nrow(liana_top_qin_d3) > 0) {
  p_qin_d3_f <- ggplot(liana_top_qin_d3, aes(x = reorder(paste0(source, " -> ", target), aggregate_rank), y = aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "Top 15 L-R Pairs - LIANA (Qin Day 3)", x = "Source -> Target", y = "Aggregate Rank") + theme_minimal() + PLOT_TITLE_THEME
  vals_cons <- numeric(0)
  if (exists("liana_top_lee_d3") && is.data.frame(liana_top_lee_d3) && nrow(liana_top_lee_d3) > 0) vals_cons <- c(vals_cons, liana_top_lee_d3$aggregate_rank)
  if (exists("liana_top_wang_d3") && is.data.frame(liana_top_wang_d3) && nrow(liana_top_wang_d3) > 0) vals_cons <- c(vals_cons, liana_top_wang_d3$aggregate_rank)
  if (exists("liana_top_qin_d1") && is.data.frame(liana_top_qin_d1) && nrow(liana_top_qin_d1) > 0) vals_cons <- c(vals_cons, liana_top_qin_d1$aggregate_rank)
  if (exists("liana_top_qin_d3") && is.data.frame(liana_top_qin_d3) && nrow(liana_top_qin_d3) > 0) vals_cons <- c(vals_cons, liana_top_qin_d3$aggregate_rank)
  vals_cons <- vals_cons[is.finite(vals_cons)]
  UNIFY_CONS_AR <- if (length(vals_cons) > 0) range(vals_cons) else c(0, 1)
  if (UNIFY_CONS_AR[1] == UNIFY_CONS_AR[2]) UNIFY_CONS_AR[2] <- UNIFY_CONS_AR[1] + 1e-6
  print(p_qin_d3_f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
}
unique_sources_qin_d3 <- unique(qin_day3_liana_result_df$source)
neutrophil_targets_qin_d3 <- NEUTROPHIL_STATES[NEUTROPHIL_STATES %in% unique(qin_day3_liana_result_df$target)]
dot_sources_qin_d3 <- setdiff(unique_sources_qin_d3, neutrophil_targets_qin_d3)
if (length(dot_sources_qin_d3) == 0) dot_sources_qin_d3 <- unique_sources_qin_d3
p_qin_d3_dot <- NULL
if (length(neutrophil_targets_qin_d3) > 0 && length(dot_sources_qin_d3) > 0) p_qin_d3_dot <- tryCatch(liana::liana_dotplot(qin_day3_liana_result_df, source_groups = dot_sources_qin_d3, target_groups = neutrophil_targets_qin_d3, ntop = 20, size_range = LIANA_DOTPLOT_SIZE_RANGE), error = function(e) NULL)
if (!is.null(p_qin_d3_dot)) { p_qin_d3_dot <- p_qin_d3_dot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)); print(p_qin_d3_dot) }
receptor_freq_qin_d3 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(qin_day3_liana_consensus, receptor.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
p_qin_d3_h <- tryCatch(ggplot(receptor_freq_qin_d3, aes(x = reorder(receptor.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Receptors - Qin Day 3", x = "Receptor", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_qin_d3_h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
ligand_freq_qin_d3 <- tryCatch(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(qin_day3_liana_consensus, ligand.complex), count = dplyr::n(), mean_rank = mean(aggregate_rank)), desc(count)), n = 15), error = function(e) data.frame())
p_qin_d3_i <- tryCatch(ggplot(ligand_freq_qin_d3, aes(x = reorder(ligand.complex, -count), y = count, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top 15 Ligands - Qin Day 3", x = "Ligand", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
print(p_qin_d3_i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
interaction_matrix_qin_d3 <- tryCatch(tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(qin_day3_liana_consensus, source, target), interaction_count = dplyr::n(), .groups = "drop"), names_from = target, values_from = interaction_count, values_fill = 0), error = function(e) data.frame())
interaction_matrix_qin_d3_mat <- tryCatch({ cols_num <- setdiff(colnames(interaction_matrix_qin_d3), "source"); mat <- as.matrix(interaction_matrix_qin_d3[, cols_num, drop = FALSE]); rownames(mat) <- interaction_matrix_qin_d3$source; mat }, error = function(e) matrix(0, nrow = 0, ncol = 0))
p_qin_d3_j <- tryCatch(pheatmap::pheatmap(interaction_matrix_qin_d3_mat, color = colorRampPalette(c("white", "yellow", "orange", "red"))(100), main = "Cell Type Interaction Frequency - Qin Day 3", display_numbers = TRUE), error = function(e) NULL)
if (!is.null(p_qin_d3_j)) { print(p_qin_d3_j); grid::grid.newpage(); grid::grid.draw(p_qin_d3_j$gtable) }
liana_trunc_qin_d3 <- dplyr::filter(qin_day3_liana_result_df, aggregate_rank <= 0.01)
if (nrow(liana_trunc_qin_d3) > 0) {
  liana_trunc_qin_d3$source <- as.character(liana_trunc_qin_d3$source)
  liana_trunc_qin_d3$target <- as.character(liana_trunc_qin_d3$target)
  p_heat_freq_qin_d3 <- tryCatch(liana::heat_freq(liana_trunc_qin_d3), error = function(e) { message("heat_freq Qin D3: ", conditionMessage(e)); NULL })
  if (!is.null(p_heat_freq_qin_d3)) print(p_heat_freq_qin_d3)
  unique_sources_chord_qin_d3 <- unique(liana_trunc_qin_d3$source)
  unique_targets_chord_qin_d3 <- unique(liana_trunc_qin_d3$target)
  grDevices::png(file.path(OUTPUT_DIR, "QinDay3_LIANA_ChordFreq.png"), width = 1400, height = 1400, res = 150)
  tryCatch(liana::chord_freq(liana_trunc_qin_d3, source_groups = unique_sources_chord_qin_d3, target_groups = unique_targets_chord_qin_d3), error = function(e) message("chord_freq Qin D3 (PNG): ", conditionMessage(e)))
  grDevices::dev.off()
  p_chord_freq_qin_d3 <- tryCatch(liana::chord_freq(liana_trunc_qin_d3, source_groups = unique_sources_chord_qin_d3, target_groups = unique_targets_chord_qin_d3), error = function(e) { message("chord_freq Qin D3: ", conditionMessage(e)); NULL })
  if (!is.null(p_chord_freq_qin_d3)) print(p_chord_freq_qin_d3)
  liana_mat_qin_d3 <- as.matrix(table(liana_trunc_qin_d3$source, liana_trunc_qin_d3$target))
  p_liana_heatmap_qin_d3 <- NULL
  if (nrow(liana_mat_qin_d3) > 0 && ncol(liana_mat_qin_d3) > 0) p_liana_heatmap_qin_d3 <- tryCatch(liana::liana_heatmap(liana_mat_qin_d3), error = function(e) { message("liana_heatmap Qin D3: ", conditionMessage(e)); NULL })
  if (!is.null(p_liana_heatmap_qin_d3)) ComplexHeatmap::draw(p_liana_heatmap_qin_d3)
}
if (nrow(liana_trunc_qin_d3) == 0) print("No interactions with aggregate_rank <= 0.01 in Qin Day 3; skipping heat_freq, chord_freq (PNG+screen), and liana_heatmap on truncated table.")
source_importance_qin_d3 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(qin_day3_liana_consensus, source), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), importance_score = dplyr::n() * (1 - mean(aggregate_rank))), desc(importance_score)), error = function(e) data.frame())
p_qin_d3_k <- tryCatch(ggplot(source_importance_qin_d3, aes(x = reorder(source, importance_score), y = importance_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkred") + labs(title = "Source Cell Importance - Qin Day 3", x = "Cell Type", y = "Importance Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_iy <- numeric(0)
vals_imr <- numeric(0)
imp_tabnames <- c("source_importance_lee_d1", "source_importance_lee_d3", "source_importance_wang_d3", "source_importance_qin_d1", "source_importance_qin_d3")
for (fn in imp_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("importance_score" %in% names(d)) vals_iy <- c(vals_iy, d$importance_score)
  if ("mean_rank" %in% names(d)) vals_imr <- c(vals_imr, d$mean_rank)
}
vals_iy <- vals_iy[is.finite(vals_iy)]
vals_imr <- vals_imr[is.finite(vals_imr)]
UNIFY_SRC_IMP_Y <- if (length(vals_iy) > 0) range(vals_iy) else c(0, 1)
UNIFY_SRC_IMP_MR <- if (length(vals_imr) > 0) range(vals_imr) else c(0, 1)
if (UNIFY_SRC_IMP_Y[1] == UNIFY_SRC_IMP_Y[2]) UNIFY_SRC_IMP_Y[2] <- UNIFY_SRC_IMP_Y[1] + 1e-6
if (UNIFY_SRC_IMP_MR[1] == UNIFY_SRC_IMP_MR[2]) UNIFY_SRC_IMP_MR[2] <- UNIFY_SRC_IMP_MR[1] + 1e-6
print(p_qin_d3_k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
p_qin_d3_n <- tryCatch(ggplot(head(qin_day3_liana_arg1pos_topranked, 15), aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), aggregate_rank), y = -aggregate_rank, fill = aggregate_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "red", high = "green") + labs(title = "Top Arg1pos Interactions - Qin Day 3", x = "L-R Pair", y = "Rank Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_a1 <- numeric(0)
if (exists("lee_day1_liana_arg1pos_topranked") && is.data.frame(lee_day1_liana_arg1pos_topranked) && nrow(lee_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("lee_day3_liana_arg1pos_topranked") && is.data.frame(lee_day3_liana_arg1pos_topranked) && nrow(lee_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("wang_day3_liana_arg1pos_topranked") && is.data.frame(wang_day3_liana_arg1pos_topranked) && nrow(wang_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(wang_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day1_liana_arg1pos_topranked") && is.data.frame(qin_day1_liana_arg1pos_topranked) && nrow(qin_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day3_liana_arg1pos_topranked") && is.data.frame(qin_day3_liana_arg1pos_topranked) && nrow(qin_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day3_liana_arg1pos_topranked$aggregate_rank, 15))
vals_a1 <- vals_a1[is.finite(vals_a1)]
UNIFY_ARG1_AR <- if (length(vals_a1) > 0) range(vals_a1) else c(0, 1)
if (UNIFY_ARG1_AR[1] == UNIFY_ARG1_AR[2]) UNIFY_ARG1_AR[2] <- UNIFY_ARG1_AR[1] + 1e-6
UNIFY_ARG1_NEGY <- c(-UNIFY_ARG1_AR[2], -UNIFY_ARG1_AR[1])
print(p_qin_d3_n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
qin_day3_arg1pos_specific_plot <- head(qin_day3_arg1pos_specific, 15)
p_qin_d3_o2 <- tryCatch(ggplot(qin_day3_arg1pos_specific_plot, aes(x = reorder(paste0(ligand.complex, " -> ", receptor.complex), specificity_index), y = specificity_index, fill = source)) + geom_bar(stat = "identity") + coord_flip() + labs(title = "Arg1+-Specific Signals (Qin Day 3)", x = "L-R Pair", y = "Specificity Index") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
print(p_qin_d3_o2)
target_specificity_qin_d3 <- tryCatch(dplyr::arrange(dplyr::summarise(dplyr::group_by(dplyr::filter(qin_day3_liana_result_df, target %in% c("NeutrophilArg1pos", "NeutrophilArg1neg")), target), interaction_count = dplyr::n(), mean_rank = mean(aggregate_rank), specificity_score = dplyr::n() * (1 - mean(aggregate_rank)), .groups = "drop"), desc(specificity_score)), error = function(e) data.frame())
p_qin_d3_l <- tryCatch(ggplot(target_specificity_qin_d3, aes(x = reorder(target, specificity_score), y = specificity_score, fill = mean_rank)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_gradient(low = "lightblue", high = "darkgreen") + labs(title = "Target Cell Specificity - Qin Day 3", x = "Cell Type", y = "Specificity Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_sy <- numeric(0)
vals_smr <- numeric(0)
spec_tabnames <- c("target_specificity_lee_d1", "target_specificity_lee_d3", "target_specificity_wang_d3", "target_specificity_qin_d1", "target_specificity_qin_d3")
for (fn in spec_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("specificity_score" %in% names(d)) vals_sy <- c(vals_sy, d$specificity_score)
  if ("mean_rank" %in% names(d)) vals_smr <- c(vals_smr, d$mean_rank)
}
vals_sy <- vals_sy[is.finite(vals_sy)]
vals_smr <- vals_smr[is.finite(vals_smr)]
UNIFY_TGT_SPEC_Y <- if (length(vals_sy) > 0) range(vals_sy) else c(0, 1)
UNIFY_TGT_SPEC_MR <- if (length(vals_smr) > 0) range(vals_smr) else c(0, 1)
if (UNIFY_TGT_SPEC_Y[1] == UNIFY_TGT_SPEC_Y[2]) UNIFY_TGT_SPEC_Y[2] <- UNIFY_TGT_SPEC_Y[1] + 1e-6
if (UNIFY_TGT_SPEC_MR[1] == UNIFY_TGT_SPEC_MR[2]) UNIFY_TGT_SPEC_MR[2] <- UNIFY_TGT_SPEC_MR[1] + 1e-6
print(p_qin_d3_l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
p_qin_d3_m <- tryCatch(ggplot(qin_day3_liana_consensus, aes(x = aggregate_rank)) + geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) + labs(title = "Distribution of L-R Pair Aggregate Ranks - Qin Day 3", x = "Aggregate Rank Score", y = "Frequency") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void())
vals_hist <- numeric(0)
if (exists("lee_day1_liana_consensus") && is.data.frame(lee_day1_liana_consensus) && nrow(lee_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day1_liana_consensus)) vals_hist <- c(vals_hist, lee_day1_liana_consensus$aggregate_rank)
if (exists("lee_day3_liana_consensus") && is.data.frame(lee_day3_liana_consensus) && nrow(lee_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day3_liana_consensus)) vals_hist <- c(vals_hist, lee_day3_liana_consensus$aggregate_rank)
if (exists("wang_day3_liana_consensus") && is.data.frame(wang_day3_liana_consensus) && nrow(wang_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(wang_day3_liana_consensus)) vals_hist <- c(vals_hist, wang_day3_liana_consensus$aggregate_rank)
if (exists("qin_day1_liana_consensus") && is.data.frame(qin_day1_liana_consensus) && nrow(qin_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day1_liana_consensus)) vals_hist <- c(vals_hist, qin_day1_liana_consensus$aggregate_rank)
if (exists("qin_day3_liana_consensus") && is.data.frame(qin_day3_liana_consensus) && nrow(qin_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day3_liana_consensus)) vals_hist <- c(vals_hist, qin_day3_liana_consensus$aggregate_rank)
vals_hist <- vals_hist[is.finite(vals_hist)]
UNIFY_HIST_AR <- if (length(vals_hist) > 0) range(vals_hist) else c(0, 1)
if (UNIFY_HIST_AR[1] == UNIFY_HIST_AR[2]) UNIFY_HIST_AR[2] <- UNIFY_HIST_AR[1] + 1e-6
print(p_qin_d3_m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
method_cols_qin_d3 <- colnames(qin_day3_liana_result_df)[grep("pval_", colnames(qin_day3_liana_result_df))]
consensus_qin_d3_step1 <- dplyr::slice_head(qin_day3_liana_consensus, n = 15)
consensus_qin_d3 <- NULL
consensus_qin_d3_has_cols <- length(method_cols_qin_d3) > 0 && nrow(consensus_qin_d3_step1) > 0
if (consensus_qin_d3_has_cols) consensus_qin_d3 <- as.data.frame(dplyr::select(consensus_qin_d3_step1, dplyr::all_of(method_cols_qin_d3)))
if (consensus_qin_d3_has_cols) consensus_qin_d3_rownames <- paste0(consensus_qin_d3_step1$source, " | ", consensus_qin_d3_step1$ligand.complex, " -> ", consensus_qin_d3_step1$receptor.complex)
if (consensus_qin_d3_has_cols) rownames(consensus_qin_d3) <- make.unique(as.character(consensus_qin_d3_rownames))
consensus_qin_d3_ready <- !is.null(consensus_qin_d3) && nrow(consensus_qin_d3) > 0
if (consensus_qin_d3_ready) { p_qin_d3_o <- pheatmap::pheatmap(consensus_qin_d3, color = colorRampPalette(c("red", "white", "blue"))(100), main = "Method Consensus (p-values) - Qin Day 3"); print(p_qin_d3_o) }


qin_day3_neut_cells_vln <- colnames(qin_day3)[c(qin_day3_pos, qin_day3_neg)]
top_receptors_qin_d3 <- head(unique(qin_day3_arg1pos_received$receptor.complex), LIANA_TOP_RECEPTOR_VLN)
top_receptors_qin_d3_single <- top_receptors_qin_d3[!grepl("[_+]", top_receptors_qin_d3)]
top_receptors_qin_d3_in_data <- top_receptors_qin_d3_single[top_receptors_qin_d3_single %in% rownames(qin_day3)]
if (length(top_receptors_qin_d3_in_data) > 0) {
  qin_day3_neut_obj_vln <- subset(qin_day3, cells = qin_day3_neut_cells_vln)
  p_receptor_vln_qin_d3 <- Seurat::VlnPlot(qin_day3_neut_obj_vln, features = top_receptors_qin_d3_in_data, group.by = "arg1_status", pt.size = 0.1, ncol = min(3L, length(top_receptors_qin_d3_in_data)))
  print(p_receptor_vln_qin_d3)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "QinDay3_TopLIANA_Receptors_VlnPlot.png"), p_receptor_vln_qin_d3, width = 10, height = 6, dpi = 150)
}
print("✓ Qin Day 3 LIANA analysis complete")


# Checkpoint: resume from here to run Qin Day 3 PROGENy/BARTsc or STEP 0. load(file.path(OUTPUT_DIR, "Workspace_AfterQinDay3_LIANA.RData"))
save.image(file.path(OUTPUT_DIR, "Workspace_AfterQinDay3_LIANA.RData"))
saveRDS(list(qin_day3_liana_result_df = qin_day3_liana_result_df, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay3_LIANA.rds"))
# checkpoint_qd3 <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay3_LIANA.rds")); list2env(checkpoint_qd3, envir = .GlobalEnv)
print("✓ Checkpoint saved: Workspace_AfterQinDay3_LIANA.RData, Checkpoint_AfterQinDay3_LIANA.rds")

# -------- Qin Day 3 PROGENy Pathway Analysis (Exploratory; downstream pathway support only) --------
print("--- PROGENy Pathway Analysis: Qin Day 3 (Exploratory) ---")
# Secondary support only: PROGENy summarizes downstream pathway state in Arg1pos vs Arg1neg neutrophils.
progeny_model_mouse_qin_d3 <- progeny::model_mouse_full
colnames(progeny_model_mouse_qin_d3) <- tolower(colnames(progeny_model_mouse_qin_d3))
colnames(progeny_model_mouse_qin_d3)[colnames(progeny_model_mouse_qin_d3) == "p.value"] <- "p_value"
stopifnot(all(c("gene", "pathway", "weight") %in% colnames(progeny_model_mouse_qin_d3)))
progeny_network_qin_d3 <- data.frame(source = progeny_model_mouse_qin_d3$pathway, target = progeny_model_mouse_qin_d3$gene, weight = progeny_model_mouse_qin_d3$weight, stringsAsFactors = FALSE)
progeny_network_qin_d3 <- progeny_network_qin_d3[progeny_network_qin_d3$weight != 0, ]
qin_day3_neut_cells <- colnames(qin_day3)[c(qin_day3_pos, qin_day3_neg)]
qin_day3_neut_expr_mat <- tryCatch(Seurat::GetAssayData(qin_day3, layer = "data")[, qin_day3_neut_cells, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day3_progeny_result <- tryCatch(decoupleR::run_wmean(mat = qin_day3_neut_expr_mat, network = progeny_network_qin_d3, .source = "source", .target = "target", .mor = "weight", minsize = 5), error = function(e) data.frame())
qin_day3_progeny_for_wide <- if (nrow(qin_day3_progeny_result) > 0 && "statistic" %in% colnames(qin_day3_progeny_result)) dplyr::filter(qin_day3_progeny_result, statistic == "norm_wmean") else qin_day3_progeny_result
if (nrow(qin_day3_progeny_for_wide) == 0 && nrow(qin_day3_progeny_result) > 0 && "statistic" %in% colnames(qin_day3_progeny_result)) qin_day3_progeny_for_wide <- dplyr::filter(qin_day3_progeny_result, statistic == "wmean")
qin_day3_pw_wide <- tryCatch(tidyr::pivot_wider(qin_day3_progeny_for_wide, names_from = "condition", values_from = "score", id_cols = "source"), error = function(e) data.frame())
qin_day3_pw_cols_num <- tryCatch(sapply(qin_day3_pw_wide[, -1, drop = FALSE], function(x) as.numeric(unlist(x))), error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day3_progeny_scores_mat <- tryCatch(as.matrix(qin_day3_pw_cols_num), error = function(e) matrix(0, nrow = 0, ncol = 0))
rownames(qin_day3_progeny_scores_mat) <- tryCatch(as.character(qin_day3_pw_wide$source), error = function(e) character(0))
colnames(qin_day3_progeny_scores_mat) <- tryCatch(colnames(qin_day3_pw_wide)[-1], error = function(e) character(0))
qin_day3_arg1pos_cells <- colnames(qin_day3)[qin_day3_pos]
qin_day3_arg1neg_cells <- colnames(qin_day3)[qin_day3_neg]
qin_day3_arg1pos_cells_in_mat <- qin_day3_arg1pos_cells[qin_day3_arg1pos_cells %in% colnames(qin_day3_progeny_scores_mat)]
qin_day3_arg1neg_cells_in_mat <- qin_day3_arg1neg_cells[qin_day3_arg1neg_cells %in% colnames(qin_day3_progeny_scores_mat)]
qin_day3_progeny_arg1pos_scores <- tryCatch(qin_day3_progeny_scores_mat[, qin_day3_arg1pos_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day3_progeny_arg1neg_scores <- tryCatch(qin_day3_progeny_scores_mat[, qin_day3_arg1neg_cells_in_mat, drop = FALSE], error = function(e) matrix(0, nrow = 0, ncol = 0))
qin_day3_pathway_means_arg1pos <- tryCatch(rowMeans(qin_day3_progeny_arg1pos_scores, na.rm = TRUE), error = function(e) numeric(0))
qin_day3_pathway_means_arg1neg <- tryCatch(rowMeans(qin_day3_progeny_arg1neg_scores, na.rm = TRUE), error = function(e) numeric(0))
qin_day3_pathway_medians_arg1pos <- tryCatch(apply(qin_day3_progeny_arg1pos_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
qin_day3_pathway_medians_arg1neg <- tryCatch(apply(qin_day3_progeny_arg1neg_scores, 1, function(x) stats::median(x, na.rm = TRUE)), error = function(e) numeric(0))
qin_day3_pathway_shift <- pmax(0, -pmin(qin_day3_pathway_means_arg1pos, qin_day3_pathway_means_arg1neg, na.rm = TRUE)) + 1e-10
qin_day3_pathway_log2fc <- log2((qin_day3_pathway_means_arg1pos + qin_day3_pathway_shift) / (qin_day3_pathway_means_arg1neg + qin_day3_pathway_shift))
qin_day3_pathway_comparison <- tryCatch(data.frame(
  pathway = names(qin_day3_pathway_means_arg1pos),
  Arg1pos_mean = qin_day3_pathway_means_arg1pos,
  Arg1neg_mean = qin_day3_pathway_means_arg1neg,
  Arg1pos_median = qin_day3_pathway_medians_arg1pos,
  Arg1neg_median = qin_day3_pathway_medians_arg1neg,
  mean_difference = qin_day3_pathway_means_arg1pos - qin_day3_pathway_means_arg1neg,
  median_difference = qin_day3_pathway_medians_arg1pos - qin_day3_pathway_medians_arg1neg,
  log2FC = qin_day3_pathway_log2fc,
  stringsAsFactors = FALSE
), error = function(e) data.frame())
qin_day3_all_pathways <- tryCatch(unique(progeny_network_qin_d3$source), error = function(e) character(0))
qin_day3_pathway_results <- qin_day3_pathway_comparison
qin_day3_arg1pos_count <- length(qin_day3_pos)
qin_day3_arg1neg_count <- length(qin_day3_neg)
print(paste("Qin Day 3: Arg1pos neutrophils =", qin_day3_arg1pos_count, ", Arg1neg neutrophils =", qin_day3_arg1neg_count))
qin_day3_pathway_pvals <- numeric(length(qin_day3_all_pathways))
names(qin_day3_pathway_pvals) <- qin_day3_all_pathways
qin_day3_pathway_effect_sizes <- numeric(length(qin_day3_all_pathways))
names(qin_day3_pathway_effect_sizes) <- qin_day3_all_pathways
qin_day3_pathway_ci_lower <- numeric(length(qin_day3_all_pathways))
names(qin_day3_pathway_ci_lower) <- qin_day3_all_pathways
qin_day3_pathway_ci_upper <- numeric(length(qin_day3_all_pathways))
names(qin_day3_pathway_ci_upper) <- qin_day3_all_pathways
qin_day3_pathway_adequate_n <- logical(length(qin_day3_all_pathways))
names(qin_day3_pathway_adequate_n) <- qin_day3_all_pathways
for (pw_idx in seq_along(qin_day3_all_pathways)) {
  pw <- qin_day3_all_pathways[pw_idx]
  arg1pos_scores_pw <- tryCatch(as.numeric(qin_day3_progeny_arg1pos_scores[pw, ]), error = function(e) numeric(0))
  arg1neg_scores_pw <- tryCatch(as.numeric(qin_day3_progeny_arg1neg_scores[pw, ]), error = function(e) numeric(0))
  arg1pos_scores_pw <- arg1pos_scores_pw[!is.na(arg1pos_scores_pw)]
  arg1neg_scores_pw <- arg1neg_scores_pw[!is.na(arg1neg_scores_pw)]
  n_arg1pos <- length(arg1pos_scores_pw)
  n_arg1neg <- length(arg1neg_scores_pw)
  qin_day3_pathway_adequate_n[pw] <- (n_arg1pos >= 5) & (n_arg1neg >= 5)
  median_arg1pos <- tryCatch(stats::median(arg1pos_scores_pw, na.rm = TRUE), error = function(e) 0)
  median_arg1neg <- tryCatch(stats::median(arg1neg_scores_pw, na.rm = TRUE), error = function(e) 0)
  qin_day3_pathway_effect_sizes[pw] <- median_arg1pos - median_arg1neg
  qin_day3_pathway_pvals[pw] <- 1.0
  qin_day3_pathway_ci_lower[pw] <- NA_real_
  qin_day3_pathway_ci_upper[pw] <- NA_real_
  wilcox_result <- tryCatch(stats::wilcox.test(arg1pos_scores_pw, arg1neg_scores_pw, conf.int = TRUE, conf.level = 0.95), error = function(e) NULL)
  wilcox_pvalue <- tryCatch(if (!is.null(wilcox_result)) wilcox_result$p.value else 1.0, error = function(e) 1.0)
  wilcox_pvalue_length <- tryCatch(length(wilcox_pvalue), error = function(e) 0)
  wilcox_pvalue_final <- tryCatch(if (wilcox_pvalue_length > 0) wilcox_pvalue[1] else 1.0, error = function(e) 1.0)
  qin_day3_pathway_pvals[pw] <- wilcox_pvalue_final
  wilcox_ci_lower <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[1] else NA_real_, error = function(e) NA_real_)
  wilcox_ci_upper <- tryCatch(if (!is.null(wilcox_result) && !is.null(wilcox_result$conf.int)) wilcox_result$conf.int[2] else NA_real_, error = function(e) NA_real_)
  qin_day3_pathway_ci_lower[pw] <- wilcox_ci_lower
  qin_day3_pathway_ci_upper[pw] <- wilcox_ci_upper
  qin_day3_pathway_warning_msg <- tryCatch(paste("Warning: Pathway", pw, "has insufficient sample size (Arg1pos n =", n_arg1pos, ", Arg1neg n =", n_arg1neg, "). Skipping statistical test."), error = function(e) "")
  qin_day3_pathway_warning_vector <- c("", qin_day3_pathway_warning_msg)
  qin_day3_pathway_warning_index <- tryCatch(as.numeric(!qin_day3_pathway_adequate_n[pw]) + 1, error = function(e) 1)
  print(qin_day3_pathway_warning_vector[qin_day3_pathway_warning_index])
}
qin_day3_pathway_pvals_adj <- p.adjust(qin_day3_pathway_pvals, method = "BH")
qin_day3_pathway_results$p_value <- tryCatch(qin_day3_pathway_pvals[qin_day3_pathway_results$pathway], error = function(e) rep(1.0, nrow(qin_day3_pathway_results)))
qin_day3_pathway_results$p_adj <- tryCatch(qin_day3_pathway_pvals_adj[qin_day3_pathway_results$pathway], error = function(e) rep(1.0, nrow(qin_day3_pathway_results)))
qin_day3_pathway_results$effect_size_median_diff <- tryCatch(qin_day3_pathway_effect_sizes[qin_day3_pathway_results$pathway], error = function(e) rep(0.0, nrow(qin_day3_pathway_results)))
qin_day3_pathway_results$ci_lower_95 <- tryCatch(qin_day3_pathway_ci_lower[qin_day3_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(qin_day3_pathway_results)))
qin_day3_pathway_results$ci_upper_95 <- tryCatch(qin_day3_pathway_ci_upper[qin_day3_pathway_results$pathway], error = function(e) rep(NA_real_, nrow(qin_day3_pathway_results)))
qin_day3_pathway_results$adequate_sample_size <- tryCatch(qin_day3_pathway_adequate_n[qin_day3_pathway_results$pathway], error = function(e) rep(FALSE, nrow(qin_day3_pathway_results)))
qin_day3_pathway_results$significant <- tryCatch((qin_day3_pathway_results$p_adj < 0.05) & qin_day3_pathway_results$adequate_sample_size, error = function(e) rep(FALSE, nrow(qin_day3_pathway_results)))
write.csv(qin_day3_pathway_comparison, file.path(OUTPUT_DIR, "QinDay3_PROGENy_PathwayComparison.csv"), row.names = FALSE)
write.csv(qin_day3_pathway_results, file.path(OUTPUT_DIR, "QinDay3_PROGENy_PathwayResults.csv"), row.names = FALSE)
saveRDS(list(pathway_comparison = qin_day3_pathway_comparison, pathway_results = qin_day3_pathway_results, pathway_means_arg1pos = qin_day3_pathway_means_arg1pos, pathway_means_arg1neg = qin_day3_pathway_means_arg1neg, progeny_network = progeny_network_qin_d3), file.path(OUTPUT_DIR, "QinDay3_PROGENy_Results.rds"))
qin_day3_identified_ligands <- tryCatch(unique(c(qin_day3_liana_arg1pos_topranked$ligand.complex, qin_day3_liana_arg1neg_topranked$ligand.complex)), error = function(e) character(0))
qin_day3_ligand_pathway_map <- data.frame(ligand = character(0), pathway = character(0), weight = numeric(0), stringsAsFactors = FALSE)
for (lig in qin_day3_identified_ligands) { lig_genes <- trimws(unlist(strsplit(lig, "[_+]"))); lig_pathways <- dplyr::filter(progeny_network_qin_d3, target %in% lig_genes); if (nrow(lig_pathways) > 0) qin_day3_ligand_pathway_map <- rbind(qin_day3_ligand_pathway_map, data.frame(ligand = lig, pathway = lig_pathways$source, weight = lig_pathways$weight, stringsAsFactors = FALSE)) }
qin_day3_pathway_long <- tryCatch(tidyr::pivot_longer(qin_day3_pathway_results, cols = c("Arg1pos_mean", "Arg1neg_mean"), names_to = "Group", values_to = "Pathway_Score"), error = function(e) data.frame())
p_qin_d3_progeny1 <- tryCatch(ggplot(qin_day3_pathway_long, aes(x = pathway, y = Pathway_Score, fill = Group)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("Arg1pos_mean" = "red", "Arg1neg_mean" = "lightblue"), labels = c("Arg1pos", "Arg1neg")) + labs(title = "PROGENy Pathway Activity - Qin Day 3", x = "Pathway", y = "Pathway Activity Score") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
vals_pw <- numeric(0)
if (exists("lee_day1_pathway_long") && is.data.frame(lee_day1_pathway_long) && nrow(lee_day1_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day1_pathway_long)) vals_pw <- c(vals_pw, lee_day1_pathway_long$Pathway_Score)
if (exists("lee_day3_pathway_long") && is.data.frame(lee_day3_pathway_long) && nrow(lee_day3_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day3_pathway_long)) vals_pw <- c(vals_pw, lee_day3_pathway_long$Pathway_Score)
if (exists("wang_day3_pathway_long") && is.data.frame(wang_day3_pathway_long) && nrow(wang_day3_pathway_long) > 0 && "Pathway_Score" %in% names(wang_day3_pathway_long)) vals_pw <- c(vals_pw, wang_day3_pathway_long$Pathway_Score)
if (exists("qin_day1_pathway_long") && is.data.frame(qin_day1_pathway_long) && nrow(qin_day1_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day1_pathway_long)) vals_pw <- c(vals_pw, qin_day1_pathway_long$Pathway_Score)
if (exists("qin_day3_pathway_long") && is.data.frame(qin_day3_pathway_long) && nrow(qin_day3_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day3_pathway_long)) vals_pw <- c(vals_pw, qin_day3_pathway_long$Pathway_Score)
vals_pw <- vals_pw[is.finite(vals_pw)]
UNIFY_PW_Y <- if (length(vals_pw) > 0) range(vals_pw) else c(-1, 1)
if (UNIFY_PW_Y[1] == UNIFY_PW_Y[2]) UNIFY_PW_Y[2] <- UNIFY_PW_Y[1] + 1e-6
print(p_qin_d3_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
top_pw_qin_d3 <- character(0)
if (nrow(qin_day3_pathway_comparison) > 0 && "mean_difference" %in% colnames(qin_day3_pathway_comparison)) {
  ord_pw_qin_d3 <- order(-abs(qin_day3_pathway_comparison$mean_difference))
  take_pw_qin_d3 <- min(6L, length(ord_pw_qin_d3))
  top_pw_qin_d3 <- as.character(qin_day3_pathway_comparison$pathway[ord_pw_qin_d3[seq_len(take_pw_qin_d3)]])
}
if (length(top_pw_qin_d3) > 0 && nrow(qin_day3_progeny_scores_mat) > 0) {
  n_cells_pw_q3 <- ncol(qin_day3_progeny_scores_mat)
  n_pw_sel_q3 <- length(top_pw_qin_d3)
  cell_vec_pw_q3 <- rep(colnames(qin_day3_progeny_scores_mat), each = n_pw_sel_q3)
  pathway_vec_pw_q3 <- rep(top_pw_qin_d3, times = n_cells_pw_q3)
  score_vec_pw_q3 <- as.vector(t(qin_day3_progeny_scores_mat[top_pw_qin_d3, , drop = FALSE]))
  pw_scores_long_qin_d3 <- data.frame(cell = cell_vec_pw_q3, pathway = pathway_vec_pw_q3, score = score_vec_pw_q3, stringsAsFactors = FALSE)
  pw_scores_long_qin_d3$arg1_status <- ifelse(pw_scores_long_qin_d3$cell %in% qin_day3_arg1pos_cells, "Arg1pos", "Arg1neg")
  p_progeny_vln_qin_d3 <- ggplot(pw_scores_long_qin_d3, aes(x = arg1_status, y = score, fill = arg1_status)) + geom_violin(trim = FALSE) + geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) + facet_wrap(~pathway, scales = "free_y") + labs(title = "PROGENy pathway scores: Arg1+ vs Arg1- neutrophils (Qin Day 3)", x = "ARG1 status", y = "Activity score (norm_wmean)") + theme_minimal() + PLOT_TITLE_THEME
  print(p_progeny_vln_qin_d3)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "QinDay3_PROGENy_TopPathways_Violin.png"), p_progeny_vln_qin_d3, width = 12, height = 8, dpi = 150)
}
progeny_heatmap_mat_qin_d3 <- tryCatch(rbind(Arg1pos = qin_day3_pathway_means_arg1pos, Arg1neg = qin_day3_pathway_means_arg1neg), error = function(e) matrix(0, nrow = 0, ncol = 0))
if (nrow(progeny_heatmap_mat_qin_d3) > 0 && ncol(progeny_heatmap_mat_qin_d3) > 0) {
  colors_progeny_qin_d3 <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors_use_progeny_qin_d3 <- grDevices::colorRampPalette(colors = colors_progeny_qin_d3)(100)
  p_qin_d3_progeny_heat <- pheatmap::pheatmap(progeny_heatmap_mat_qin_d3, color = colors_use_progeny_qin_d3, border_color = "white", cellwidth = 20, cellheight = 20, main = "PROGENy Pathway Activity: Arg1+ vs Arg1- Neutrophils (Qin Day 3)")
  print(p_qin_d3_progeny_heat)
}
qin_day3_ligand_map_plot <- if (nrow(qin_day3_ligand_pathway_map) > 0 && "ligand" %in% names(qin_day3_ligand_pathway_map)) { top15_lig_qin_d3 <- head(unique(qin_day3_ligand_pathway_map$ligand), 15); qin_day3_ligand_pathway_map[qin_day3_ligand_pathway_map$ligand %in% top15_lig_qin_d3, ] } else data.frame()
p_qin_d3_progeny2 <- if (nrow(qin_day3_ligand_map_plot) > 0) tryCatch(ggplot(qin_day3_ligand_map_plot, aes(x = ligand, y = pathway, size = abs(weight), color = weight)) + geom_point(alpha = 0.7) + scale_color_gradient2(low = "blue", mid = "white", high = "red") + labs(title = "Top 15 Ligands Linked to PROGENy Pathways - Qin Day 3", x = "Ligand", y = "Pathway") + theme_minimal() + PLOT_TITLE_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1)), error = function(e) ggplot() + theme_void() + labs(title = "No data available")) else ggplot() + theme_void() + labs(title = "No ligands from LIANA (Qin Day 3)")
vals_w <- numeric(0)
map_tabnames <- c("lee_day1_ligand_pathway_map", "lee_day3_ligand_pathway_map", "wang_day3_ligand_pathway_map", "qin_day1_ligand_pathway_map", "qin_day3_ligand_pathway_map")
for (fn in map_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("weight" %in% names(d)) vals_w <- c(vals_w, d$weight)
}
vals_w <- vals_w[is.finite(vals_w)]
UNIFY_LIGW_ABS <- if (length(vals_w) > 0) range(abs(vals_w)) else c(0, 1)
if (UNIFY_LIGW_ABS[1] == UNIFY_LIGW_ABS[2]) UNIFY_LIGW_ABS[2] <- UNIFY_LIGW_ABS[1] + 1e-6
mxw <- if (length(vals_w) > 0) max(abs(vals_w)) else 1
if (!is.finite(mxw) || mxw <= 0) mxw <- 1
UNIFY_LIGW_COL <- c(-mxw, mxw)
print(p_qin_d3_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))
print("✓ Qin Day 3 PROGENy pathway analysis complete")


saveRDS(list(qin_day3 = qin_day3, qin_day3_neut_cells = qin_day3_neut_cells, qin_day3_pos = qin_day3_pos, qin_day3_neg = qin_day3_neg, qin_day3_pathway_results = qin_day3_pathway_results, progeny_network_qin_d3 = progeny_network_qin_d3, OUTPUT_DIR = OUTPUT_DIR), file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay3_PROGENy.rds"))
# ck_qin_d3_p <- readRDS(file.path(OUTPUT_DIR, "Checkpoint_AfterQinDay3_PROGENy.rds")); list2env(ck_qin_d3_p, envir = .GlobalEnv)

# -------- Qin Day 3 BARTsc TF Analysis (Arg1pos vs Arg1neg; downstream TF support only) --------
if (bartsc_initialized) {
  bartsc_initialized <- suppressWarnings(tryCatch({ BARTsc::load_bart2(); TRUE }, error = function(e) FALSE))
  if (!bartsc_initialized && exists("bart2", envir = .GlobalEnv)) bartsc_initialized <- TRUE
}
print("--- BARTsc TF Analysis: Qin Day 3 (Arg1pos vs Arg1neg neutrophils) ---")
if (bartsc_initialized) {
  options("mc.cores" = min(6L, parallel::detectCores()))
  qin_day3_neut_subset <- subset(qin_day3, cells = qin_day3_neut_cells)
  qin_day3_bart_label <- setNames(ifelse(colnames(qin_day3_neut_subset) %in% colnames(qin_day3)[qin_day3_pos], "Arg1pos", "Arg1neg"), colnames(qin_day3_neut_subset))
  qin_day3_bart_label <- factor(qin_day3_bart_label, levels = c("Arg1pos", "Arg1neg"))
  qin_day3_bart_proj <- BARTsc::bartsc(name = "QinDay3_Arg1", genome = "mm10", label = qin_day3_bart_label, cell_types_used = c("Arg1pos", "Arg1neg"), RNA_cnt_matrix = Seurat::GetAssayData(qin_day3_neut_subset, layer = "counts"))
  qin_day3_bart_proj <- BARTsc::normalize_RNA(qin_day3_bart_proj)
  qin_day3_bart_proj <- BARTsc::find_signature_genes(qin_day3_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  qin_day3_bart_proj <- BARTsc::find_pairwise_deg(qin_day3_bart_proj, min.pct = BART_MIN_PCT, min.diff.pct = BART_MIN_DIFF_PCT, log2fc.thr = BART_LOG2FC_THR, pval.thr = NULL, padj.thr = BART_PADJ_THR, auc.thr = BART_AUC_THR, max.cells.per.ident = Inf)
  n_qin_d3_posneg <- 0L
  n_qin_d3_negpos <- 0L
  qin_d3_pw <- qin_day3_bart_proj@data$pairwise_DEG
  if (!is.null(qin_d3_pw) && "Arg1pos::Arg1neg" %in% names(qin_d3_pw)) n_qin_d3_posneg <- { x <- qin_d3_pw[["Arg1pos::Arg1neg"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  if (!is.null(qin_d3_pw) && "Arg1neg::Arg1pos" %in% names(qin_d3_pw)) n_qin_d3_negpos <- { x <- qin_d3_pw[["Arg1neg::Arg1pos"]]; if (is.data.frame(x)) as.integer(nrow(x)) else as.integer(length(x)) }
  print(paste0("Qin Day 3 BARTsc pairwise DEG count (@data$pairwise_DEG): Arg1pos::Arg1neg=", n_qin_d3_posneg, ", Arg1neg::Arg1pos=", n_qin_d3_negpos))
  bart_deg_ok_qin_d3 <- if (BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT > 0) n_qin_d3_posneg >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT && n_qin_d3_negpos >= BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT else TRUE
  bart_nt_qin_d3 <- table(qin_day3_bart_label)
  bart_unequal_qin_d3 <- length(bart_nt_qin_d3) == 2L && length(unique(as.vector(bart_nt_qin_d3))) > 1L
  bart_crossct_allowed_qin_d3 <- !isTRUE(BARTSC_SKIP_CROSSCT_TEST) && length(bart_nt_qin_d3) >= 2L && all(as.integer(bart_nt_qin_d3) >= BARTSC_MIN_CELLS_CROSSCT) && !(isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_qin_d3) && bart_deg_ok_qin_d3
  qin_day3_bart_proj <- BARTsc::run_signature_RNA(qin_day3_bart_proj)
  qin_day3_bart_sig <- BARTsc::get_result(qin_day3_bart_proj, analysis = "cell type signature", mod = "RNA")
  qin_day3_bart_proj <- BARTsc::calc_crossCT_auc_RNA(qin_day3_bart_proj)
  print(bart_nt_qin_d3)
  if (!bart_crossct_allowed_qin_d3) message(paste0(
    "BARTsc Qin Day 3: crossCT_test / find_key_regulators skipped (BARTSC_SKIP_CROSSCT_TEST=", isTRUE(BARTSC_SKIP_CROSSCT_TEST),
    "; min_cells_ok=", all(as.integer(bart_nt_qin_d3) >= BARTSC_MIN_CELLS_CROSSCT),
    "; skip_if_unequal_n=", isTRUE(BARTSC_SKIP_CROSSCT_IF_UNEQUAL_ARG1_N) && bart_unequal_qin_d3,
    "; deg_ok=", bart_deg_ok_qin_d3,
    " [Arg1pos::Arg1neg=", n_qin_d3_posneg,
    ", Arg1neg::Arg1pos=", n_qin_d3_negpos,
    ", min_each_dir=", BART_MIN_DEG_EACH_DIRECTION_FOR_CROSSCT,
    "]). Vignette scRNA-seq.md §6: dot_plot/deviation_heatmap follow crossCT_test."
  ))
  if (bart_crossct_allowed_qin_d3) qin_day3_bart_proj <- BARTsc::crossCT_test(qin_day3_bart_proj, mod = "RNA")
  qin_day3_bart_cross <- if (bart_crossct_allowed_qin_d3) BARTsc::get_result(qin_day3_bart_proj, analysis = "cross-cell-type", mod = "RNA") else list()
  qin_day3_bart_cross_dev <- if (is.null(qin_day3_bart_cross)) list() else if (is.list(qin_day3_bart_cross) && "deviation" %in% names(qin_day3_bart_cross) && is.list(qin_day3_bart_cross$deviation)) qin_day3_bart_cross$deviation else qin_day3_bart_cross
  bart_outdir_qin_d3 <- file.path(OUTPUT_DIR, "BARTsc_QinDay3")
  dir.create(bart_outdir_qin_d3, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(qin_day3_bart_sig) && !is.null(qin_day3_bart_sig[["Arg1pos"]])) write.csv(qin_day3_bart_sig[["Arg1pos"]], file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_Signature_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(qin_day3_bart_sig) && !is.null(qin_day3_bart_sig[["Arg1neg"]])) write.csv(qin_day3_bart_sig[["Arg1neg"]], file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_Signature_Arg1neg.csv"), row.names = FALSE)
  qin_day3_bart_key <- NULL
  if (bart_crossct_allowed_qin_d3) {
    qin_day3_bart_proj <- BARTsc::find_key_regulators(qin_day3_bart_proj, mod = "RNA", min.N.profile = 3)
    qin_day3_bart_key <- BARTsc::get_result(qin_day3_bart_proj, analysis = "Key regs ident", mod = "RNA")
  }
  if (!is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1pos"]])) write.csv(qin_day3_bart_key[["Arg1pos"]], file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_KeyRegulators_Arg1pos.csv"), row.names = FALSE)
  if (!is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1neg"]])) write.csv(qin_day3_bart_key[["Arg1neg"]], file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_KeyRegulators_Arg1neg.csv"), row.names = FALSE)
  saveRDS(qin_day3_bart_proj, file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_Object.rds"))
  print("✓ Qin Day 3 BARTsc TF analysis complete — plots: Qin Day 3 BARTsc visualizations section below")
}

# -------- Qin Day 3: BARTsc visualizations (reload QinDay3_BARTsc_Object.rds; no re-run of bartsc pipeline) --------
bart_viz_qd3_rds <- file.path(OUTPUT_DIR, "BARTsc_QinDay3", "QinDay3_BARTsc_Object.rds")
if (!requireNamespace("BARTsc", quietly = TRUE)) message("Qin Day 3 BARTsc visualizations skipped: package BARTsc not installed.")
if (requireNamespace("BARTsc", quietly = TRUE) && !file.exists(bart_viz_qd3_rds)) message(paste0("Qin Day 3 BARTsc visualizations skipped: missing ", bart_viz_qd3_rds))
if (requireNamespace("BARTsc", quietly = TRUE) && file.exists(bart_viz_qd3_rds)) {
  if (!exists("BARTSC_N_TF_DOTHEAT")) BARTSC_N_TF_DOTHEAT <- 5L
  if (!exists("BARTSC_N_LABELED_TFS")) BARTSC_N_LABELED_TFS <- 5L
  if (!exists("BARTSC_KEYREG_SCATTER_W_IN")) BARTSC_KEYREG_SCATTER_W_IN <- 11
  if (!exists("BARTSC_KEYREG_SCATTER_H_IN")) BARTSC_KEYREG_SCATTER_H_IN <- 7
  if (!exists("BARTSC_KEYREG_XLIM")) BARTSC_KEYREG_XLIM <- c(-1, 1)
  if (!exists("BARTSC_KEYREG_YLIM_SIG")) BARTSC_KEYREG_YLIM_SIG <- c(0, 7)
  if (!exists("BARTSC_KEYREG_ZLIM_DMDR")) BARTSC_KEYREG_ZLIM_DMDR <- c(-2, 2)
  if (!exists("BARTSC_KEYREG_RANK_MAX")) BARTSC_KEYREG_RANK_MAX <- 60L
  suppressWarnings(tryCatch({ BARTsc::load_bart2(); NULL }, error = function(e) NULL))
  if (exists("bart2", envir = .GlobalEnv)) {
    if (!exists("types", envir = .GlobalEnv)) types <<- reticulate::import("types")
    bart2_mods_viz_qd3 <- reticulate::py_to_r(reticulate::py_get_attr(bart2, "__all__"))
    for (m_viz_qd3 in bart2_mods_viz_qd3) {
      if (!exists(m_viz_qd3, envir = .GlobalEnv)) assign(m_viz_qd3, reticulate::import(paste0("bart2.", m_viz_qd3), delay_load = TRUE), envir = .GlobalEnv)
    }
  }
  if (!exists("qin_day3_bart_proj", envir = .GlobalEnv, inherits = FALSE)) qin_day3_bart_proj <- readRDS(bart_viz_qd3_rds)
  bart_outdir_qin_d3 <- file.path(OUTPUT_DIR, "BARTsc_QinDay3")
  dir.create(bart_outdir_qin_d3, showWarnings = FALSE, recursive = TRUE)
  if (interactive() && grDevices::dev.cur() == 1L) grDevices::dev.new()
  qin_day3_bart_cross_viz <- BARTsc::get_result(qin_day3_bart_proj, analysis = "cross-cell-type", mod = "RNA")
  qin_day3_bart_cross_dev <- if (is.null(qin_day3_bart_cross_viz)) list() else if (is.list(qin_day3_bart_cross_viz) && "deviation" %in% names(qin_day3_bart_cross_viz) && is.list(qin_day3_bart_cross_viz$deviation)) qin_day3_bart_cross_viz$deviation else qin_day3_bart_cross_viz
  qin_day3_bart_key <- BARTsc::get_result(qin_day3_bart_proj, analysis = "Key regs ident", mod = "RNA")
  bart_tf_names_qin_d3 <- names(qin_day3_bart_cross_dev)
  tf_example_qin_d3 <- if (length(bart_tf_names_qin_d3) > 0) bart_tf_names_qin_d3[1] else NULL
  if (length(qin_day3_bart_cross_dev) > 0 && !is.null(tf_example_qin_d3)) {
    p_qin_d3_bart_dot <- BARTsc::dot_plot(qin_day3_bart_proj, mod = "RNA", tf = tf_example_qin_d3, max_dot_size = 22)
    print(p_qin_d3_bart_dot)
    ggplot2::ggsave(file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_DotPlot.png"), p_qin_d3_bart_dot, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_DotPlot.pdf"), p_qin_d3_bart_dot, width = 8, height = 6)
    p_qin_d3_bart_heat <- BARTsc::deviation_heatmap(qin_day3_bart_proj, mod = "RNA", tf = tf_example_qin_d3, tile_fontsize = 6)
    print(p_qin_d3_bart_heat)
    ggplot2::ggsave(file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_DeviationHeatmap.png"), p_qin_d3_bart_heat, width = 8, height = 6)
    ggplot2::ggsave(file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_DeviationHeatmap.pdf"), p_qin_d3_bart_heat, width = 8, height = 6)
  }
  top_tf_qin_d3_pos <- if (length(qin_day3_bart_cross_dev) > 0 && !is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1pos"]]) && "TF" %in% colnames(qin_day3_bart_key[["Arg1pos"]])) head(qin_day3_bart_key[["Arg1pos"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_tf_qin_d3_neg <- if (length(qin_day3_bart_cross_dev) > 0 && !is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1neg"]]) && "TF" %in% colnames(qin_day3_bart_key[["Arg1neg"]])) head(qin_day3_bart_key[["Arg1neg"]]$TF, BARTSC_N_TF_DOTHEAT) else character(0)
  top_key_tfs_qin_d3 <- intersect(unique(c(top_tf_qin_d3_pos, top_tf_qin_d3_neg)), names(qin_day3_bart_cross_dev))
  if (length(qin_day3_bart_cross_dev) > 0 && length(top_key_tfs_qin_d3) > 0) {
    i_tf_qd3 <- 1L
    while (i_tf_qd3 <= length(top_key_tfs_qin_d3)) {
      tf_cur_qin_d3 <- top_key_tfs_qin_d3[i_tf_qd3]
      p_dot_qd3 <- BARTsc::dot_plot(qin_day3_bart_proj, mod = "RNA", tf = tf_cur_qin_d3, max_dot_size = 22)
      print(p_dot_qd3)
      ggplot2::ggsave(file.path(bart_outdir_qin_d3, paste0("QinDay3_BARTsc_DotPlot_", tf_cur_qin_d3, ".png")), p_dot_qd3, width = 8, height = 6)
      p_heat_qd3 <- BARTsc::deviation_heatmap(qin_day3_bart_proj, mod = "RNA", tf = tf_cur_qin_d3, tile_fontsize = 6)
      print(p_heat_qd3)
      ggplot2::ggsave(file.path(bart_outdir_qin_d3, paste0("QinDay3_BARTsc_DeviationHeatmap_", tf_cur_qin_d3, ".png")), p_heat_qd3, width = 8, height = 6)
      i_tf_qd3 <- i_tf_qd3 + 1L
    }
  }
  tfs_qin_d3_pos <- character(0)
  if (!is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1pos"]])) {
    df_q3p <- qin_day3_bart_key[["Arg1pos"]]
    if (nrow(df_q3p) > 0) {
      tf_c_q3p <- intersect(colnames(df_q3p), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_q3p <- if (length(tf_c_q3p) == 0) colnames(df_q3p)[1] else tf_c_q3p[1]
      ord_q3p <- if ("final_rank" %in% colnames(df_q3p)) order(df_q3p$final_rank, na.last = TRUE) else seq_len(nrow(df_q3p))
      tfs_qin_d3_pos <- as.character(head(df_q3p[ord_q3p, , drop = FALSE][[tf_col_q3p]], BARTSC_N_LABELED_TFS))
    }
  }
  tfs_qin_d3_neg <- character(0)
  if (!is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1neg"]])) {
    df_q3n <- qin_day3_bart_key[["Arg1neg"]]
    if (nrow(df_q3n) > 0) {
      tf_c_q3n <- intersect(colnames(df_q3n), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
      tf_col_q3n <- if (length(tf_c_q3n) == 0) colnames(df_q3n)[1] else tf_c_q3n[1]
      ord_q3n <- if ("final_rank" %in% colnames(df_q3n)) order(df_q3n$final_rank, na.last = TRUE) else seq_len(nrow(df_q3n))
      tfs_qin_d3_neg <- as.character(head(df_q3n[ord_q3n, , drop = FALSE][[tf_col_q3n]], BARTSC_N_LABELED_TFS))
    }
  }
  if (length(tfs_qin_d3_pos) > 0) {
    grDevices::png(file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_KeyRegScatter_Official_Arg1pos.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(qin_day3_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_qin_d3_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_qin_d3_pos) > 0) message("Displaying official key_regulator_scatter for Arg1pos (Qin Day 3)...")
  if (interactive() && length(tfs_qin_d3_pos) > 0) key_regulator_scatter_unified(qin_day3_bart_proj, mod = "RNA", cell_type = "Arg1pos", tfs_labeled = tfs_qin_d3_pos, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1pos", subtitle = "Anti-inflammatory hypothesis (exploratory)")
  if (length(tfs_qin_d3_neg) > 0) {
    grDevices::png(file.path(bart_outdir_qin_d3, "QinDay3_BARTsc_KeyRegScatter_Official_Arg1neg.png"), width = BARTSC_KEYREG_SCATTER_W_IN, height = BARTSC_KEYREG_SCATTER_H_IN, units = "in", res = 150)
    key_regulator_scatter_unified(qin_day3_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_qin_d3_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
    grDevices::dev.off()
  }
  if (interactive() && length(tfs_qin_d3_neg) > 0) message("Displaying official key_regulator_scatter for Arg1neg (Qin Day 3)...")
  if (interactive() && length(tfs_qin_d3_neg) > 0) key_regulator_scatter_unified(qin_day3_bart_proj, mod = "RNA", cell_type = "Arg1neg", tfs_labeled = tfs_qin_d3_neg, xlim = BARTSC_KEYREG_XLIM, ylim_sig = BARTSC_KEYREG_YLIM_SIG, zlim_dmdr = BARTSC_KEYREG_ZLIM_DMDR, rank_color_max = BARTSC_KEYREG_RANK_MAX, main = "Arg1neg", subtitle = "Pro-inflammatory hypothesis (exploratory)")
  print("✓ Qin Day 3 BARTsc visualizations complete")
}

# ============================================================================
# INTEGRATION: LIANA + PROGENy footprint co-membership (exploratory) + OmniPath receptor→TF (Qin Day 3)
# ============================================================================
# Linear: footprint co-membership CSV → OmniPath receptor→TF → full-chain CSV → sender×pathway heatmap (BARTsc plots: section above).
print("--- Integration: LIANA + PROGENy footprint overlap (exploratory) + BARTsc/OmniPath (Qin Day 3) ---")
# Exploratory overlay only: LIANA provides incoming-signal candidates; PROGENy/BARTsc/OmniPath add downstream support, not causal proof.
if (!exists("progeny_network_qin_d3") || !all(c("source", "target") %in% colnames(progeny_network_qin_d3)) || !exists("qin_day3_arg1pos_received") || nrow(qin_day3_arg1pos_received) == 0) {
  message("Qin Day 3 integration skipped (need PROGENy network + LIANA Arg1pos received with rows).")
} else {
liana_receptors_qin_d3 <- unique(qin_day3_arg1pos_received$receptor.complex)
liana_receptors_qin_d3_split <- unique(trimws(unlist(strsplit(liana_receptors_qin_d3, "[_+]"))))
print(paste("Unique receptors received by Arg1+ neutrophils (LIANA, Qin Day 3):", length(liana_receptors_qin_d3_split)))
liana_expanded_qin_d3 <- qin_day3_arg1pos_received
liana_expanded_qin_d3$receptor_subunits <- lapply(strsplit(liana_expanded_qin_d3$receptor.complex, "[_+]"), function(x) trimws(x))
liana_expanded_qin_d3 <- tidyr::unnest_longer(liana_expanded_qin_d3, receptor_subunits)
pathways_active_qin_d3 <- names(qin_day3_pathway_means_arg1pos)[qin_day3_pathway_means_arg1pos > qin_day3_pathway_means_arg1neg]
# Exploratory only: PROGENy links pathway -> gene targets (footprints), not receptor -> pathway activation.
receptor_pathway_links_qin_d3 <- progeny_network_qin_d3[progeny_network_qin_d3$target %in% liana_receptors_qin_d3_split & progeny_network_qin_d3$source %in% pathways_active_qin_d3, ]
receptor_pathway_links_qin_d3 <- receptor_pathway_links_qin_d3[, c("target", "source")]
colnames(receptor_pathway_links_qin_d3) <- c("receptor_gene", "pathway")
pathway_diff_qin_d3 <- qin_day3_pathway_means_arg1pos - qin_day3_pathway_means_arg1neg
receptor_pathway_links_qin_d3$pathway_activity_diff <- pathway_diff_qin_d3[receptor_pathway_links_qin_d3$pathway]
integration_receptor_pathway_qin_d3 <- merge(liana_expanded_qin_d3, receptor_pathway_links_qin_d3, by.x = "receptor_subunits", by.y = "receptor_gene", all.x = TRUE)
integration_receptor_pathway_qin_d3 <- integration_receptor_pathway_qin_d3[!is.na(integration_receptor_pathway_qin_d3$pathway), ]
write.csv(integration_receptor_pathway_qin_d3, file.path(OUTPUT_DIR, "QinDay3_LIANA_PROGENy_FootprintCoMembership_Exploratory.csv"), row.names = FALSE)
print(paste("LIANA x PROGENy footprint co-membership (exploratory, Qin Day 3):", nrow(integration_receptor_pathway_qin_d3), "rows (not causal receptor->pathway)"))
receptor_tf_links_qin_d3 <- data.frame(receptor_gene = character(0), tf = character(0), omnipath_via = character(0), omnipath_hops = integer(0), stringsAsFactors = FALSE)
active_tfs_qin_d3 <- character(0)
has_bart_key_qin_d3 <- exists("qin_day3_bart_key") && !is.null(qin_day3_bart_key) && !is.null(qin_day3_bart_key[["Arg1pos"]]) && nrow(qin_day3_bart_key[["Arg1pos"]]) > 0
if (has_bart_key_qin_d3) {
  bart_tfs_qin_d3 <- qin_day3_bart_key[["Arg1pos"]]
  tf_col_qin_d3 <- intersect(colnames(bart_tfs_qin_d3), c("tf", "TF", "regulator", "Regulator", "gene", "Gene"))
  if (length(tf_col_qin_d3) == 0) tf_col_qin_d3 <- colnames(bart_tfs_qin_d3)[1] else tf_col_qin_d3 <- tf_col_qin_d3[1]
  active_tfs_qin_d3 <- unique(as.character(bart_tfs_qin_d3[[tf_col_qin_d3]]))
}
omnipath_ok_qin_d3 <- requireNamespace("OmnipathR", quietly = TRUE) && length(active_tfs_qin_d3) > 0

pathway_interactions_qin_d3 <- data.frame()
omni_qin_d3_first <- if (!omnipath_ok_qin_d3) list(dat = NULL, err1 = NA_character_) else tryCatch(list(dat = OmnipathR::import_omnipath_interactions(datasets = c("omnipath", "pathwayextra"), organism = 10090, genesymbols = TRUE), err1 = NA_character_), error = function(e1) list(dat = NULL, err1 = conditionMessage(e1)))
if (omnipath_ok_qin_d3) pathway_interactions_qin_d3 <- omni_qin_d3_first$dat
omnipath_qin_d3_err1 <- omni_qin_d3_first$err1
if (omnipath_ok_qin_d3 && is.null(pathway_interactions_qin_d3)) pathway_interactions_qin_d3 <- tryCatch(OmnipathR::import_pathwayextra_interactions(organism = 10090, genesymbols = TRUE), error = function(e2) { message("OmniPath curated+pathwayextra import failed (Qin Day 3): ", omnipath_qin_d3_err1, "; fallback: ", conditionMessage(e2)); data.frame() })

omni_qin_d3_cols_ok <- omnipath_ok_qin_d3 && nrow(pathway_interactions_qin_d3) > 0 && all(c("source_genesymbol", "target_genesymbol") %in% colnames(pathway_interactions_qin_d3))
omni_direct_qin_d3 <- if (omni_qin_d3_cols_ok) subset(pathway_interactions_qin_d3, source_genesymbol %in% liana_receptors_qin_d3_split & target_genesymbol %in% active_tfs_qin_d3, select = c("source_genesymbol", "target_genesymbol")) else data.frame()
if (omni_qin_d3_cols_ok && nrow(omni_direct_qin_d3) > 0) receptor_tf_links_qin_d3 <- rbind(receptor_tf_links_qin_d3, data.frame(receptor_gene = omni_direct_qin_d3$source_genesymbol, tf = omni_direct_qin_d3$target_genesymbol, omnipath_via = NA_character_, omnipath_hops = 1L, stringsAsFactors = FALSE))

hop1_qin_d3 <- if (omni_qin_d3_cols_ok) pathway_interactions_qin_d3[pathway_interactions_qin_d3$source_genesymbol %in% liana_receptors_qin_d3_split, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_qin_d3_cols_ok) colnames(hop1_qin_d3) <- c("receptor_gene", "via")
hop2_qin_d3 <- if (omni_qin_d3_cols_ok) pathway_interactions_qin_d3[pathway_interactions_qin_d3$target_genesymbol %in% active_tfs_qin_d3, c("source_genesymbol", "target_genesymbol"), drop = FALSE] else data.frame()
if (omni_qin_d3_cols_ok) colnames(hop2_qin_d3) <- c("via", "tf")
omni_2hop_qin_d3 <- if (omni_qin_d3_cols_ok) merge(hop1_qin_d3, hop2_qin_d3, by = "via") else data.frame()
if (omni_qin_d3_cols_ok && nrow(omni_2hop_qin_d3) > 0) omni_2hop_qin_d3 <- unique(omni_2hop_qin_d3[, c("receptor_gene", "tf", "via")])
if (omni_qin_d3_cols_ok && nrow(omni_2hop_qin_d3) > 0) receptor_tf_links_qin_d3 <- rbind(receptor_tf_links_qin_d3, data.frame(receptor_gene = omni_2hop_qin_d3$receptor_gene, tf = omni_2hop_qin_d3$tf, omnipath_via = omni_2hop_qin_d3$via, omnipath_hops = 2L, stringsAsFactors = FALSE))
if (omni_qin_d3_cols_ok) receptor_tf_links_qin_d3 <- receptor_tf_links_qin_d3[!duplicated(paste(receptor_tf_links_qin_d3$receptor_gene, receptor_tf_links_qin_d3$tf)), ]

integration_receptor_tf_qin_d3 <- data.frame()
if (nrow(receptor_tf_links_qin_d3) > 0) {
  integration_receptor_tf_qin_d3 <- merge(liana_expanded_qin_d3, receptor_tf_links_qin_d3, by.x = "receptor_subunits", by.y = "receptor_gene", all.x = TRUE)
  integration_receptor_tf_qin_d3 <- integration_receptor_tf_qin_d3[!is.na(integration_receptor_tf_qin_d3$tf), ]
}
if (nrow(integration_receptor_tf_qin_d3) > 0) {
  write.csv(integration_receptor_tf_qin_d3, file.path(OUTPUT_DIR, "QinDay3_LIANA_BARTsc_Integration_DataDriven.csv"), row.names = FALSE)
  n1_qd3 <- sum(integration_receptor_tf_qin_d3$omnipath_hops == 1L, na.rm = TRUE)
  n2_qd3 <- sum(integration_receptor_tf_qin_d3$omnipath_hops == 2L, na.rm = TRUE)
  print(paste0("LIANA -> BARTsc integration (Qin Day 3): ", nrow(integration_receptor_tf_qin_d3), " receptor-TF rows (OmniPath curated+pathwayextra; direct=", n1_qd3, ", two-hop=", n2_qd3, ")"))
}
if (nrow(receptor_tf_links_qin_d3) == 0 && length(active_tfs_qin_d3) > 0) {
  write.csv(data.frame(receptor = liana_receptors_qin_d3_split, note = "Active TFs in Arg1pos (no OmniPath direct/two-hop link):", active_tfs = paste(active_tfs_qin_d3, collapse = "; "), stringsAsFactors = FALSE), file.path(OUTPUT_DIR, "QinDay3_LIANA_BARTsc_Receptors_and_TFs.csv"), row.names = FALSE)
  print("LIANA receptors and BARTsc TFs saved separately (Qin Day 3, no OmniPath link)")
}
has_pathway_qin_d3 <- nrow(integration_receptor_pathway_qin_d3) > 0
has_tf_qin_d3 <- nrow(integration_receptor_tf_qin_d3) > 0
full_chain_qin_d3 <- data.frame()
if (has_pathway_qin_d3 && has_tf_qin_d3) full_chain_qin_d3 <- merge(integration_receptor_pathway_qin_d3, integration_receptor_tf_qin_d3[, c("source", "ligand.complex", "receptor_subunits", "aggregate_rank", "tf", "omnipath_via", "omnipath_hops")], by = c("source", "ligand.complex", "receptor_subunits", "aggregate_rank"), all = TRUE)
if (has_pathway_qin_d3 && !has_tf_qin_d3) { full_chain_qin_d3 <- integration_receptor_pathway_qin_d3; full_chain_qin_d3$tf <- NA_character_ }
if (!has_pathway_qin_d3 && has_tf_qin_d3) { full_chain_qin_d3 <- integration_receptor_tf_qin_d3; full_chain_qin_d3$pathway <- NA_character_; full_chain_qin_d3$pathway_activity_diff <- NA_real_ }
if (nrow(full_chain_qin_d3) > 0) full_chain_qin_d3$evidence_level <- ifelse(!is.na(full_chain_qin_d3$pathway) & !is.na(full_chain_qin_d3$tf), "STRONG (pathway + TF linked)", ifelse(!is.na(full_chain_qin_d3$pathway) | !is.na(full_chain_qin_d3$tf), "MODERATE (pathway or TF linked)", "WEAK (no link)"))
if (nrow(full_chain_qin_d3) > 0) full_chain_qin_d3 <- full_chain_qin_d3[order(full_chain_qin_d3$evidence_level, full_chain_qin_d3$aggregate_rank), ]
# Legacy filename retained for compatibility; contents are an exploratory overlay, not a causal signal chain.
if (nrow(full_chain_qin_d3) > 0) write.csv(full_chain_qin_d3, file.path(OUTPUT_DIR, "QinDay3_FullSignalChain_DataDriven.csv"), row.names = FALSE)
if (nrow(full_chain_qin_d3) > 0) {
  cols_show_qin_d3 <- intersect(c("source", "receptor_subunits", "pathway", "tf", "omnipath_via", "omnipath_hops", "evidence_level"), colnames(full_chain_qin_d3))
  print(head(full_chain_qin_d3[, cols_show_qin_d3, drop = FALSE], 20))
}
heatmap_data_qin_d3 <- data.frame()
if (nrow(full_chain_qin_d3) > 0 && has_pathway_qin_d3) heatmap_data_qin_d3 <- as.data.frame.matrix(table(integration_receptor_pathway_qin_d3$source, integration_receptor_pathway_qin_d3$pathway))
if (nrow(heatmap_data_qin_d3) > 0 && ncol(heatmap_data_qin_d3) > 0) pheatmap::pheatmap(heatmap_data_qin_d3, cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette(c("white", "blue", "red"))(50), main = "Qin Day 3: Sender x Pathway (LIANA receptor x PROGENy footprint co-membership; exploratory)")
if (nrow(full_chain_qin_d3) == 0) print("No integration rows Qin Day 3 (PROGENy/BARTsc may not overlap LIANA receptors)")
print("✓ Qin Day 3 LIANA <-> BARTsc <-> PROGENy integration complete (exploratory overlay)")
}

print("✓✓✓ STEP 0 COMPLETE: LIANA CONSENSUS ANALYSIS ✓✓✓")
print("")

# ============================================================================
# SECTION 5: STEP 1 - LEE DAY 1 DATA OVERVIEW & QC PLOTS
# ============================================================================
# LIANA plots (1F-1O): in Lee Day 1 LIANA section above
# PROGENy plots (1X-1Y): in Lee Day 1 PROGENy section above
# BARTsc plots: per-cohort "BARTsc visualizations" sections (after each BARTsc analysis block; can rerun from *BARTsc_Object.rds)
# Below: Arg1 validation, UMAP, DEG, population, feature, pathway overview

print(">>> STEP 1: LEE DAY 1 DATA OVERVIEW & QC PLOTS <<<")

# Note: lee_day1 and lee_day1_arg1_status are already defined above (before LIANA section)
# They are reused here for visualization

# -------- PLOT 1A: Arg1 Validation Violin --------
print("Generating: Plot 1A - Arg1 Validation Violin (Lee Day 1)")
p_lee_d1_1a <- VlnPlot(
  lee_day1, 
  features = "Arg1", 
  group.by = "arg1_status",
  pt.size = 0.1, 
  cols = c("Arg1neg" = "grey", "Arg1pos" = "red")
) + 
  labs(title = "Validation: Arg1 Expression - Lee Day 1", y = "Log Normalized Expression") +
  PLOT_TITLE_THEME

print(p_lee_d1_1a)

# -------- PLOT 1B: UMAP showing Arg1 distribution --------
print("Generating: Plot 1B - UMAP Arg1 Distribution (Lee Day 1)")
p_lee_d1_1b <- DimPlot(
  lee_day1, 
  group.by = "arg1_status", 
  cols = c("Arg1neg" = "lightblue", "Arg1pos" = "red"),
  pt.size = 2
) + 
  labs(title = "Arg1pos vs Arg1neg Distribution - Lee Day 1") +
  PLOT_TITLE_THEME

print(p_lee_d1_1b)

# -------- PLOT 1C: Cell type composition --------
print("Generating: Plot 1C - Cell Type Composition (Lee Day 1)")
comp_lee_d1 <- as.data.frame(table(lee_day1_pruned_labels))
colnames(comp_lee_d1) <- c("CellType", "Count")
comp_lee_d1 <- comp_lee_d1[order(comp_lee_d1$Count, decreasing = TRUE), ]

p_lee_d1_1c <- ggplot(comp_lee_d1, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Cell Type Composition - Lee Day 1", x = "Cell Type", y = "Number of Cells") +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_lee_d1_1c)

# -------- PLOT 1D: DEG Volcano (Arg1pos vs Arg1neg neutrophils) --------
print("Generating: Plot 1D - DEG Volcano (Lee Day 1)")
lee_day1_neutrophils <- subset(lee_day1, cells = lee_day1_neut_cells)
Idents(lee_day1_neutrophils) <- lee_day1_neutrophils$arg1_status
lee_day1_deg <- FindMarkers(lee_day1_neutrophils, ident.1 = "Arg1pos", ident.2 = "Arg1neg", logfc.threshold = 0.25, min.pct = 0.1)
lee_day1_deg$gene <- rownames(lee_day1_deg)
# FindMarkers already returns p_val_adj - do NOT double-adjust
lee_day1_deg$significant <- (lee_day1_deg$p_val_adj < 0.05) & (abs(lee_day1_deg$avg_log2FC) > 0.5)

p_lee_d1_1d <- ggplot(lee_day1_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  geom_text_repel(
    data = head(lee_day1_deg[lee_day1_deg$significant, ], 10),
    aes(label = gene),
    size = 3
  ) +
  labs(
    title = "DEG Volcano Plot - Lee Day 1 (Arg1pos vs Arg1neg)",
    x = "Log2 Fold Change",
    y = "-log10 (FDR-adjusted p-value)"
  ) +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_lee_d1_1d)
write.csv(lee_day1_deg, file.path(OUTPUT_DIR, "01d_LeeDay1_DEG_Results.csv"))

# -------- PLOT 1E: Heatmap top DEGs --------
print("Generating: Plot 1E - DEG Heatmap (Lee Day 1)")
top_deg_lee_d1 <- head(lee_day1_deg[lee_day1_deg$significant, ], 20)
expr_matrix_lee_d1 <- Seurat::GetAssayData(lee_day1_neutrophils, layer = "data")[top_deg_lee_d1$gene, , drop = FALSE]
# Scale genes (rows) before plotting.
expr_matrix_lee_d1_scaled <- t(scale(t(expr_matrix_lee_d1)))

p_lee_d1_1e <- pheatmap::pheatmap(
  expr_matrix_lee_d1_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 20 DEGs - Lee Day 1 (Arg1pos vs Arg1neg)",
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 10
)

print(p_lee_d1_1e)

# Lee Day 1 LIANA plots 1F-1O: generated in Lee Day 1 LIANA section above

# -------- PLOT 1P: Lee Day 1 Data Overview - Violin plots marker genes --------
print("Generating: Plot 1P - Marker Gene Violin Plot (Lee Day 1)")
marker_genes_lee_d1 <- c("Cd4", "Cd8a", "Iba1", "Gfap", "Pdgfra")
marker_genes_present <- marker_genes_lee_d1[marker_genes_lee_d1 %in% rownames(lee_day1)]
p_lee_d1_1p <- tryCatch(VlnPlot(lee_day1, features = marker_genes_present[1:min(3, length(marker_genes_present))], ncol = 1) + labs(title = "Marker Gene Expression - Lee Day 1") + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1p)

# -------- PLOT 1Q: Feature plots for selected genes --------
print("Generating: Plot 1Q - Feature Plots (Lee Day 1)")
feature_genes_lee_d1 <- c("Arg1", "Il10", "Tnf")
feature_genes_present <- feature_genes_lee_d1[feature_genes_lee_d1 %in% rownames(lee_day1)]
p_lee_d1_1q <- tryCatch(FeaturePlot(lee_day1, features = feature_genes_present[1:min(3, length(feature_genes_present))], ncol = min(3, length(feature_genes_present)), pt.size = 2), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1q)

# -------- PLOT 1R: Arg1pos population size comparison --------
print("Generating: Plot 1R - Arg1pos Population Size (Lee Day 1)")
lee_day1_neut_total <- length(lee_day1_pos) + length(lee_day1_neg)
pop_size_lee_d1 <- data.frame(
  Status = c("Arg1pos", "Arg1neg"),
  Count = c(length(lee_day1_pos), length(lee_day1_neg)),
  Percentage = c(
    length(lee_day1_pos) / lee_day1_neut_total * 100,
    length(lee_day1_neg) / lee_day1_neut_total * 100
  )
)

p_lee_d1_1r <- ggplot(pop_size_lee_d1, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), vjust = -0.5) +
  scale_fill_manual(values = c("Arg1pos" = "red", "Arg1neg" = "lightblue")) +
  labs(title = "Arg1pos Population - Lee Day 1", x = "Neutrophil Status", y = "Cell Count") +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_lee_d1_1r)

# -------- PLOT 1S: Top 20 DEGs expression profile --------
print("Generating: Plot 1S - DEG Expression Profile (Lee Day 1)")
top_deg_genes_lee_d1 <- tryCatch(head(lee_day1_deg[lee_day1_deg$significant, ], 5)$gene, error = function(e) character(0))
expr_profile_lee_d1_step1 <- tryCatch(data.frame(gene = top_deg_genes_lee_d1, Arg1pos = rowMeans(Seurat::GetAssayData(lee_day1, layer = "data")[top_deg_genes_lee_d1, lee_day1_arg1pos_cells, drop = FALSE]), Arg1neg = rowMeans(Seurat::GetAssayData(lee_day1, layer = "data")[top_deg_genes_lee_d1, lee_day1_arg1neg_cells, drop = FALSE])), error = function(e) data.frame())
expr_profile_lee_d1 <- tryCatch(tidyr::pivot_longer(expr_profile_lee_d1_step1, cols = -gene, names_to = "Status", values_to = "Expression"), error = function(e) data.frame())
p_lee_d1_1s <- tryCatch(ggplot(expr_profile_lee_d1, aes(x = gene, y = Expression, fill = Status)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("Arg1pos" = "red", "Arg1neg" = "lightblue")) + coord_flip() + labs(title = "Top 5 DEG Expression Comparison - Lee Day 1", x = "Gene", y = "Mean Expression") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1s)

# -------- PLOT 1T: Cell cycle status --------
print("Generating: Plot 1T - Cell Cycle Analysis (Lee Day 1)")
phase_comp_lee_d1 <- tryCatch({phase_df <- as.data.frame(table(lee_day1$Phase[c(lee_day1_pos, lee_day1_neg)], lee_day1_arg1_status[c(lee_day1_pos, lee_day1_neg)])); colnames(phase_df) <- c("Phase", "Arg1Status", "Count"); phase_df}, error = function(e) data.frame())
p_lee_d1_1t <- tryCatch(ggplot(phase_comp_lee_d1, aes(x = Arg1Status, y = Count, fill = Phase)) + geom_bar(stat = "identity", position = "fill") + labs(title = "Cell Cycle Phase Distribution - Lee Day 1", x = "Arg1 Status", y = "Proportion") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1t)

# -------- PLOT 1U: Pathway analysis signatures --------
print("Generating: Plot 1U - Pathway Signature Analysis (Lee Day 1)")
pathway_list <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", collection = "H"), error = function(e) data.frame())
pathway_names <- tryCatch(unique(pathway_list$gs_name)[1:5], error = function(e) character(0))
score_list <- list()
for (pw in pathway_names) {
  genes <- tryCatch(pathway_list$gene_symbol[pathway_list$gs_name == pw], error = function(e) character(0))
  genes_present <- genes[genes %in% rownames(lee_day1)]
  score_list[[pw]] <- tryCatch(colMeans(Seurat::GetAssayData(lee_day1, layer = "data")[genes_present, ]), error = function(e) numeric(0))
}
score_df <- tryCatch({df <- as.data.frame(do.call(rbind, score_list)); df$pathway <- rownames(df); df}, error = function(e) data.frame())
score_df_long <- tryCatch(tidyr::pivot_longer(score_df, cols = -pathway, names_to = "cell", values_to = "score"), error = function(e) data.frame())
p_lee_d1_1u <- tryCatch(ggplot(score_df_long, aes(x = pathway, y = score, fill = pathway)) + geom_boxplot() + coord_flip() + labs(title = "Hallmark Pathway Signatures - Lee Day 1", x = "Pathway", y = "Enrichment Score") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1u)

# -------- PLOT 1V: Source-target interaction alluvial plot --------
print("Generating: Plot 1V - Interaction Alluvial Plot (Lee Day 1)")
lee_day1_liana_consensus_nrow_1v <- nrow(lee_day1_liana_consensus)
alluvial_lee_d1_step1 <- slice_head(lee_day1_liana_consensus, n = 10)
alluvial_lee_d1_step2 <- tryCatch(dplyr::select(alluvial_lee_d1_step1, source, target, aggregate_rank), error = function(e) data.frame())
alluvial_lee_d1 <- tryCatch(mutate(alluvial_lee_d1_step2, lr_pair = paste0(source, " -> ", target)), error = function(e) data.frame())
p_lee_d1_1v <- tryCatch(ggplot(alluvial_lee_d1, aes(axis1 = source, axis2 = target, y = -aggregate_rank)) + geom_alluvium(aes(fill = source), alpha = 0.6) + geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + scale_x_discrete(limits = c("Source", "Target")) + labs(title = "Top 10 L-R Communication Flow - Lee Day 1", y = "Consensus support (-aggregate rank)") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1v)

# -------- PLOT 1W: Score distribution comparison --------
print("Generating: Plot 1W - Rank Score Distribution (Lee Day 1)")
p_lee_d1_1w <- tryCatch(ggplot(lee_day1_liana_result_df, aes(x = aggregate_rank, fill = aggregate_rank)) + geom_density(alpha = 0.7) + scale_fill_gradient(low = "lightblue", high = "darkblue") + labs(title = "L-R Pair Rank Score Distribution - Lee Day 1", x = "Aggregate Rank Score", y = "Density") + theme_minimal() + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d1_1w)

# Lee Day 1 PROGENy plots 1X-1Y: generated in Lee Day 1 PROGENy section above

print("✓ Lee Day 1: Data overview plots complete (1A-1E, 1P-1W). LIANA 1F-1O and PROGENy 1X-1Y in their sections above.")

# ============================================================================
# SECTION 6: STEP 2 - LEE DAY 3 DATA OVERVIEW & QC PLOTS
# ============================================================================
# LIANA plots (2E-2O2): in Lee Day 3 LIANA section above
# PROGENy plots (2P-2Q): in Lee Day 3 PROGENy section above
# Below: Arg1, UMAP, DEG, population, feature, marker violin

if (exists("lee_day3") && exists("lee_day3_pos") && exists("lee_day3_neg") && (!exists("lee_day3_arg1_status") || length(lee_day3_arg1_status) != ncol(lee_day3))) {
  lee_day3_arg1_status <- rep("Arg1neg", ncol(lee_day3))
  lee_day3_arg1_status[lee_day3_pos] <- "Arg1pos"
  lee_day3$arg1_status <- lee_day3_arg1_status
}
print(">>> STEP 2: LEE DAY 3 DATA OVERVIEW & QC PLOTS <<<")

# Note: lee_day3 and lee_day3_arg1_status are already defined above (before LIANA section)
# They are reused here for visualization

# -------- PLOT 2A: Arg1 Validation Violin --------
print("Generating: Plot 2A - Arg1 Validation Violin (Lee Day 3)")
p_lee_d3_2a <- VlnPlot(
  lee_day3, 
  features = "Arg1", 
  group.by = "arg1_status",
  pt.size = 0.1, 
  cols = c("Arg1neg" = "grey", "Arg1pos" = "red")
) + 
  labs(title = "Validation: Arg1 Expression - Lee Day 3", y = "Log Normalized Expression") +
  PLOT_TITLE_THEME

print(p_lee_d3_2a)

# -------- PLOT 2B: UMAP showing Arg1 distribution --------
print("Generating: Plot 2B - UMAP Arg1 Distribution (Lee Day 3)")
p_lee_d3_2b <- DimPlot(
  lee_day3, 
  group.by = "arg1_status", 
  cols = c("Arg1neg" = "lightblue", "Arg1pos" = "red"),
  pt.size = 2
) + 
  labs(title = "Arg1pos vs Arg1neg Distribution - Lee Day 3") +
  PLOT_TITLE_THEME

print(p_lee_d3_2b)

# -------- PLOT 2C: DEG Volcano (Arg1pos vs Arg1neg) --------
print("Generating: Plot 2C - DEG Volcano (Lee Day 3)")
lee_day3_neutrophils <- subset(lee_day3, cells = lee_day3_neut_cells)
Idents(lee_day3_neutrophils) <- lee_day3_neutrophils$arg1_status
lee_day3_deg <- FindMarkers(lee_day3_neutrophils, ident.1 = "Arg1pos", ident.2 = "Arg1neg", logfc.threshold = 0.25, min.pct = 0.1)
lee_day3_deg$gene <- rownames(lee_day3_deg)
# FindMarkers already returns p_val_adj - do NOT double-adjust
lee_day3_deg$significant <- (lee_day3_deg$p_val_adj < 0.05) & (abs(lee_day3_deg$avg_log2FC) > 0.5)

p_lee_d3_2c <- ggplot(lee_day3_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  geom_text_repel(
    data = head(lee_day3_deg[lee_day3_deg$significant, ], 10),
    aes(label = gene),
    size = 3
  ) +
  labs(
    title = "DEG Volcano Plot - Lee Day 3 (Arg1pos vs Arg1neg)",
    x = "Log2 Fold Change",
    y = "-log10 (FDR-adjusted p-value)"
  ) +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_lee_d3_2c)
write.csv(lee_day3_deg, file.path(OUTPUT_DIR, "02c_LeeDay3_DEG_Results.csv"))

# -------- PLOT 2D: Heatmap top DEGs --------
print("Generating: Plot 2D - DEG Heatmap (Lee Day 3)")
top_deg_lee_d3 <- head(lee_day3_deg[lee_day3_deg$significant, ], 20)
expr_matrix_lee_d3 <- Seurat::GetAssayData(lee_day3_neutrophils, layer = "data")[top_deg_lee_d3$gene, , drop = FALSE]
# Scale genes (rows) before plotting.
expr_matrix_lee_d3_scaled <- t(scale(t(expr_matrix_lee_d3)))

p_lee_d3_2d <- pheatmap::pheatmap(
  expr_matrix_lee_d3_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 20 DEGs - Lee Day 3 (Arg1pos vs Arg1neg)",
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 10
)

print(p_lee_d3_2d)

# Lee Day 3 LIANA plots 2E-2O2: generated in Lee Day 3 LIANA section above

# -------- PLOT 2J: Lee Day 3 Data Overview - Arg1pos population size --------
print("Generating: Plot 2J - Arg1pos Population Size (Lee Day 3)")
pop_size_lee_d3 <- data.frame(
  Status = c("Arg1pos", "Arg1neg"),
  Count = c(sum(lee_day3_arg1_status == "Arg1pos"), sum(lee_day3_arg1_status == "Arg1neg")),
  Percentage = c(
    sum(lee_day3_arg1_status == "Arg1pos") / length(lee_day3_arg1_status) * 100,
    sum(lee_day3_arg1_status == "Arg1neg") / length(lee_day3_arg1_status) * 100
  )
)

p_lee_d3_2j <- ggplot(pop_size_lee_d3, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), vjust = -0.5) +
  scale_fill_manual(values = c("Arg1pos" = "red", "Arg1neg" = "lightblue")) +
  labs(title = "Arg1pos Population - Lee Day 3", x = "Neutrophil Status", y = "Cell Count") +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_lee_d3_2j)

# -------- PLOT 2K: Feature plots --------
print("Generating: Plot 2K - Feature Plots (Lee Day 3)")
feature_genes_lee_d3 <- c("Arg1", "Il10", "Tnf")
feature_genes_present <- feature_genes_lee_d3[feature_genes_lee_d3 %in% rownames(lee_day3)]
p_lee_d3_2k <- tryCatch(FeaturePlot(lee_day3, features = feature_genes_present[1:min(3, length(feature_genes_present))], ncol = min(3, length(feature_genes_present)), pt.size = 2), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d3_2k)

# -------- PLOT 2N: Lee Day 3 Data Overview - Violin plots for marker genes --------
print("Generating: Plot 2N - Marker Gene Violin Plot (Lee Day 3)")
marker_genes_lee_d3 <- c("Cd4", "Cd8a", "Iba1", "Gfap", "Pdgfra")
marker_genes_present <- marker_genes_lee_d3[marker_genes_lee_d3 %in% rownames(lee_day3)]
p_lee_d3_2n <- tryCatch(VlnPlot(lee_day3, features = marker_genes_present[1:min(3, length(marker_genes_present))], ncol = 1) + labs(title = "Marker Gene Expression - Lee Day 3") + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_lee_d3_2n)

# Lee Day 3 PROGENy plots 2P-2Q: generated in Lee Day 3 PROGENy section above

print("✓ Lee Day 3: Data overview plots complete. LIANA and PROGENy in their sections above.")

# ============================================================================
# SECTION 7: STEP 3 - WANG DAY 3 DATA OVERVIEW & QC PLOTS
# ============================================================================
# LIANA plots (3E-3G3, 3O2): in Wang Day 3 LIANA section above
# PROGENy plots (3M-3N): in Wang Day 3 PROGENy section above
# Below: Arg1, UMAP, DEG, population, feature, marker violin

if (exists("wang_day3") && exists("wang_day3_pos") && exists("wang_day3_neg") && (!exists("wang_day3_arg1_status") || length(wang_day3_arg1_status) != ncol(wang_day3))) {
  wang_day3_arg1_status <- rep("Arg1neg", ncol(wang_day3))
  wang_day3_arg1_status[wang_day3_pos] <- "Arg1pos"
  wang_day3$arg1_status <- wang_day3_arg1_status
}
print(">>> STEP 3: WANG DAY 3 DATA OVERVIEW & QC PLOTS <<<")

# Note: wang_day3 and wang_day3_arg1_status are already defined above (before LIANA section)
# They are reused here for visualization

# -------- PLOT 3A: Arg1 Validation Violin --------
print("Generating: Plot 3A - Arg1 Validation Violin (Wang Day 3)")
p_wang_d3_3a <- VlnPlot(
  wang_day3, 
  features = "Arg1", 
  group.by = "arg1_status",
  pt.size = 0.1, 
  cols = c("Arg1neg" = "grey", "Arg1pos" = "red")
) + 
  labs(title = "Validation: Arg1 Expression - Wang Day 3", y = "Log Normalized Expression") +
  PLOT_TITLE_THEME

print(p_wang_d3_3a)

# -------- PLOT 3B: UMAP showing Arg1 distribution --------
print("Generating: Plot 3B - UMAP Arg1 Distribution (Wang Day 3)")
p_wang_d3_3b <- DimPlot(
  wang_day3, 
  group.by = "arg1_status", 
  cols = c("Arg1neg" = "lightblue", "Arg1pos" = "red"),
  pt.size = 2
) + 
  labs(title = "Arg1pos vs Arg1neg Distribution - Wang Day 3") +
  PLOT_TITLE_THEME

print(p_wang_d3_3b)

# -------- PLOT 3C: DEG Volcano (Arg1pos vs Arg1neg) --------
print("Generating: Plot 3C - DEG Volcano (Wang Day 3)")
wang_day3_neutrophils <- subset(wang_day3, cells = wang_day3_neut_cells)
Idents(wang_day3_neutrophils) <- wang_day3_neutrophils$arg1_status
wang_day3_deg <- FindMarkers(wang_day3_neutrophils, ident.1 = "Arg1pos", ident.2 = "Arg1neg", logfc.threshold = 0.25, min.pct = 0.1)
wang_day3_deg$gene <- rownames(wang_day3_deg)
# FindMarkers already returns p_val_adj - do NOT double-adjust
wang_day3_deg$significant <- (wang_day3_deg$p_val_adj < 0.05) & (abs(wang_day3_deg$avg_log2FC) > 0.5)

p_wang_d3_3c <- ggplot(wang_day3_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  geom_text_repel(
    data = head(wang_day3_deg[wang_day3_deg$significant, ], 10),
    aes(label = gene),
    size = 3
  ) +
  labs(
    title = "DEG Volcano Plot - Wang Day 3 (Arg1pos vs Arg1neg)",
    x = "Log2 Fold Change",
    y = "-log10 (FDR-adjusted p-value)"
  ) +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_wang_d3_3c)
write.csv(wang_day3_deg, file.path(OUTPUT_DIR, "03c_WangDay3_DEG_Results.csv"))

# -------- PLOT 3D: Heatmap top DEGs --------
print("Generating: Plot 3D - DEG Heatmap (Wang Day 3)")
top_deg_wang_d3 <- head(wang_day3_deg[wang_day3_deg$significant, ], 20)
expr_matrix_wang_d3 <- Seurat::GetAssayData(wang_day3_neutrophils, layer = "data")[top_deg_wang_d3$gene, , drop = FALSE]
# Scale genes (rows) before plotting.
expr_matrix_wang_d3_scaled <- t(scale(t(expr_matrix_wang_d3)))

p_wang_d3_3d <- pheatmap::pheatmap(
  expr_matrix_wang_d3_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 20 DEGs - Wang Day 3 (Arg1pos vs Arg1neg)",
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 10
)

print(p_wang_d3_3d)

# Wang Day 3 LIANA plots 3E-3G3: generated in Wang Day 3 LIANA section above

# -------- PLOT 3H: Wang Day 3 Data Overview - Population size --------
print("Generating: Plot 3H - Arg1pos Population Size (Wang Day 3)")
pop_size_wang_d3 <- data.frame(
  Status = c("Arg1pos", "Arg1neg"),
  Count = c(sum(wang_day3_arg1_status == "Arg1pos"), sum(wang_day3_arg1_status == "Arg1neg")),
  Percentage = c(
    sum(wang_day3_arg1_status == "Arg1pos") / length(wang_day3_arg1_status) * 100,
    sum(wang_day3_arg1_status == "Arg1neg") / length(wang_day3_arg1_status) * 100
  )
)

p_wang_d3_3h <- ggplot(pop_size_wang_d3, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), vjust = -0.5) +
  scale_fill_manual(values = c("Arg1pos" = "red", "Arg1neg" = "lightblue")) +
  labs(title = "Arg1pos Population - Wang Day 3", x = "Neutrophil Status", y = "Cell Count") +
  theme_minimal() +
  PLOT_TITLE_THEME

print(p_wang_d3_3h)

# -------- PLOT 3I: Feature plots --------
print("Generating: Plot 3I - Feature Plots (Wang Day 3)")
feature_genes_wang_d3 <- c("Arg1", "Il10", "Tnf")
feature_genes_present <- feature_genes_wang_d3[feature_genes_wang_d3 %in% rownames(wang_day3)]
p_wang_d3_3i <- tryCatch(FeaturePlot(wang_day3, features = feature_genes_present[1:min(3, length(feature_genes_present))], ncol = min(3, length(feature_genes_present)), pt.size = 2), error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_wang_d3_3i)

# -------- PLOT 3L: Wang Day 3 Data Overview - Marker gene violin plots --------
print("Generating: Plot 3L - Marker Gene Violin Plot (Wang Day 3)")
marker_genes_wang_d3 <- c("Cd4", "Cd8a", "Iba1", "Gfap", "Pdgfra")
marker_genes_present <- marker_genes_wang_d3[marker_genes_wang_d3 %in% rownames(wang_day3)]
p_wang_d3_3l <- tryCatch(VlnPlot(wang_day3, features = marker_genes_present[1:min(3, length(marker_genes_present))], ncol = 1) + labs(title = "Marker Gene Expression - Wang Day 3") + PLOT_TITLE_THEME, error = function(e) ggplot() + theme_void() + labs(title = "No data available"))
print(p_wang_d3_3l)

# Wang Day 3 PROGENy plots 3M-3N: generated in Wang Day 3 PROGENy section above

print("✓ Wang Day 3: Data overview plots complete. LIANA and PROGENy in their sections above.")
# Unified-scale LIANA/PROGENy ggplots: printed inline where each p_* is created (recompute UNIFY_* from all existing cohort tables before each print).


# ============================================================================
# Cross-dataset Arg1+ received L-R consensus (positional rank in each cohort top-N; exploratory)
# ============================================================================
CROSS_LR_TOP_N <- 20L
CROSS_LR_STRONG_MAX <- 10L
CROSS_LR_MIN_COHORTS <- 2L
cross_lr_parts <- list()
if (exists("lee_day1_arg1pos_top_signals") && is.data.frame(lee_day1_arg1pos_top_signals) && nrow(lee_day1_arg1pos_top_signals) > 0) {
  d_ld1 <- lee_day1_arg1pos_top_signals[seq_len(min(CROSS_LR_TOP_N, nrow(lee_day1_arg1pos_top_signals))), , drop = FALSE]
  cross_lr_parts[[length(cross_lr_parts) + 1L]] <- data.frame(lr = paste0(d_ld1$ligand.complex, "->", d_ld1$receptor.complex), rank = seq_len(nrow(d_ld1)), dataset = "Lee_D1", stringsAsFactors = FALSE)
}
if (exists("lee_day3_arg1pos_top_signals") && is.data.frame(lee_day3_arg1pos_top_signals) && nrow(lee_day3_arg1pos_top_signals) > 0) {
  d_ld3 <- lee_day3_arg1pos_top_signals[seq_len(min(CROSS_LR_TOP_N, nrow(lee_day3_arg1pos_top_signals))), , drop = FALSE]
  cross_lr_parts[[length(cross_lr_parts) + 1L]] <- data.frame(lr = paste0(d_ld3$ligand.complex, "->", d_ld3$receptor.complex), rank = seq_len(nrow(d_ld3)), dataset = "Lee_D3", stringsAsFactors = FALSE)
}
if (exists("wang_day3_arg1pos_top_signals") && is.data.frame(wang_day3_arg1pos_top_signals) && nrow(wang_day3_arg1pos_top_signals) > 0) {
  d_wd3 <- wang_day3_arg1pos_top_signals[seq_len(min(CROSS_LR_TOP_N, nrow(wang_day3_arg1pos_top_signals))), , drop = FALSE]
  cross_lr_parts[[length(cross_lr_parts) + 1L]] <- data.frame(lr = paste0(d_wd3$ligand.complex, "->", d_wd3$receptor.complex), rank = seq_len(nrow(d_wd3)), dataset = "Wang_D3", stringsAsFactors = FALSE)
}
if (exists("qin_day1_arg1pos_top_signals") && is.data.frame(qin_day1_arg1pos_top_signals) && nrow(qin_day1_arg1pos_top_signals) > 0) {
  d_qd1 <- qin_day1_arg1pos_top_signals[seq_len(min(CROSS_LR_TOP_N, nrow(qin_day1_arg1pos_top_signals))), , drop = FALSE]
  cross_lr_parts[[length(cross_lr_parts) + 1L]] <- data.frame(lr = paste0(d_qd1$ligand.complex, "->", d_qd1$receptor.complex), rank = seq_len(nrow(d_qd1)), dataset = "Qin_D1", stringsAsFactors = FALSE)
}
if (exists("qin_day3_arg1pos_top_signals") && is.data.frame(qin_day3_arg1pos_top_signals) && nrow(qin_day3_arg1pos_top_signals) > 0) {
  d_qd3 <- qin_day3_arg1pos_top_signals[seq_len(min(CROSS_LR_TOP_N, nrow(qin_day3_arg1pos_top_signals))), , drop = FALSE]
  cross_lr_parts[[length(cross_lr_parts) + 1L]] <- data.frame(lr = paste0(d_qd3$ligand.complex, "->", d_qd3$receptor.complex), rank = seq_len(nrow(d_qd3)), dataset = "Qin_D3", stringsAsFactors = FALSE)
}
if (length(cross_lr_parts) >= CROSS_LR_MIN_COHORTS) {
  all_ds_lr <- do.call(rbind, cross_lr_parts)
  all_ds_lr <- stats::aggregate(rank ~ lr + dataset, all_ds_lr, min)
  cross_lr_wide <- tidyr::pivot_wider(all_ds_lr, names_from = dataset, values_from = rank, values_fill = CROSS_LR_TOP_N + 1L)
  cross_lr_mat <- as.matrix(cross_lr_wide[, -1L, drop = FALSE])
  rownames(cross_lr_mat) <- cross_lr_wide$lr
  keep_lr <- apply(cross_lr_mat, 1L, function(x) sum(x <= CROSS_LR_STRONG_MAX, na.rm = TRUE)) >= CROSS_LR_MIN_COHORTS
  cross_lr_mat_f <- cross_lr_mat[keep_lr, , drop = FALSE]
  if (nrow(cross_lr_mat_f) > 0L && ncol(cross_lr_mat_f) > 0L) {
    p_cross_lr <- pheatmap::pheatmap(cross_lr_mat_f, color = grDevices::colorRampPalette(c("red", "white", "blue"))(100L), main = "Cross-dataset Arg1+ neutrophil: L-R consensus\n(red = stronger rank in that cohort's top list)", cluster_cols = FALSE)
    print(p_cross_lr)
    write.csv(data.frame(lr = rownames(cross_lr_mat_f), cross_lr_mat_f, check.names = FALSE), file.path(OUTPUT_DIR, "CrossDataset_Arg1pos_TopLR_Consensus.csv"), row.names = FALSE)
  } else {
    message("Cross-dataset L-R heatmap: no L-R pairs ranked in top ", CROSS_LR_STRONG_MAX, " in at least ", CROSS_LR_MIN_COHORTS, " cohorts.")
  }
} else {
  message("Cross-dataset L-R heatmap skipped (need top-signal tables from at least ", CROSS_LR_MIN_COHORTS, " cohorts).")
}

# ============================================================================
# POST-HOC: Cross-cohort unified ggplot replot (all five cohorts in memory)
# Inline prints use progressive pooling (limits widen as later cohorts run). This
# block recomputes every UNIFY_* once from all cohort tables and reprints the
# same p_* objects so Lee D1 and Qin D3 share identical scale limits.
# ============================================================================
print("")
print(">>> POST-HOC cross-cohort unified scales — replot LIANA + PROGENy comparison ggplots <<<")
if (!exists("LIANA_TOP_SIGNAL_POINT_SIZE_RANGE")) LIANA_TOP_SIGNAL_POINT_SIZE_RANGE <- c(2, 8)
vals_nlr <- numeric(0)
if (exists("liana_network_lee_d1_both") && is.data.frame(liana_network_lee_d1_both) && nrow(liana_network_lee_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d1_both$aggregate_rank + 1e-10))
if (exists("liana_network_lee_d3_both") && is.data.frame(liana_network_lee_d3_both) && nrow(liana_network_lee_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_lee_d3_both$aggregate_rank + 1e-10))
if (exists("liana_network_wang_d3_both") && is.data.frame(liana_network_wang_d3_both) && nrow(liana_network_wang_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_wang_d3_both$aggregate_rank + 1e-10))
if (exists("liana_network_qin_d1_both") && is.data.frame(liana_network_qin_d1_both) && nrow(liana_network_qin_d1_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d1_both$aggregate_rank + 1e-10))
if (exists("liana_network_qin_d3_both") && is.data.frame(liana_network_qin_d3_both) && nrow(liana_network_qin_d3_both) > 0) vals_nlr <- c(vals_nlr, -log10(liana_network_qin_d3_both$aggregate_rank + 1e-10))
vals_nlr <- vals_nlr[is.finite(vals_nlr)]
UNIFY_LIANA_NEGLOG10 <- if (length(vals_nlr) > 0) range(vals_nlr) else c(0, 10)
if (length(UNIFY_LIANA_NEGLOG10) != 2 || !all(is.finite(UNIFY_LIANA_NEGLOG10))) UNIFY_LIANA_NEGLOG10 <- c(0, 10)
if (UNIFY_LIANA_NEGLOG10[1] == UNIFY_LIANA_NEGLOG10[2]) UNIFY_LIANA_NEGLOG10[2] <- UNIFY_LIANA_NEGLOG10[1] + 1e-6
vals_ext <- numeric(0)
if (exists("lee_day1_external_to_neutrophils") && is.data.frame(lee_day1_external_to_neutrophils) && nrow(lee_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day1_external_to_neutrophils$aggregate_rank)
if (exists("lee_day3_external_to_neutrophils") && is.data.frame(lee_day3_external_to_neutrophils) && nrow(lee_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(lee_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, lee_day3_external_to_neutrophils$aggregate_rank)
if (exists("wang_day3_external_to_neutrophils") && is.data.frame(wang_day3_external_to_neutrophils) && nrow(wang_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(wang_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, wang_day3_external_to_neutrophils$aggregate_rank)
if (exists("qin_day1_external_to_neutrophils") && is.data.frame(qin_day1_external_to_neutrophils) && nrow(qin_day1_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day1_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day1_external_to_neutrophils$aggregate_rank)
if (exists("qin_day3_external_to_neutrophils") && is.data.frame(qin_day3_external_to_neutrophils) && nrow(qin_day3_external_to_neutrophils) > 0 && "aggregate_rank" %in% names(qin_day3_external_to_neutrophils)) vals_ext <- c(vals_ext, qin_day3_external_to_neutrophils$aggregate_rank)
if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0 && "aggregate_rank" %in% names(liana_top_lee_d1)) vals_ext <- c(vals_ext, liana_top_lee_d1$aggregate_rank)
if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_lee_d3_external)) vals_ext <- c(vals_ext, liana_top_lee_d3_external$aggregate_rank)
if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_wang_d3_external)) vals_ext <- c(vals_ext, liana_top_wang_d3_external$aggregate_rank)
if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d1_external)) vals_ext <- c(vals_ext, liana_top_qin_d1_external$aggregate_rank)
if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0 && "aggregate_rank" %in% names(liana_top_qin_d3_external)) vals_ext <- c(vals_ext, liana_top_qin_d3_external$aggregate_rank)
vals_ext <- suppressWarnings(as.numeric(vals_ext))
vals_ext <- vals_ext[is.finite(vals_ext)]
plot_ext_ar <- numeric(0)
if (exists("liana_top_lee_d1") && is.data.frame(liana_top_lee_d1) && nrow(liana_top_lee_d1) > 0) plot_ext_ar <- c(plot_ext_ar, suppressWarnings(as.numeric(liana_top_lee_d1[["aggregate_rank"]])))
if (exists("liana_top_lee_d3_external") && is.data.frame(liana_top_lee_d3_external) && nrow(liana_top_lee_d3_external) > 0) plot_ext_ar <- c(plot_ext_ar, suppressWarnings(as.numeric(liana_top_lee_d3_external[["aggregate_rank"]])))
if (exists("liana_top_wang_d3_external") && is.data.frame(liana_top_wang_d3_external) && nrow(liana_top_wang_d3_external) > 0) plot_ext_ar <- c(plot_ext_ar, suppressWarnings(as.numeric(liana_top_wang_d3_external[["aggregate_rank"]])))
if (exists("liana_top_qin_d1_external") && is.data.frame(liana_top_qin_d1_external) && nrow(liana_top_qin_d1_external) > 0) plot_ext_ar <- c(plot_ext_ar, suppressWarnings(as.numeric(liana_top_qin_d1_external[["aggregate_rank"]])))
if (exists("liana_top_qin_d3_external") && is.data.frame(liana_top_qin_d3_external) && nrow(liana_top_qin_d3_external) > 0) plot_ext_ar <- c(plot_ext_ar, suppressWarnings(as.numeric(liana_top_qin_d3_external[["aggregate_rank"]])))
plot_ext_ar <- plot_ext_ar[is.finite(plot_ext_ar)]
UNIFY_EXT_AR <- range(c(vals_ext, plot_ext_ar), na.rm = TRUE)
if (length(plot_ext_ar) == 0 && length(vals_ext) == 0) UNIFY_EXT_AR <- c(0, 1)
if (!all(is.finite(UNIFY_EXT_AR))) UNIFY_EXT_AR <- c(0, 1)
if (UNIFY_EXT_AR[1] == UNIFY_EXT_AR[2]) UNIFY_EXT_AR[2] <- UNIFY_EXT_AR[1] + max(abs(UNIFY_EXT_AR[1]) * 1e-6, 1e-12)
vals_cons <- numeric(0)
if (exists("liana_top_lee_d3") && is.data.frame(liana_top_lee_d3) && nrow(liana_top_lee_d3) > 0) vals_cons <- c(vals_cons, liana_top_lee_d3$aggregate_rank)
if (exists("liana_top_wang_d3") && is.data.frame(liana_top_wang_d3) && nrow(liana_top_wang_d3) > 0) vals_cons <- c(vals_cons, liana_top_wang_d3$aggregate_rank)
if (exists("liana_top_qin_d1") && is.data.frame(liana_top_qin_d1) && nrow(liana_top_qin_d1) > 0) vals_cons <- c(vals_cons, liana_top_qin_d1$aggregate_rank)
if (exists("liana_top_qin_d3") && is.data.frame(liana_top_qin_d3) && nrow(liana_top_qin_d3) > 0) vals_cons <- c(vals_cons, liana_top_qin_d3$aggregate_rank)
vals_cons <- vals_cons[is.finite(vals_cons)]
UNIFY_CONS_AR <- if (length(vals_cons) > 0) range(vals_cons) else c(0, 1)
if (UNIFY_CONS_AR[1] == UNIFY_CONS_AR[2]) UNIFY_CONS_AR[2] <- UNIFY_CONS_AR[1] + 1e-6
vals_fc <- numeric(0)
vals_fmr <- numeric(0)
freq_tabnames <- c("receptor_freq_lee_d1", "receptor_freq_lee_d3", "receptor_freq_wang_d3", "receptor_freq_qin_d1", "receptor_freq_qin_d3", "ligand_freq_lee_d1", "ligand_freq_lee_d3", "ligand_freq_wang_d3", "ligand_freq_qin_d1", "ligand_freq_qin_d3")
for (fn in freq_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("count" %in% names(d)) vals_fc <- c(vals_fc, d$count)
  if ("mean_rank" %in% names(d)) vals_fmr <- c(vals_fmr, d$mean_rank)
}
vals_fc <- vals_fc[is.finite(vals_fc)]
vals_fmr <- vals_fmr[is.finite(vals_fmr)]
UNIFY_LIANA_FREQ_COUNT <- if (length(vals_fc) > 0) range(vals_fc) else c(0, 1)
UNIFY_LIANA_FREQ_MR <- if (length(vals_fmr) > 0) range(vals_fmr) else c(0, 1)
if (UNIFY_LIANA_FREQ_COUNT[1] == UNIFY_LIANA_FREQ_COUNT[2]) UNIFY_LIANA_FREQ_COUNT[2] <- UNIFY_LIANA_FREQ_COUNT[1] + 1e-6
if (UNIFY_LIANA_FREQ_MR[1] == UNIFY_LIANA_FREQ_MR[2]) UNIFY_LIANA_FREQ_MR[2] <- UNIFY_LIANA_FREQ_MR[1] + 1e-6
vals_iy <- numeric(0)
vals_imr <- numeric(0)
imp_tabnames <- c("source_importance_lee_d1", "source_importance_lee_d3", "source_importance_wang_d3", "source_importance_qin_d1", "source_importance_qin_d3")
for (fn in imp_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("importance_score" %in% names(d)) vals_iy <- c(vals_iy, d$importance_score)
  if ("mean_rank" %in% names(d)) vals_imr <- c(vals_imr, d$mean_rank)
}
vals_iy <- vals_iy[is.finite(vals_iy)]
vals_imr <- vals_imr[is.finite(vals_imr)]
UNIFY_SRC_IMP_Y <- if (length(vals_iy) > 0) range(vals_iy) else c(0, 1)
UNIFY_SRC_IMP_MR <- if (length(vals_imr) > 0) range(vals_imr) else c(0, 1)
if (UNIFY_SRC_IMP_Y[1] == UNIFY_SRC_IMP_Y[2]) UNIFY_SRC_IMP_Y[2] <- UNIFY_SRC_IMP_Y[1] + 1e-6
if (UNIFY_SRC_IMP_MR[1] == UNIFY_SRC_IMP_MR[2]) UNIFY_SRC_IMP_MR[2] <- UNIFY_SRC_IMP_MR[1] + 1e-6
vals_sy <- numeric(0)
vals_smr <- numeric(0)
spec_tabnames <- c("target_specificity_lee_d1", "target_specificity_lee_d3", "target_specificity_wang_d3", "target_specificity_qin_d1", "target_specificity_qin_d3")
for (fn in spec_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("specificity_score" %in% names(d)) vals_sy <- c(vals_sy, d$specificity_score)
  if ("mean_rank" %in% names(d)) vals_smr <- c(vals_smr, d$mean_rank)
}
vals_sy <- vals_sy[is.finite(vals_sy)]
vals_smr <- vals_smr[is.finite(vals_smr)]
UNIFY_TGT_SPEC_Y <- if (length(vals_sy) > 0) range(vals_sy) else c(0, 1)
UNIFY_TGT_SPEC_MR <- if (length(vals_smr) > 0) range(vals_smr) else c(0, 1)
if (UNIFY_TGT_SPEC_Y[1] == UNIFY_TGT_SPEC_Y[2]) UNIFY_TGT_SPEC_Y[2] <- UNIFY_TGT_SPEC_Y[1] + 1e-6
if (UNIFY_TGT_SPEC_MR[1] == UNIFY_TGT_SPEC_MR[2]) UNIFY_TGT_SPEC_MR[2] <- UNIFY_TGT_SPEC_MR[1] + 1e-6
vals_a1 <- numeric(0)
if (exists("lee_day1_liana_arg1pos_topranked") && is.data.frame(lee_day1_liana_arg1pos_topranked) && nrow(lee_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("lee_day3_liana_arg1pos_topranked") && is.data.frame(lee_day3_liana_arg1pos_topranked) && nrow(lee_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(lee_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("wang_day3_liana_arg1pos_topranked") && is.data.frame(wang_day3_liana_arg1pos_topranked) && nrow(wang_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(wang_day3_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day1_liana_arg1pos_topranked") && is.data.frame(qin_day1_liana_arg1pos_topranked) && nrow(qin_day1_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day1_liana_arg1pos_topranked$aggregate_rank, 15))
if (exists("qin_day3_liana_arg1pos_topranked") && is.data.frame(qin_day3_liana_arg1pos_topranked) && nrow(qin_day3_liana_arg1pos_topranked) > 0) vals_a1 <- c(vals_a1, head(qin_day3_liana_arg1pos_topranked$aggregate_rank, 15))
vals_a1 <- vals_a1[is.finite(vals_a1)]
UNIFY_ARG1_AR <- if (length(vals_a1) > 0) range(vals_a1) else c(0, 1)
if (UNIFY_ARG1_AR[1] == UNIFY_ARG1_AR[2]) UNIFY_ARG1_AR[2] <- UNIFY_ARG1_AR[1] + 1e-6
UNIFY_ARG1_NEGY <- c(-UNIFY_ARG1_AR[2], -UNIFY_ARG1_AR[1])
vals_hist <- numeric(0)
if (exists("lee_day1_liana_consensus") && is.data.frame(lee_day1_liana_consensus) && nrow(lee_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day1_liana_consensus)) vals_hist <- c(vals_hist, lee_day1_liana_consensus$aggregate_rank)
if (exists("lee_day3_liana_consensus") && is.data.frame(lee_day3_liana_consensus) && nrow(lee_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(lee_day3_liana_consensus)) vals_hist <- c(vals_hist, lee_day3_liana_consensus$aggregate_rank)
if (exists("wang_day3_liana_consensus") && is.data.frame(wang_day3_liana_consensus) && nrow(wang_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(wang_day3_liana_consensus)) vals_hist <- c(vals_hist, wang_day3_liana_consensus$aggregate_rank)
if (exists("qin_day1_liana_consensus") && is.data.frame(qin_day1_liana_consensus) && nrow(qin_day1_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day1_liana_consensus)) vals_hist <- c(vals_hist, qin_day1_liana_consensus$aggregate_rank)
if (exists("qin_day3_liana_consensus") && is.data.frame(qin_day3_liana_consensus) && nrow(qin_day3_liana_consensus) > 0 && "aggregate_rank" %in% names(qin_day3_liana_consensus)) vals_hist <- c(vals_hist, qin_day3_liana_consensus$aggregate_rank)
vals_hist <- vals_hist[is.finite(vals_hist)]
UNIFY_HIST_AR <- if (length(vals_hist) > 0) range(vals_hist) else c(0, 1)
if (UNIFY_HIST_AR[1] == UNIFY_HIST_AR[2]) UNIFY_HIST_AR[2] <- UNIFY_HIST_AR[1] + 1e-6
vals_pw <- numeric(0)
if (exists("lee_day1_pathway_long") && is.data.frame(lee_day1_pathway_long) && nrow(lee_day1_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day1_pathway_long)) vals_pw <- c(vals_pw, lee_day1_pathway_long$Pathway_Score)
if (exists("lee_day3_pathway_long") && is.data.frame(lee_day3_pathway_long) && nrow(lee_day3_pathway_long) > 0 && "Pathway_Score" %in% names(lee_day3_pathway_long)) vals_pw <- c(vals_pw, lee_day3_pathway_long$Pathway_Score)
if (exists("wang_day3_pathway_long") && is.data.frame(wang_day3_pathway_long) && nrow(wang_day3_pathway_long) > 0 && "Pathway_Score" %in% names(wang_day3_pathway_long)) vals_pw <- c(vals_pw, wang_day3_pathway_long$Pathway_Score)
if (exists("qin_day1_pathway_long") && is.data.frame(qin_day1_pathway_long) && nrow(qin_day1_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day1_pathway_long)) vals_pw <- c(vals_pw, qin_day1_pathway_long$Pathway_Score)
if (exists("qin_day3_pathway_long") && is.data.frame(qin_day3_pathway_long) && nrow(qin_day3_pathway_long) > 0 && "Pathway_Score" %in% names(qin_day3_pathway_long)) vals_pw <- c(vals_pw, qin_day3_pathway_long$Pathway_Score)
vals_pw <- vals_pw[is.finite(vals_pw)]
UNIFY_PW_Y <- if (length(vals_pw) > 0) range(vals_pw) else c(-1, 1)
if (UNIFY_PW_Y[1] == UNIFY_PW_Y[2]) UNIFY_PW_Y[2] <- UNIFY_PW_Y[1] + 1e-6
vals_w <- numeric(0)
map_tabnames <- c("lee_day1_ligand_pathway_map", "lee_day3_ligand_pathway_map", "wang_day3_ligand_pathway_map", "qin_day1_ligand_pathway_map", "qin_day3_ligand_pathway_map")
for (fn in map_tabnames) {
  if (!exists(fn)) next
  d <- get(fn)
  if (!is.data.frame(d) || nrow(d) == 0) next
  if ("weight" %in% names(d)) vals_w <- c(vals_w, d$weight)
}
vals_w <- vals_w[is.finite(vals_w)]
UNIFY_LIGW_ABS <- if (length(vals_w) > 0) range(abs(vals_w)) else c(0, 1)
if (UNIFY_LIGW_ABS[1] == UNIFY_LIGW_ABS[2]) UNIFY_LIGW_ABS[2] <- UNIFY_LIGW_ABS[1] + 1e-6
mxw <- if (length(vals_w) > 0) max(abs(vals_w)) else 1
if (!is.finite(mxw) || mxw <= 0) mxw <- 1
UNIFY_LIGW_COL <- c(-mxw, mxw)

print("--- Lee Day 1 (final unified scales) ---")
if (exists("p_01f")) print(p_01f + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = LIANA_TOP_SIGNAL_POINT_SIZE_RANGE))
if (exists("p_01g")) print(p_01g + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
if (exists("p_01h")) print(p_01h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_01i")) print(p_01i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_01k")) print(p_01k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
if (exists("p_01n")) print(p_01n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
if (exists("p_01l")) print(p_01l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
if (exists("p_01m")) print(p_01m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
if (exists("p_lee_d1_progeny1")) print(p_lee_d1_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
if (exists("p_lee_d1_progeny2")) print(p_lee_d1_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("--- Lee Day 3 (final unified scales) ---")
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02e")) print(p_02e + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = LIANA_TOP_SIGNAL_POINT_SIZE_RANGE))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02f0")) print(p_02f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02f")) print(p_02f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02h")) print(p_02h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02i")) print(p_02i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02k")) print(p_02k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02n")) print(p_02n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02l")) print(p_02l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
if (exists("lee_d3_liana_viz_ok") && lee_d3_liana_viz_ok && exists("p_02m")) print(p_02m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
if (exists("p_lee_d3_progeny1")) print(p_lee_d3_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
if (exists("p_lee_d3_progeny2")) print(p_lee_d3_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("--- Wang Day 3 (final unified scales) ---")
if (exists("p_03e")) print(p_03e + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = LIANA_TOP_SIGNAL_POINT_SIZE_RANGE))
if (exists("p_03f0")) print(p_03f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
if (exists("p_03f")) print(p_03f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
if (exists("p_03h")) print(p_03h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_03i")) print(p_03i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_03k")) print(p_03k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
if (exists("p_03n")) print(p_03n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
if (exists("p_03l")) print(p_03l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
if (exists("p_03m_rank")) print(p_03m_rank + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
if (exists("p_wang_d3_progeny1")) print(p_wang_d3_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
if (exists("p_wang_d3_progeny2")) print(p_wang_d3_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("--- Qin Day 1 (final unified scales) ---")
if (exists("p_qin_d1_network")) print(p_qin_d1_network + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = LIANA_TOP_SIGNAL_POINT_SIZE_RANGE))
if (exists("p_qin_d1_f0")) print(p_qin_d1_f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
if (exists("p_qin_d1_f")) print(p_qin_d1_f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
if (exists("p_qin_d1_h")) print(p_qin_d1_h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_qin_d1_i")) print(p_qin_d1_i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_qin_d1_k")) print(p_qin_d1_k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
if (exists("p_qin_d1_n")) print(p_qin_d1_n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
if (exists("p_qin_d1_l")) print(p_qin_d1_l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
if (exists("p_qin_d1_m")) print(p_qin_d1_m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
if (exists("p_qin_d1_progeny1")) print(p_qin_d1_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
if (exists("p_qin_d1_progeny2")) print(p_qin_d1_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))

print("--- Qin Day 3 (final unified scales) ---")
if (exists("p_qin_d3_network")) print(p_qin_d3_network + ggplot2::scale_size_continuous(limits = UNIFY_LIANA_NEGLOG10, range = LIANA_TOP_SIGNAL_POINT_SIZE_RANGE))
if (exists("p_qin_d3_f0")) print(p_qin_d3_f0 + ggplot2::scale_y_continuous(limits = UNIFY_EXT_AR, oob = scales::squish) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_EXT_AR, oob = scales::squish))
if (exists("p_qin_d3_f")) print(p_qin_d3_f + ggplot2::scale_y_continuous(limits = UNIFY_CONS_AR) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue", limits = UNIFY_CONS_AR))
if (exists("p_qin_d3_h")) print(p_qin_d3_h + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_qin_d3_i")) print(p_qin_d3_i + ggplot2::scale_y_continuous(limits = UNIFY_LIANA_FREQ_COUNT) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_LIANA_FREQ_MR))
if (exists("p_qin_d3_k")) print(p_qin_d3_k + ggplot2::scale_y_continuous(limits = UNIFY_SRC_IMP_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred", limits = UNIFY_SRC_IMP_MR))
if (exists("p_qin_d3_n")) print(p_qin_d3_n + ggplot2::scale_y_continuous(limits = UNIFY_ARG1_NEGY) + ggplot2::scale_fill_gradient(low = "red", high = "green", limits = UNIFY_ARG1_AR))
if (exists("p_qin_d3_l")) print(p_qin_d3_l + ggplot2::scale_y_continuous(limits = UNIFY_TGT_SPEC_Y) + ggplot2::scale_fill_gradient(low = "lightblue", high = "darkgreen", limits = UNIFY_TGT_SPEC_MR))
if (exists("p_qin_d3_m")) print(p_qin_d3_m + ggplot2::scale_x_continuous(limits = UNIFY_HIST_AR))
if (exists("p_qin_d3_progeny1")) print(p_qin_d3_progeny1 + ggplot2::scale_y_continuous(limits = UNIFY_PW_Y))
if (exists("p_qin_d3_progeny2")) print(p_qin_d3_progeny2 + ggplot2::scale_size_continuous(limits = UNIFY_LIGW_ABS) + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = UNIFY_LIGW_COL, midpoint = 0))
print("✓ POST-HOC cross-cohort unified ggplot replot complete")

# ============================================================================
# SECTION 8: SUMMARY & CLEANUP
# ============================================================================

print("")
print("✓✓✓ MASTER PIPELINE COMPLETE ✓✓✓")
print("")
print("SUMMARY:")
print("  • Step 0: LIANA Consensus Analysis + PROGENy Pathway Analysis (3 datasets)")
print("  • Step 1: Lee Day 1 Analysis - 30 visualizations (including PROGENy pathway analysis)")
print("  • Step 2: Lee Day 3 Analysis - 21 visualizations (including PROGENy pathway analysis)")
print("  • Step 3: Wang Day 3 Analysis - 18 visualizations (including PROGENy pathway analysis)")
print("  ==========================================")
print("  TOTAL: 69 PUBLICATION-READY PLOTS")
print("  ==========================================")
print("")
print("PROGENy Pathway Analysis (Exploratory):")
print("  • Pathway activity scores: Arg1pos vs Arg1neg neutrophils")
print("  • Statistical comparison: Wilcoxon tests for pathway differences (all PROGENy pathways)")
print("  • Ligand-to-pathway mapping: Links identified ligands to all PROGENy pathways")
print("  • Output files: *_PROGENy_PathwayComparison.csv, *_PROGENy_PathwayResults.csv")
print("  • Output files: *_LigandToPathwayMapping.csv")
print("")
print("Output Directory: ", OUTPUT_DIR)
print("All results saved as:")
print("  • .png files (high-resolution publication-ready figures)")
print("  • .csv files (data tables)")
print("  • .rds files (R objects for downstream analysis)")
print("")

# Save session info for reproducibility
sessionInfo_output <- capture.output(sessionInfo())
writeLines(sessionInfo_output, file.path(OUTPUT_DIR, "SessionInfo.txt"))
print("✓ Session information saved to SessionInfo.txt")

# Final cleanup
gc()

print("Pipeline execution completed successfully!")
