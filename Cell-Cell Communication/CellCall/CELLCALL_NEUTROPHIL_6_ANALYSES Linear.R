#!/usr/bin/env Rscript

# ============================================================================
# CELLCALL NEUTROPHIL SIGNALING ANALYSIS - 6 SEPARATE ANALYSES
# Research Question: What signals do neutrophils receive to turn on ARG1?
# Comparing Day 1 vs Day 3 to identify temporal changes
#  - No helpers/wrappers, no long if-else
#  - Step-by-step CellCall workflow for each analysis
#  - Uses ACTUAL ARG1 gene expression (not random assignment)
#  - 6 separate analyses: Lee Day1/3 ARG1+/-, Wang Day3 ARG1+/-
# ============================================================================
#
# RUN ORDER: after CellChat_ARG1_Final.R (optional but recommended for parallel LR priors in downstream work).
# Required in getwd(): LeeDat.rds, WangDat.rds
# Rerun from: line 22 (suppressPackageStartupMessages) for a full run; minimum refresh after the Arg1+ fix:
#   from line 183 (STEP 1.7 TransCommuProfile Lee Day 1 ARG1+) through the saveRDS for that analysis.
#
# ---- Working directory ----
CELLCOMM_PROJECT_ROOT <- Sys.getenv("CELLCOMM_PROJECT_ROOT", unset = "")
if (nzchar(CELLCOMM_PROJECT_ROOT) && dir.exists(CELLCOMM_PROJECT_ROOT)) setwd(CELLCOMM_PROJECT_ROOT)

suppressPackageStartupMessages({
  library(Seurat)
  library(cellcall)
  library(ggplot2)
  library(dplyr)
  library(htmlwidgets)
})

# Memory management - Full 256GB RAM system (no limitations)
options(future.globals.maxSize = 256000 * 1024^2)  # 256GB limit for future operations
if (.Platform$OS.type == "windows" && exists("memory.limit", mode = "function")) memory.limit(size = 256000)
options(expressions = 2000000)  # Maximum expression limit for large datasets
gc()  # Clean up memory

cat("Memory limits set for full 256GB RAM system:\n")
cat("- future.globals.maxSize: 256GB\n")
cat("- memory.limit: 256GB\n")
cat("- expressions: 2000000\n")
cat("No limitations - using full physical RAM capacity.\n\n")

# ============================================================================
# STEP 1: LOAD DATASETS
# ============================================================================

cat("=== STEP 1: LOADING DATASETS ===\n")
cat("Lee dataset: timepoints 0, 1, 3, 7 (HAS both Day 1 and Day 3)\n")
cat("Wang dataset: timepoints 0, 3, 14 (HAS Day 3 only)\n")
cat("Brennan dataset: timepoints 0, 7, 28 (NO Day 1 or Day 3)\n\n")

# Load Lee dataset (for Day 1 and Day 3)
lee_data <- readRDS("LeeDat.rds")
DefaultAssay(lee_data) <- "RNA"

# Load Wang dataset (for Day 3 only)
wang_data <- readRDS("WangDat.rds")
DefaultAssay(wang_data) <- "RNA"

cat("Lee dataset:", ncol(lee_data), "cells\n")
cat("Wang dataset:", ncol(wang_data), "cells\n")

# Fix Lee dataset time values
cat("Fixing Lee dataset time values...\n")
lee_data$time <- as.character(lee_data$time)
lee_data$time[lee_data$time == "Uninjured"] <- 0
lee_data$time[lee_data$time == "1dpi"] <- 1
lee_data$time[lee_data$time == "3dpi"] <- 3
lee_data$time[lee_data$time == "7dpi"] <- 7
lee_data$time <- as.numeric(lee_data$time)

cat("Lee dataset timepoints:", unique(lee_data$time), "\n")
cat("Wang dataset timepoints:", unique(wang_data$time), "\n")

# ============================================================================
# STEP 2: DETERMINE ARG1 STATUS FOR BOTH DATASETS
# ============================================================================

cat("\n=== STEP 2: DETERMINING ARG1 STATUS ===\n")

# Lee dataset ARG1 status
cat("Determining ARG1 status for Lee dataset...\n")
DefaultAssay(lee_data) <- "RNA"
if ("Arg1" %in% rownames(lee_data)) {
  lee_arg1_expr <- GetAssayData(lee_data, assay = "RNA", layer = "data")["Arg1", ]
  lee_data$Arg1_status <- ifelse(lee_arg1_expr > 0, "Arg1pos", "Arg1neg")
  cat("Lee ARG1+ cells:", sum(lee_data$Arg1_status == "Arg1pos"), "\n")
  cat("Lee ARG1- cells:", sum(lee_data$Arg1_status == "Arg1neg"), "\n")
} else {
  cat("ARG1 gene not found in Lee dataset\n")
  lee_data$Arg1_status <- "Unknown"
}

# Wang dataset ARG1 status
cat("Determining ARG1 status for Wang dataset...\n")
DefaultAssay(wang_data) <- "RNA"
if ("Arg1" %in% rownames(wang_data)) {
  wang_arg1_expr <- GetAssayData(wang_data, assay = "RNA", layer = "data")["Arg1", ]
  wang_data$Arg1_status <- ifelse(wang_arg1_expr > 0, "Arg1pos", "Arg1neg")
  cat("Wang ARG1+ cells:", sum(wang_data$Arg1_status == "Arg1pos"), "\n")
  cat("Wang ARG1- cells:", sum(wang_data$Arg1_status == "Arg1neg"), "\n")
} else {
  cat("ARG1 gene not found in Wang dataset\n")
  wang_data$Arg1_status <- "Unknown"
}

# Clean up memory
gc()

# Remove intermediate objects to save memory
# Note: Initial data objects will be loaded fresh in each analysis block
gc()

# ============================================================================
# ANALYSIS 1: LEE DAY 1 ARG1+ SIGNALING
# ============================================================================
# PURPOSE: Analyze cell-cell communication signals received by ARG1+ neutrophils on Day 1
# DATASET: Lee dataset, timepoint 1 (1dpi)
# TARGET: Neutrophils with ARG1 expression > 0
# COMPARISON: Will be compared with ARG1- neutrophils (Analysis 2) and Day 3 (Analysis 3)
# OUTPUT: 5 CellCall visualizations + RDS file for re-analysis

cat("\n=== ANALYSIS 1: LEE DAY 1 ARG1+ SIGNALING ===\n")
cat("Research Question: What signals do ARG1+ neutrophils receive on Day 1?\n")

# --- STEP 1.1: Data Preparation ---
cat("STEP 1.1: Data Preparation\n")
lee_day1_subset <- subset(lee_data, time == 1)
cat("Lee Day 1 cells:", ncol(lee_day1_subset), "\n")

# --- STEP 1.2: Prepare ARG1 status ---
cat("STEP 1.2: Preparing cell type labels...\n")
lee_day1_subset$cell_type_clean <- as.character(lee_day1_subset$celltype)
lee_day1_subset$label_arg1 <- lee_day1_subset$cell_type_clean

# Assign ARG1 status to neutrophils
neutrophil_cells <- colnames(lee_day1_subset)[grepl("Neutrophil", lee_day1_subset$cell_type_clean, ignore.case = TRUE)]
cat("Lee Day 1 neutrophils:", length(neutrophil_cells), "\n")

for (cell in neutrophil_cells) {
  arg1_status <- lee_day1_subset$Arg1_status[colnames(lee_day1_subset) == cell]
  lee_day1_subset$label_arg1[colnames(lee_day1_subset) == cell] <- paste0("Neutrophil_", arg1_status)
}

cat("Lee Day 1 ARG1+ neutrophils:", sum(lee_day1_subset$label_arg1 == "Neutrophil_Arg1pos", na.rm=TRUE), "\n")
cat("Lee Day 1 ARG1- neutrophils:", sum(lee_day1_subset$label_arg1 == "Neutrophil_Arg1neg", na.rm=TRUE), "\n")

# --- STEP 1.3: Subset for ARG1+ analysis ---
cat("STEP 1.3: Preparing data for ARG1+ analysis...\n")
lee_day1_arg1pos_data <- subset(lee_day1_subset, label_arg1 == "Neutrophil_Arg1pos" | !grepl("Neutrophil", label_arg1, ignore.case = TRUE))

# --- STEP 1.4: Fix cell IDs and cell types for CellCall ---
cat("STEP 1.4: Fixing cell IDs and cell types for CellCall compatibility...\n")
new_cell_ids <- gsub("-", "_", colnames(lee_day1_arg1pos_data))
lee_day1_arg1pos_data <- RenameCells(lee_day1_arg1pos_data, new.names = new_cell_ids)
lee_day1_arg1pos_data$label_arg1 <- gsub("-", "", lee_day1_arg1pos_data$label_arg1)

# --- STEP 1.5: Memory optimization ---
rm(lee_day1_subset)
gc()

cat("Lee Day 1 ARG1+ analysis cells:", ncol(lee_day1_arg1pos_data), "\n")
cat("Cell types:", unique(lee_day1_arg1pos_data$label_arg1), "\n")

if (ncol(lee_day1_arg1pos_data) < 100) {
  cat("WARNING: Very few cells for analysis. Consider adjusting parameters.\n")
}

# --- STEP 1.6: Create CellCall object ---
cat("STEP 1.6: Creating CellCall object for Lee Day 1 ARG1+...\n")
gc()

cc_lee_day1_arg1pos <- CreateObject_fromSeurat(
  Seurat.object = lee_day1_arg1pos_data,
  slot = "counts",
  cell_type = "label_arg1",
  Org = "Mus musculus"
)
cat("CellCall object created successfully!\n")

rm(lee_day1_arg1pos_data)
gc()

# --- STEP 1.7: Run TransCommuProfile ---
cat("STEP 1.7: Running TransCommuProfile for Lee Day 1 ARG1+...\n")
gc()

cc_lee_day1_arg1pos <- TransCommuProfile(
  object = cc_lee_day1_arg1pos,
  Org = "Mus musculus",
  pValueCor = 0.1,
  CorValue = 0.15,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  IS_core = TRUE
)
cat("TransCommuProfile completed successfully!\n")

gc()
# --- STEP 1.8: Create Visualizations ---
cat("STEP 1.8: Creating visualizations for Lee Day 1 ARG1+...\n")

# 1.8.1. Circular Plot (Official CellCall Format) - FIXED
cat("1.8.1. Creating ViewInterCircos (FIXED)...\n")
# Identify target routes for the plot
all_routes <- colnames(cc_lee_day1_arg1pos@data$expr_l_r_log2_scale)
target_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)

if (length(target_routes) > 0) {
  # Get cell types for color mapping
  cell_types_in_obj <- unique(sapply(strsplit(target_routes, "-"), `[`, 1))
  cell_types_in_obj <- unique(c(cell_types_in_obj, sapply(strsplit(target_routes, "-"), `[`, 2)))
  cell_types_in_obj <- cell_types_in_obj[!is.na(cell_types_in_obj) & cell_types_in_obj != ""]
  
  if (length(cell_types_in_obj) >= 2) {
    colors_for_circos <- grDevices::rainbow(length(cell_types_in_obj))
    cell_color_df <- data.frame(color = colors_for_circos, row.names = cell_types_in_obj, check.names = FALSE)
    circlize::circos.clear()
    grid::grid.newpage()
    
    # Use correct parameters from official documentation
    p_circos <- ViewInterCircos(
      object = cc_lee_day1_arg1pos,
      font = 2,
      cellColor = cell_color_df,
      lrColor = c("#F16B6F", "#84B1ED"),
      arr.type = "big.arrow",
      arr.length = 0.04,
      trackhight1 = 0.05,
      slot = "expr_l_r_log2_scale",
      linkcolor.from.sender = TRUE,
      gap.degree = 0.5,  # Much smaller gap for targeted analysis
      order.vector = cell_types_in_obj,
      trackhight2 = 0.032,
      track.margin2 = c(0.01, 0.12),
      DIY = FALSE
    )
    print(p_circos)
    cat("  ✓ FIXED ViewInterCircos completed successfully!\n")
  } else {
    cat("  ✗ Not enough cell types for ViewInterCircos\n")
  }
} else {
  cat("  ✗ No neutrophil-receiving routes found for ViewInterCircos\n")
}

# 1.8.2. Heatmap (Official CellCall Format)
cat("1.8.2. Creating viewPheatmap...\n")
p_heatmap <- viewPheatmap(
  object = cc_lee_day1_arg1pos,
  slot = "expr_l_r_log2_scale",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  angle_col = 45,
  main = "Lee Day 1 ARG1+ Neutrophil Signaling"
)
print(p_heatmap)
cat("  viewPheatmap completed successfully!\n")

# 1.8.3. Pathway Analysis (Official CellCall Format - TARGETED)
cat("1.8.3. Running TARGETED pathway analysis (Neutrophils as receivers)...\n")
# Get all routes for pathway analysis
all_routes <- colnames(cc_lee_day1_arg1pos@data$expr_l_r_log2_scale)

# Filter for routes where Neutrophils are the RECEIVER
neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
cat("  Found", length(neutrophil_routes), "routes where Neutrophils are receivers\n")

if (length(neutrophil_routes) > 0) {
  # Run pathway analysis ONLY on targeted routes
  pathway.hyper.list <- lapply(neutrophil_routes, function(route){
    cat("  Processing TARGETED route:", route, "\n")
    tmp <- getHyperPathway(data = cc_lee_day1_arg1pos@data$expr_l_r_log2_scale, 
                          object = cc_lee_day1_arg1pos, cella_cellb = route, Org="Mus musculus")
    return(tmp)
  })
  
  # Filter out NULL results
  pathway.hyper.list <- pathway.hyper.list[!sapply(pathway.hyper.list, is.null)]
  pathway.hyper.list <- pathway.hyper.list[sapply(pathway.hyper.list, nrow) > 0]
  
  if (length(pathway.hyper.list) > 0) {
    myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names(pathway.hyper.list))
    p_bubble <- plotBubble(myPub.df) + 
      ggtitle("Signals Received by ARG1+ Neutrophils: Pathway Analysis")
    print(p_bubble)
    cat("  TARGETED pathway analysis completed successfully!\n")
  } else {
    cat("  No significant pathways were found for targeted analysis.\n")
  }
} else {
  cat("  No neutrophil-receiving routes found for pathway analysis.\n")
}

# 1.8.4. Sankey Plot (Official CellCall Format) - FIXED
cat("1.8.4. Creating Sankey plot (FINAL FIX)...\n")
# Find the best route where Neutrophils are RECEIVERS
comm_matrix <- cc_lee_day1_arg1pos@data$expr_l_r_log2_scale
all_routes <- colnames(comm_matrix)
neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)

if (length(neutrophil_routes) > 0) {
  route_totals <- colSums(comm_matrix[, neutrophil_routes, drop = FALSE], na.rm = TRUE)
  best_route <- names(sort(route_totals, decreasing = TRUE))[1]
  
  cat("  Using TARGETED route:", best_route, "\n")
  
  # Step 1: Split the route into sender and receiver
  sender <- sub("-.*$", "", best_route)
  receiver <- sub(".*-", "", best_route)
  
  # Step 2: Run LR2TF analysis
  cat("  Step 1: Running LR2TF analysis...\n")
  cc_lee_day1_arg1pos <- LR2TF(
    object = cc_lee_day1_arg1pos,
    sender_cell = sender,
    recevier_cell = receiver,
    slot = "expr_l_r_log2_scale",
    org = "Mus musculus"
  )
  cat("  ✓ LR2TF analysis completed!\n")
  
  # **CRUCIAL FIX: Check if the analysis produced any results**
  cat("  Step 2: Checking for significant L-R-TF connections...\n")
  lr2tf_results <- cc_lee_day1_arg1pos@tf.result$LR2TF.network[[best_route]]
  
  if (!is.null(lr2tf_results) && nrow(lr2tf_results) > 0) {
    cat("  ✓ Found", nrow(lr2tf_results), "connections. Creating plot...\n")
    
    # Now, create the plot
    p_sankey <- LRT.Dimplot(cc_lee_day1_arg1pos)
    print(p_sankey)
    cat("  ✓ Sankey plot completed successfully!\n")
    
  } else {
    # This is the new, informative message if no results are found
    cat("  ✗ Sankey plot failed: No significant L-R-TF connections were found by the LR2TF analysis for this route.\n")
    cat("  💡 TIP: The analysis worked, but the biological signal might be too weak with default filters. Try relaxing the pValue/corValue in the LR2TF function or test a different route (e.g., Macrophage-Neutrophil).\n")
  }
  
} else {
  cat("  ✗ No neutrophil-receiving routes found for Sankey plot.\n")
}

# 1.8.5. Ridge Plot (Official CellCall Format)
cat("1.8.5. Creating Ridge plot...\n")
# Get GSEA object
if (length(cc_lee_day1_arg1pos@data$gsea.list) > 0) {
  cell_type <- names(cc_lee_day1_arg1pos@data$gsea.list)[1]
  egmt <- cc_lee_day1_arg1pos@data$gsea.list[[cell_type]]
  
  # Filter TFs
  egmt.df <- data.frame(egmt)
  flag.index <- which(egmt.df$p.adjust < 0.05)
  
  if (length(flag.index) > 0) {
    p_ridge <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, 
                             core_enrichment = T, orderBy = "NES", decreasing = FALSE)
    print(p_ridge)
    cat("  Ridge plot completed successfully!\n")
  } else {
    cat("  No significant TFs found for Ridge plot\n")
  }
} else {
  cat("  No GSEA data available for Ridge plot\n")
}

# 1.8.6. TF Enrichment Plot (Official CellCall Format)
cat("1.8.6. Creating TF enrichment plot...\n")
if (length(cc_lee_day1_arg1pos@data$gsea.list) > 0) {
  cell_type <- names(cc_lee_day1_arg1pos@data$gsea.list)[1]
  tf_names <- names(cc_lee_day1_arg1pos@data$gsea.list[[cell_type]]@geneSets)
  
  if (length(tf_names) > 0) {
    selected_tfs <- head(tf_names, 3)
    p_tf <- getGSEAplot(gsea.list=cc_lee_day1_arg1pos@data$gsea.list, 
                        geneSetID=selected_tfs, 
                        myCelltype=cell_type, 
                        fc.list=cc_lee_day1_arg1pos@data$fc.list)
    print(p_tf)
    cat("  TF enrichment plot completed successfully!\n")
  } else {
    cat("  No TFs available for enrichment plot\n")
  }
} else {
  cat("  No GSEA data available for TF enrichment plot\n")
}

# --- STEP 1.9: Save Results ---
cat("STEP 1.9: Saving results...\n")
saveRDS(cc_lee_day1_arg1pos, "CellCall_Lee_Day1_Arg1pos.rds")
cat("  Analysis saved!\n")

# --- STEP 1.10: Cleanup ---
cat("STEP 1.10: Cleaning up memory after Analysis 1...\n")
rm(lee_day1_subset, lee_day1_arg1pos_data, pw_lee_day1_arg1pos)
gc()
cat("Memory cleaned. cc_lee_day1_arg1pos kept for visualization.\n")

cat("=== ANALYSIS 1 COMPLETED ===\n\n")

# ============================================================================
# ANALYSIS 2: LEE DAY 1 ARG1- SIGNALING
# ============================================================================
# PURPOSE: Analyze cell-cell communication signals received by ARG1- neutrophils on Day 1
# DATASET: Lee dataset, timepoint 1 (1dpi)
# TARGET: Neutrophils with ARG1 expression = 0
# COMPARISON: Direct comparison with ARG1+ neutrophils (Analysis 1) to identify ARG1-specific signals
# OUTPUT: 5 CellCall visualizations + RDS file for re-analysis

cat("\n=== ANALYSIS 2: LEE DAY 1 ARG1- SIGNALING ===\n")
cat("Research Question: What signals do ARG1- neutrophils receive on Day 1?\n")

# --- STEP 2.1: Data Preparation ---
cat("STEP 2.1: Data Preparation\n")
lee_day1_subset <- subset(lee_data, time == 1)
cat("Lee Day 1 cells:", ncol(lee_day1_subset), "\n")

# --- STEP 2.2: Prepare ARG1 status ---
cat("STEP 2.2: Preparing cell type labels...\n")
lee_day1_subset$cell_type_clean <- as.character(lee_day1_subset$celltype)
lee_day1_subset$label_arg1 <- lee_day1_subset$cell_type_clean

# Assign ARG1 status to neutrophils
neutrophil_cells <- colnames(lee_day1_subset)[grepl("Neutrophil", lee_day1_subset$cell_type_clean, ignore.case = TRUE)]
cat("Lee Day 1 neutrophils:", length(neutrophil_cells), "\n")

for (cell in neutrophil_cells) {
  arg1_status <- lee_day1_subset$Arg1_status[colnames(lee_day1_subset) == cell]
  lee_day1_subset$label_arg1[colnames(lee_day1_subset) == cell] <- paste0("Neutrophil_", arg1_status)
}

cat("Lee Day 1 ARG1+ neutrophils:", sum(lee_day1_subset$label_arg1 == "Neutrophil_Arg1pos", na.rm=TRUE), "\n")
cat("Lee Day 1 ARG1- neutrophils:", sum(lee_day1_subset$label_arg1 == "Neutrophil_Arg1neg", na.rm=TRUE), "\n")

# --- STEP 2.3: Subset for ARG1- analysis ---
cat("STEP 2.3: Preparing data for ARG1- analysis...\n")
lee_day1_arg1neg_data <- subset(lee_day1_subset, label_arg1 == "Neutrophil_Arg1neg" | !grepl("Neutrophil", label_arg1, ignore.case = TRUE))

# --- STEP 2.4: Fix cell IDs and cell types for CellCall ---
cat("STEP 2.4: Fixing cell IDs and cell types for CellCall compatibility...\n")
new_cell_ids <- gsub("-", "_", colnames(lee_day1_arg1neg_data))
lee_day1_arg1neg_data <- RenameCells(lee_day1_arg1neg_data, new.names = new_cell_ids)
lee_day1_arg1neg_data$label_arg1 <- gsub("-", "", lee_day1_arg1neg_data$label_arg1)

# --- STEP 2.5: Memory optimization ---
rm(lee_day1_subset)
gc()

cat("Lee Day 1 ARG1- analysis cells:", ncol(lee_day1_arg1neg_data), "\n")
cat("Cell types:", unique(lee_day1_arg1neg_data$label_arg1), "\n")

if (ncol(lee_day1_arg1neg_data) < 100) {
  cat("WARNING: Very few cells for analysis. Consider adjusting parameters.\n")
}

# --- STEP 2.6: Create CellCall object ---
cat("STEP 2.6: Creating CellCall object for Lee Day 1 ARG1-...\n")
gc()

cc_lee_day1_arg1neg <- CreateObject_fromSeurat(
  Seurat.object = lee_day1_arg1neg_data,
  slot = "counts",
  cell_type = "label_arg1",
  Org = "Mus musculus"
)
cat("CellCall object created successfully!\n")

rm(lee_day1_arg1neg_data)
gc()

# --- STEP 2.7: Run TransCommuProfile ---
cat("STEP 2.7: Running TransCommuProfile for Lee Day 1 ARG1-...\n")
gc()

cc_lee_day1_arg1neg <- TransCommuProfile(
  object = cc_lee_day1_arg1neg,
  Org = "Mus musculus",
  pValueCor = 0.1,
  CorValue = 0.15,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  IS_core = TRUE
)
cat("TransCommuProfile completed successfully!\n")

gc()

# --- STEP 2.8: Create Visualizations ---
cat("STEP 2.8: Creating visualizations for Lee Day 1 ARG1-...\n")

# 2.8.1. Circular Plot (Official CellCall Format) - FIXED
cat("2.8.1. Creating ViewInterCircos (FIXED)...\n")
# Identify target routes for the plot
all_routes <- colnames(cc_lee_day1_arg1neg@data$expr_l_r_log2_scale)
target_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)

if (length(target_routes) > 0) {
  # Get cell types for color mapping
  cell_types_in_obj <- unique(sapply(strsplit(target_routes, "-"), `[`, 1))
  cell_types_in_obj <- unique(c(cell_types_in_obj, sapply(strsplit(target_routes, "-"), `[`, 2)))
  cell_types_in_obj <- cell_types_in_obj[!is.na(cell_types_in_obj) & cell_types_in_obj != ""]
  
  if (length(cell_types_in_obj) >= 2) {
    colors_for_circos <- grDevices::rainbow(length(cell_types_in_obj))
    cell_color_df <- data.frame(color = colors_for_circos, row.names = cell_types_in_obj, check.names = FALSE)
    circlize::circos.clear()
    grid::grid.newpage()
    
    # Use correct parameters from official documentation
    p_circos <- ViewInterCircos(
      object = cc_lee_day1_arg1neg,
      font = 2,
      cellColor = cell_color_df,
      lrColor = c("#F16B6F", "#84B1ED"),
      arr.type = "big.arrow",
      arr.length = 0.04,
      trackhight1 = 0.05,
      slot = "expr_l_r_log2_scale",
      linkcolor.from.sender = TRUE,
      gap.degree = 0.5,  # Much smaller gap for targeted analysis
      order.vector = cell_types_in_obj,
      trackhight2 = 0.032,
      track.margin2 = c(0.01, 0.12),
      DIY = FALSE
    )
    print(p_circos)
    cat("  ✓ FIXED ViewInterCircos completed successfully!\n")
  } else {
    cat("  ✗ Not enough cell types for ViewInterCircos\n")
  }
} else {
  cat("  ✗ No neutrophil-receiving routes found for ViewInterCircos\n")
}

# 2.8.2. Heatmap (Official CellCall Format)
cat("2.8.2. Creating viewPheatmap...\n")
p_heatmap <- viewPheatmap(
  object = cc_lee_day1_arg1neg,
  slot = "expr_l_r_log2_scale",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  angle_col = 45,
  main = "Lee Day 1 ARG1- Neutrophil Signaling"
)
print(p_heatmap)
cat("  viewPheatmap completed successfully!\n")

# 2.8.3. Pathway Analysis (Official CellCall Format - TARGETED)
cat("2.8.3. Running TARGETED pathway analysis (Neutrophils as receivers)...\n")
# Get all routes for pathway analysis
all_routes <- colnames(cc_lee_day1_arg1neg@data$expr_l_r_log2_scale)

# Filter for routes where Neutrophils are the RECEIVER
neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
cat("  Found", length(neutrophil_routes), "routes where Neutrophils are receivers\n")

if (length(neutrophil_routes) > 0) {
  # Run pathway analysis ONLY on targeted routes
  pathway.hyper.list <- lapply(neutrophil_routes, function(route){
    cat("  Processing TARGETED route:", route, "\n")
    tmp <- getHyperPathway(data = cc_lee_day1_arg1neg@data$expr_l_r_log2_scale, 
                          object = cc_lee_day1_arg1neg, cella_cellb = route, Org="Mus musculus")
    return(tmp)
  })
  
  # Filter out NULL results
  pathway.hyper.list <- pathway.hyper.list[!sapply(pathway.hyper.list, is.null)]
  pathway.hyper.list <- pathway.hyper.list[sapply(pathway.hyper.list, nrow) > 0]
  
  if (length(pathway.hyper.list) > 0) {
    myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names(pathway.hyper.list))
    p_bubble <- plotBubble(myPub.df) + 
      ggtitle("Signals Received by ARG1- Neutrophils: Pathway Analysis")
    print(p_bubble)
    cat("  TARGETED pathway analysis completed successfully!\n")
  } else {
    cat("  No significant pathways were found for targeted analysis.\n")
  }
} else {
  cat("  No neutrophil-receiving routes found for pathway analysis.\n")
}

# 2.8.4. Sankey Plot (Official CellCall Format) - FIXED
cat("2.8.4. Creating Sankey plot (FINAL FIX)...\n")

  # Find the best route where Neutrophils are RECEIVERS
  comm_matrix <- cc_lee_day1_arg1neg@data$expr_l_r_log2_scale
  all_routes <- colnames(comm_matrix)
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(neutrophil_routes) > 0) {
    route_totals <- colSums(comm_matrix[, neutrophil_routes, drop = FALSE], na.rm = TRUE)
    best_route <- names(sort(route_totals, decreasing = TRUE))[1]
    
    cat("  Using TARGETED route:", best_route, "\n")
    
    # Step 1: Split the route into sender and receiver
    sender <- sub("-.*$", "", best_route)
    receiver <- sub(".*-", "", best_route)
    
    # Step 2: Run LR2TF analysis
    cat("  Step 1: Running LR2TF analysis...\n")
    cc_lee_day1_arg1neg <- LR2TF(
      object = cc_lee_day1_arg1neg,
      sender_cell = sender,
      recevier_cell = receiver,
      slot = "expr_l_r_log2_scale",
      org = "Mus musculus"
    )
    cat("  ✓ LR2TF analysis completed!\n")
    
    # **CRUCIAL FIX: Check if the analysis produced any results**
    cat("  Step 2: Checking for significant L-R-TF connections...\n")
    lr2tf_results <- cc_lee_day1_arg1neg@tf.result$LR2TF.network[[best_route]]
    
    if (!is.null(lr2tf_results) && nrow(lr2tf_results) > 0) {
      cat("  ✓ Found", nrow(lr2tf_results), "connections. Creating plot...\n")
      
      # Now, create the plot
      p_sankey <- LRT.Dimplot(cc_lee_day1_arg1neg)
      print(p_sankey)
      cat("  ✓ Sankey plot completed successfully!\n")
      
    } else {
      # This is the new, informative message if no results are found
      cat("  ✗ Sankey plot failed: No significant L-R-TF connections were found by the LR2TF analysis for this route.\n")
      cat("  💡 TIP: The analysis worked, but the biological signal might be too weak with default filters. Try relaxing the pValue/corValue in the LR2TF function or test a different route (e.g., Macrophage-Neutrophil).\n")
    }
    
  } else {
    cat("  ✗ No neutrophil-receiving routes found for Sankey plot.\n")
  }
  

# 2.8.5. Ridge Plot (Official CellCall Format)
cat("2.8.5. Creating Ridge plot...\n")

  # Get GSEA object
  if (length(cc_lee_day1_arg1neg@data$gsea.list) > 0) {
    cell_type <- names(cc_lee_day1_arg1neg@data$gsea.list)[1]
    egmt <- cc_lee_day1_arg1neg@data$gsea.list[[cell_type]]
    
    # Filter TFs
    egmt.df <- data.frame(egmt)
    flag.index <- which(egmt.df$p.adjust < 0.05)
    
    if (length(flag.index) > 0) {
      p_ridge <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, 
                               core_enrichment = T, orderBy = "NES", decreasing = FALSE)
      print(p_ridge)
      cat("  Ridge plot completed successfully!\n")
    } else {
      cat("  No significant TFs found for Ridge plot\n")
    }
  } else {
    cat("  No GSEA data available for Ridge plot\n")
  }

# 2.8.6. TF Enrichment Plot (Official CellCall Format)
cat("2.8.6. Creating TF enrichment plot...\n")

  if (length(cc_lee_day1_arg1neg@data$gsea.list) > 0) {
    cell_type <- names(cc_lee_day1_arg1neg@data$gsea.list)[1]
    tf_names <- names(cc_lee_day1_arg1neg@data$gsea.list[[cell_type]]@geneSets)
    
    if (length(tf_names) > 0) {
      selected_tfs <- head(tf_names, 3)
      p_tf <- getGSEAplot(gsea.list=cc_lee_day1_arg1neg@data$gsea.list, 
                          geneSetID=selected_tfs, 
                          myCelltype=cell_type, 
                          fc.list=cc_lee_day1_arg1neg@data$fc.list)
      print(p_tf)
      cat("  TF enrichment plot completed successfully!\n")
    } else {
      cat("  No TFs available for enrichment plot\n")
    }
  } else {
    cat("  No GSEA data available for TF enrichment plot\n")
  }

# --- STEP 2.9: Save Results ---
cat("STEP 2.9: Saving results...\n")
saveRDS(cc_lee_day1_arg1neg, "CellCall_Lee_Day1_Arg1neg.rds")
cat("  Analysis saved!\n")

# --- STEP 2.10: Cleanup ---
cat("STEP 2.10: Cleaning up memory after Analysis 2...\n")
rm(lee_day1_subset, lee_day1_arg1neg_data, pw_lee_day1_arg1neg)
gc()
cat("Memory cleaned. cc_lee_day1_arg1neg kept for visualization.\n")

cat("=== ANALYSIS 2 COMPLETED ===\n\n")

# ============================================================================
# ANALYSIS 3: LEE DAY 3 ARG1+ SIGNALING
# ============================================================================
# PURPOSE: Analyze cell-cell communication signals received by ARG1+ neutrophils on Day 3
# DATASET: Lee dataset, timepoint 3 (3dpi)
# TARGET: Neutrophils with ARG1 expression > 0
# COMPARISON: Temporal comparison with Day 1 ARG1+ (Analysis 1) to identify time-dependent changes
# OUTPUT: 5 CellCall visualizations + RDS file for re-analysis

cat("\n=== ANALYSIS 3: LEE DAY 3 ARG1+ SIGNALING ===\n")
cat("Research Question: What signals do ARG1+ neutrophils receive on Day 3?\n")

# --- STEP 3.1: Data Preparation ---
cat("STEP 3.1: Data Preparation\n")
lee_day3_subset <- subset(lee_data, time == 3)
cat("Lee Day 3 cells:", ncol(lee_day3_subset), "\n")

# --- STEP 3.2: Prepare ARG1 status ---
cat("STEP 3.2: Preparing cell type labels...\n")
lee_day3_subset$cell_type_clean <- as.character(lee_day3_subset$celltype)
lee_day3_subset$label_arg1 <- lee_day3_subset$cell_type_clean

# Assign ARG1 status to neutrophils
neutrophil_cells <- colnames(lee_day3_subset)[grepl("Neutrophil", lee_day3_subset$cell_type_clean, ignore.case = TRUE)]
cat("Lee Day 3 neutrophils:", length(neutrophil_cells), "\n")

for (cell in neutrophil_cells) {
  arg1_status <- lee_day3_subset$Arg1_status[colnames(lee_day3_subset) == cell]
  lee_day3_subset$label_arg1[colnames(lee_day3_subset) == cell] <- paste0("Neutrophil_", arg1_status)
}

cat("Lee Day 3 ARG1+ neutrophils:", sum(lee_day3_subset$label_arg1 == "Neutrophil_Arg1pos", na.rm=TRUE), "\n")
cat("Lee Day 3 ARG1- neutrophils:", sum(lee_day3_subset$label_arg1 == "Neutrophil_Arg1neg", na.rm=TRUE), "\n")

# --- STEP 3.3: Subset for ARG1+ analysis ---
cat("STEP 3.3: Preparing data for ARG1+ analysis...\n")
lee_day3_arg1pos_data <- subset(lee_day3_subset, label_arg1 == "Neutrophil_Arg1pos" | !grepl("Neutrophil", label_arg1, ignore.case = TRUE))

# --- STEP 3.4: Fix cell IDs and cell types for CellCall ---
cat("STEP 3.4: Fixing cell IDs and cell types for CellCall compatibility...\n")
new_cell_ids <- gsub("-", "_", colnames(lee_day3_arg1pos_data))
lee_day3_arg1pos_data <- RenameCells(lee_day3_arg1pos_data, new.names = new_cell_ids)
lee_day3_arg1pos_data$label_arg1 <- gsub("-", "", lee_day3_arg1pos_data$label_arg1)
# --- STEP 3.5: Memory optimization ---
rm(lee_day3_subset)
gc()

cat("Lee Day 3 ARG1+ analysis cells:", ncol(lee_day3_arg1pos_data), "\n")
cat("Cell types:", unique(lee_day3_arg1pos_data$label_arg1), "\n")

if (ncol(lee_day3_arg1pos_data) < 100) {
  cat("WARNING: Very few cells for analysis. Consider adjusting parameters.\n")
}

# --- STEP 3.6: Create CellCall object ---
cat("STEP 3.6: Creating CellCall object for Lee Day 3 ARG1+...\n")
gc()

cc_lee_day3_arg1pos <- CreateObject_fromSeurat(
  Seurat.object = lee_day3_arg1pos_data,
  slot = "counts",
  cell_type = "label_arg1",
  Org = "Mus musculus"
)
cat("CellCall object created successfully!\n")

rm(lee_day3_arg1pos_data)
gc()

# --- STEP 3.7: Run TransCommuProfile ---
cat("STEP 3.7: Running TransCommuProfile for Lee Day 3 ARG1+...\n")
gc()

cc_lee_day3_arg1pos <- TransCommuProfile(
  object = cc_lee_day3_arg1pos,
  Org = "Mus musculus",
  pValueCor = 0.1,
  CorValue = 0.15,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  IS_core = TRUE
)
cat("TransCommuProfile completed successfully!\n")

gc()

# --- STEP 3.8: Create Visualizations ---
cat("STEP 3.8: Creating visualizations for Lee Day 3 ARG1+...\n")

# 3.8.1. Circular Plot (Official CellCall Format) - FIXED
cat("3.8.1. Creating ViewInterCircos (FIXED)...\n")

  # Identify target routes for the plot
  all_routes <- colnames(cc_lee_day3_arg1pos@data$expr_l_r_log2_scale)
  target_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(target_routes) > 0) {
    # Get cell types for color mapping
    cell_types_in_obj <- unique(sapply(strsplit(target_routes, "-"), `[`, 1))
    cell_types_in_obj <- unique(c(cell_types_in_obj, sapply(strsplit(target_routes, "-"), `[`, 2)))
    cell_types_in_obj <- cell_types_in_obj[!is.na(cell_types_in_obj) & cell_types_in_obj != ""]
    
    if (length(cell_types_in_obj) >= 2) {
      colors_for_circos <- grDevices::rainbow(length(cell_types_in_obj))
      cell_color_df <- data.frame(color = colors_for_circos, row.names = cell_types_in_obj, check.names = FALSE)
      circlize::circos.clear()
      grid::grid.newpage()
      
      # Use correct parameters from official documentation
      p_circos <- ViewInterCircos(
        object = cc_lee_day3_arg1pos,
        font = 2,
        cellColor = cell_color_df,
        lrColor = c("#F16B6F", "#84B1ED"),
        arr.type = "big.arrow",
        arr.length = 0.04,
        trackhight1 = 0.05,
        slot = "expr_l_r_log2_scale",
        linkcolor.from.sender = TRUE,
        gap.degree = 0.5,  # Much smaller gap for targeted analysis
        order.vector = cell_types_in_obj,
        trackhight2 = 0.032,
        track.margin2 = c(0.01, 0.12),
        DIY = FALSE
      )
      print(p_circos)
      cat("  ✓ FIXED ViewInterCircos completed successfully!\n")
    } else {
      cat("  ✗ Not enough cell types for ViewInterCircos\n")
    }
  } else {
    cat("  ✗ No neutrophil-receiving routes found for ViewInterCircos\n")
  }

# 3.8.2. Heatmap (Official CellCall Format)
cat("3.8.2. Creating viewPheatmap...\n")

  p_heatmap <- viewPheatmap(
    object = cc_lee_day3_arg1pos,
    slot = "expr_l_r_log2_scale",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    angle_col = 45,
    main = "Lee Day 3 ARG1+ Neutrophil Signaling"
  )
  print(p_heatmap)
  cat("  viewPheatmap completed successfully!\n")

# 3.8.3. Pathway Analysis (Official CellCall Format - TARGETED)
cat("3.8.3. Running TARGETED pathway analysis (Neutrophils as receivers)...\n")

  # Get all routes for pathway analysis
  all_routes <- colnames(cc_lee_day3_arg1pos@data$expr_l_r_log2_scale)
  
  # Filter for routes where Neutrophils are the RECEIVER
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  cat("  Found", length(neutrophil_routes), "routes where Neutrophils are receivers\n")
  
  if (length(neutrophil_routes) > 0) {
    # Run pathway analysis ONLY on targeted routes
    pathway.hyper.list <- lapply(neutrophil_routes, function(route){
      cat("  Processing TARGETED route:", route, "\n")
      tmp <- getHyperPathway(data = cc_lee_day3_arg1pos@data$expr_l_r_log2_scale, 
                            object = cc_lee_day3_arg1pos, cella_cellb = route, Org="Mus musculus")
      return(tmp)
    })
    
    # Filter out NULL results
    pathway.hyper.list <- pathway.hyper.list[!sapply(pathway.hyper.list, is.null)]
    pathway.hyper.list <- pathway.hyper.list[sapply(pathway.hyper.list, nrow) > 0]
    
    if (length(pathway.hyper.list) > 0) {
      myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names(pathway.hyper.list))
      p_bubble <- plotBubble(myPub.df) + 
        ggtitle("Signals Received by ARG1+ Neutrophils: Pathway Analysis")
      print(p_bubble)
      cat("  TARGETED pathway analysis completed successfully!\n")
    } else {
      cat("  No significant pathways were found for targeted analysis.\n")
    }
  } else {
    cat("  No neutrophil-receiving routes found for pathway analysis.\n")
  }

# 3.8.4. Sankey Plot (Official CellCall Format) - FIXED
cat("3.8.4. Creating Sankey plot (FINAL FIX)...\n")

  # Find the best route where Neutrophils are RECEIVERS
  comm_matrix <- cc_lee_day3_arg1pos@data$expr_l_r_log2_scale
  all_routes <- colnames(comm_matrix)
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(neutrophil_routes) > 0) {
    route_totals <- colSums(comm_matrix[, neutrophil_routes, drop = FALSE], na.rm = TRUE)
    best_route <- names(sort(route_totals, decreasing = TRUE))[1]
    
    cat("  Using TARGETED route:", best_route, "\n")
    
    # Step 1: Split the route into sender and receiver
    sender <- sub("-.*$", "", best_route)
    receiver <- sub(".*-", "", best_route)
    
    # Step 2: Run LR2TF analysis
    cat("  Step 1: Running LR2TF analysis...\n")
    cc_lee_day3_arg1pos <- LR2TF(
      object = cc_lee_day3_arg1pos,
      sender_cell = sender,
      recevier_cell = receiver,
      slot = "expr_l_r_log2_scale",
      org = "Mus musculus"
    )
    cat("  ✓ LR2TF analysis completed!\n")
    
    # **CRUCIAL FIX: Check if the analysis produced any results**
    cat("  Step 2: Checking for significant L-R-TF connections...\n")
    lr2tf_results <- cc_lee_day3_arg1pos@tf.result$LR2TF.network[[best_route]]
    
    if (!is.null(lr2tf_results) && nrow(lr2tf_results) > 0) {
      cat("  ✓ Found", nrow(lr2tf_results), "connections. Creating plot...\n")
      
      # Now, create the plot
      p_sankey <- LRT.Dimplot(cc_lee_day3_arg1pos)
      print(p_sankey)
      cat("  ✓ Sankey plot completed successfully!\n")
      
    } else {
      # This is the new, informative message if no results are found
      cat("  ✗ Sankey plot failed: No significant L-R-TF connections were found by the LR2TF analysis for this route.\n")
      cat("  💡 TIP: The analysis worked, but the biological signal might be too weak with default filters. Try relaxing the pValue/corValue in the LR2TF function or test a different route (e.g., Macrophage-Neutrophil).\n")
    }
    
  } else {
    cat("  ✗ No neutrophil-receiving routes found for Sankey plot.\n")
  }
  

# 3.8.5. Ridge Plot (Official CellCall Format)
cat("3.8.5. Creating Ridge plot...\n")

  # Get GSEA object
  if (length(cc_lee_day3_arg1pos@data$gsea.list) > 0) {
    cell_type <- names(cc_lee_day3_arg1pos@data$gsea.list)[1]
    egmt <- cc_lee_day3_arg1pos@data$gsea.list[[cell_type]]
    
    # Filter TFs
    egmt.df <- data.frame(egmt)
    flag.index <- which(egmt.df$p.adjust < 0.05)
    
    if (length(flag.index) > 0) {
      p_ridge <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, 
                               core_enrichment = T, orderBy = "NES", decreasing = FALSE)
      print(p_ridge)
      cat("  Ridge plot completed successfully!\n")
    } else {
      cat("  No significant TFs found for Ridge plot\n")
    }
  } else {
    cat("  No GSEA data available for Ridge plot\n")
  }

# 3.8.6. TF Enrichment Plot (Official CellCall Format)
cat("3.8.6. Creating TF enrichment plot...\n")

  if (length(cc_lee_day3_arg1pos@data$gsea.list) > 0) {
    cell_type <- names(cc_lee_day3_arg1pos@data$gsea.list)[1]
    tf_names <- names(cc_lee_day3_arg1pos@data$gsea.list[[cell_type]]@geneSets)
    
    if (length(tf_names) > 0) {
      selected_tfs <- head(tf_names, 3)
      p_tf <- getGSEAplot(gsea.list=cc_lee_day3_arg1pos@data$gsea.list, 
                          geneSetID=selected_tfs, 
                          myCelltype=cell_type, 
                          fc.list=cc_lee_day3_arg1pos@data$fc.list)
      print(p_tf)
      cat("  TF enrichment plot completed successfully!\n")
    } else {
      cat("  No TFs available for enrichment plot\n")
    }
  } else {
    cat("  No GSEA data available for TF enrichment plot\n")
  }

# --- STEP 3.9: Save Results ---
cat("STEP 3.9: Saving results...\n")
saveRDS(cc_lee_day3_arg1pos, "CellCall_Lee_Day3_Arg1pos.rds")
cat("  Analysis saved!\n")

# --- STEP 3.10: Cleanup ---
cat("STEP 3.10: Cleaning up memory after Analysis 3...\n")
rm(lee_day3_subset, lee_day3_arg1pos_data, pw_lee_day3_arg1pos)
gc()
cat("Memory cleaned. cc_lee_day3_arg1pos kept for visualization.\n")

cat("=== ANALYSIS 3 COMPLETED ===\n\n")

# ============================================================================
# ANALYSIS 4: LEE DAY 3 ARG1- SIGNALING
# ============================================================================
# PURPOSE: Analyze cell-cell communication signals received by ARG1- neutrophils on Day 3
# DATASET: Lee dataset, timepoint 3 (3dpi)
# TARGET: Neutrophils with ARG1 expression = 0
# COMPARISON: Direct comparison with ARG1+ neutrophils (Analysis 3) and temporal with Day 1 ARG1- (Analysis 2)
# OUTPUT: 5 CellCall visualizations + RDS file for re-analysis

cat("\n=== ANALYSIS 4: LEE DAY 3 ARG1- SIGNALING ===\n")
cat("Research Question: What signals do ARG1- neutrophils receive on Day 3?\n")

# --- STEP 4.1: Data Preparation ---
cat("STEP 4.1: Data Preparation\n")
lee_day3_subset <- subset(lee_data, time == 3)
cat("Lee Day 3 cells:", ncol(lee_day3_subset), "\n")

# --- STEP 4.2: Prepare ARG1 status ---
cat("STEP 4.2: Preparing cell type labels...\n")
lee_day3_subset$cell_type_clean <- as.character(lee_day3_subset$celltype)
lee_day3_subset$label_arg1 <- lee_day3_subset$cell_type_clean

# Assign ARG1 status to neutrophils
neutrophil_cells <- colnames(lee_day3_subset)[grepl("Neutrophil", lee_day3_subset$cell_type_clean, ignore.case = TRUE)]
cat("Lee Day 3 neutrophils:", length(neutrophil_cells), "\n")

for (cell in neutrophil_cells) {
  arg1_status <- lee_day3_subset$Arg1_status[colnames(lee_day3_subset) == cell]
  lee_day3_subset$label_arg1[colnames(lee_day3_subset) == cell] <- paste0("Neutrophil_", arg1_status)
}

cat("Lee Day 3 ARG1+ neutrophils:", sum(lee_day3_subset$label_arg1 == "Neutrophil_Arg1pos", na.rm=TRUE), "\n")
cat("Lee Day 3 ARG1- neutrophils:", sum(lee_day3_subset$label_arg1 == "Neutrophil_Arg1neg", na.rm=TRUE), "\n")

# --- STEP 4.3: Subset for ARG1- analysis ---
cat("STEP 4.3: Preparing data for ARG1- analysis...\n")
lee_day3_arg1neg_data <- subset(lee_day3_subset, label_arg1 == "Neutrophil_Arg1neg" | !grepl("Neutrophil", label_arg1, ignore.case = TRUE))

# --- STEP 4.4: Fix cell IDs and cell types for CellCall ---
cat("STEP 4.4: Fixing cell IDs and cell types for CellCall compatibility...\n")
new_cell_ids <- gsub("-", "_", colnames(lee_day3_arg1neg_data))
lee_day3_arg1neg_data <- RenameCells(lee_day3_arg1neg_data, new.names = new_cell_ids)
lee_day3_arg1neg_data$label_arg1 <- gsub("-", "", lee_day3_arg1neg_data$label_arg1)

# --- STEP 4.5: Memory optimization ---
rm(lee_day3_subset)
gc()

cat("Lee Day 3 ARG1- analysis cells:", ncol(lee_day3_arg1neg_data), "\n")
cat("Cell types:", unique(lee_day3_arg1neg_data$label_arg1), "\n")

if (ncol(lee_day3_arg1neg_data) < 100) {
  cat("WARNING: Very few cells for analysis. Consider adjusting parameters.\n")
}

# --- STEP 4.6: Create CellCall object ---
cat("STEP 4.6: Creating CellCall object for Lee Day 3 ARG1-...\n")
gc()

cc_lee_day3_arg1neg <- CreateObject_fromSeurat(
  Seurat.object = lee_day3_arg1neg_data,
  slot = "counts",
  cell_type = "label_arg1",
  Org = "Mus musculus"
)
cat("CellCall object created successfully!\n")

rm(lee_day3_arg1neg_data)
gc()

# --- STEP 4.7: Run TransCommuProfile ---
cat("STEP 4.7: Running TransCommuProfile for Lee Day 3 ARG1-...\n")
gc()

cc_lee_day3_arg1neg <- TransCommuProfile(
  object = cc_lee_day3_arg1neg,
  Org = "Mus musculus",
  pValueCor = 0.1,
  CorValue = 0.15,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  IS_core = TRUE
)
cat("TransCommuProfile completed successfully!\n")

gc()

# --- STEP 4.8: Create Visualizations ---
cat("STEP 4.8: Creating visualizations for Lee Day 3 ARG1-...\n")

# 4.8.1. Circular Plot (Official CellCall Format) - FIXED
cat("4.8.1. Creating ViewInterCircos (FIXED)...\n")

  # Identify target routes for the plot
  all_routes <- colnames(cc_lee_day3_arg1neg@data$expr_l_r_log2_scale)
  target_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(target_routes) > 0) {
    # Get cell types for color mapping
    cell_types_in_obj <- unique(sapply(strsplit(target_routes, "-"), `[`, 1))
    cell_types_in_obj <- unique(c(cell_types_in_obj, sapply(strsplit(target_routes, "-"), `[`, 2)))
    cell_types_in_obj <- cell_types_in_obj[!is.na(cell_types_in_obj) & cell_types_in_obj != ""]
    
    if (length(cell_types_in_obj) >= 2) {
      colors_for_circos <- grDevices::rainbow(length(cell_types_in_obj))
      cell_color_df <- data.frame(color = colors_for_circos, row.names = cell_types_in_obj, check.names = FALSE)
      circlize::circos.clear()
      grid::grid.newpage()
      
      # Use correct parameters from official documentation
      p_circos <- ViewInterCircos(
        object = cc_lee_day3_arg1neg,
        font = 2,
        cellColor = cell_color_df,
        lrColor = c("#F16B6F", "#84B1ED"),
        arr.type = "big.arrow",
        arr.length = 0.04,
        trackhight1 = 0.05,
        slot = "expr_l_r_log2_scale",
        linkcolor.from.sender = TRUE,
        gap.degree = 0.5,  # Much smaller gap for targeted analysis
        order.vector = cell_types_in_obj,
        trackhight2 = 0.032,
        track.margin2 = c(0.01, 0.12),
        DIY = FALSE
      )
      print(p_circos)
      cat("  ✓ FIXED ViewInterCircos completed successfully!\n")
    } else {
      cat("  ✗ Not enough cell types for ViewInterCircos\n")
    }
  } else {
    cat("  ✗ No neutrophil-receiving routes found for ViewInterCircos\n")
  }

# 4.8.2. Heatmap (Official CellCall Format)
cat("4.8.2. Creating viewPheatmap...\n")

  p_heatmap <- viewPheatmap(
    object = cc_lee_day3_arg1neg,
    slot = "expr_l_r_log2_scale",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    angle_col = 45,
    main = "Lee Day 3 ARG1- Neutrophil Signaling"
  )
  print(p_heatmap)
  cat("  viewPheatmap completed successfully!\n")

# 4.8.3. Pathway Analysis (Official CellCall Format - TARGETED)
cat("4.8.3. Running TARGETED pathway analysis (Neutrophils as receivers)...\n")

  # Get all routes for pathway analysis
  all_routes <- colnames(cc_lee_day3_arg1neg@data$expr_l_r_log2_scale)
  
  # Filter for routes where Neutrophils are the RECEIVER
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  cat("  Found", length(neutrophil_routes), "routes where Neutrophils are receivers\n")
  
  if (length(neutrophil_routes) > 0) {
    # Run pathway analysis ONLY on targeted routes
    pathway.hyper.list <- lapply(neutrophil_routes, function(route){
      cat("  Processing TARGETED route:", route, "\n")
      tmp <- getHyperPathway(data = cc_lee_day3_arg1neg@data$expr_l_r_log2_scale, 
                            object = cc_lee_day3_arg1neg, cella_cellb = route, Org="Mus musculus")
      return(tmp)
    })
    
    # Filter out NULL results
    pathway.hyper.list <- pathway.hyper.list[!sapply(pathway.hyper.list, is.null)]
    pathway.hyper.list <- pathway.hyper.list[sapply(pathway.hyper.list, nrow) > 0]
    
    if (length(pathway.hyper.list) > 0) {
      myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names(pathway.hyper.list))
      p_bubble <- plotBubble(myPub.df) + 
        ggtitle("Signals Received by ARG1- Neutrophils: Pathway Analysis")
      print(p_bubble)
      cat("  TARGETED pathway analysis completed successfully!\n")
    } else {
      cat("  No significant pathways were found for targeted analysis.\n")
    }
  } else {
    cat("  No neutrophil-receiving routes found for pathway analysis.\n")
  }

# 4.8.4. Sankey Plot (Official CellCall Format) - FIXED
cat("4.8.4. Creating Sankey plot (FINAL FIX)...\n")

  # Find the best route where Neutrophils are RECEIVERS
  comm_matrix <- cc_lee_day3_arg1neg@data$expr_l_r_log2_scale
  all_routes <- colnames(comm_matrix)
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(neutrophil_routes) > 0) {
    route_totals <- colSums(comm_matrix[, neutrophil_routes, drop = FALSE], na.rm = TRUE)
    best_route <- names(sort(route_totals, decreasing = TRUE))[1]
    
    cat("  Using TARGETED route:", best_route, "\n")
    
    # Step 1: Split the route into sender and receiver
    sender <- sub("-.*$", "", best_route)
    receiver <- sub(".*-", "", best_route)
    
    # Step 2: Run LR2TF analysis
    cat("  Step 1: Running LR2TF analysis...\n")
    cc_lee_day3_arg1neg <- LR2TF(
      object = cc_lee_day3_arg1neg,
      sender_cell = sender,
      recevier_cell = receiver,
      slot = "expr_l_r_log2_scale",
      org = "Mus musculus"
    )
    cat("  ✓ LR2TF analysis completed!\n")
    
    # **CRUCIAL FIX: Check if the analysis produced any results**
    cat("  Step 2: Checking for significant L-R-TF connections...\n")
    lr2tf_results <- cc_lee_day3_arg1neg@tf.result$LR2TF.network[[best_route]]
    
    if (!is.null(lr2tf_results) && nrow(lr2tf_results) > 0) {
      cat("  ✓ Found", nrow(lr2tf_results), "connections. Creating plot...\n")
      
      # Now, create the plot
      p_sankey <- LRT.Dimplot(cc_lee_day3_arg1neg)
      print(p_sankey)
      cat("  ✓ Sankey plot completed successfully!\n")
      
    } else {
      # This is the new, informative message if no results are found
      cat("  ✗ Sankey plot failed: No significant L-R-TF connections were found by the LR2TF analysis for this route.\n")
      cat("  💡 TIP: The analysis worked, but the biological signal might be too weak with default filters. Try relaxing the pValue/corValue in the LR2TF function or test a different route (e.g., Macrophage-Neutrophil).\n")
    }
    
  } else {
    cat("  ✗ No neutrophil-receiving routes found for Sankey plot.\n")
  }
  

# 4.8.5. Ridge Plot (Official CellCall Format)
cat("4.8.5. Creating Ridge plot...\n")

  # Get GSEA object
  if (length(cc_lee_day3_arg1neg@data$gsea.list) > 0) {
    cell_type <- names(cc_lee_day3_arg1neg@data$gsea.list)[1]
    egmt <- cc_lee_day3_arg1neg@data$gsea.list[[cell_type]]
    
    # Filter TFs
    egmt.df <- data.frame(egmt)
    flag.index <- which(egmt.df$p.adjust < 0.05)
    
    if (length(flag.index) > 0) {
      p_ridge <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, 
                               core_enrichment = T, orderBy = "NES", decreasing = FALSE)
      print(p_ridge)
      cat("  Ridge plot completed successfully!\n")
    } else {
      cat("  No significant TFs found for Ridge plot\n")
    }
  } else {
    cat("  No GSEA data available for Ridge plot\n")
  }

# 4.8.6. TF Enrichment Plot (Official CellCall Format)
cat("4.8.6. Creating TF enrichment plot...\n")

  if (length(cc_lee_day3_arg1neg@data$gsea.list) > 0) {
    cell_type <- names(cc_lee_day3_arg1neg@data$gsea.list)[1]
    tf_names <- names(cc_lee_day3_arg1neg@data$gsea.list[[cell_type]]@geneSets)
    
    if (length(tf_names) > 0) {
      selected_tfs <- head(tf_names, 3)
      p_tf <- getGSEAplot(gsea.list=cc_lee_day3_arg1neg@data$gsea.list, 
                          geneSetID=selected_tfs, 
                          myCelltype=cell_type, 
                          fc.list=cc_lee_day3_arg1neg@data$fc.list)
      print(p_tf)
      cat("  TF enrichment plot completed successfully!\n")
    } else {
      cat("  No TFs available for enrichment plot\n")
    }
  } else {
    cat("  No GSEA data available for TF enrichment plot\n")
  }

# --- STEP 4.9: Save Results ---
cat("STEP 4.9: Saving results...\n")
saveRDS(cc_lee_day3_arg1neg, "CellCall_Lee_Day3_Arg1neg.rds")
cat("  Analysis saved!\n")

# --- STEP 4.10: Cleanup ---
cat("STEP 4.10: Cleaning up memory after Analysis 4...\n")
rm(lee_day3_subset, lee_day3_arg1neg_data, pw_lee_day3_arg1neg)
gc()
cat("Memory cleaned. cc_lee_day3_arg1neg kept for visualization.\n")

cat("=== ANALYSIS 4 COMPLETED ===\n\n")

# ============================================================================
# ANALYSIS 5: WANG DAY 3 ARG1+ SIGNALING
# ============================================================================
# PURPOSE: Analyze cell-cell communication signals received by ARG1+ neutrophils on Day 3
# DATASET: Wang dataset, timepoint 3 (3dpi)
# TARGET: Neutrophils with ARG1 expression > 0
# COMPARISON: Cross-dataset validation with Lee Day 3 ARG1+ (Analysis 3) to confirm findings
# OUTPUT: 5 CellCall visualizations + RDS file for re-analysis

cat("\n=== ANALYSIS 5: WANG DAY 3 ARG1+ SIGNALING ===\n")
cat("Research Question: What signals do ARG1+ neutrophils receive on Day 3?\n")

# --- STEP 5.1: Data Preparation ---
cat("STEP 5.1: Data Preparation\n")
wang_day3_subset <- subset(wang_data, time == 3)
cat("Wang Day 3 cells:", ncol(wang_day3_subset), "\n")

# --- STEP 5.2: Prepare ARG1 status ---
cat("STEP 5.2: Preparing cell type labels...\n")
wang_day3_subset$cell_type_clean <- gsub("[-_ ]", "", wang_day3_subset$pruned_labels)
wang_day3_subset$label_arg1 <- wang_day3_subset$cell_type_clean

# Get ARG1 expression for Wang Day 3
wang_day3_arg1_expr <- GetAssayData(wang_day3_subset, assay = "RNA", layer = "data")["Arg1", ]

# Assign ARG1 status to neutrophils
neutrophil_cells <- colnames(wang_day3_subset)[grepl("Neutrophil", wang_day3_subset$cell_type_clean, ignore.case = TRUE)]
cat("Wang Day 3 neutrophils:", length(neutrophil_cells), "\n")

for (cell in neutrophil_cells) {
  arg1_expr <- wang_day3_arg1_expr[cell]
  arg1_status <- ifelse(arg1_expr > 0, "Arg1pos", "Arg1neg")
  wang_day3_subset$label_arg1[colnames(wang_day3_subset) == cell] <- paste0("Neutrophil_", arg1_status)
}

cat("Wang Day 3 ARG1+ neutrophils:", sum(wang_day3_subset$label_arg1 == "Neutrophil_Arg1pos", na.rm=TRUE), "\n")
cat("Wang Day 3 ARG1- neutrophils:", sum(wang_day3_subset$label_arg1 == "Neutrophil_Arg1neg", na.rm=TRUE), "\n")

# --- STEP 5.3: Subset for ARG1+ analysis ---
cat("STEP 5.3: Preparing data for ARG1+ analysis...\n")
wang_day3_arg1pos_data <- subset(wang_day3_subset, label_arg1 == "Neutrophil_Arg1pos" | !grepl("Neutrophil", label_arg1, ignore.case = TRUE))

# --- STEP 5.4: Fix cell IDs and cell types for CellCall ---
cat("STEP 5.4: Fixing cell IDs and cell types for CellCall compatibility...\n")
new_cell_ids <- gsub("-", "_", colnames(wang_day3_arg1pos_data))
wang_day3_arg1pos_data <- RenameCells(wang_day3_arg1pos_data, new.names = new_cell_ids)
wang_day3_arg1pos_data$label_arg1 <- gsub("-", "", wang_day3_arg1pos_data$label_arg1)

# --- STEP 5.5: Memory optimization ---
rm(wang_day3_subset)
gc()

cat("Wang Day 3 ARG1+ analysis cells:", ncol(wang_day3_arg1pos_data), "\n")
cat("Cell types:", unique(wang_day3_arg1pos_data$label_arg1), "\n")

if (ncol(wang_day3_arg1pos_data) < 100) {
  cat("WARNING: Very few cells for analysis. Consider adjusting parameters.\n")
}

# --- STEP 5.6: Create CellCall object ---
cat("STEP 5.6: Creating CellCall object for Wang Day 3 ARG1+...\n")
gc()

cc_wang_day3_arg1pos <- CreateObject_fromSeurat(
  Seurat.object = wang_day3_arg1pos_data,
  slot = "counts",
  cell_type = "label_arg1",
  Org = "Mus musculus"
)
cat("CellCall object created successfully!\n")

rm(wang_day3_arg1pos_data)
gc()

# --- STEP 5.7: Run TransCommuProfile ---
cat("STEP 5.7: Running TransCommuProfile for Wang Day 3 ARG1+...\n")
gc()

cc_wang_day3_arg1pos <- TransCommuProfile(
  object = cc_wang_day3_arg1pos,
  Org = "Mus musculus",
  pValueCor = 0.1,
  CorValue = 0.15,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  IS_core = TRUE
)
cat("TransCommuProfile completed successfully!\n")

gc()

# --- STEP 5.8: Create Visualizations ---
cat("STEP 5.8: Creating visualizations for Wang Day 3 ARG1+...\n")

# 5.8.1. Circular Plot (Official CellCall Format) - FIXED
cat("5.8.1. Creating ViewInterCircos (FIXED)...\n")

  # Identify target routes for the plot
  all_routes <- colnames(cc_wang_day3_arg1pos@data$expr_l_r_log2_scale)
  target_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(target_routes) > 0) {
    # Get cell types for color mapping
    cell_types_in_obj <- unique(sapply(strsplit(target_routes, "-"), `[`, 1))
    cell_types_in_obj <- unique(c(cell_types_in_obj, sapply(strsplit(target_routes, "-"), `[`, 2)))
    cell_types_in_obj <- cell_types_in_obj[!is.na(cell_types_in_obj) & cell_types_in_obj != ""]
    
    if (length(cell_types_in_obj) >= 2) {
      colors_for_circos <- grDevices::rainbow(length(cell_types_in_obj))
      cell_color_df <- data.frame(color = colors_for_circos, row.names = cell_types_in_obj, check.names = FALSE)
      circlize::circos.clear()
      grid::grid.newpage()
      
      # Use correct parameters from official documentation
      p_circos <- ViewInterCircos(
        object = cc_wang_day3_arg1pos,
        font = 2,
        cellColor = cell_color_df,
        lrColor = c("#F16B6F", "#84B1ED"),
        arr.type = "big.arrow",
        arr.length = 0.04,
        trackhight1 = 0.05,
        slot = "expr_l_r_log2_scale",
        linkcolor.from.sender = TRUE,
        gap.degree = 0.5,  # Much smaller gap for targeted analysis
        order.vector = cell_types_in_obj,
        trackhight2 = 0.032,
        track.margin2 = c(0.01, 0.12),
        DIY = FALSE
      )
      print(p_circos)
      cat("  ✓ FIXED ViewInterCircos completed successfully!\n")
    } else {
      cat("  ✗ Not enough cell types for ViewInterCircos\n")
    }
  } else {
    cat("  ✗ No neutrophil-receiving routes found for ViewInterCircos\n")
  }

# 5.8.2. Heatmap (Official CellCall Format)
cat("5.8.2. Creating viewPheatmap...\n")

  p_heatmap <- viewPheatmap(
    object = cc_wang_day3_arg1pos,
    slot = "expr_l_r_log2_scale",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    angle_col = 45,
    main = "Wang Day 3 ARG1+ Neutrophil Signaling"
  )
  print(p_heatmap)
  cat("  viewPheatmap completed successfully!\n")

# 5.8.3. Pathway Analysis (Official CellCall Format - TARGETED)
cat("5.8.3. Running TARGETED pathway analysis (Neutrophils as receivers)...\n")

  # Get all routes for pathway analysis
  all_routes <- colnames(cc_wang_day3_arg1pos@data$expr_l_r_log2_scale)
  
  # Filter for routes where Neutrophils are the RECEIVER
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  cat("  Found", length(neutrophil_routes), "routes where Neutrophils are receivers\n")
  
  if (length(neutrophil_routes) > 0) {
    # Run pathway analysis ONLY on targeted routes
    pathway.hyper.list <- lapply(neutrophil_routes, function(route){
      cat("  Processing TARGETED route:", route, "\n")
      tmp <- getHyperPathway(data = cc_wang_day3_arg1pos@data$expr_l_r_log2_scale, 
                            object = cc_wang_day3_arg1pos, cella_cellb = route, Org="Mus musculus")
      return(tmp)
    })
    
    # Filter out NULL results
    pathway.hyper.list <- pathway.hyper.list[!sapply(pathway.hyper.list, is.null)]
    pathway.hyper.list <- pathway.hyper.list[sapply(pathway.hyper.list, nrow) > 0]
    
    if (length(pathway.hyper.list) > 0) {
      myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names(pathway.hyper.list))
      p_bubble <- plotBubble(myPub.df) + 
        ggtitle("Signals Received by ARG1+ Neutrophils: Pathway Analysis")
      print(p_bubble)
      cat("  TARGETED pathway analysis completed successfully!\n")
    } else {
      cat("  No significant pathways were found for targeted analysis.\n")
    }
  } else {
    cat("  No neutrophil-receiving routes found for pathway analysis.\n")
  }

# 5.8.4. Sankey Plot (Official CellCall Format) - FIXED
cat("5.8.4. Creating Sankey plot (FINAL FIX)...\n")

  # Find the best route where Neutrophils are RECEIVERS
  comm_matrix <- cc_wang_day3_arg1pos@data$expr_l_r_log2_scale
  all_routes <- colnames(comm_matrix)
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(neutrophil_routes) > 0) {
    route_totals <- colSums(comm_matrix[, neutrophil_routes, drop = FALSE], na.rm = TRUE)
    best_route <- names(sort(route_totals, decreasing = TRUE))[1]
    
    cat("  Using TARGETED route:", best_route, "\n")
    
    # Step 1: Split the route into sender and receiver
    sender <- sub("-.*$", "", best_route)
    receiver <- sub(".*-", "", best_route)
    
    # Step 2: Run LR2TF analysis
    cat("  Step 1: Running LR2TF analysis...\n")
    cc_wang_day3_arg1pos <- LR2TF(
      object = cc_wang_day3_arg1pos,
      sender_cell = sender,
      recevier_cell = receiver,
      slot = "expr_l_r_log2_scale",
      org = "Mus musculus"
    )
    cat("  ✓ LR2TF analysis completed!\n")
    
    # **CRUCIAL FIX: Check if the analysis produced any results**
    cat("  Step 2: Checking for significant L-R-TF connections...\n")
    lr2tf_results <- cc_wang_day3_arg1pos@tf.result$LR2TF.network[[best_route]]
    
    if (!is.null(lr2tf_results) && nrow(lr2tf_results) > 0) {
      cat("  ✓ Found", nrow(lr2tf_results), "connections. Creating plot...\n")
      
      # Now, create the plot
      p_sankey <- LRT.Dimplot(cc_wang_day3_arg1pos)
      print(p_sankey)
      cat("  ✓ Sankey plot completed successfully!\n")
      
    } else {
      # This is the new, informative message if no results are found
      cat("  ✗ Sankey plot failed: No significant L-R-TF connections were found by the LR2TF analysis for this route.\n")
      cat("  💡 TIP: The analysis worked, but the biological signal might be too weak with default filters. Try relaxing the pValue/corValue in the LR2TF function or test a different route (e.g., Macrophage-Neutrophil).\n")
    }
    
  } else {
    cat("  ✗ No neutrophil-receiving routes found for Sankey plot.\n")
  }
  

# 5.8.5. Ridge Plot (Official CellCall Format)
cat("5.8.5. Creating Ridge plot...\n")

  # Get GSEA object
  if (length(cc_wang_day3_arg1pos@data$gsea.list) > 0) {
    cell_type <- names(cc_wang_day3_arg1pos@data$gsea.list)[1]
    egmt <- cc_wang_day3_arg1pos@data$gsea.list[[cell_type]]
    
    # Filter TFs
    egmt.df <- data.frame(egmt)
    flag.index <- which(egmt.df$p.adjust < 0.05)
    
    if (length(flag.index) > 0) {
      p_ridge <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, 
                               core_enrichment = T, orderBy = "NES", decreasing = FALSE)
      print(p_ridge)
      cat("  Ridge plot completed successfully!\n")
    } else {
      cat("  No significant TFs found for Ridge plot\n")
    }
  } else {
    cat("  No GSEA data available for Ridge plot\n")
  }

# 5.8.6. TF Enrichment Plot (Official CellCall Format)
cat("5.8.6. Creating TF enrichment plot...\n")

  if (length(cc_wang_day3_arg1pos@data$gsea.list) > 0) {
    cell_type <- names(cc_wang_day3_arg1pos@data$gsea.list)[1]
    tf_names <- names(cc_wang_day3_arg1pos@data$gsea.list[[cell_type]]@geneSets)
    
    if (length(tf_names) > 0) {
      selected_tfs <- head(tf_names, 3)
      p_tf <- getGSEAplot(gsea.list=cc_wang_day3_arg1pos@data$gsea.list, 
                          geneSetID=selected_tfs, 
                          myCelltype=cell_type, 
                          fc.list=cc_wang_day3_arg1pos@data$fc.list)
      print(p_tf)
      cat("  TF enrichment plot completed successfully!\n")
    } else {
      cat("  No TFs available for enrichment plot\n")
    }
  } else {
    cat("  No GSEA data available for TF enrichment plot\n")
  }

# --- STEP 5.9: Save Results ---
cat("STEP 5.9: Saving results...\n")
saveRDS(cc_wang_day3_arg1pos, "CellCall_Wang_Day3_Arg1pos.rds")
cat("  Analysis saved!\n")

# --- STEP 5.10: Cleanup ---
cat("STEP 5.10: Cleaning up memory after Analysis 5...\n")
rm(wang_day3_subset, wang_day3_arg1pos_data, pw_wang_day3_arg1pos)
gc()
cat("Memory cleaned. cc_wang_day3_arg1pos kept for visualization.\n")

cat("=== ANALYSIS 5 COMPLETED ===\n\n")

# ============================================================================
# ANALYSIS 6: WANG DAY 3 ARG1- SIGNALING
# ============================================================================
# PURPOSE: Analyze cell-cell communication signals received by ARG1- neutrophils on Day 3
# DATASET: Wang dataset, timepoint 3 (3dpi)
# TARGET: Neutrophils with ARG1 expression = 0
# COMPARISON: Cross-dataset validation with Lee Day 3 ARG1- (Analysis 4) and direct with Wang ARG1+ (Analysis 5)
# OUTPUT: 5 CellCall visualizations + RDS file for re-analysis

cat("\n=== ANALYSIS 6: WANG DAY 3 ARG1- SIGNALING ===\n")
cat("Research Question: What signals do ARG1- neutrophils receive on Day 3?\n")

# --- STEP 6.1: Data Preparation ---
cat("STEP 6.1: Data Preparation\n")
wang_day3_subset <- subset(wang_data, time == 3)
cat("Wang Day 3 cells:", ncol(wang_day3_subset), "\n")

# --- STEP 6.2: Prepare ARG1 status ---
cat("STEP 6.2: Preparing cell type labels...\n")
wang_day3_subset$cell_type_clean <- gsub("[-_ ]", "", wang_day3_subset$pruned_labels)
wang_day3_subset$label_arg1 <- wang_day3_subset$cell_type_clean

# Get ARG1 expression for Wang Day 3
wang_day3_arg1_expr <- GetAssayData(wang_day3_subset, assay = "RNA", layer = "data")["Arg1", ]

# Assign ARG1 status to neutrophils
neutrophil_cells <- colnames(wang_day3_subset)[grepl("Neutrophil", wang_day3_subset$cell_type_clean, ignore.case = TRUE)]
cat("Wang Day 3 neutrophils:", length(neutrophil_cells), "\n")

for (cell in neutrophil_cells) {
  arg1_expr <- wang_day3_arg1_expr[cell]
  arg1_status <- ifelse(arg1_expr > 0, "Arg1pos", "Arg1neg")
  wang_day3_subset$label_arg1[colnames(wang_day3_subset) == cell] <- paste0("Neutrophil_", arg1_status)
}

cat("Wang Day 3 ARG1+ neutrophils:", sum(wang_day3_subset$label_arg1 == "Neutrophil_Arg1pos", na.rm=TRUE), "\n")
cat("Wang Day 3 ARG1- neutrophils:", sum(wang_day3_subset$label_arg1 == "Neutrophil_Arg1neg", na.rm=TRUE), "\n")

# --- STEP 6.3: Subset for ARG1- analysis ---
cat("STEP 6.3: Preparing data for ARG1- analysis...\n")
wang_day3_arg1neg_data <- subset(wang_day3_subset, label_arg1 == "Neutrophil_Arg1neg" | !grepl("Neutrophil", label_arg1, ignore.case = TRUE))

# --- STEP 6.4: Fix cell IDs and cell types for CellCall ---
cat("STEP 6.4: Fixing cell IDs and cell types for CellCall compatibility...\n")
new_cell_ids <- gsub("-", "_", colnames(wang_day3_arg1neg_data))
wang_day3_arg1neg_data <- RenameCells(wang_day3_arg1neg_data, new.names = new_cell_ids)
wang_day3_arg1neg_data$label_arg1 <- gsub("-", "", wang_day3_arg1neg_data$label_arg1)

# --- STEP 6.5: Memory optimization ---
rm(wang_day3_subset)
gc()

cat("Wang Day 3 ARG1- analysis cells:", ncol(wang_day3_arg1neg_data), "\n")
cat("Cell types:", unique(wang_day3_arg1neg_data$label_arg1), "\n")

if (ncol(wang_day3_arg1neg_data) < 100) {
  cat("WARNING: Very few cells for analysis. Consider adjusting parameters.\n")
}

# --- STEP 6.6: Create CellCall object ---
cat("STEP 6.6: Creating CellCall object for Wang Day 3 ARG1-...\n")
gc()

cc_wang_day3_arg1neg <- CreateObject_fromSeurat(
  Seurat.object = wang_day3_arg1neg_data,
  slot = "counts",
  cell_type = "label_arg1",
  Org = "Mus musculus"
)
cat("CellCall object created successfully!\n")

rm(wang_day3_arg1neg_data)
gc()

# --- STEP 6.7: Run TransCommuProfile ---
cat("STEP 6.7: Running TransCommuProfile for Wang Day 3 ARG1-...\n")
gc()

cc_wang_day3_arg1neg <- TransCommuProfile(
  object = cc_wang_day3_arg1neg,
  Org = "Mus musculus",
  pValueCor = 0.1,
  CorValue = 0.15,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  IS_core = TRUE
)
cat("TransCommuProfile completed successfully!\n")

gc()

# --- STEP 6.8: Create Visualizations ---
cat("STEP 6.8: Creating visualizations for Wang Day 3 ARG1-...\n")

# 6.8.1. Circular Plot (Official CellCall Format) - FIXED
cat("6.8.1. Creating ViewInterCircos (FIXED)...\n")

  # Identify target routes for the plot
  all_routes <- colnames(cc_wang_day3_arg1neg@data$expr_l_r_log2_scale)
  target_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(target_routes) > 0) {
    # Get cell types for color mapping
    cell_types_in_obj <- unique(sapply(strsplit(target_routes, "-"), `[`, 1))
    cell_types_in_obj <- unique(c(cell_types_in_obj, sapply(strsplit(target_routes, "-"), `[`, 2)))
    cell_types_in_obj <- cell_types_in_obj[!is.na(cell_types_in_obj) & cell_types_in_obj != ""]
    
    if (length(cell_types_in_obj) >= 2) {
      colors_for_circos <- grDevices::rainbow(length(cell_types_in_obj))
      cell_color_df <- data.frame(color = colors_for_circos, row.names = cell_types_in_obj, check.names = FALSE)
      circlize::circos.clear()
      grid::grid.newpage()
      
      # Use correct parameters from official documentation
      p_circos <- ViewInterCircos(
        object = cc_wang_day3_arg1neg,
        font = 2,
        cellColor = cell_color_df,
        lrColor = c("#F16B6F", "#84B1ED"),
        arr.type = "big.arrow",
        arr.length = 0.04,
        trackhight1 = 0.05,
        slot = "expr_l_r_log2_scale",
        linkcolor.from.sender = TRUE,
        gap.degree = 0.5,  # Much smaller gap for targeted analysis
        order.vector = cell_types_in_obj,
        trackhight2 = 0.032,
        track.margin2 = c(0.01, 0.12),
        DIY = FALSE
      )
      print(p_circos)
      cat("  ✓ FIXED ViewInterCircos completed successfully!\n")
    } else {
      cat("  ✗ Not enough cell types for ViewInterCircos\n")
    }
  } else {
    cat("  ✗ No neutrophil-receiving routes found for ViewInterCircos\n")
  }

# 6.8.2. Heatmap (Official CellCall Format)
cat("6.8.2. Creating viewPheatmap...\n")

  p_heatmap <- viewPheatmap(
    object = cc_wang_day3_arg1neg,
    slot = "expr_l_r_log2_scale",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    angle_col = 45,
    main = "Wang Day 3 ARG1- Neutrophil Signaling"
  )
  print(p_heatmap)
  cat("  viewPheatmap completed successfully!\n")

# 6.8.3. Pathway Analysis (Official CellCall Format - TARGETED)
cat("6.8.3. Running TARGETED pathway analysis (Neutrophils as receivers)...\n")

  # Get all routes for pathway analysis
  all_routes <- colnames(cc_wang_day3_arg1neg@data$expr_l_r_log2_scale)
  
  # Filter for routes where Neutrophils are the RECEIVER
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  cat("  Found", length(neutrophil_routes), "routes where Neutrophils are receivers\n")
  
  if (length(neutrophil_routes) > 0) {
    # Run pathway analysis ONLY on targeted routes
    pathway.hyper.list <- lapply(neutrophil_routes, function(route){
      cat("  Processing TARGETED route:", route, "\n")
      tmp <- getHyperPathway(data = cc_wang_day3_arg1neg@data$expr_l_r_log2_scale, 
                            object = cc_wang_day3_arg1neg, cella_cellb = route, Org="Mus musculus")
      return(tmp)
    })
    
    # Filter out NULL results
    pathway.hyper.list <- pathway.hyper.list[!sapply(pathway.hyper.list, is.null)]
    pathway.hyper.list <- pathway.hyper.list[sapply(pathway.hyper.list, nrow) > 0]
    
    if (length(pathway.hyper.list) > 0) {
      myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names(pathway.hyper.list))
      p_bubble <- plotBubble(myPub.df) + 
        ggtitle("Signals Received by ARG1- Neutrophils: Pathway Analysis")
      print(p_bubble)
      cat("  TARGETED pathway analysis completed successfully!\n")
    } else {
      cat("  No significant pathways were found for targeted analysis.\n")
    }
  } else {
    cat("  No neutrophil-receiving routes found for pathway analysis.\n")
  }

# 6.8.4. Sankey Plot (Official CellCall Format) - FIXED
cat("6.8.4. Creating Sankey plot (FINAL FIX)...\n")

  # Find the best route where Neutrophils are RECEIVERS
  comm_matrix <- cc_wang_day3_arg1neg@data$expr_l_r_log2_scale
  all_routes <- colnames(comm_matrix)
  neutrophil_routes <- grep(".*-Neutrophil.*", all_routes, value = TRUE)
  
  if (length(neutrophil_routes) > 0) {
    route_totals <- colSums(comm_matrix[, neutrophil_routes, drop = FALSE], na.rm = TRUE)
    best_route <- names(sort(route_totals, decreasing = TRUE))[1]
    
    cat("  Using TARGETED route:", best_route, "\n")
    
    # Step 1: Split the route into sender and receiver
    sender <- sub("-.*$", "", best_route)
    receiver <- sub(".*-", "", best_route)
    
    # Step 2: Run LR2TF analysis
    cat("  Step 1: Running LR2TF analysis...\n")
    cc_wang_day3_arg1neg <- LR2TF(
      object = cc_wang_day3_arg1neg,
      sender_cell = sender,
      recevier_cell = receiver,
      slot = "expr_l_r_log2_scale",
      org = "Mus musculus"
    )
    cat("  ✓ LR2TF analysis completed!\n")
    
    # **CRUCIAL FIX: Check if the analysis produced any results**
    cat("  Step 2: Checking for significant L-R-TF connections...\n")
    lr2tf_results <- cc_wang_day3_arg1neg@tf.result$LR2TF.network[[best_route]]
    
    if (!is.null(lr2tf_results) && nrow(lr2tf_results) > 0) {
      cat("  ✓ Found", nrow(lr2tf_results), "connections. Creating plot...\n")
      
      # Now, create the plot
      p_sankey <- LRT.Dimplot(cc_wang_day3_arg1neg)
      print(p_sankey)
      cat("  ✓ Sankey plot completed successfully!\n")
      
    } else {
      # This is the new, informative message if no results are found
      cat("  ✗ Sankey plot failed: No significant L-R-TF connections were found by the LR2TF analysis for this route.\n")
      cat("  💡 TIP: The analysis worked, but the biological signal might be too weak with default filters. Try relaxing the pValue/corValue in the LR2TF function or test a different route (e.g., Macrophage-Neutrophil).\n")
    }
    
  } else {
    cat("  ✗ No neutrophil-receiving routes found for Sankey plot.\n")
  }
  

# 6.8.5. Ridge Plot (Official CellCall Format)
cat("6.8.5. Creating Ridge plot...\n")

  # Get GSEA object
  if (length(cc_wang_day3_arg1neg@data$gsea.list) > 0) {
    cell_type <- names(cc_wang_day3_arg1neg@data$gsea.list)[1]
    egmt <- cc_wang_day3_arg1neg@data$gsea.list[[cell_type]]
    
    # Filter TFs
    egmt.df <- data.frame(egmt)
    flag.index <- which(egmt.df$p.adjust < 0.05)
    
    if (length(flag.index) > 0) {
      p_ridge <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, 
                               core_enrichment = T, orderBy = "NES", decreasing = FALSE)
      print(p_ridge)
      cat("  Ridge plot completed successfully!\n")
    } else {
      cat("  No significant TFs found for Ridge plot\n")
    }
  } else {
    cat("  No GSEA data available for Ridge plot\n")
  }

# 6.8.6. TF Enrichment Plot (Official CellCall Format)
cat("6.8.6. Creating TF enrichment plot...\n")

  if (length(cc_wang_day3_arg1neg@data$gsea.list) > 0) {
    cell_type <- names(cc_wang_day3_arg1neg@data$gsea.list)[1]
    tf_names <- names(cc_wang_day3_arg1neg@data$gsea.list[[cell_type]]@geneSets)
    
    if (length(tf_names) > 0) {
      selected_tfs <- head(tf_names, 3)
      p_tf <- getGSEAplot(gsea.list=cc_wang_day3_arg1neg@data$gsea.list, 
                          geneSetID=selected_tfs, 
                          myCelltype=cell_type, 
                          fc.list=cc_wang_day3_arg1neg@data$fc.list)
      print(p_tf)
      cat("  TF enrichment plot completed successfully!\n")
    } else {
      cat("  No TFs available for enrichment plot\n")
    }
  } else {
    cat("  No GSEA data available for TF enrichment plot\n")
  }

# --- STEP 6.9: Save Results ---
cat("STEP 6.9: Saving results...\n")
saveRDS(cc_wang_day3_arg1neg, "CellCall_Wang_Day3_Arg1neg.rds")
cat("  Analysis saved!\n")

# --- STEP 6.10: Cleanup ---
cat("STEP 6.10: Cleaning up memory after Analysis 6...\n")
rm(wang_day3_subset, wang_day3_arg1neg_data, pw_wang_day3_arg1neg)
gc()
cat("Memory cleaned. cc_wang_day3_arg1neg kept for visualization.\n")

cat("=== ANALYSIS 6 COMPLETED ===\n\n")

# All analyses completed successfully!
# Individual RDS files have been saved for each analysis:
# - CellCall_Lee_Day1_Arg1pos.rds
# - CellCall_Lee_Day1_Arg1neg.rds  
# - CellCall_Lee_Day3_Arg1pos.rds
# - CellCall_Lee_Day3_Arg1neg.rds
# - CellCall_Wang_Day3_Arg1pos.rds
# - CellCall_Wang_Day3_Arg1neg.rds

# Final cleanup
gc()
