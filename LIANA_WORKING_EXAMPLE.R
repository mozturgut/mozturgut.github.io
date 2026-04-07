#!/usr/bin/env Rscript

# ============================================================================
# LIANA SCE FIX - Complete Working Example
# ============================================================================
# This script demonstrates the correct way to create SCE objects for LIANA
# that survive LIANA's internal filtering
#
# Based on actual error from Wang Day 3 analysis showing:
#   Error 1: 'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'
#   Error 2: 'nrow' of 'int_colData' not equal to 'ncol(object)'
# ============================================================================

# ============================================================================
# APPROACH 1: SIMPLE (RECOMMENDED) - Let LIANA Do All Filtering
# ============================================================================

approach_1_simple <- function() {
  cat("\n=== APPROACH 1: Simple - Let LIANA Filter ===\n")
  
  # Assuming wang_day3 is your Seurat object and wang_day3_liana_labels are your labels
  
  # Step 1: Extract matrices (no filtering!)
  SeuratObject::DefaultAssay(wang_day3) <- "RNA"
  if (inherits(wang_day3[["RNA"]], "Assay5")) {
    wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
  }
  
  wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
  wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
  
  # Step 2: Verify dimensions match between counts and data
  stopifnot(
    "Counts and data must have same dimensions" = 
      nrow(wang_day3_cnt) == nrow(wang_day3_dat) && 
      ncol(wang_day3_cnt) == ncol(wang_day3_dat)
  )
  
  # Step 3: Verify labels length matches cell count
  stopifnot(
    "Labels must match number of cells" = 
      ncol(wang_day3_cnt) == length(wang_day3_liana_labels)
  )
  
  cat("  Dimensions verified:\n")
  cat("    Genes:", nrow(wang_day3_cnt), "\n")
  cat("    Cells:", ncol(wang_day3_cnt), "\n")
  cat("    Labels:", length(wang_day3_liana_labels), "\n")
  
  # Step 4: Create SCE with colData passed during construction
  wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = wang_day3_cnt,
      logcounts = wang_day3_dat
    ),
    colData = S4Vectors::DataFrame(
      label = factor(wang_day3_liana_labels)
    )
  )
  
  # Step 5: Set colLabels (this ensures proper internal structure)
  SingleCellExperiment::colLabels(wang_day3_sce) <- factor(wang_day3_liana_labels)
  
  # Step 6: Verify SCE is valid
  stopifnot("SCE must be valid" = validObject(wang_day3_sce))
  
  cat("  ✓ SCE created and validated\n")
  cat("    SCE genes:", nrow(wang_day3_sce), "\n")
  cat("    SCE cells:", ncol(wang_day3_sce), "\n")
  cat("    colData rows:", nrow(SingleCellExperiment::colData(wang_day3_sce)), "\n")
  cat("    colLabels length:", length(SingleCellExperiment::colLabels(wang_day3_sce)), "\n")
  
  # Step 7: Run LIANA (it will do the filtering internally)
  wang_day3_liana_args <- list(
    sce = wang_day3_sce,
    resource = "Consensus",  # or wang_day3_liana_resource
    expr_prop = 0.05,
    verbose = TRUE,
    min_cells = 0
  )
  
  # Add external_resource if using custom
  if (!is.null(wang_day3_liana_external)) {
    wang_day3_liana_args$external_resource <- wang_day3_liana_external
    wang_day3_liana_args$resource <- "custom"
  }
  
  cat("  Running LIANA...\n")
  wang_day3_liana_result <- suppressWarnings(
    do.call(liana::liana_wrap, wang_day3_liana_args)
  )
  
  cat("  ✓ LIANA completed successfully\n")
  
  return(wang_day3_liana_result)
}

# ============================================================================
# APPROACH 2: COMPLEX - Pre-filter Correctly with Index Tracking
# ============================================================================

approach_2_prefilter <- function() {
  cat("\n=== APPROACH 2: Pre-filter with Index Tracking ===\n")
  
  # Step 1: Extract matrices
  SeuratObject::DefaultAssay(wang_day3) <- "RNA"
  if (inherits(wang_day3[["RNA"]], "Assay5")) {
    wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
  }
  
  wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
  wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
  
  # Ensure sparse matrix format
  if (!inherits(wang_day3_cnt, "matrix") && !inherits(wang_day3_cnt, "Matrix")) {
    wang_day3_cnt <- as(wang_day3_cnt, "sparseMatrix")
  }
  if (!inherits(wang_day3_dat, "matrix") && !inherits(wang_day3_dat, "Matrix")) {
    wang_day3_dat <- as(wang_day3_dat, "sparseMatrix")
  }
  
  cat("  Starting dimensions:\n")
  cat("    Genes:", nrow(wang_day3_cnt), "\n")
  cat("    Cells:", ncol(wang_day3_cnt), "\n")
  
  # Step 2: CRITICAL - Track cell indices from the start
  wang_day3_original_ncells <- ncol(wang_day3_cnt)
  wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)
  
  cat("  Cell indices initialized: 1 to", wang_day3_original_ncells, "\n")
  
  # Step 3: Filter genes (optional - to LR entity genes)
  if (exists("wang_day3_entity_genes") && length(wang_day3_entity_genes) > 0) {
    wang_day3_genes_in_data <- rownames(wang_day3_cnt)
    wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
    
    if (sum(wang_day3_keep_genes) < 3) {
      cat("  Warning: Less than 3 entity genes found, keeping all genes\n")
      wang_day3_keep_genes <- rep(TRUE, length(wang_day3_genes_in_data))
    }
    
    wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
    wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
    
    cat("  After gene filtering:", nrow(wang_day3_cnt), "genes\n")
  }
  
  # Step 4: Filter cells with zero counts
  wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
  wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
  wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]
  
  # Step 5: CRITICAL - Update cell indices after filtering
  wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]
  
  cat("  After cell filtering:", ncol(wang_day3_cnt), "cells\n")
  cat("  Cell indices now:", length(wang_day3_cell_indices), "tracked\n")
  
  # Step 6: Filter genes with zero counts
  wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
  wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
  wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
  
  cat("  After gene re-filtering:", nrow(wang_day3_cnt), "genes\n")
  
  # Step 7: CRITICAL - Create labels using tracked indices (not boolean filter!)
  wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])
  
  cat("  Labels created:", length(wang_day3_lab), "\n")
  
  # Step 8: Verify dimensions match
  stopifnot(
    "Counts and data must have same dimensions" = 
      nrow(wang_day3_cnt) == nrow(wang_day3_dat) && 
      ncol(wang_day3_cnt) == ncol(wang_day3_dat)
  )
  stopifnot(
    "Labels must match number of cells" = 
      ncol(wang_day3_cnt) == length(wang_day3_lab)
  )
  
  cat("  ✓ Dimension checks passed\n")
  
  # Step 9: Create SCE with properly aligned colData
  wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = wang_day3_cnt,
      logcounts = wang_day3_dat
    ),
    colData = S4Vectors::DataFrame(
      label = wang_day3_lab
    )
  )
  
  # Step 10: Set colLabels
  SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab
  
  # Step 11: Verify SCE is valid
  stopifnot("SCE must be valid" = validObject(wang_day3_sce))
  
  cat("  ✓ SCE created and validated\n")
  cat("    Final SCE dimensions:\n")
  cat("      Genes:", nrow(wang_day3_sce), "\n")
  cat("      Cells:", ncol(wang_day3_sce), "\n")
  cat("      colData rows:", nrow(SingleCellExperiment::colData(wang_day3_sce)), "\n")
  cat("      colLabels length:", length(SingleCellExperiment::colLabels(wang_day3_sce)), "\n")
  
  # Step 12: Run LIANA
  wang_day3_liana_args <- list(
    sce = wang_day3_sce,
    resource = "Consensus",
    expr_prop = 0.05,
    verbose = TRUE,
    min_cells = 0
  )
  
  if (!is.null(wang_day3_liana_external)) {
    wang_day3_liana_args$external_resource <- wang_day3_liana_external
    wang_day3_liana_args$resource <- "custom"
  }
  
  cat("  Running LIANA...\n")
  wang_day3_liana_result <- suppressWarnings(
    do.call(liana::liana_wrap, wang_day3_liana_args)
  )
  
  cat("  ✓ LIANA completed successfully\n")
  
  return(wang_day3_liana_result)
}

# ============================================================================
# USAGE INSTRUCTIONS
# ============================================================================

cat("\n=== LIANA SCE Fix - Working Example ===\n")
cat("\nThis script provides two approaches:\n")
cat("\n1. SIMPLE (Recommended):")
cat("\n   - Let LIANA do all filtering internally")
cat("\n   - Just create SCE with all cells/genes and proper labels")
cat("\n   - Minimal code, less error-prone\n")
cat("\n2. COMPLEX (Advanced):")
cat("\n   - Pre-filter genes/cells before LIANA")
cat("\n   - Track cell indices through all filtering steps")
cat("\n   - More control but requires careful implementation\n")
cat("\nBoth approaches will work if implemented correctly.\n")
cat("\nKey principle: ALWAYS ensure ncol(matrix) == length(labels)\n")

# ============================================================================
# COMPLETE EXAMPLE FOR WANG DAY 3 (Copy this into your script)
# ============================================================================

# Uncomment and use this in your actual script:

# # Wang Day 3 LIANA Analysis (SIMPLE APPROACH - RECOMMENDED)
# cat("--- LIANA Analysis: Wang Day 3 ---\n")
# 
# # Create LIANA labels (Arg1+ vs Arg1- for neutrophils)
# wang_day3_liana_labels <- wang_day3_pruned_labels
# wang_day3_liana_labels[wang_day3_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
# wang_day3_liana_labels[wang_day3_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
# Seurat::Idents(wang_day3) <- factor(wang_day3_liana_labels)
# 
# # Extract matrices
# SeuratObject::DefaultAssay(wang_day3) <- "RNA"
# if (inherits(wang_day3[["RNA"]], "Assay5")) {
#   wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
# }
# wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
# wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
# 
# # Verify dimensions
# stopifnot(ncol(wang_day3_cnt) == length(wang_day3_liana_labels))
# 
# # Create SCE
# wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
#   assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
#   colData = S4Vectors::DataFrame(label = factor(wang_day3_liana_labels))
# )
# SingleCellExperiment::colLabels(wang_day3_sce) <- factor(wang_day3_liana_labels)
# 
# # Verify SCE is valid
# stopifnot(validObject(wang_day3_sce))
# 
# # Run LIANA
# wang_day3_liana_result <- liana::liana_wrap(
#   sce = wang_day3_sce,
#   resource = wang_day3_liana_resource,
#   expr_prop = 0.05,
#   verbose = TRUE,
#   min_cells = 0
# )
