#!/usr/bin/env Rscript

# ============================================================================
# LIANA ANALYSIS - FIXED VERSION
# Fix for SingleCellExperiment dimension mismatch error
# 
# ERROR FIXED: 'nrow' of 'int_colData' not equal to 'ncol(object)'
# ROOT CAUSE: Metadata labels not properly aligned with filtered expression matrices
# SOLUTION: Track cell indices through all filtering steps
# ============================================================================

# This script provides the corrected SCE building pattern for LIANA analysis
# Apply this pattern to all dataset sections: Lee Day 1/3, Wang Day 3, Qin Day 1/3

# ============================================================================
# CORRECTED PATTERN: Build SingleCellExperiment for LIANA
# ============================================================================

build_liana_sce_fixed <- function(seurat_obj, liana_labels, entity_genes, dataset_name = "Dataset") {
  # This is a reference implementation showing the correct pattern
  # The actual code in your script should follow this pattern inline (no helper functions)
  
  cat("Building SCE for", dataset_name, "...\n")
  
  # Ensure RNA assay is active and layers are joined
  SeuratObject::DefaultAssay(seurat_obj) <- "RNA"
  if (inherits(seurat_obj[["RNA"]], "Assay5")) {
    seurat_obj[["RNA"]] <- SeuratObject::JoinLayers(seurat_obj[["RNA"]])
  }
  
  # Extract count and data matrices
  cnt <- SeuratObject::GetAssayData(seurat_obj, layer = "counts")
  dat <- SeuratObject::GetAssayData(seurat_obj, layer = "data")
  if (is.null(cnt) || nrow(cnt) == 0 || ncol(cnt) == 0) cnt <- dat
  if (!inherits(cnt, "matrix") && !inherits(cnt, "Matrix")) cnt <- as(cnt, "sparseMatrix")
  if (!inherits(dat, "matrix") && !inherits(dat, "Matrix")) dat <- as(dat, "sparseMatrix")
  
  # CRITICAL: Store original cell indices before any filtering
  original_ncells <- ncol(cnt)
  cell_indices <- seq_len(original_ncells)
  
  cat("  Original cells:", original_ncells, "\n")
  
  # Step 1: Filter genes to keep only LR entity genes
  genes_in_data <- rownames(cnt)
  keep_genes <- genes_in_data %in% entity_genes
  if (sum(keep_genes) < 3) {
    cat("  Warning: Less than 3 entity genes found, keeping all genes\n")
    keep_genes <- rep(TRUE, length(genes_in_data))
  }
  cnt <- cnt[keep_genes, , drop = FALSE]
  dat <- dat[keep_genes, , drop = FALSE]
  
  cat("  Genes after entity filtering:", nrow(cnt), "\n")
  
  # Step 2: Filter cells based on nonzero counts
  nonzero_cells <- Matrix::colSums(cnt) > 0
  cnt <- cnt[, nonzero_cells, drop = FALSE]
  dat <- dat[, nonzero_cells, drop = FALSE]
  
  # CRITICAL FIX: Update cell_indices to track which cells remain
  cell_indices <- cell_indices[nonzero_cells]
  
  cat("  Cells after nonzero filtering:", ncol(cnt), "\n")
  
  # Step 3: Filter genes again to remove any that became all-zero
  nonzero_genes <- Matrix::rowSums(cnt) > 0
  cnt <- cnt[nonzero_genes, , drop = FALSE]
  dat <- dat[nonzero_genes, , drop = FALSE]
  
  cat("  Genes after second filtering:", nrow(cnt), "\n")
  
  # CRITICAL FIX: Create labels using tracked cell indices
  # This ensures alignment between filtered matrices and metadata
  lab <- factor(liana_labels[cell_indices])
  
  # Verify dimensions match (safety check)
  if (ncol(cnt) != ncol(dat)) {
    stop("Dimension mismatch: cnt and dat have different number of columns")
  }
  if (ncol(cnt) != length(lab)) {
    stop("Dimension mismatch: matrices have ", ncol(cnt), " cells but labels have ", length(lab))
  }
  
  cat("  Final dimensions - Genes:", nrow(cnt), ", Cells:", ncol(cnt), ", Labels:", length(lab), "\n")
  cat("  Label distribution:\n")
  print(table(lab))
  
  # Create SingleCellExperiment with properly aligned data
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = cnt, logcounts = dat),
    colData = S4Vectors::DataFrame(label = lab)
  )
  SingleCellExperiment::colLabels(sce) <- lab
  
  cat("  ✓ SCE created successfully\n\n")
  
  return(sce)
}

# ============================================================================
# INLINE PATTERN FOR ACTUAL SCRIPT (NO HELPER FUNCTIONS)
# ============================================================================

# For use in the actual analysis script, copy the pattern below and adapt variable names
# This shows the pattern for Wang Day 3, adapt for other datasets

# --- EXAMPLE: Wang Day 3 SCE Building (Corrected) ---
# 
# SeuratObject::DefaultAssay(wang_day3) <- "RNA"
# if (inherits(wang_day3[["RNA"]], "Assay5")) wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
# wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
# wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
# if (is.null(wang_day3_cnt) || nrow(wang_day3_cnt) == 0 || ncol(wang_day3_cnt) == 0) wang_day3_cnt <- wang_day3_dat
# if (!inherits(wang_day3_cnt, "matrix") && !inherits(wang_day3_cnt, "Matrix")) wang_day3_cnt <- as(wang_day3_cnt, "sparseMatrix")
# if (!inherits(wang_day3_dat, "matrix") && !inherits(wang_day3_dat, "Matrix")) wang_day3_dat <- as(wang_day3_dat, "sparseMatrix")
# 
# # CRITICAL: Store original cell indices before filtering
# wang_day3_original_ncells <- ncol(wang_day3_cnt)
# wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)
# 
# # Filter genes to keep only LR entity genes
# wang_day3_genes_in_data <- rownames(wang_day3_cnt)
# wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
# if (sum(wang_day3_keep_genes) < 3) wang_day3_keep_genes <- rep(TRUE, length(wang_day3_genes_in_data))
# wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
# wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
# 
# # Filter cells based on nonzero counts
# wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
# wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
# wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]
# 
# # CRITICAL FIX: Update cell indices
# wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]
# 
# # Filter genes again to remove zeros
# wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
# wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
# wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
# 
# # CRITICAL FIX: Create labels using tracked indices
# wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])
# 
# # Verify dimensions
# stopifnot(ncol(wang_day3_cnt) == ncol(wang_day3_dat))
# stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))
# 
# # Create SCE with aligned data
# wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
#   assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
#   colData = S4Vectors::DataFrame(label = wang_day3_lab)
# )
# SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab

# ============================================================================
# SUMMARY OF CHANGES
# ============================================================================

# BEFORE (BUGGY):
# 1. Filter cells: nonzero_cells <- colSums(cnt) > 0
# 2. Apply filter: cnt <- cnt[, nonzero_cells]
# 3. Create labels: lab <- factor(liana_labels[nonzero_cells])  # BUG! Wrong indices!
# 
# AFTER (FIXED):
# 1. Track indices: cell_indices <- seq_len(ncol(cnt))
# 2. Filter cells: nonzero_cells <- colSums(cnt) > 0
# 3. Apply filter: cnt <- cnt[, nonzero_cells]
# 4. Update indices: cell_indices <- cell_indices[nonzero_cells]  # KEY FIX!
# 5. Create labels: lab <- factor(liana_labels[cell_indices])  # Correct alignment!

cat("LIANA SCE fix patterns loaded.\n")
cat("Apply the inline pattern to all dataset sections in your analysis script.\n")
cat("See LIANA_SCE_FIX.md for detailed explanation.\n")
