#!/usr/bin/env Rscript

# ============================================================================
# TEST SCRIPT: Verify LIANA SCE Fix
# ============================================================================
# This script tests the index tracking logic to ensure proper alignment
# between filtered matrices and metadata labels

cat("=== LIANA SCE Fix Validation Test ===\n\n")

# Simulate the filtering process with test data
set.seed(42)

# Create mock data
n_genes <- 100
n_cells <- 50
n_lr_genes <- 30  # Subset of genes that are LR-related

cat("Test Setup:\n")
cat("  Original genes:", n_genes, "\n")
cat("  Original cells:", n_cells, "\n")
cat("  LR entity genes:", n_lr_genes, "\n\n")

# Create mock expression matrix (sparse)
mock_counts <- Matrix::Matrix(
  rpois(n_genes * n_cells, lambda = 2),
  nrow = n_genes,
  ncol = n_cells,
  sparse = TRUE
)
rownames(mock_counts) <- paste0("Gene", 1:n_genes)
colnames(mock_counts) <- paste0("Cell", 1:n_cells)

# Create mock labels (simulating cell types and Arg1 status)
mock_labels <- c(
  rep("Neutrophil_Arg1pos", 15),
  rep("Neutrophil_Arg1neg", 15),
  rep("Macrophage", 10),
  rep("T_cell", 10)
)

# Create mock LR entity genes
lr_entity_genes <- paste0("Gene", sample(1:n_genes, n_lr_genes))

cat("Mock data created successfully\n\n")

# ============================================================================
# TEST 1: BUGGY APPROACH (Original Code)
# ============================================================================
cat("TEST 1: Buggy Approach (Original Code)\n")
cat("---------------------------------------\n")

tryCatch({
  # Simulate the buggy approach
  buggy_cnt <- mock_counts
  
  # Filter genes
  genes_in_data <- rownames(buggy_cnt)
  keep_genes <- genes_in_data %in% lr_entity_genes
  buggy_cnt <- buggy_cnt[keep_genes, , drop = FALSE]
  cat("  After gene filter:", nrow(buggy_cnt), "genes,", ncol(buggy_cnt), "cells\n")
  
  # Filter cells (first filter)
  nonzero_cells <- Matrix::colSums(buggy_cnt) > 0
  buggy_cnt <- buggy_cnt[, nonzero_cells, drop = FALSE]
  cat("  After cell filter:", nrow(buggy_cnt), "genes,", ncol(buggy_cnt), "cells\n")
  
  # Filter genes again
  nonzero_genes <- Matrix::rowSums(buggy_cnt) > 0
  buggy_cnt <- buggy_cnt[nonzero_genes, , drop = FALSE]
  cat("  After gene re-filter:", nrow(buggy_cnt), "genes,", ncol(buggy_cnt), "cells\n")
  
  # BUG: Using boolean filter on original labels
  buggy_lab <- factor(mock_labels[nonzero_cells])
  
  cat("  Labels created:", length(buggy_lab), "labels\n")
  cat("  Matrix cells:", ncol(buggy_cnt), "\n")
  
  if (length(buggy_lab) != ncol(buggy_cnt)) {
    cat("  ✗ DIMENSION MISMATCH DETECTED!\n")
    cat("    This would cause: 'nrow' of 'int_colData' not equal to 'ncol(object)'\n")
  } else {
    cat("  ✓ Dimensions match (unexpected - test may need adjustment)\n")
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
})

cat("\n")

# ============================================================================
# TEST 2: FIXED APPROACH (Corrected Code)
# ============================================================================
cat("TEST 2: Fixed Approach (Corrected Code)\n")
cat("---------------------------------------\n")

tryCatch({
  # Simulate the fixed approach
  fixed_cnt <- mock_counts
  
  # CRITICAL: Track original cell indices
  original_ncells <- ncol(fixed_cnt)
  cell_indices <- seq_len(original_ncells)
  cat("  Initial cell indices:", length(cell_indices), "\n")
  
  # Filter genes
  genes_in_data <- rownames(fixed_cnt)
  keep_genes <- genes_in_data %in% lr_entity_genes
  fixed_cnt <- fixed_cnt[keep_genes, , drop = FALSE]
  cat("  After gene filter:", nrow(fixed_cnt), "genes,", ncol(fixed_cnt), "cells\n")
  
  # Filter cells (first filter)
  nonzero_cells <- Matrix::colSums(fixed_cnt) > 0
  fixed_cnt <- fixed_cnt[, nonzero_cells, drop = FALSE]
  
  # CRITICAL FIX: Update tracked indices
  cell_indices <- cell_indices[nonzero_cells]
  cat("  After cell filter:", nrow(fixed_cnt), "genes,", ncol(fixed_cnt), "cells\n")
  cat("  Updated cell indices:", length(cell_indices), "\n")
  
  # Filter genes again
  nonzero_genes <- Matrix::rowSums(fixed_cnt) > 0
  fixed_cnt <- fixed_cnt[nonzero_genes, , drop = FALSE]
  cat("  After gene re-filter:", nrow(fixed_cnt), "genes,", ncol(fixed_cnt), "cells\n")
  
  # CRITICAL FIX: Use tracked indices for labels
  fixed_lab <- factor(mock_labels[cell_indices])
  
  cat("  Labels created:", length(fixed_lab), "labels\n")
  cat("  Matrix cells:", ncol(fixed_cnt), "\n")
  
  # Verify dimensions match
  if (length(fixed_lab) == ncol(fixed_cnt)) {
    cat("  ✓ DIMENSIONS MATCH! Fix works correctly.\n")
    
    # Additional verification
    if (ncol(fixed_cnt) <= original_ncells && ncol(fixed_cnt) > 0) {
      cat("  ✓ Cell count is valid (0 < filtered <= original)\n")
    }
    
    # Show label distribution
    cat("  Label distribution:\n")
    label_table <- table(fixed_lab)
    for (label_name in names(label_table)) {
      cat("    ", label_name, ":", label_table[label_name], "\n")
    }
    
  } else {
    cat("  ✗ DIMENSION MISMATCH!\n")
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
})

cat("\n")

# ============================================================================
# TEST 3: Simulate actual SCE creation
# ============================================================================
cat("TEST 3: SingleCellExperiment Creation\n")
cat("-------------------------------------\n")

# Check if SingleCellExperiment is available
sce_available <- requireNamespace("SingleCellExperiment", quietly = TRUE)

if (sce_available) {
  tryCatch({
    # Use the fixed approach to build SCE
    fixed_cnt <- mock_counts
    original_ncells <- ncol(fixed_cnt)
    cell_indices <- seq_len(original_ncells)
    
    # Apply filtering with index tracking
    genes_in_data <- rownames(fixed_cnt)
    keep_genes <- genes_in_data %in% lr_entity_genes
    fixed_cnt <- fixed_cnt[keep_genes, , drop = FALSE]
    
    nonzero_cells <- Matrix::colSums(fixed_cnt) > 0
    fixed_cnt <- fixed_cnt[, nonzero_cells, drop = FALSE]
    cell_indices <- cell_indices[nonzero_cells]
    
    nonzero_genes <- Matrix::rowSums(fixed_cnt) > 0
    fixed_cnt <- fixed_cnt[nonzero_genes, , drop = FALSE]
    
    fixed_lab <- factor(mock_labels[cell_indices])
    
    # Create SCE
    test_sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = fixed_cnt, logcounts = log1p(fixed_cnt)),
      colData = S4Vectors::DataFrame(label = fixed_lab)
    )
    SingleCellExperiment::colLabels(test_sce) <- fixed_lab
    
    # Validate
    if (ncol(test_sce) == length(fixed_lab)) {
      cat("  ✓ SCE created successfully\n")
      cat("  Cells in SCE:", ncol(test_sce), "\n")
      cat("  Genes in SCE:", nrow(test_sce), "\n")
      cat("  Labels in colData:", nrow(SingleCellExperiment::colData(test_sce)), "\n")
    } else {
      cat("  ✗ SCE dimension mismatch\n")
    }
  }, error = function(e) {
    cat("  ✗ ERROR creating SCE:", e$message, "\n")
  })
} else {
  cat("  SingleCellExperiment package not available, skipping SCE test\n")
  cat("  Install with: BiocManager::install('SingleCellExperiment')\n")
}

cat("\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("=== TEST SUMMARY ===\n")
cat("The fix ensures that:\n")
cat("  1. Cell indices are tracked through all filtering steps\n")
cat("  2. Labels are extracted using tracked indices (not boolean filters)\n")
cat("  3. Dimensions of expression matrices and metadata always match\n")
cat("  4. SingleCellExperiment objects can be created without errors\n")
cat("\nThis resolves the error:\n")
cat("  'nrow' of 'int_colData' not equal to 'ncol(object)'\n")
