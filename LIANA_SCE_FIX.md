# Fix for LIANA SingleCellExperiment Dimension Mismatch Error

## Problem
Error: `'nrow' of 'int_colData' not equal to 'ncol(object)'`

This occurs when creating SingleCellExperiment objects for LIANA analysis because the metadata (colData) doesn't align with the filtered expression matrices.

## Root Cause
In the Wang Day 3 section (and similar sections), the code:
1. Filters cells based on `wang_day3_nonzero_cells`
2. Creates labels using `wang_day3_liana_labels[wang_day3_nonzero_cells]`
3. But the indices may not align properly after multiple filtering steps

The issue is that the filtering operations on the count and data matrices don't maintain proper alignment with the metadata labels.

## Solution

Replace the SCE building section for Wang Day 3 (starting around line with "Wang Day 3: Build SCE") with the following corrected code:

```r
# Wang Day 3: Build SCE (pre-filter LR genes + LR-expressing cells so LIANA subset is no-op)
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
if (inherits(wang_day3[["RNA"]], "Assay5")) wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
if (is.null(wang_day3_cnt) || nrow(wang_day3_cnt) == 0 || ncol(wang_day3_cnt) == 0) wang_day3_cnt <- wang_day3_dat
if (!inherits(wang_day3_cnt, "matrix") && !inherits(wang_day3_cnt, "Matrix")) wang_day3_cnt <- as(wang_day3_cnt, "sparseMatrix")
if (!inherits(wang_day3_dat, "matrix") && !inherits(wang_day3_dat, "Matrix")) wang_day3_dat <- as(wang_day3_dat, "sparseMatrix")

# Store original cell indices before filtering
wang_day3_original_ncells <- ncol(wang_day3_cnt)
wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)

# Filter genes to keep only LR entity genes
wang_day3_genes_in_data <- rownames(wang_day3_cnt)
wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
if (sum(wang_day3_keep_genes) < 3) wang_day3_keep_genes <- rep(TRUE, length(wang_day3_genes_in_data))
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]

# Filter cells based on nonzero counts
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]

# Update cell indices to track which cells remain after filtering
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]

# Filter genes again to remove any that became all-zero after cell filtering
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]

# Create labels using the tracked cell indices (CRITICAL FIX)
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])

# Verify dimensions match before creating SCE
stopifnot(ncol(wang_day3_cnt) == ncol(wang_day3_dat))
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))

# Create SingleCellExperiment with properly aligned data
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = wang_day3_lab)
)
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab
```

## Key Changes

1. **Track Cell Indices**: Added `wang_day3_cell_indices` to track which cells remain after each filtering step
2. **Proper Label Alignment**: Use `wang_day3_liana_labels[wang_day3_cell_indices]` instead of `wang_day3_liana_labels[wang_day3_nonzero_cells]`
3. **Dimension Verification**: Added `stopifnot` checks to verify dimensions match before creating SCE
4. **Clear Comments**: Added comments to explain each filtering step

## Apply to All Similar Sections

The same fix pattern should be applied to:
- Lee Day 1 LIANA section
- Lee Day 3 LIANA section
- Qin Day 1 LIANA section
- Qin Day 3 LIANA section

Each section follows the same pattern and has the same bug.

## Why This Works

The original code applied boolean indexing (`wang_day3_nonzero_cells`) after the data had already been filtered, causing the indices to become misaligned. By tracking the original cell indices and updating them through each filtering step, we maintain proper alignment between the expression matrices and the metadata.
