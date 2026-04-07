# How to Apply the LIANA SCE Fix

## Quick Summary

The error `'nrow' of 'int_colData' not equal to 'ncol(object)'` occurs because cell filtering creates misalignment between expression matrices and metadata labels.

**The Fix:** Track cell indices through all filtering steps.

## Step-by-Step Application

### For Wang Day 3 Section

Find the section that builds the SCE object (around the comment `# Wang Day 3: Build SCE`). Replace the existing code with:

```r
# Wang Day 3: Build SCE (pre-filter LR genes + LR-expressing cells so LIANA subset is no-op)
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
if (inherits(wang_day3[["RNA"]], "Assay5")) wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
if (is.null(wang_day3_cnt) || nrow(wang_day3_cnt) == 0 || ncol(wang_day3_cnt) == 0) wang_day3_cnt <- wang_day3_dat
if (!inherits(wang_day3_cnt, "matrix") && !inherits(wang_day3_cnt, "Matrix")) wang_day3_cnt <- as(wang_day3_cnt, "sparseMatrix")
if (!inherits(wang_day3_dat, "matrix") && !inherits(wang_day3_dat, "Matrix")) wang_day3_dat <- as(wang_day3_dat, "sparseMatrix")

# CRITICAL: Store original cell indices before filtering
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

# CRITICAL FIX: Update cell indices to track which cells remain
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]

# Filter genes again to remove any that became all-zero
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]

# CRITICAL FIX: Create labels using tracked cell indices (not boolean filter)
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

### For Lee Day 1 Section

Apply the same pattern with `lee_day1` prefixes:

**Key changes:**
1. Add: `lee_day1_original_ncells <- ncol(lee_day1_cnt)`
2. Add: `lee_day1_cell_indices <- seq_len(lee_day1_original_ncells)`
3. After cell filtering, add: `lee_day1_cell_indices <- lee_day1_cell_indices[lee_day1_nonzero_cells]`
4. Change: `lee_day1_lab <- factor(lee_day1_liana_labels[lee_day1_nonzero_cells])` 
   to: `lee_day1_lab <- factor(lee_day1_liana_labels[lee_day1_cell_indices])`

### For Lee Day 3 Section

Apply the same pattern with `lee_day3` prefixes (same as Lee Day 1, just different prefix).

### For Qin Day 1 Section

Apply the same pattern with `qin_day1` prefixes (same as Lee Day 1, just different prefix).

### For Qin Day 3 Section

Apply the same pattern with `qin_day3` prefixes (same as Lee Day 1, just different prefix).

## What Changed and Why

### Before (Buggy)
```r
# Filter cells
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells]

# BUG: nonzero_cells is a boolean vector that was used AFTER filtering
# So it doesn't correspond to the original indices anymore
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_nonzero_cells])
```

### After (Fixed)
```r
# Track original indices
wang_day3_cell_indices <- seq_len(ncol(wang_day3_cnt))

# Filter cells
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells]

# FIX: Update indices to track which cells survived filtering
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]

# FIX: Use tracked indices to get the right labels
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])
```

## Verification

After applying the fix, you should see:
- No `validObject()` errors
- LIANA analysis completes successfully
- The `stopifnot` checks pass, confirming dimensions match

## Additional Safety Checks

The fixed code includes these safety checks:
```r
stopifnot(ncol(wang_day3_cnt) == ncol(wang_day3_dat))
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))
```

If these fail, it means there's still a dimension mismatch and the filtering logic needs review.

## Common Mistakes to Avoid

1. **Don't** use the boolean filter directly for indexing after matrices are filtered:
   ```r
   # WRONG
   labels <- liana_labels[nonzero_cells]  # After matrices already filtered
   ```

2. **Do** track indices through all filtering steps:
   ```r
   # CORRECT
   cell_indices <- seq_len(original_ncells)
   cell_indices <- cell_indices[nonzero_cells]
   labels <- liana_labels[cell_indices]
   ```

3. **Don't** reuse the same boolean vector name after filtering:
   ```r
   # CONFUSING (avoid)
   cells <- colSums(mat) > 0
   mat <- mat[, cells]
   mat <- mat[rowSums(mat) > 0, ]
   cells <- colSums(mat) > 0  # Reusing same name - confusing!
   ```

4. **Do** use descriptive names that indicate the filtering stage:
   ```r
   # CLEAR
   nonzero_cells_stage1 <- colSums(mat) > 0
   mat <- mat[, nonzero_cells_stage1]
   # or better yet, use index tracking as shown above
   ```
