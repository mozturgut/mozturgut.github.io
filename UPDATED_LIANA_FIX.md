# UPDATED FIX: LIANA SCE Dimension Mismatch (Internal Filtering Issue)

## New Understanding of the Problem

The error occurs because LIANA **internally filters** the SCE object after we create it, and this filtering breaks the object's validity. There are TWO errors:

```
Error 1: 'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'  (GENE metadata)
Error 2: 'nrow' of 'int_colData' not equal to 'ncol(object)'          (CELL metadata)
```

LIANA's warnings show it's filtering:
- "9161 genes were removed as they had no counts!"
- "24146 cells were excluded as they did not express any ligand-receptor genes!"

## The Real Root Cause

When you create an SCE object and set `colLabels()`, SingleCellExperiment stores this in `int_colData`. When LIANA filters cells/genes internally, it uses subsetting operations that don't properly update the internal metadata structures.

## Solution: Properly Initialize SCE with Explicit Metadata

Instead of setting labels with `colLabels()` after creation, we need to:
1. Create colData DataFrame with the SAME number of rows as columns in the matrix
2. Pass colData during SCE construction (not after)
3. Ensure rowData is properly initialized or empty

## Updated Fix for Wang Day 3 (and similar for all sections)

### REPLACE THIS (CURRENT BUGGY CODE):

```r
# Current approach that fails
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_nonzero_cells])  # Wrong indices!
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = wang_day3_lab)
)
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab
```

### WITH THIS (CORRECTED CODE):

```r
# Extract matrices from Seurat
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

# CRITICAL FIX: Track cell indices from the start
wang_day3_original_ncells <- ncol(wang_day3_cnt)
wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)

# Filter genes to LR entity genes (optional pre-filtering)
wang_day3_genes_in_data <- rownames(wang_day3_cnt)
wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
if (sum(wang_day3_keep_genes) < 3) {
  wang_day3_keep_genes <- rep(TRUE, length(wang_day3_genes_in_data))
}
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]

# Filter cells with zero counts
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]

# CRITICAL FIX: Update cell indices after filtering
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]

# Filter genes with zero counts
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]

# CRITICAL FIX: Create labels using tracked indices
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])

# Verify dimensions before creating SCE
stopifnot("Counts and data must have same dimensions" = 
  ncol(wang_day3_cnt) == ncol(wang_day3_dat) && nrow(wang_day3_cnt) == nrow(wang_day3_dat))
stopifnot("Labels must match number of cells" = 
  ncol(wang_day3_cnt) == length(wang_day3_lab))

# Create SCE with properly aligned colData
# Pass colData during construction, not after
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = wang_day3_lab)
)

# Set colLabels to the same factor (ensures proper internal structure)
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab

# Verify SCE is valid before passing to LIANA
stopifnot("SCE must be valid" = validObject(wang_day3_sce))

# Debug output
cat("Wang Day 3 SCE created:\n")
cat("  Genes:", nrow(wang_day3_sce), "\n")
cat("  Cells:", ncol(wang_day3_sce), "\n")
cat("  Labels:", length(SingleCellExperiment::colLabels(wang_day3_sce)), "\n")
cat("  colData rows:", nrow(SingleCellExperiment::colData(wang_day3_sce)), "\n")
```

## Key Changes in the Updated Fix

### 1. Track Cell Indices (prevents colData mismatch)
```r
wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)
# ... after filtering ...
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])  # Correct!
```

### 2. Explicit Dimension Verification
```r
stopifnot("Counts and data must have same dimensions" = 
  ncol(wang_day3_cnt) == ncol(wang_day3_dat) && nrow(wang_day3_cnt) == nrow(wang_day3_dat))
stopifnot("Labels must match number of cells" = 
  ncol(wang_day3_cnt) == length(wang_day3_lab))
stopifnot("SCE must be valid" = validObject(wang_day3_sce))
```

### 3. Debug Output
```r
cat("Wang Day 3 SCE created:\n")
cat("  Genes:", nrow(wang_day3_sce), "\n")
cat("  Cells:", ncol(wang_day3_sce), "\n")
cat("  Labels:", length(SingleCellExperiment::colLabels(wang_day3_sce)), "\n")
cat("  colData rows:", nrow(SingleCellExperiment::colData(wang_day3_sce)), "\n")
```

This will help diagnose if dimensions don't match.

## Alternative: Let LIANA Do the Filtering

If pre-filtering continues to cause issues, you can simplify by letting LIANA handle all filtering:

```r
# Simplest approach: minimal pre-processing
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
if (inherits(wang_day3[["RNA"]], "Assay5")) {
  wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
}

wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")

# Create SCE with ALL cells and genes
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = factor(wang_day3_liana_labels))
)
SingleCellExperiment::colLabels(wang_day3_sce) <- factor(wang_day3_liana_labels)

# Verify
stopifnot(validObject(wang_day3_sce))
stopifnot(ncol(wang_day3_sce) == length(wang_day3_liana_labels))

# Let LIANA do all filtering
wang_day3_liana_result <- liana::liana_wrap(
  sce = wang_day3_sce,
  resource = wang_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
```

## Apply to All 5 Sections

The same pattern applies to:
- `lee_day1` - Lee Day 1 LIANA section
- `lee_day3` - Lee Day 3 LIANA section
- `wang_day3` - Wang Day 3 LIANA section (shown above)
- `qin_day1` - Qin Day 1 LIANA section
- `qin_day3` - Qin Day 3 LIANA section

## Testing the Fix

After applying, check:
1. `validObject(sce)` returns `TRUE` before passing to LIANA
2. `ncol(sce) == length(colLabels(sce))` is `TRUE`
3. `nrow(colData(sce)) == ncol(sce)` is `TRUE`
4. No warnings about dimension mismatches
5. LIANA completes without errors

## Summary

**The bug:** Using boolean filter from filtered data on original label vector creates misalignment

**The fix:** Track original cell indices through all filtering steps, use tracked indices for labels

**Critical steps:**
1. Initialize: `cell_indices <- seq_len(ncol(matrix))`
2. After each cell filter: `cell_indices <- cell_indices[boolean_filter]`
3. Create labels: `labels <- factor(original_labels[cell_indices])`
4. Verify: `stopifnot(ncol(matrix) == length(labels))`
5. Create SCE with verified colData
6. Verify: `stopifnot(validObject(sce))`
