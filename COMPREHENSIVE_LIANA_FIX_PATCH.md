# COMPREHENSIVE FIX FOR LIANA SCE DIMENSION MISMATCH
# ============================================================================
# This document provides the exact code replacement for all affected sections
# Apply these changes to fix: 'nrow' of 'int_colData' not equal to 'ncol(object)'
# ============================================================================

## SECTION: Wang Day 3 LIANA

### ORIGINAL (BUGGY) CODE:
```r
# Wang Day 3: Build SCE (pre-filter LR genes + LR-expressing cells so LIANA subset is no-op)
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
if (inherits(wang_day3[["RNA"]], "Assay5")) wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
if (is.null(wang_day3_cnt) || nrow(wang_day3_cnt) == 0 || ncol(wang_day3_cnt) == 0) wang_day3_cnt <- wang_day3_dat
if (!inherits(wang_day3_cnt, "matrix") && !inherits(wang_day3_cnt, "Matrix")) wang_day3_cnt <- as(wang_day3_cnt, "sparseMatrix")
if (!inherits(wang_day3_dat, "matrix") && !inherits(wang_day3_dat, "Matrix")) wang_day3_dat <- as(wang_day3_dat, "sparseMatrix")
wang_day3_genes_in_data <- rownames(wang_day3_cnt)
wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
if (sum(wang_day3_keep_genes) < 3) wang_day3_keep_genes <- rep(TRUE, length(wang_day3_genes_in_data))
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells]
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_nonzero_cells])  # BUG HERE!
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat), colData = S4Vectors::DataFrame(label = wang_day3_lab))
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab
```

### REPLACEMENT (FIXED) CODE:
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
# Filter genes
wang_day3_genes_in_data <- rownames(wang_day3_cnt)
wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
if (sum(wang_day3_keep_genes) < 3) wang_day3_keep_genes <- rep(TRUE, length(wang_day3_genes_in_data))
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
# Filter cells
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]
# CRITICAL FIX: Update cell indices
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]
# Filter genes again
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
# CRITICAL FIX: Use tracked indices for labels
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])
# Verify dimensions
stopifnot(ncol(wang_day3_cnt) == ncol(wang_day3_dat))
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))
# Create SCE
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat), colData = S4Vectors::DataFrame(label = wang_day3_lab))
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab
```

## SECTION: Lee Day 1 LIANA

Apply the same fix pattern. Replace the SCE building section with:

```r
# Lee Day 1: Build SCE
SeuratObject::DefaultAssay(lee_day1) <- "RNA"
if (inherits(lee_day1[["RNA"]], "Assay5")) lee_day1[["RNA"]] <- SeuratObject::JoinLayers(lee_day1[["RNA"]])
lee_day1_cnt <- SeuratObject::GetAssayData(lee_day1, layer = "counts")
lee_day1_dat <- SeuratObject::GetAssayData(lee_day1, layer = "data")
if (is.null(lee_day1_cnt) || nrow(lee_day1_cnt) == 0 || ncol(lee_day1_cnt) == 0) lee_day1_cnt <- lee_day1_dat
if (!inherits(lee_day1_cnt, "matrix") && !inherits(lee_day1_cnt, "Matrix")) lee_day1_cnt <- as(lee_day1_cnt, "sparseMatrix")
if (!inherits(lee_day1_dat, "matrix") && !inherits(lee_day1_dat, "Matrix")) lee_day1_dat <- as(lee_day1_dat, "sparseMatrix")
lee_day1_original_ncells <- ncol(lee_day1_cnt)
lee_day1_cell_indices <- seq_len(lee_day1_original_ncells)
lee_day1_genes_in_data <- rownames(lee_day1_cnt)
lee_day1_keep_genes <- lee_day1_genes_in_data %in% lee_day1_entity_genes
if (sum(lee_day1_keep_genes) < 3) lee_day1_keep_genes <- rep(TRUE, length(lee_day1_genes_in_data))
lee_day1_cnt <- lee_day1_cnt[lee_day1_keep_genes, , drop = FALSE]
lee_day1_dat <- lee_day1_dat[lee_day1_keep_genes, , drop = FALSE]
lee_day1_nonzero_cells <- Matrix::colSums(lee_day1_cnt) > 0
lee_day1_cnt <- lee_day1_cnt[, lee_day1_nonzero_cells, drop = FALSE]
lee_day1_dat <- lee_day1_dat[, lee_day1_nonzero_cells, drop = FALSE]
lee_day1_cell_indices <- lee_day1_cell_indices[lee_day1_nonzero_cells]
lee_day1_nonzero_genes <- Matrix::rowSums(lee_day1_cnt) > 0
lee_day1_cnt <- lee_day1_cnt[lee_day1_nonzero_genes, , drop = FALSE]
lee_day1_dat <- lee_day1_dat[lee_day1_nonzero_genes, , drop = FALSE]
lee_day1_lab <- factor(lee_day1_liana_labels[lee_day1_cell_indices])
stopifnot(ncol(lee_day1_cnt) == ncol(lee_day1_dat))
stopifnot(ncol(lee_day1_cnt) == length(lee_day1_lab))
lee_day1_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = lee_day1_cnt, logcounts = lee_day1_dat), colData = S4Vectors::DataFrame(label = lee_day1_lab))
SingleCellExperiment::colLabels(lee_day1_sce) <- lee_day1_lab
```

## SECTION: Lee Day 3 LIANA

```r
# Lee Day 3: Build SCE
SeuratObject::DefaultAssay(lee_day3) <- "RNA"
if (inherits(lee_day3[["RNA"]], "Assay5")) lee_day3[["RNA"]] <- SeuratObject::JoinLayers(lee_day3[["RNA"]])
lee_day3_cnt <- SeuratObject::GetAssayData(lee_day3, layer = "counts")
lee_day3_dat <- SeuratObject::GetAssayData(lee_day3, layer = "data")
if (is.null(lee_day3_cnt) || nrow(lee_day3_cnt) == 0 || ncol(lee_day3_cnt) == 0) lee_day3_cnt <- lee_day3_dat
if (!inherits(lee_day3_cnt, "matrix") && !inherits(lee_day3_cnt, "Matrix")) lee_day3_cnt <- as(lee_day3_cnt, "sparseMatrix")
if (!inherits(lee_day3_dat, "matrix") && !inherits(lee_day3_dat, "Matrix")) lee_day3_dat <- as(lee_day3_dat, "sparseMatrix")
lee_day3_original_ncells <- ncol(lee_day3_cnt)
lee_day3_cell_indices <- seq_len(lee_day3_original_ncells)
lee_day3_genes_in_data <- rownames(lee_day3_cnt)
lee_day3_keep_genes <- lee_day3_genes_in_data %in% lee_day3_entity_genes
if (sum(lee_day3_keep_genes) < 3) lee_day3_keep_genes <- rep(TRUE, length(lee_day3_genes_in_data))
lee_day3_cnt <- lee_day3_cnt[lee_day3_keep_genes, , drop = FALSE]
lee_day3_dat <- lee_day3_dat[lee_day3_keep_genes, , drop = FALSE]
lee_day3_nonzero_cells <- Matrix::colSums(lee_day3_cnt) > 0
lee_day3_cnt <- lee_day3_cnt[, lee_day3_nonzero_cells, drop = FALSE]
lee_day3_dat <- lee_day3_dat[, lee_day3_nonzero_cells, drop = FALSE]
lee_day3_cell_indices <- lee_day3_cell_indices[lee_day3_nonzero_cells]
lee_day3_nonzero_genes <- Matrix::rowSums(lee_day3_cnt) > 0
lee_day3_cnt <- lee_day3_cnt[lee_day3_nonzero_genes, , drop = FALSE]
lee_day3_dat <- lee_day3_dat[lee_day3_nonzero_genes, , drop = FALSE]
lee_day3_lab <- factor(lee_day3_liana_labels[lee_day3_cell_indices])
stopifnot(ncol(lee_day3_cnt) == ncol(lee_day3_dat))
stopifnot(ncol(lee_day3_cnt) == length(lee_day3_lab))
lee_day3_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = lee_day3_cnt, logcounts = lee_day3_dat), colData = S4Vectors::DataFrame(label = lee_day3_lab))
SingleCellExperiment::colLabels(lee_day3_sce) <- lee_day3_lab
```

## SECTION: Qin Day 1 LIANA

```r
# Qin Day 1: Build SCE
SeuratObject::DefaultAssay(qin_day1) <- "RNA"
if (inherits(qin_day1[["RNA"]], "Assay5")) qin_day1[["RNA"]] <- SeuratObject::JoinLayers(qin_day1[["RNA"]])
qin_day1_cnt <- SeuratObject::GetAssayData(qin_day1, layer = "counts")
qin_day1_dat <- SeuratObject::GetAssayData(qin_day1, layer = "data")
if (is.null(qin_day1_cnt) || nrow(qin_day1_cnt) == 0 || ncol(qin_day1_cnt) == 0) qin_day1_cnt <- qin_day1_dat
if (!inherits(qin_day1_cnt, "matrix") && !inherits(qin_day1_cnt, "Matrix")) qin_day1_cnt <- as(qin_day1_cnt, "sparseMatrix")
if (!inherits(qin_day1_dat, "matrix") && !inherits(qin_day1_dat, "Matrix")) qin_day1_dat <- as(qin_day1_dat, "sparseMatrix")
qin_day1_original_ncells <- ncol(qin_day1_cnt)
qin_day1_cell_indices <- seq_len(qin_day1_original_ncells)
qin_day1_genes_in_data <- rownames(qin_day1_cnt)
qin_day1_keep_genes <- qin_day1_genes_in_data %in% qin_day1_entity_genes
if (sum(qin_day1_keep_genes) < 3) qin_day1_keep_genes <- rep(TRUE, length(qin_day1_genes_in_data))
qin_day1_cnt <- qin_day1_cnt[qin_day1_keep_genes, , drop = FALSE]
qin_day1_dat <- qin_day1_dat[qin_day1_keep_genes, , drop = FALSE]
qin_day1_nonzero_cells <- Matrix::colSums(qin_day1_cnt) > 0
qin_day1_cnt <- qin_day1_cnt[, qin_day1_nonzero_cells, drop = FALSE]
qin_day1_dat <- qin_day1_dat[, qin_day1_nonzero_cells, drop = FALSE]
qin_day1_cell_indices <- qin_day1_cell_indices[qin_day1_nonzero_cells]
qin_day1_nonzero_genes <- Matrix::rowSums(qin_day1_cnt) > 0
qin_day1_cnt <- qin_day1_cnt[qin_day1_nonzero_genes, , drop = FALSE]
qin_day1_dat <- qin_day1_dat[qin_day1_nonzero_genes, , drop = FALSE]
qin_day1_lab <- factor(qin_day1_liana_labels[qin_day1_cell_indices])
stopifnot(ncol(qin_day1_cnt) == ncol(qin_day1_dat))
stopifnot(ncol(qin_day1_cnt) == length(qin_day1_lab))
qin_day1_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = qin_day1_cnt, logcounts = qin_day1_dat), colData = S4Vectors::DataFrame(label = qin_day1_lab))
SingleCellExperiment::colLabels(qin_day1_sce) <- qin_day1_lab
```

## SECTION: Qin Day 3 LIANA

```r
# Qin Day 3: Build SCE
SeuratObject::DefaultAssay(qin_day3) <- "RNA"
if (inherits(qin_day3[["RNA"]], "Assay5")) qin_day3[["RNA"]] <- SeuratObject::JoinLayers(qin_day3[["RNA"]])
qin_day3_cnt <- SeuratObject::GetAssayData(qin_day3, layer = "counts")
qin_day3_dat <- SeuratObject::GetAssayData(qin_day3, layer = "data")
if (is.null(qin_day3_cnt) || nrow(qin_day3_cnt) == 0 || ncol(qin_day3_cnt) == 0) qin_day3_cnt <- qin_day3_dat
if (!inherits(qin_day3_cnt, "matrix") && !inherits(qin_day3_cnt, "Matrix")) qin_day3_cnt <- as(qin_day3_cnt, "sparseMatrix")
if (!inherits(qin_day3_dat, "matrix") && !inherits(qin_day3_dat, "Matrix")) qin_day3_dat <- as(qin_day3_dat, "sparseMatrix")
qin_day3_original_ncells <- ncol(qin_day3_cnt)
qin_day3_cell_indices <- seq_len(qin_day3_original_ncells)
qin_day3_genes_in_data <- rownames(qin_day3_cnt)
qin_day3_keep_genes <- qin_day3_genes_in_data %in% qin_day3_entity_genes
if (sum(qin_day3_keep_genes) < 3) qin_day3_keep_genes <- rep(TRUE, length(qin_day3_genes_in_data))
qin_day3_cnt <- qin_day3_cnt[qin_day3_keep_genes, , drop = FALSE]
qin_day3_dat <- qin_day3_dat[qin_day3_keep_genes, , drop = FALSE]
qin_day3_nonzero_cells <- Matrix::colSums(qin_day3_cnt) > 0
qin_day3_cnt <- qin_day3_cnt[, qin_day3_nonzero_cells, drop = FALSE]
qin_day3_dat <- qin_day3_dat[, qin_day3_nonzero_cells, drop = FALSE]
qin_day3_cell_indices <- qin_day3_cell_indices[qin_day3_nonzero_cells]
qin_day3_nonzero_genes <- Matrix::rowSums(qin_day3_cnt) > 0
qin_day3_cnt <- qin_day3_cnt[qin_day3_nonzero_genes, , drop = FALSE]
qin_day3_dat <- qin_day3_dat[qin_day3_nonzero_genes, , drop = FALSE]
qin_day3_lab <- factor(qin_day3_liana_labels[qin_day3_cell_indices])
stopifnot(ncol(qin_day3_cnt) == ncol(qin_day3_dat))
stopifnot(ncol(qin_day3_cnt) == length(qin_day3_lab))
qin_day3_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = qin_day3_cnt, logcounts = qin_day3_dat), colData = S4Vectors::DataFrame(label = qin_day3_lab))
SingleCellExperiment::colLabels(qin_day3_sce) <- qin_day3_lab
```

## KEY CHANGES IN ALL SECTIONS

1. **Added before filtering:**
   - `{dataset}_original_ncells <- ncol({dataset}_cnt)`
   - `{dataset}_cell_indices <- seq_len({dataset}_original_ncells)`

2. **Added after cell filtering:**
   - `{dataset}_cell_indices <- {dataset}_cell_indices[{dataset}_nonzero_cells]`

3. **Changed label creation from:**
   - `{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_nonzero_cells])`
   
   **To:**
   - `{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_cell_indices])`

4. **Added dimension verification:**
   - `stopifnot(ncol({dataset}_cnt) == ncol({dataset}_dat))`
   - `stopifnot(ncol({dataset}_cnt) == length({dataset}_lab))`

These changes ensure proper alignment between filtered expression matrices and metadata labels, preventing the `'nrow' of 'int_colData' not equal to 'ncol(object)'` error.
