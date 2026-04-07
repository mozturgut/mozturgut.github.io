# LIANA SCE Fix - Quick Reference Card

## The Error
```
Error in `validObject()`:
! invalid class "SingleCellExperiment" object: 
    'nrow' of 'int_colData' not equal to 'ncol(object)'
```

## The Fix (3 Steps per Section)

### STEP 1: Initialize Index Tracking
**Add BEFORE any filtering:**
```r
{dataset}_original_ncells <- ncol({dataset}_cnt)
{dataset}_cell_indices <- seq_len({dataset}_original_ncells)
```

### STEP 2: Update Indices After Cell Filtering
**Add AFTER this line: `{dataset}_cnt <- {dataset}_cnt[, {dataset}_nonzero_cells, drop = FALSE]`**
```r
{dataset}_cell_indices <- {dataset}_cell_indices[{dataset}_nonzero_cells]
```

### STEP 3: Use Tracked Indices for Labels
**Replace:**
```r
{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_nonzero_cells])
```
**With:**
```r
{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_cell_indices])
```

### STEP 4: Add Safety Checks
**Add BEFORE creating SCE:**
```r
stopifnot(ncol({dataset}_cnt) == ncol({dataset}_dat))
stopifnot(ncol({dataset}_cnt) == length({dataset}_lab))
```

## Apply to These 5 Sections

1. **Lee Day 1**: Replace `{dataset}` with `lee_day1`
2. **Lee Day 3**: Replace `{dataset}` with `lee_day3`
3. **Wang Day 3**: Replace `{dataset}` with `wang_day3`
4. **Qin Day 1**: Replace `{dataset}` with `qin_day1`
5. **Qin Day 3**: Replace `{dataset}` with `qin_day3`

## Complete Example (Wang Day 3)

### BEFORE (Buggy)
```r
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
# ... matrix conversions ...
wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells]
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_nonzero_cells])  # ❌ BUG!
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(...)
```

### AFTER (Fixed)
```r
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")
# ... matrix conversions ...
wang_day3_original_ncells <- ncol(wang_day3_cnt)              # ✓ STEP 1
wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)  # ✓ STEP 1
wang_day3_keep_genes <- wang_day3_genes_in_data %in% wang_day3_entity_genes
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]  # ✓ STEP 2
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])  # ✓ STEP 3
stopifnot(ncol(wang_day3_cnt) == ncol(wang_day3_dat))                    # ✓ STEP 4
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))                  # ✓ STEP 4
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(...)
```

## Checklist

For each of the 5 sections, verify:

- [ ] Added `{dataset}_original_ncells <- ncol({dataset}_cnt)`
- [ ] Added `{dataset}_cell_indices <- seq_len({dataset}_original_ncells)`
- [ ] Added `{dataset}_cell_indices <- {dataset}_cell_indices[{dataset}_nonzero_cells]` after cell filtering
- [ ] Changed `{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_nonzero_cells])` to use `[{dataset}_cell_indices]`
- [ ] Added `stopifnot(ncol({dataset}_cnt) == ncol({dataset}_dat))`
- [ ] Added `stopifnot(ncol({dataset}_cnt) == length({dataset}_lab))`

## Why This Works

**Before:** Boolean filter from filtered data → applied to original labels → ❌ dimension mismatch

**After:** Track original indices → update through all filters → always aligned → ✓ dimensions match

## Need More Help?

- Full code for all sections: **COMPREHENSIVE_LIANA_FIX_PATCH.md**
- Detailed guide: **LIANA_FIX_HOWTO.md**
- Visual explanation: **VISUAL_EXPLANATION.md**
- Root cause analysis: **LIANA_SCE_FIX.md**
- Complete README: **README_LIANA_FIX.md**
