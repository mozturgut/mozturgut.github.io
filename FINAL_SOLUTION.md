# LIANA Fix - Final Solution Summary

## The Problem (Actual Error from Your Code)

```
Error in `validObject()`:
! invalid class "SingleCellExperiment" object: 1: 
    'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'
invalid class "SingleCellExperiment" object: 2: 
    'nrow' of 'int_colData' not equal to 'ncol(object)'
```

With warnings:
- "9161 genes were removed as they had no counts!"
- "24146 cells were excluded as they did not express any ligand-receptor genes!"

## Root Cause

When filtering cells/genes before creating the SCE object, the code used boolean filters from AFTER filtering to index into labels from BEFORE filtering, causing dimension misalignment.

## Solution: Two Approaches

### APPROACH 1: SIMPLE ✅ (RECOMMENDED)

**Don't pre-filter. Let LIANA do all filtering.**

```r
# Wang Day 3 Example
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
if (inherits(wang_day3[["RNA"]], "Assay5")) {
  wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
}

wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")

# Verify dimensions
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_liana_labels))

# Create SCE (no filtering!)
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = factor(wang_day3_liana_labels))
)
SingleCellExperiment::colLabels(wang_day3_sce) <- factor(wang_day3_liana_labels)

# Verify
stopifnot(validObject(wang_day3_sce))

# Run LIANA (it handles filtering)
wang_day3_liana_result <- liana::liana_wrap(
  sce = wang_day3_sce,
  resource = wang_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
```

**Why this works:** LIANA's internal filtering is designed to work on a complete SCE object. Let it do its job.

### APPROACH 2: COMPLEX (Advanced Users)

**If you must pre-filter, track cell indices correctly.**

```r
# Extract matrices
wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")

# CRITICAL: Track cell indices
wang_day3_cell_indices <- seq_len(ncol(wang_day3_cnt))

# Filter genes
wang_day3_keep_genes <- rownames(wang_day3_cnt) %in% wang_day3_entity_genes
wang_day3_cnt <- wang_day3_cnt[wang_day3_keep_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_keep_genes, , drop = FALSE]

# Filter cells
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]

# CRITICAL: Update tracked indices
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]

# Filter genes again
wang_day3_nonzero_genes <- Matrix::rowSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[wang_day3_nonzero_genes, , drop = FALSE]
wang_day3_dat <- wang_day3_dat[wang_day3_nonzero_genes, , drop = FALSE]

# CRITICAL: Use tracked indices for labels
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])

# Verify
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))

# Create SCE
wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = wang_day3_lab)
)
SingleCellExperiment::colLabels(wang_day3_sce) <- wang_day3_lab

# Verify
stopifnot(validObject(wang_day3_sce))

# Run LIANA
wang_day3_liana_result <- liana::liana_wrap(
  sce = wang_day3_sce,
  resource = wang_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
```

**Why this works:** By tracking which original cell indices survive filtering, we maintain perfect alignment between matrices and labels.

## Quick Decision Guide

**Use APPROACH 1 (Simple) if:**
- You just want it to work
- You don't need to control pre-filtering
- You're okay with LIANA doing the filtering

**Use APPROACH 2 (Complex) if:**
- You need specific control over which genes/cells to include
- You understand index tracking
- You want to minimize LIANA's internal processing

## Apply to All 5 Sections

Both approaches work for:
1. Lee Day 1 (`lee_day1`)
2. Lee Day 3 (`lee_day3`)
3. Wang Day 3 (`wang_day3`)
4. Qin Day 1 (`qin_day1`)
5. Qin Day 3 (`qin_day3`)

Just replace `wang_day3` with the appropriate dataset name.

## Verification Checklist

Before calling `liana_wrap()`, verify:

```r
# 1. Matrices have same dimensions
stopifnot(nrow(cnt) == nrow(dat))
stopifnot(ncol(cnt) == ncol(dat))

# 2. Labels match cell count
stopifnot(ncol(cnt) == length(labels))

# 3. SCE is valid
stopifnot(validObject(sce))

# 4. colData has correct rows
stopifnot(nrow(SingleCellExperiment::colData(sce)) == ncol(sce))
```

## Files in This Fix Package

| File | Purpose | Read This If... |
|------|---------|-----------------|
| **FINAL_SOLUTION.md** | **This file** | You want the complete solution |
| **LIANA_WORKING_EXAMPLE.R** | **Complete working code** | You want copy-paste solution |
| **UPDATED_LIANA_FIX.md** | Updated fix based on real error | You want detailed explanation |
| QUICK_REFERENCE.md | 1-page cheat sheet | You want quick lookup |
| README_LIANA_FIX.md | Complete guide | You want overview |
| VISUAL_EXPLANATION.md | Visual diagrams | You want to understand why |
| COMPREHENSIVE_LIANA_FIX_PATCH.md | Code for all 5 sections | You want all sections at once |
| INDEX_LIANA_FIX.md | Master index | You want navigation |

## Recommended Next Steps

1. **Choose your approach:** Simple (recommended) or Complex (advanced)
2. **Copy the code:** From `LIANA_WORKING_EXAMPLE.R` or this file
3. **Apply to each section:** Lee Day 1, Lee Day 3, Wang Day 3, Qin Day 1, Qin Day 3
4. **Test with one section first:** Verify it works before applying to all
5. **Check the output:** Ensure LIANA completes without errors

## Common Mistakes to Avoid

### ❌ DON'T:
```r
# Filter cells
nonzero_cells <- colSums(filtered_matrix) > 0
filtered_matrix <- filtered_matrix[, nonzero_cells]

# BUG: Using boolean from filtered data on original labels
labels <- original_labels[nonzero_cells]  # WRONG!
```

### ✅ DO (Simple):
```r
# No filtering - just create SCE
sce <- SingleCellExperiment(
  assays = list(counts = cnt, logcounts = dat),
  colData = S4Vectors::DataFrame(label = factor(labels))
)
```

### ✅ DO (Complex):
```r
# Track indices
cell_indices <- seq_len(ncol(matrix))

# Filter
nonzero_cells <- colSums(matrix) > 0
matrix <- matrix[, nonzero_cells]

# Update indices
cell_indices <- cell_indices[nonzero_cells]

# Use tracked indices
labels <- original_labels[cell_indices]  # CORRECT!
```

## Support

If issues persist:
1. Check you applied the fix to ALL 5 sections
2. Verify `validObject(sce)` returns `TRUE` before LIANA
3. Check `ncol(sce) == length(colLabels(sce))`
4. Review `LIANA_WORKING_EXAMPLE.R` for complete working code
5. Try the SIMPLE approach first - it's more reliable

## Success Criteria

You'll know it works when:
- ✅ No `validObject()` errors
- ✅ No dimension mismatch warnings
- ✅ LIANA completes without errors
- ✅ Results are generated for all methods

## Summary

**Problem:** Dimension mismatch when creating SCE for LIANA  
**Cause:** Boolean filters from filtered data applied to original labels  
**Solution 1:** Let LIANA do all filtering (simple, recommended)  
**Solution 2:** Track cell indices through all filtering (complex, advanced)  
**Result:** SCE object with perfectly aligned data and metadata

Choose the approach that fits your needs and coding comfort level. Both will work!
