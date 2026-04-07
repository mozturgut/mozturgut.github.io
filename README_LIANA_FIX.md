# LIANA SCE Fix - Complete Implementation Guide

## Overview

This repository contains a comprehensive fix for the LIANA analysis error:
```
Error in `validObject()`:
! invalid class "SingleCellExperiment" object: 
    'nrow' of 'int_colData' not equal to 'ncol(object)'
```

## What's Included

1. **LIANA_SCE_FIX.md** - Detailed explanation of the root cause and solution
2. **LIANA_analysis_fix.R** - Reference implementation showing the correct pattern
3. **LIANA_FIX_HOWTO.md** - Step-by-step application instructions
4. **COMPREHENSIVE_LIANA_FIX_PATCH.md** - Ready-to-use code for all 5 sections
5. **test_liana_fix.R** - Validation script to verify the fix works
6. **README_LIANA_FIX.md** - This file

## Quick Start

### Step 1: Understand the Problem

The error occurs when building SingleCellExperiment objects for LIANA analysis. After filtering cells and genes, the metadata labels don't align with the filtered expression matrices.

**Root Cause:** Using boolean indices after matrices are already filtered creates misalignment.

### Step 2: Identify Where to Apply the Fix

Look for code sections that build SCE objects, specifically:
- Lee Day 1 LIANA section
- Lee Day 3 LIANA section  
- Wang Day 3 LIANA section
- Qin Day 1 LIANA section
- Qin Day 3 LIANA section

Each section follows this pattern:
```r
# Build SCE
cnt <- GetAssayData(obj, layer = "counts")
# ... filter genes ...
# ... filter cells ...
nonzero_cells <- colSums(cnt) > 0
cnt <- cnt[, nonzero_cells]
lab <- factor(labels[nonzero_cells])  # <-- BUG HERE!
sce <- SingleCellExperiment(...)
```

### Step 3: Apply the Fix

For each affected section, add cell index tracking:

**BEFORE (Buggy):**
```r
# Filter cells
nonzero_cells <- Matrix::colSums(cnt) > 0
cnt <- cnt[, nonzero_cells, drop = FALSE]
dat <- dat[, nonzero_cells, drop = FALSE]

# BUG: Using boolean filter on original labels
lab <- factor(liana_labels[nonzero_cells])
```

**AFTER (Fixed):**
```r
# Track original indices
original_ncells <- ncol(cnt)
cell_indices <- seq_len(original_ncells)

# Filter cells
nonzero_cells <- Matrix::colSums(cnt) > 0
cnt <- cnt[, nonzero_cells, drop = FALSE]
dat <- dat[, nonzero_cells, drop = FALSE]

# FIX: Update tracked indices
cell_indices <- cell_indices[nonzero_cells]

# FIX: Use tracked indices for labels
lab <- factor(liana_labels[cell_indices])

# Verify dimensions match
stopifnot(ncol(cnt) == length(lab))
```

### Step 4: Use the Comprehensive Patch

Open **COMPREHENSIVE_LIANA_FIX_PATCH.md** and copy the corrected code for each section. This file provides complete, ready-to-use code for all 5 affected sections.

### Step 5: Verify the Fix

After applying the fix:

1. **Dimension checks pass:**
   ```r
   stopifnot(ncol(cnt) == ncol(dat))
   stopifnot(ncol(cnt) == length(lab))
   ```

2. **No validObject() errors** when creating SCE

3. **LIANA analysis completes** successfully

4. **(Optional) Run test script:**
   ```r
   Rscript test_liana_fix.R
   ```

## The Three Key Changes

For **each** of the 5 dataset sections (Lee Day 1, Lee Day 3, Wang Day 3, Qin Day 1, Qin Day 3):

### Change 1: Add Index Tracking (BEFORE filtering)
```r
# Add these two lines BEFORE any filtering starts
{dataset}_original_ncells <- ncol({dataset}_cnt)
{dataset}_cell_indices <- seq_len({dataset}_original_ncells)
```

### Change 2: Update Indices (AFTER cell filtering)
```r
# After filtering cells, add this line
{dataset}_cell_indices <- {dataset}_cell_indices[{dataset}_nonzero_cells]
```

### Change 3: Use Tracked Indices (for labels)
```r
# Change FROM:
{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_nonzero_cells])

# Change TO:
{dataset}_lab <- factor({dataset}_liana_labels[{dataset}_cell_indices])
```

Replace `{dataset}` with: `lee_day1`, `lee_day3`, `wang_day3`, `qin_day1`, or `qin_day3`

## Example: Wang Day 3 Complete Fix

See **COMPREHENSIVE_LIANA_FIX_PATCH.md** for the complete code, but here's the key difference:

### Before (causes error):
```r
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells]
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_nonzero_cells])  # BUG!
```

### After (works correctly):
```r
wang_day3_original_ncells <- ncol(wang_day3_cnt)
wang_day3_cell_indices <- seq_len(wang_day3_original_ncells)
wang_day3_nonzero_cells <- Matrix::colSums(wang_day3_cnt) > 0
wang_day3_cnt <- wang_day3_cnt[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_dat <- wang_day3_dat[, wang_day3_nonzero_cells, drop = FALSE]
wang_day3_cell_indices <- wang_day3_cell_indices[wang_day3_nonzero_cells]  # TRACK!
wang_day3_lab <- factor(wang_day3_liana_labels[wang_day3_cell_indices])  # FIXED!
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_lab))  # VERIFY!
```

## Why This Works

The original code used a boolean vector (`nonzero_cells`) to index into the original label vector, but this boolean vector was created from an already-filtered matrix. This created a length mismatch.

The fix tracks which **original cell indices** remain after each filtering step, ensuring perfect alignment between the filtered expression data and the metadata labels.

## Files in This Fix Package

| File | Purpose |
|------|---------|
| LIANA_SCE_FIX.md | Root cause analysis and solution overview |
| LIANA_analysis_fix.R | Reference implementation with patterns |
| LIANA_FIX_HOWTO.md | Step-by-step application guide |
| COMPREHENSIVE_LIANA_FIX_PATCH.md | Complete code for all 5 sections |
| test_liana_fix.R | Validation script |
| README_LIANA_FIX.md | This summary guide |

## Testing

To test the fix works:

```bash
# If you have R installed
Rscript test_liana_fix.R
```

The test script will:
1. Demonstrate the buggy approach
2. Show the fixed approach
3. Validate dimension alignment
4. (If available) Create a test SCE object

## Support

If you encounter issues:

1. Check that you've added ALL three changes (index tracking, index updating, using tracked indices)
2. Verify the change is applied to ALL 5 sections
3. Ensure `stopifnot` dimension checks are included
4. Review **COMPREHENSIVE_LIANA_FIX_PATCH.md** for the exact code

## Summary

The fix is simple but must be applied consistently:
- **Track** cell indices from the start
- **Update** indices after each filtering step
- **Use** tracked indices (not boolean filters) for labels
- **Verify** dimensions match before creating SCE

This ensures the expression matrices and metadata remain perfectly aligned throughout all filtering operations.
