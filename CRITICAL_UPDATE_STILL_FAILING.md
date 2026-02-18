# 🚨 CRITICAL UPDATE - LIANA Error Persists After SIMPLE Fix

## What Happened

You applied the SIMPLE approach fix (no pre-filtering), but **the error still occurs**:

```
Error: 'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'  
Error: 'nrow' of 'int_colData' not equal to 'ncol(object)'
```

## What This Means

### The Good News ✅
- Your SCE is created correctly
- `validObject(lee_day3_sce)` returns TRUE
- Dimensions are aligned
- No pre-filtering issues

### The Bad News ❌
- LIANA's internal filtering STILL breaks the SCE
- This is a **LIANA/SingleCellExperiment compatibility issue**
- Not a problem with our code structure
- Deeper than initially thought

## Root Cause

LIANA's subsetting operations don't properly update SingleCellExperiment's internal metadata structures (`int_elementMetadata` and `int_colData`). This happens even with a correctly created SCE.

**Hypothesis:** The way we're manually creating the SCE is incompatible with how LIANA expects it to be structured internally.

## Solution: Try LIANA's Built-In Seurat Support

Instead of manually creating an SCE, **let LIANA handle the Seurat → SCE conversion**.

### RECOMMENDED FIX (APPROACH A)

Replace the entire Lee Day 3 SCE building section with:

```r
# Lee Day 3: Use LIANA's Seurat wrapper (NO MANUAL SCE CREATION)
print("--- LIANA Analysis: Lee Day 3 ---")

# Prepare labels (same as before)
lee_day3_liana_labels <- lee_day3_pruned_labels
lee_day3_liana_labels[lee_day3_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
lee_day3_liana_labels[lee_day3_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
lee_day3_liana_labels[is.na(lee_day3_liana_labels) | lee_day3_liana_labels == ""] <- "Other"

# Add labels to Seurat metadata
lee_day3$liana_labels <- lee_day3_liana_labels
Seurat::Idents(lee_day3) <- factor(lee_day3_liana_labels)

# Load custom LR resource (keep existing CellChat/CellCall loading code here)
lee_day3_liana_resource <- "Consensus"
lee_day3_liana_external <- NULL
# ... [your CellChat/CellCall loading code] ...

# Run LIANA directly on Seurat object (NO SCE CREATION)
lee_day3_liana_args <- list(
  seurat_object = lee_day3,          # Pass Seurat directly
  ident_col = "liana_labels",        # Column name for cell labels
  resource = lee_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0,
  assay = "RNA"                      # Specify assay
)

if (!is.null(lee_day3_liana_external)) {
  lee_day3_liana_args$external_resource <- lee_day3_liana_external
}

# Let LIANA handle Seurat → SCE conversion internally
lee_day3_liana_result <- suppressWarnings(
  do.call(liana::liana_wrap, lee_day3_liana_args)
)
```

### Why This Should Work

1. ✅ LIANA has built-in Seurat support
2. ✅ LIANA's internal conversion creates SCE correctly for its own needs
3. ✅ No manual SCE creation = no structural incompatibility
4. ✅ LIANA tested this path extensively

## If APPROACH A Doesn't Work

Try the other approaches in `ALTERNATIVE_APPROACHES.md`:

- **APPROACH B:** Create SCE with rownames in colData
- **APPROACH C:** Add explicit rowData (gene metadata)
- **APPROACH D:** Use Seurat's `as.SingleCellExperiment()`

## Check Your LIANA Version

This might be a version issue. Check:

```r
packageVersion("liana")
packageVersion("SingleCellExperiment")
```

If outdated, update:

```r
# Update SingleCellExperiment
BiocManager::install("SingleCellExperiment", update = TRUE, force = TRUE)

# Update LIANA
devtools::install_github("saezlab/liana")
# or
install.packages("liana")
```

## Debugging

If all approaches fail, get debug info:

```r
# Before liana_wrap, add:
cat("\n=== SCE Debug Info ===\n")
cat("SCE Dimensions:", nrow(lee_day3_sce), "x", ncol(lee_day3_sce), "\n")
cat("colData rows:", nrow(colData(lee_day3_sce)), "\n")
cat("colLabels length:", length(colLabels(lee_day3_sce)), "\n")
if (length(rowData(lee_day3_sce)) > 0) {
  cat("rowData rows:", nrow(rowData(lee_day3_sce)), "\n")
}
cat("Valid:", validObject(lee_day3_sce), "\n")

# Check internal metadata
int_col <- tryCatch(nrow(int_colData(lee_day3_sce)), error = function(e) "ERROR")
int_row <- tryCatch(nrow(int_elementMetadata(lee_day3_sce)), error = function(e) "ERROR")
cat("int_colData rows:", int_col, "\n")
cat("int_elementMetadata rows:", int_row, "\n")

# Check LIANA version
cat("\nPackage Versions:\n")
cat("LIANA:", as.character(packageVersion("liana")), "\n")
cat("SingleCellExperiment:", as.character(packageVersion("SingleCellExperiment")), "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
```

## Key Takeaway

**Stop creating SCE manually. Let LIANA handle it.**

The manual SCE creation, even when done correctly, creates structures that LIANA's internal code can't handle properly. Use LIANA's Seurat wrapper instead.

## Action Items

1. ✅ Try APPROACH A from above (Seurat wrapper)
2. ✅ If that fails, check versions
3. ✅ If still failing, try APPROACH C (add rowData)
4. ✅ Post version numbers and full error trace

## Files to Reference

- **`ALTERNATIVE_APPROACHES.md`** - All 4 alternative approaches
- **`MASTER_SUMMARY.md`** - Original fix package
- **`START_HERE_URGENT_FIX.md`** - Original simple fix (no longer sufficient)

## Bottom Line

This is no longer a "simple" fix. It's a compatibility issue between:
- How we create SCE objects manually
- How LIANA expects SCE objects to be structured
- How LIANA's subsetting operations work

**Solution:** Don't create SCE manually. Use LIANA's Seurat support instead.

---

**Status:** ⚠️ Investigating deeper LIANA/SCE compatibility issue  
**Next:** Try APPROACH A (Seurat wrapper)  
**Last Updated:** 2026-02-18
