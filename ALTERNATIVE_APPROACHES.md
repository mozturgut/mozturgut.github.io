# ALTERNATIVE FIX - LIANA Still Failing with SIMPLE Approach

## Problem Update

The SIMPLE approach was applied but the error STILL occurs:
```
Error: 'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'
Error: 'nrow' of 'int_colData' not equal to 'ncol(object)'
```

## What We Learned

1. ✅ SCE is created correctly (`validObject` returns TRUE)
2. ✅ Dimensions are aligned when we create it
3. ❌ LIANA's internal filtering STILL breaks it
4. Root cause: **LIANA's subsetting code doesn't properly update SCE metadata**

## New Solution: Try Different Approaches

### APPROACH A: Use LIANA's Seurat Wrapper Directly

LIANA has built-in support for Seurat objects. Try this:

```r
# Lee Day 3: Use LIANA's built-in Seurat support (EASIEST)
print("--- LIANA Analysis: Lee Day 3 ---")

# Prepare labels (same as before)
lee_day3_liana_labels <- lee_day3_pruned_labels
lee_day3_liana_labels[lee_day3_pos] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1pos")
lee_day3_liana_labels[lee_day3_neg] <- paste0(NEUTROPHIL_BASE_LABEL, "Arg1neg")
lee_day3_liana_labels[is.na(lee_day3_liana_labels) | lee_day3_liana_labels == ""] <- "Other"

# Add labels to Seurat metadata
lee_day3$liana_labels <- lee_day3_liana_labels
Seurat::Idents(lee_day3) <- factor(lee_day3_liana_labels)

# Load custom LR resource (same as before)
# ... [keep the CellChat/CellCall loading code] ...

# Run LIANA directly on Seurat object
lee_day3_liana_args <- list(
  seurat_object = lee_day3,
  ident_col = "liana_labels",  # Use the column we just added
  resource = lee_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)

if (!is.null(lee_day3_liana_external)) {
  lee_day3_liana_args$external_resource <- lee_day3_liana_external
}

# Let LIANA handle the Seurat -> SCE conversion
lee_day3_liana_result <- suppressWarnings(
  do.call(liana::liana_wrap, lee_day3_liana_args)
)
```

### APPROACH B: Different SCE Creation Method

If Approach A doesn't work, try creating SCE without setting colLabels separately:

```r
# Lee Day 3: Build SCE (Alternative method)
SeuratObject::DefaultAssay(lee_day3) <- "RNA"
if (inherits(lee_day3[["RNA"]], "Assay5")) {
  lee_day3[["RNA"]] <- SeuratObject::JoinLayers(lee_day3[["RNA"]])
}

lee_day3_cnt <- SeuratObject::GetAssayData(lee_day3, layer = "counts")
lee_day3_dat <- SeuratObject::GetAssayData(lee_day3, layer = "data")

# Create colData with additional metadata
lee_day3_coldata <- S4Vectors::DataFrame(
  label = factor(lee_day3_liana_labels),
  cell_name = colnames(lee_day3_cnt),
  row.names = colnames(lee_day3_cnt)
)

# Create SCE
lee_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = lee_day3_cnt, logcounts = lee_day3_dat),
  colData = lee_day3_coldata
)

# DON'T set colLabels separately - let it use colData$label
# SingleCellExperiment::colLabels(lee_day3_sce) <- factor(lee_day3_liana_labels)  # REMOVE THIS

# Verify
stopifnot(validObject(lee_day3_sce))
stopifnot(ncol(lee_day3_sce) == nrow(lee_day3_coldata))

# Run LIANA
lee_day3_liana_args <- list(
  sce = lee_day3_sce,
  resource = lee_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(lee_day3_liana_external)) {
  lee_day3_liana_args$external_resource <- lee_day3_liana_external
}

lee_day3_liana_result <- suppressWarnings(
  do.call(liana::liana_wrap, lee_day3_liana_args)
)
```

### APPROACH C: Add rowData to SCE

The error mentions `int_elementMetadata` (gene metadata). Try adding rowData:

```r
# Lee Day 3: Build SCE with rowData
SeuratObject::DefaultAssay(lee_day3) <- "RNA"
if (inherits(lee_day3[["RNA"]], "Assay5")) {
  lee_day3[["RNA"]] <- SeuratObject::JoinLayers(lee_day3[["RNA"]])
}

lee_day3_cnt <- SeuratObject::GetAssayData(lee_day3, layer = "counts")
lee_day3_dat <- SeuratObject::GetAssayData(lee_day3, layer = "data")

# Create rowData for genes
lee_day3_rowdata <- S4Vectors::DataFrame(
  gene_name = rownames(lee_day3_cnt),
  row.names = rownames(lee_day3_cnt)
)

# Create colData for cells
lee_day3_coldata <- S4Vectors::DataFrame(
  label = factor(lee_day3_liana_labels),
  row.names = colnames(lee_day3_cnt)
)

# Create SCE with both rowData and colData
lee_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = lee_day3_cnt, logcounts = lee_day3_dat),
  rowData = lee_day3_rowdata,
  colData = lee_day3_coldata
)

# Set colLabels to match colData$label
SingleCellExperiment::colLabels(lee_day3_sce) <- lee_day3_coldata$label

# Verify
stopifnot(validObject(lee_day3_sce))

# Run LIANA
lee_day3_liana_args <- list(
  sce = lee_day3_sce,
  resource = lee_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(lee_day3_liana_external)) {
  lee_day3_liana_args$external_resource <- lee_day3_liana_external
}

lee_day3_liana_result <- suppressWarnings(
  do.call(liana::liana_wrap, lee_day3_liana_args)
)
```

### APPROACH D: Use as.SingleCellExperiment from Seurat

Use Seurat's built-in conversion:

```r
# Lee Day 3: Use Seurat's SCE conversion
SeuratObject::DefaultAssay(lee_day3) <- "RNA"

# Add labels to metadata first
lee_day3$liana_labels <- lee_day3_liana_labels
Seurat::Idents(lee_day3) <- factor(lee_day3_liana_labels)

# Convert using Seurat's method
lee_day3_sce <- tryCatch({
  # Try Seurat's conversion if available
  if ("as.SingleCellExperiment" %in% methods(class = "Seurat")) {
    as.SingleCellExperiment(lee_day3)
  } else {
    # Fallback: manual creation
    cnt <- GetAssayData(lee_day3, layer = "counts")
    dat <- GetAssayData(lee_day3, layer = "data")
    
    SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = cnt, logcounts = dat),
      colData = lee_day3@meta.data
    )
  }
}, error = function(e) {
  # If conversion fails, create manually
  cnt <- GetAssayData(lee_day3, layer = "counts")
  dat <- GetAssayData(lee_day3, layer = "data")
  
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = cnt, logcounts = dat),
    colData = S4Vectors::DataFrame(lee_day3@meta.data)
  )
})

# Set colLabels from the metadata
SingleCellExperiment::colLabels(lee_day3_sce) <- factor(lee_day3$liana_labels)

# Verify
stopifnot(validObject(lee_day3_sce))

# Run LIANA
lee_day3_liana_args <- list(
  sce = lee_day3_sce,
  resource = lee_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(lee_day3_liana_external)) {
  lee_day3_liana_args$external_resource <- lee_day3_liana_external
}

lee_day3_liana_result <- suppressWarnings(
  do.call(liana::liana_wrap, lee_day3_liana_args)
)
```

## Recommended Order to Try

1. **APPROACH A** (Easiest) - Let LIANA handle Seurat directly
2. **APPROACH C** - Add rowData to SCE
3. **APPROACH B** - Different colData structure
4. **APPROACH D** - Use Seurat's SCE conversion

## Debugging Steps

If all approaches fail, add this debugging code before `liana_wrap`:

```r
# Debug SCE structure
cat("SCE Debug Info:\n")
cat("  Dimensions:", nrow(lee_day3_sce), "genes x", ncol(lee_day3_sce), "cells\n")
cat("  colData rows:", nrow(colData(lee_day3_sce)), "\n")
cat("  colLabels length:", length(colLabels(lee_day3_sce)), "\n")
cat("  rowData rows:", nrow(rowData(lee_day3_sce)), "\n")
cat("  int_colData rows:", nrow(int_colData(lee_day3_sce)), "\n")
cat("  int_elementMetadata rows:", nrow(int_elementMetadata(lee_day3_sce)), "\n")
cat("  Valid:", validObject(lee_day3_sce), "\n")
```

## Alternative: Different LIANA Version

This might be a version compatibility issue. Check:

```r
packageVersion("liana")
packageVersion("SingleCellExperiment")
packageVersion("Seurat")
```

If using old versions, consider updating:
```r
BiocManager::install("SingleCellExperiment", update = TRUE)
devtools::install_github("saezlab/liana")
```

## Bottom Line

The issue is deeper than expected. LIANA's internal operations are incompatible with the SCE we're creating. Try **APPROACH A first** - let LIANA handle the Seurat object directly without manual SCE creation.
