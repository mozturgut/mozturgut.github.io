# COPY-PASTE FIX TEMPLATE FOR ALL 5 LIANA SECTIONS

## Instructions

For each of the 5 sections (Lee Day 1, Lee Day 3, Wang Day 3, Qin Day 1, Qin Day 3):

1. Find the comment: `# {Dataset} Day {N}: Build SCE`
2. Delete everything from that comment until just before `lee_day3_liana_result <- ...` (or equivalent for other datasets)
3. Replace with the template below, changing `{dataset}` to the appropriate prefix

## Universal Template (Works for All 5 Sections)

```r
# {Dataset} Day {N}: Build SCE (SIMPLE APPROACH - Let LIANA filter)
SeuratObject::DefaultAssay({dataset}_day{N}) <- "RNA"
if (inherits({dataset}_day{N}[["RNA"]], "Assay5")) {
  {dataset}_day{N}[["RNA"]] <- SeuratObject::JoinLayers({dataset}_day{N}[["RNA"]])
}

# Extract matrices (NO FILTERING)
{dataset}_day{N}_cnt <- SeuratObject::GetAssayData({dataset}_day{N}, layer = "counts")
{dataset}_day{N}_dat <- SeuratObject::GetAssayData({dataset}_day{N}, layer = "data")

# Verify dimensions
stopifnot(ncol({dataset}_day{N}_cnt) == ncol({dataset}_day{N}_dat))
stopifnot(ncol({dataset}_day{N}_cnt) == length({dataset}_day{N}_liana_labels))

# Create SCE with ALL data (no pre-filtering)
{dataset}_day{N}_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = {dataset}_day{N}_cnt, logcounts = {dataset}_day{N}_dat),
  colData = S4Vectors::DataFrame(label = factor({dataset}_day{N}_liana_labels))
)
SingleCellExperiment::colLabels({dataset}_day{N}_sce) <- factor({dataset}_day{N}_liana_labels)

# Verify SCE is valid
stopifnot(validObject({dataset}_day{N}_sce))

# Run LIANA (it handles all filtering)
{dataset}_day{N}_liana_args <- list(
  sce = {dataset}_day{N}_sce,
  resource = {dataset}_day{N}_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null({dataset}_day{N}_liana_external)) {
  {dataset}_day{N}_liana_args$external_resource <- {dataset}_day{N}_liana_external
}
```

## Specific Replacements

### Lee Day 1
Replace `{dataset}` with `lee` and `{N}` with `1`:
```r
# Lee Day 1: Build SCE (SIMPLE APPROACH - Let LIANA filter)
SeuratObject::DefaultAssay(lee_day1) <- "RNA"
if (inherits(lee_day1[["RNA"]], "Assay5")) {
  lee_day1[["RNA"]] <- SeuratObject::JoinLayers(lee_day1[["RNA"]])
}

lee_day1_cnt <- SeuratObject::GetAssayData(lee_day1, layer = "counts")
lee_day1_dat <- SeuratObject::GetAssayData(lee_day1, layer = "data")

stopifnot(ncol(lee_day1_cnt) == ncol(lee_day1_dat))
stopifnot(ncol(lee_day1_cnt) == length(lee_day1_liana_labels))

lee_day1_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = lee_day1_cnt, logcounts = lee_day1_dat),
  colData = S4Vectors::DataFrame(label = factor(lee_day1_liana_labels))
)
SingleCellExperiment::colLabels(lee_day1_sce) <- factor(lee_day1_liana_labels)

stopifnot(validObject(lee_day1_sce))

lee_day1_liana_args <- list(
  sce = lee_day1_sce,
  resource = lee_day1_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(lee_day1_liana_external)) {
  lee_day1_liana_args$external_resource <- lee_day1_liana_external
}
```

### Lee Day 3
Replace `{dataset}` with `lee` and `{N}` with `3`:
```r
# Lee Day 3: Build SCE (SIMPLE APPROACH - Let LIANA filter)
SeuratObject::DefaultAssay(lee_day3) <- "RNA"
if (inherits(lee_day3[["RNA"]], "Assay5")) {
  lee_day3[["RNA"]] <- SeuratObject::JoinLayers(lee_day3[["RNA"]])
}

lee_day3_cnt <- SeuratObject::GetAssayData(lee_day3, layer = "counts")
lee_day3_dat <- SeuratObject::GetAssayData(lee_day3, layer = "data")

stopifnot(ncol(lee_day3_cnt) == ncol(lee_day3_dat))
stopifnot(ncol(lee_day3_cnt) == length(lee_day3_liana_labels))

lee_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = lee_day3_cnt, logcounts = lee_day3_dat),
  colData = S4Vectors::DataFrame(label = factor(lee_day3_liana_labels))
)
SingleCellExperiment::colLabels(lee_day3_sce) <- factor(lee_day3_liana_labels)

stopifnot(validObject(lee_day3_sce))

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
```

### Wang Day 3
Replace `{dataset}` with `wang` and `{N}` with `3`:
```r
# Wang Day 3: Build SCE (SIMPLE APPROACH - Let LIANA filter)
SeuratObject::DefaultAssay(wang_day3) <- "RNA"
if (inherits(wang_day3[["RNA"]], "Assay5")) {
  wang_day3[["RNA"]] <- SeuratObject::JoinLayers(wang_day3[["RNA"]])
}

wang_day3_cnt <- SeuratObject::GetAssayData(wang_day3, layer = "counts")
wang_day3_dat <- SeuratObject::GetAssayData(wang_day3, layer = "data")

stopifnot(ncol(wang_day3_cnt) == ncol(wang_day3_dat))
stopifnot(ncol(wang_day3_cnt) == length(wang_day3_liana_labels))

wang_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = wang_day3_cnt, logcounts = wang_day3_dat),
  colData = S4Vectors::DataFrame(label = factor(wang_day3_liana_labels))
)
SingleCellExperiment::colLabels(wang_day3_sce) <- factor(wang_day3_liana_labels)

stopifnot(validObject(wang_day3_sce))

wang_day3_liana_args <- list(
  sce = wang_day3_sce,
  resource = wang_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(wang_day3_liana_external)) {
  wang_day3_liana_args$external_resource <- wang_day3_liana_external
}
```

### Qin Day 1
Replace `{dataset}` with `qin` and `{N}` with `1`:
```r
# Qin Day 1: Build SCE (SIMPLE APPROACH - Let LIANA filter)
SeuratObject::DefaultAssay(qin_day1) <- "RNA"
if (inherits(qin_day1[["RNA"]], "Assay5")) {
  qin_day1[["RNA"]] <- SeuratObject::JoinLayers(qin_day1[["RNA"]])
}

qin_day1_cnt <- SeuratObject::GetAssayData(qin_day1, layer = "counts")
qin_day1_dat <- SeuratObject::GetAssayData(qin_day1, layer = "data")

stopifnot(ncol(qin_day1_cnt) == ncol(qin_day1_dat))
stopifnot(ncol(qin_day1_cnt) == length(qin_day1_liana_labels))

qin_day1_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = qin_day1_cnt, logcounts = qin_day1_dat),
  colData = S4Vectors::DataFrame(label = factor(qin_day1_liana_labels))
)
SingleCellExperiment::colLabels(qin_day1_sce) <- factor(qin_day1_liana_labels)

stopifnot(validObject(qin_day1_sce))

qin_day1_liana_args <- list(
  sce = qin_day1_sce,
  resource = qin_day1_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(qin_day1_liana_external)) {
  qin_day1_liana_args$external_resource <- qin_day1_liana_external
}
```

### Qin Day 3
Replace `{dataset}` with `qin` and `{N}` with `3`:
```r
# Qin Day 3: Build SCE (SIMPLE APPROACH - Let LIANA filter)
SeuratObject::DefaultAssay(qin_day3) <- "RNA"
if (inherits(qin_day3[["RNA"]], "Assay5")) {
  qin_day3[["RNA"]] <- SeuratObject::JoinLayers(qin_day3[["RNA"]])
}

qin_day3_cnt <- SeuratObject::GetAssayData(qin_day3, layer = "counts")
qin_day3_dat <- SeuratObject::GetAssayData(qin_day3, layer = "data")

stopifnot(ncol(qin_day3_cnt) == ncol(qin_day3_dat))
stopifnot(ncol(qin_day3_cnt) == length(qin_day3_liana_labels))

qin_day3_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = qin_day3_cnt, logcounts = qin_day3_dat),
  colData = S4Vectors::DataFrame(label = factor(qin_day3_liana_labels))
)
SingleCellExperiment::colLabels(qin_day3_sce) <- factor(qin_day3_liana_labels)

stopifnot(validObject(qin_day3_sce))

qin_day3_liana_args <- list(
  sce = qin_day3_sce,
  resource = qin_day3_liana_resource,
  expr_prop = 0.05,
  verbose = TRUE,
  min_cells = 0
)
if (!is.null(qin_day3_liana_external)) {
  qin_day3_liana_args$external_resource <- qin_day3_liana_external
}
```

## What Gets Deleted

In each section, delete these lines:
- `if (is.null({dataset}_day{N}_cnt) || ...` matrix validation
- `if (!inherits({dataset}_day{N}_cnt, ...` sparse matrix conversion
- `{dataset}_day{N}_original_ncells <- ...` index tracking
- `{dataset}_day{N}_cell_indices <- ...` index initialization
- `{dataset}_day{N}_genes_in_data <- ...` gene filtering
- `{dataset}_day{N}_keep_genes <- ...` gene keep filter
- `{dataset}_day{N}_cnt <- {dataset}_day{N}_cnt[{dataset}_day{N}_keep_genes, ...` gene filtering
- `{dataset}_day{N}_dat <- {dataset}_day{N}_dat[{dataset}_day{N}_keep_genes, ...` gene filtering
- `{dataset}_day{N}_nonzero_cells <- ...` cell zero filter
- `{dataset}_day{N}_cnt <- {dataset}_day{N}_cnt[, {dataset}_day{N}_nonzero_cells, ...` cell filtering
- `{dataset}_day{N}_dat <- {dataset}_day{N}_dat[, {dataset}_day{N}_nonzero_cells, ...` cell filtering
- `{dataset}_day{N}_cell_indices <- {dataset}_day{N}_cell_indices[{dataset}_day{N}_nonzero_cells]` index update
- `{dataset}_day{N}_nonzero_genes <- ...` gene zero filter
- `{dataset}_day{N}_cnt <- {dataset}_day{N}_cnt[{dataset}_day{N}_nonzero_genes, ...` gene re-filtering
- `{dataset}_day{N}_dat <- {dataset}_day{N}_dat[{dataset}_day{N}_nonzero_genes, ...` gene re-filtering
- `{dataset}_day{N}_lab <- factor({dataset}_day{N}_liana_labels[{dataset}_day{N}_cell_indices])` tracked labels

## Summary

**Delete ~20 lines of filtering code per section.**
**Replace with ~15 lines of simple SCE creation.**
**Total: 5 sections × 20 lines = 100 lines deleted, 75 lines added.**
**Net: 25 fewer lines, much more reliable.**

## After Applying

Your code will:
1. ✅ Be simpler and easier to maintain
2. ✅ Work reliably with LIANA's internal filtering
3. ✅ Have no dimension mismatch errors
4. ✅ Be consistent across all 5 sections
5. ✅ Follow the recommended LIANA usage pattern
