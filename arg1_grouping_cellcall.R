# Simple CellCall analysis helper for Arg1-positive vs Arg1-negative neutrophils
# or for timepoint comparison.
#
# Usage:
#   result <- run_arg1_cellcall(Neu_integrated, mode = "status", timepoint = "3")
#   result_time <- run_arg1_cellcall(Neu_integrated, mode = "time")

library(Seurat)
library(cellcall)

run_arg1_cellcall <- function(seurat_obj, mode = c("status", "time"), timepoint = "3", arg1_gene = "Arg1") {
  mode <- match.arg(mode)
  DefaultAssay(seurat_obj) <- "RNA"

  if (mode == "status") {
    subset_obj <- subset(seurat_obj, subset = time == timepoint)
    arg_expr <- FetchData(subset_obj, vars = arg1_gene, layer = "data")[, 1]
    subset_obj$Arg1_status <- ifelse(arg_expr > 0, "Arg1_positive", "Arg1_negative")
    Idents(subset_obj) <- "Arg1_status"
  } else {
    subset_obj <- seurat_obj
    Idents(subset_obj) <- "time"
  }

  expr <- GetAssayData(subset_obj, slot = "data")
  group <- as.character(Idents(subset_obj))
  colnames(expr) <- paste0(colnames(expr), "_", group)

  cc_obj <- CreateNichConObject(
    data = as.data.frame(as.matrix(expr)),
    min.feature = 3,
    names.field = 2,
    names.delim = "_",
    source = "TPM",
    scale.factor = 1e6,
    Org = "Mus musculus",
    project = "Arg1_Grouping"
  )

  cc_obj <- TransCommuProfile(
    object = cc_obj,
    pValueCor = 0.05,
    CorValue = 0.1,
    topTargetCor = 1,
    p.adjust = 0.05,
    use.type = "median",
    probs = 0.9,
    method = "weighted",
    IS_core = TRUE,
    Org = "Mus musculus"
  )

  cc_obj
}
