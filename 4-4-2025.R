###############################################################################
#                               Load libraries
###############################################################################
library(dplyr)
library(Seurat)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(celldex)
library(R.matlab)
library(scater)
library(scran)
library(gt)
library(topGO)
library(org.Mm.eg.db)
library(SoupX)
library(DoubletFinder)
library(patchwork)
###############################################################################
#                             Load Directories
###############################################################################
# Get file names in directory for analysis
setwd("G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples")
tenXDat <- list.dirs(recursive = FALSE)

data_dirs_Hou <- c(
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/Uninjured",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/0.5d",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/1d",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/3d",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/7d",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/14d",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/60d",
  "G:/Mustafa/singlecell-2/SRA to seurat/Hou/Samples/90d"
)

sample_names <- c(
  "Uninjured",
  "0.5d",
  "1d",
  "3d",
  "7d",
  "14d",
  "60d",
  "90d"
)

# Read in all files and turn into seurat objects
matrices_Hou <- lapply(data_dirs_Hou, Read10X)
seuObj_Hou <- lapply(matrices_Hou, CreateSeuratObject)

# Immediately add the sample name to each object's metadata.
for (i in seq_along(seuObj_Hou)) {
  if (!is.null(seuObj_Hou[[i]])) {
    seuObj_Hou[[i]]@meta.data$sample <- sample_names[i]
  }
}
###############################################################################
#                                   Doublet Detection
###############################################################################
# Backup the original Seurat objects before any doublet detection modifications.
seuObj_Hou_preDoublet <- lapply(seuObj_Hou, identity)

## 4a. Check and Recalculate Doublet Detection (if needed)

# Check for doublet detection column ("DF.classifications") and recalculate if needed
seuObj_Hou <- lapply(seuObj_Hou, function(x) {
  if (is.null(x)) return(NULL)
  if (!("DF.classifications" %in% colnames(x@meta.data))) {
    cat(sprintf("Object %s does not have doublet annotations. Restoring from backup and recalculating...\n", 
                unique(x@meta.data$sample)))
    
    # Find corresponding backup object based on sample name
    backup_index <- which(sapply(seuObj_Hou_preDoublet, function(y) {
      !is.null(y) && (unique(y@meta.data$sample) == unique(x@meta.data$sample))
    }))
    if(length(backup_index) == 0){
      warning("No valid backup found; skipping doublet recalculation.")
      return(x)
    }
    x <- seuObj_Hou_preDoublet[[backup_index[1]]]
    
    # Set default assay and compute necessary counts if missing
    DefaultAssay(x) <- "RNA"
    if (!("nCount_RNA" %in% colnames(x@meta.data))) {
      x[["nCount_RNA"]] <- Matrix::colSums(GetAssayData(x, slot = "counts"))
    }
    if (!("nFeature_RNA" %in% colnames(x@meta.data))) {
      x[["nFeature_RNA"]] <- Matrix::colSums(GetAssayData(x, slot = "counts") > 0)
    }
    
    # More explicit normalization/scaling steps with verbose = FALSE.
    x <- NormalizeData(x, verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(features = VariableFeatures(x), verbose = FALSE)
    
    # Run parameter sweep for DoubletFinder
    sweep.res <- paramSweep(x, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    nExp <- round(0.075 * ncol(x))
    
    # Re-run doubletFinder to add doublet classification column
    x <- doubletFinder(
      x,
      PCs = 1:10,
      pN = 0.25,
      pK = optimal_pK,
      nExp = nExp,
      reuse.pANN = FALSE,
      sct = FALSE
    )
  }
  return(x)
})

# Loop through each Seurat object in seuObj_Hou_with_doublets
for (i in seq_along(seuObj_Hou_with_doublets)) {
  obj <- seuObj_Hou_with_doublets[[i]]
  
  if (!is.null(obj) && inherits(obj, "Seurat")) {
    # Check if the DF.classifications column exists and is non-empty.
    if ("DF.classifications" %in% colnames(obj@meta.data) && nrow(obj@meta.data) > 0) {
      # Create the doublet_clarification column based on DF.classifications value manually.
      obj@meta.data$doublet_clarification <- ifelse(
        as.character(obj@meta.data$DF.classifications) == "Doublet",
        "Doublet Detected",
        "Singlet"
      )
    } else {
      # If DF.classifications is missing or empty, assume all cells are singlets.
      obj@meta.data$doublet_clarification <- rep("Singlet", nrow(obj@meta.data))
    }
    
    # Update the object in the list.
    seuObj_Hou_with_doublets[[i]] <- obj
  } else {
    warning(sprintf("Skipping object %s: It is either NULL or not a Seurat object.", i))
  }
}

# Optionally, verify by printing a frequency table for the new column in each Seurat object.
for (i in seq_along(seuObj_Hou_with_doublets)) {
  obj <- seuObj_Hou_with_doublets[[i]]
  if (!is.null(obj) && inherits(obj, "Seurat")) {
    cat(sprintf("Seurat Object %s (Sample: %s) doublet_clarification summary:\n",
                i, sample_names[i]))
    print(table(obj@meta.data$doublet_clarification))
    cat("\n")
  }
}
###############################################################################
# Report the doublet detection results without altering the backup object.
seuObj_Hou_with_doublets <- lapply(seq_along(seuObj_Hou), function(i) {
  if (!is.null(seuObj_Hou[[i]])) {
    cat("Seurat Object", i, "Doublet Detection Results:\n")
    current_table <- table(as.character(seuObj_Hou[[i]]@meta.data$DF.classifications))
    print(current_table)
    # Return the original Seurat object (not the printed table)
    return(seuObj_Hou[[i]])
  } else {
    return(NULL)
  }
})

# Save the backup object for later use.
saveRDS(seuObj_Hou_with_doublets, "seuObj_Hou_with_doublets.rds")
cat("Completed reporting of singlets and doublets.\n\n")
#################################################################Check metadata

# Load the RDS file containing the Seurat objects with doublet information.
seuObj_Hou_with_doublets <- readRDS("seuObj_Hou_with_doublets.rds")

cat("Starting metadata verification...\n\n")
for (i in seq_along(seuObj_Hou_with_doublets)) {
  obj <- seuObj_Hou_with_doublets[[i]]
  
  cat("======================================================\n")
  # Use the corresponding sample name if available.
  sample_label <- if (i <= length(sample_names)) sample_names[i] else paste("Element", i)
  cat("Seurat Object:", sample_label, "\n")
  
  if (is.null(obj)) {
    cat("This Seurat object is NULL.\n")
    next
  }
  
  if (!inherits(obj, "Seurat")) {
    cat("Warning: This element is not a Seurat object. Its class is:", class(obj), "\n")
    next
  }
  
  cat("Metadata Columns:\n")
  print(colnames(obj@meta.data))
  
  cat("Metadata Preview:\n")
  print(head(obj@meta.data))
  
  if ("sample" %in% colnames(obj@meta.data)) {
    cat("Unique sample values:\n")
    print(unique(obj@meta.data$sample))
  } else {
    cat("Warning: 'sample' column not found in metadata!\n")
  }
  
  doublet_cols <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(doublet_cols) > 0) {
    cat("Doublet classification column(s) found:\n")
    print(doublet_cols)
    for (col in doublet_cols) {
      cat("Frequency table for", col, ":\n")
      print(table(as.character(obj@meta.data[[col]])))
    }
  } else {
    cat("Warning: No doublet classification column found in metadata!\n")
  }
  cat("======================================================\n\n")
}
cat("Metadata verification completed.\n")
############################################################Doublet backup
# Perform doublet detection backup: print the doublet detection results and save the objects.
# This stores the original Seurat objects with the metadata column 'DF.classifications'.
seuObj_Hou_with_doublets <- lapply(seq_along(seuObj_Hou), function(i) {
  if (!is.null(seuObj_Hou[[i]])) {
    cat("Seurat Object", i, "Doublet classification table:\n")
    print(table(as.character(seuObj_Hou[[i]]@meta.data$DF.classifications)))
    return(seuObj_Hou[[i]])
  } else {
    return(NULL)
  }
})

# Save the backup object to an RDS file for later use.
saveRDS(seuObj_Hou_with_doublets, "seuObj_Hou_with_doublets.rds")
cat("Completed backup of singlets and doublets.\n\n")


###############################################################################
#                             Summary Table
###############################################################################
# Load the backup Seurat objects with doublet information
seuObj_Hou_with_doublets <- readRDS("seuObj_Hou_with_doublets.rds")


# Function to produce a summary table of singlet and doublet counts.
# This function searches for a doublet annotation column whose name starts with "DF.classifications".
reportSingletDoubletCounts <- function(object_list, sample_names) {
  summary_list <- list()
  
  for (i in seq_along(object_list)) {
    seurat_obj <- object_list[[i]]
    if (is.null(seurat_obj)) next
    
    # Use the sample name stored in metadata if available; otherwise, use the predefined sample name.
    if ("sample" %in% colnames(seurat_obj@meta.data)) {
      sample_val <- unique(seurat_obj@meta.data$sample)
      sample_name <- if (length(sample_val) == 0 || all(is.na(sample_val))) sample_names[i] else as.character(sample_val[1])
    } else {
      sample_name <- sample_names[i]
    }
    
    total_cells <- ncol(seurat_obj)
    
    # Detect the doublet annotation column using a regex to allow for varied column endings.
    df_cols <- grep("^DF\\.classifications", colnames(seurat_obj@meta.data), value = TRUE)
    if (length(df_cols) > 0) {
      # Use the first matching column by default
      df_col <- df_cols[1]
      freq <- table(as.character(seurat_obj@meta.data[[df_col]]))
      singlet_count <- if ("Singlet" %in% names(freq)) as.integer(freq["Singlet"]) else 0L
      doublet_count <- if ("Doublet" %in% names(freq)) as.integer(freq["Doublet"]) else 0L
    } else if ("doublet_class" %in% colnames(seurat_obj@meta.data)) {
      freq <- table(as.character(seurat_obj@meta.data$doublet_class))
      singlet_count <- if ("Singlet" %in% names(freq)) as.integer(freq["Singlet"]) else 0L
      doublet_count <- if ("Doublet" %in% names(freq)) as.integer(freq["Doublet"]) else 0L
    } else {
      # If no doublet annotation is present, assume all cells are singlets.
      singlet_count <- total_cells
      doublet_count <- 0L
    }
    
    summary_list[[i]] <- data.frame(
      Object_Index = i,
      Sample_Name = sample_name,
      Total_Cells = total_cells,
      Singlets = singlet_count,
      Doublets = doublet_count,
      stringsAsFactors = FALSE
    )
  }
  
  summary_table <- do.call(rbind, summary_list)
  # Order sample names according to the predefined order.
  summary_table$Sample_Name <- factor(summary_table$Sample_Name, levels = sample_names)
  summary_table <- summary_table[order(summary_table$Sample_Name), ]
  return(summary_table)
}

# Generate and print the summary table.
summary_table <- reportSingletDoubletCounts(seuObj_Hou_with_doublets, sample_names)
print(summary_table)

###############################################################################
#                                 Remove Doublets
###############################################################################

# Remove doublets from each Seurat object in the backup.
seuObj_Hou <- lapply(seuObj_Hou_with_doublets, function(x) {
  if (is.null(x)) return(NULL)
  
  # Find the doublet classification column by matching the pattern.
  df_cols <- grep("^DF\\.classifications", colnames(x@meta.data), value = TRUE)
  if (length(df_cols) > 0) {
    # Use the first matching column; update the object's metadata for convenience.
    df_col <- df_cols[1]
    x$doublet_class <- x@meta.data[[df_col]]
    # Subset the object to keep only cells labeled as singlets.
    x <- subset(x, subset = doublet_class == "Singlet")
  } else {
    # If no doublet column is found, assume all cells are singlets.
    warning("No DF.classifications column found; assuming all cells are singlets.")
  }
  
  return(x)
})

# Print the first few entries to verify changes.
head(seuObj_Hou)
print(seuObj_Hou)

# Create a variable for singlets to be used in downstream analyses.
seuObj_Hou_singlets <- seuObj_Hou

# Generate and print the summary table after doublet removal.
summary_table <- reportSingletDoubletCounts(seuObj_Hou_singlets, sample_names)
print(summary_table)

# Save the processed Seurat objects with singlets only.
saveRDS(seuObj_Hou_singlets, "seuObj_Hou_singlets.rds")
###############################################################################
#                                   QC
###############################################################################

# Iterate over each Seurat object in the list and set the default assay to "RNA"
seuObj_Hou <- lapply(seuObj_Hou, function(x) {
  DefaultAssay(x) <- "RNA"
  x
})

# Apply NormalizeData to each Seurat object in the list seuObj_Hou
seuObj_Hou <- lapply(seuObj_Hou, function(x) {
  NormalizeData(x)
})

# Assign mitochondrial gene percentages to cells
seuObj_Hou <- lapply(X = seuObj_Hou, FUN = function(x) {
  x[["percent_mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x
})

VlnPlot(seuObj_Hou[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3, pt.size = 0.01)

# Subset objects to expected expression levels
seuObj_Hou <- lapply(X = seuObj_Hou, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 10)
})

VlnPlot(seuObj_Hou[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3, pt.size = 0.01)



#####################################QC my code

# Generate QC violin plots and save files for each object
for (i in seq_along(seuObj_Hou)) {
  age <- seuObj_Hou[[i]]@meta.data$age[1]
  sample <- seuObj_Hou[[i]]@meta.data$sample[1]
  title <- paste0(age, " ", sample, " percent_mt")
  
  p1 <- VlnPlot(
    object = seuObj_Hou[[i]], 
    features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
    group.by = "sample",
    ncol = 3, 
    pt.size = 0.1
  ) + ggtitle(title)
  
  filename <- paste0("QC_Violin_", age, "_", sample, ".png")
  print(p1)
  ggsave(filename, plot = p1)
}

# Generate FeatureScatter plots and save as PNG files
for (i in seq_along(seuObj_Hou)) {
  age <- seuObj_Hou[[i]]@meta.data$age[1]
  sample <- seuObj_Hou[[i]]@meta.data$sample[1]
  title <- paste0(age, " ", sample, " Feature Scatter")
  
  p_feature1 <- FeatureScatter(
    object = seuObj_Hou[[i]], 
    feature1 = "nCount_RNA", 
    feature2 = "nFeature_RNA"
  ) + ggtitle(title)
  
  p_feature2 <- FeatureScatter(
    object = seuObj_Hou[[i]], 
    feature1 = "nCount_RNA", 
    feature2 = "percent_mt"
  ) + ggtitle(title)
  
  combinedfeatures <- p_feature1 + p_feature2
  
  filename <- paste0("FeatureScatter_", age, "_", sample, ".png")
  print(combinedfeatures)
  ggsave(filename, plot = combinedfeatures, width = 26.67, height = 13.33, dpi = 300)
}


###############################################################################
#                                     Integration
###############################################################################

# Normalize cells and find variable genes for each sample
seuObj_Hou <- lapply(seuObj_Hou, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x
})

seuObj_Hou <- lapply(seq_along(seuObj_Hou), function(i) {
  obj <- seuObj_Hou[[i]]
  obj <- RenameCells(obj, add.cell.id = sample_names[i])
  obj
})

# (Optional) Re-run normalization/feature detection
seuObj_Hou <- lapply(seuObj_Hou, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x
})

# --- Using consistent features when scaling to avoid mismatches ---
features <- SelectIntegrationFeatures(object.list = seuObj_Hou, nfeatures = 2000)
seuObj_Hou <- lapply(seuObj_Hou, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  x
})

seuObj_Hou <- Filter(Negate(is.null), seuObj_Hou)
anchors <- FindIntegrationAnchors(object.list = seuObj_Hou, anchor.features = features)
seuObj.integrated_Hou <- IntegrateData(anchorset = anchors)
DefaultAssay(seuObj.integrated_Hou) <- "integrated"
seuObj.integrated_Hou <- ScaleData(seuObj.integrated_Hou, features = features, verbose = FALSE)
seuObj.integrated_Hou <- RunPCA(seuObj.integrated_Hou, features = features, verbose = FALSE)

ElbowPlot(seuObj.integrated_Hou, ndims = 40)
DimHeatmap(seuObj.integrated_Hou, dims = 1:5, cells = 500, balanced = TRUE)
seuObj.integrated_Hou <- RunUMAP(seuObj.integrated_Hou, dims = 1:30)
DimPlot(seuObj.integrated_Hou)
seuObj.integrated_Hou <- FindNeighbors(seuObj.integrated_Hou, dims = 1:30)
seuObj.integrated_Hou <- FindClusters(seuObj.integrated_Hou)
DimPlot(seuObj.integrated_Hou, reduction = "umap", label = TRUE) + NoLegend()

FeaturePlot(seuObj.integrated_Hou, features = c("S100a9"))
FeaturePlot(seuObj.integrated_Hou, features = c("S100a8"))
FeaturePlot(seuObj.integrated_Hou, features = c("Ly6g"))
FeaturePlot(seuObj.integrated_Hou, features = c("Mmp9"))
FeaturePlot(seuObj.integrated_Hou, features = c("Ngp"))

unique_samples <- unique(seuObj.integrated_Hou@meta.data$sample)
cat("Unique sample values before setting factor levels:\n")
print(unique_samples)

saveRDS(seuObj.integrated_Hou, "seuObj_Hou_integrated.rds")
readRDS("seuObj_Hou_integrated.rds")
# Ensure sample column respects the specified order
seuObj.integrated_Hou@meta.data$sample <- factor(
  seuObj.integrated_Hou@meta.data$sample,
  levels = sample_names
)
###############################################################################
#                               Combined QC
###############################################################################
p2 <- VlnPlot(
  object = seuObj.integrated_Hou,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
  group.by = "sample",
  pt.size = 0.1,
  ncol = 3
) + ggtitle("QC Metrics Across All Samples")
print(p2)
ggsave("QC_Violin_AllSamples.png", plot = p2)

# Extract metadata from the integrated Seurat object
data_meta <- seuObj.integrated_Hou@meta.data

# Ensure the required columns exist
if (!all(c("nCount_RNA", "nFeature_RNA", "sample") %in% colnames(data_meta))) {
  stop("The required columns 'nCount_RNA', 'nFeature_RNA', and 'sample' are missing in metadata.")
}

# Proceed with plotting
p_combined <- ggplot(data_meta, aes(x = nCount_RNA, y = nFeature_RNA, color = sample)) +
  geom_point(alpha = 0.5) +
  ggtitle("Feature Scatter Plot: All Samples") +
  theme_bw() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
p_combined
# Save the plot
ggsave("FeatureScatter_AllSamples.png", plot = p_combined, width = 26.67, height = 13.33, dpi = 300)

###############################################################################
#                                 SingleR (Annotate)
###############################################################################
ref.se <- ImmGenData()
sce <- GetAssayData(seuObj.integrated_Hou, assay = "integrated", layer = "data")
pred <- SingleR(test = sce, ref = ref.se, labels = ref.se$label.main)

seuObj.integrated_Hou$cell_type <- pred$pruned.labels

table(subset(seuObj.integrated_Hou, cell_type == "Neutrophils")$seurat_clusters)

DimPlot(seuObj.integrated_Hou, group.by = "cell_type", label = TRUE, 
        label.size = 3)

##############################################################################
#                                 Neutrophils Isolated
##############################################################################
# Instead of passing the entire 'neutrophils' object, extract cell names from it.
Neutrophils_Hou <- subset(seuObj.integrated_Hou, subset = cell_type == "Neutrophils")
neutrophil_cells <- Cells(Neutrophils_Hou)
saveRDS(Neutrophils_Hou, "Neutrophils_Hou.rds")
# Create a UMAP plot highlighting the neutrophils cells in seuObj.integrated_Hou.
p_umap_specific <- DimPlot(
  object = seuObj.integrated_Hou,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  cells.highlight = neutrophil_cells  # Use a vector of cell names for highlighting
) + 
  ggtitle("UMAP with Neutrophils Highlighted") +
  # Manually set legend labels by specifying breaks and labels.
  scale_color_manual(
    breaks = c("unselected", "Group_1"),
    values = c("grey", "red"),
    labels = c("other", "Neutrophils")
  ) +
  labs(color = "Cell Group")

# Display and save the UMAP plot.
print(p_umap_specific)
ggsave("UMAP_Neutrophils_Highlighted.png", plot = p_umap_specific)

###############################################################
# Convert the sample column to a factor with the specified order
Neutrophils_Hou$sample <- factor(Neutrophils_Hou$sample, levels = sample_names)

table(Neutrophils_Hou$sample)

Neutrophils_Hou <- FindVariableFeatures(Neutrophils_Hou, selection.method = "vst", assay = "RNA")

VlnPlot(Neutrophils_Hou, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3, pt.size = 0.01, group.by = "sample")
###################################
# Get the unique sample names from the metadata.
samples <- unique(Neutrophils_Hou$sample)

# Loop through each sample, subset the object, create a scatter plot, and add a title.
plot_list <- lapply(samples, function(s) {
  obj_subset <- subset(Neutrophils_Hou, subset = sample == s)
  p <- FeatureScatter(obj_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle(paste("Sample:", s))
  return(p)
})

# Combine all scatter plots into one figure.
combined_plot <- wrap_plots(plot_list)
print(combined_plot)
########################
# first extract the metadata into a data frame.
data_meta_N <- Neutrophils_Hou@meta.data

# Proceed with plotting using ggplot2.
# The plot will display a feature scatter plot (nCount_RNA vs nFeature_RNA)
# and color the points by the 'sample' column.
library(ggplot2)

p_combined_N <- ggplot(data_meta_N, aes(x = nCount_RNA, y = nFeature_RNA, color = sample)) +
  geom_point(alpha = 0.5) +
  ggtitle("Feature Scatter Plot: All Samples") +
  theme_bw() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Display the plot.
p_combined_N
#######################################################3

top10 <- head(VariableFeatures(Neutrophils_Hou), 10)

plot1 <- VariableFeaturePlot(Neutrophils_Hou, assay = "RNA")
plot2 <- LabelPoints(plot = plot1, points = top10)
plot2

all.genes <- rownames(Neutrophils_Hou)
Neutrophils_Hou <- ScaleData(Neutrophils_Hou, features = all.genes)
Neutrophils_Hou <- RunPCA(Neutrophils_Hou, features = VariableFeatures(object = Neutrophils_Hou))
ElbowPlot(Neutrophils_Hou, ndims = 40)

Neutrophils_Hou <- RunUMAP(Neutrophils_Hou, dims = 1:15)

Neutrophils_Hou <- FindNeighbors(Neutrophils_Hou, dims = 1:15)
Neutrophils_Hou <- FindClusters(Neutrophils_Hou)

Idents(Neutrophils_Hou) <- 'sample'
DimPlot(Neutrophils_Hou, group.by = "sample")

###############################################################################
#                      Hou Split-By-Sex and Marker Analysis
###############################################################################

# Define gene lists for module scores.
bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", "Retnlg", "Lcn2", "Olfm4", "Chil3", "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", "Lyz1"))

earlyGenes_list <- list(c("Lyz2", "Ngp", "Wfdc21", "Camp", "Ifitm6", "Lcn2", "Ltf", "Arhgdib", "Mmp8", 
                          "Anxa1", "Prdx5", "Pglyrp1", "Cybb", "Ly6c2", "Chil3", "Dstn", "Cd177", 
                          "Pfn1", "Lgals3", "Ly6g", "Retnlg", "Adpgk", "Mgst1", "mt-Co1", "Serpinb1a", 
                          "Tkt", "Aldh2", "Capg", "Lamtor4", "Hmgn2"))

lateGenes_list  <- list(c("Rpl41", "Cebpb", "Fgl2", "Lst1", "Junb", "Ifitm1", "Hbb-bt", "Jund", 
                          "Rps9", "Ccl6", "Dusp1", "Csf3r", "H2-D1", "Hba-a2", "Ftl1", "Wfdc17", 
                          "Fau", "Fth1", "Hba-a1", "Tyrobp", "Srgn", "Msrb1", "Btg1", "Ifitm2", 
                          "Malat1", "Rps27", "Fxyd5", "Hbb-bs"))

library(Seurat)
library(ggplot2)
library(dplyr)

# Define a function to determine cell sex based on gene expression.
defineSex <- function(data, num) {
  curDat <- data[, num]
  Xist <- curDat[1]
  Ddx  <- curDat[2]
  
  sex <- dplyr::case_when(
    Xist > 0 & Ddx == 0 ~ "Female",
    Xist > 0 & Ddx > 0 ~ "Both",
    Xist == 0 & Ddx > 0 ~ "Male",
    Xist == 0 & Ddx == 0 ~ "Neither"
  )
  return(data.frame(number = num, sex = sex))
}

# Get unique sample names from the integrated object.
# (Assumes 'seuObj.integrated_Hou' is loaded in the environment.)
sample_list <- unique(seuObj.integrated_Hou$sample)

# Process each sample.
for (s in sample_list) {
  cat("Processing sample:", s, "\n")
  
  # Subset the integrated object for this sample.
  sample_obj <- subset(seuObj.integrated_Hou, subset = sample == s)
  
  # Optionally add a cell name field.
  sample_obj[["cellName"]] <- colnames(sample_obj)
  
  # Get expression values for Xist and Ddx3y.
  values <- GetAssayData(object = sample_obj, assay = "RNA", slot = "data")[c("Xist", "Ddx3y"), ]
  
  # Determine sex for each cell.
  sexes <- lapply(as.list(1:ncol(values)), FUN = defineSex, data = values)
  sexes <- bind_rows(sexes)
  
  # Add sex information to metadata; ensure it's a factor with levels "Male" and "Female".
  sample_obj$sex <- factor(sexes$sex, levels = c("Male", "Female"))
  
  # Subset to only Male and Female cells.
  sample_sex_labeled <- subset(sample_obj, subset = sex %in% c("Male", "Female"))
  
  # Re-run variable features, scaling, and PCA.
  sample_sex_labeled <- FindVariableFeatures(sample_sex_labeled, selection.method = "vst", nfeatures = 2000)
  sample_sex_labeled <- ScaleData(sample_sex_labeled, verbose = FALSE)
  sample_sex_labeled <- RunPCA(sample_sex_labeled, verbose = FALSE)
  ElbowPlot(sample_sex_labeled, ndims = 40)
  
  # Clustering and UMAP.
  sample_sex_labeled <- FindNeighbors(sample_sex_labeled, dims = 1:30)
  sample_sex_labeled <- FindClusters(sample_sex_labeled)
  sample_sex_labeled <- RunUMAP(sample_sex_labeled, dims = 1:30)
  
  # Set up output directory.
  base_path <- file.path("G:/Mustafa/singlecell-2/SRA to seurat/", s, "Samples", "Plots")
  dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
  
  # Save UMAP plots.
  tiff(filename = file.path(base_path, "timeUMAP.tiff"))
  DimPlot(sample_sex_labeled, group.by = "sample")
  dev.off()
  
  tiff(filename = file.path(base_path, "sexUMAP.tiff"))
  DimPlot(sample_sex_labeled, group.by = "sex")
  dev.off()
  
  # Create a heatmap of sex counts.
  sexesAtsamples <- data.frame(table(sample_sex_labeled$sample, sample_sex_labeled$sex))
  colnames(sexesAtsamples) <- c("sample", "Sex", "Count")
  tiff(filename = file.path(base_path, "sexsampleHeatmap.tiff"))
  sexesAtsamples %>% 
    filter(Count > 0) %>% 
    ggplot(aes(x = Sex, y = sample, fill = Count)) +
    geom_tile() +
    geom_text(aes(label = Count), colour = "white")
  dev.off()
  
  # Add module scores.
  sample_sex_labeled <- AddModuleScore(sample_sex_labeled, features = earlyGenes_list, search = TRUE, assay = "RNA", name = "neutrotimeEarly") 
  sample_sex_labeled <- AddModuleScore(sample_sex_labeled, features = lateGenes_list, search = TRUE, assay = "RNA", name = "neutrotimeLate")
  
  # Plot violin plots.
  tiff(filename = file.path(base_path, "sexLabeledEarlyNeutrotimeSplit.tiff"),
       units = "in", width = 5, height = 5, res = 300)
  print(
    VlnPlot(sample_sex_labeled, features = "neutrotimeEarly1", group.by = "sample", split.by = "sex",
            pt.size = 0, assay = "RNA") + 
      ggtitle("Early Neutrotime Scores") +
      xlab("Sample") +
      ylim(-1, 3.5) +
      geom_jitter(size = 0.00001, alpha = 0.2)
  )
  dev.off()
  
  tiff(filename = file.path(base_path, "sexLabeledLateNeutrotimeSplit.tiff"),
       units = "in", width = 5, height = 5, res = 300)
  print(
    VlnPlot(sample_sex_labeled, features = "neutrotimeLate1", group.by = "sample", split.by = "sex",
            pt.size = 0, assay = "RNA") + 
      ggtitle("Late Neutrotime Scores") +
      xlab("Sample") +
      ylim(-0.5, 3.5) +
      geom_jitter(size = 0.00001, alpha = 0.2)
  )
  dev.off()
  
  # Sex-specific marker analysis (only if both identities are present).
  Idents(sample_sex_labeled) <- sample_sex_labeled$sex
  current_idents <- levels(Idents(sample_sex_labeled))
  if(all(c("Male", "Female") %in% current_idents)) {
    sexMarkerMale <- FindMarkers(sample_sex_labeled, ident.1 = "Male", only.pos = TRUE)
    sexMarkerMale$pctDiff <- sexMarkerMale$pct.1 - sexMarkerMale$pct.2
    sexMarkerMale <- sexMarkerMale %>% arrange(desc(avg_log2FC))
    VlnPlot(sample_sex_labeled, features = rownames(sexMarkerMale)[4])
    
    geneListMale <- unlist(split(sexMarkerMale$p_val_adj, rownames(sexMarkerMale)))
    geneListMale <- geneListMale[geneListMale < 0.01]
    geneListMale <- names(geneListMale)
    write.csv(x = geneListMale, file = file.path(base_path, "maleGenes.csv"), 
              quote = FALSE, row.names = FALSE)
    
    sexMarkerFemale <- FindMarkers(sample_sex_labeled, ident.1 = "Female", only.pos = TRUE)
    sexMarkerFemale$pctDiff <- sexMarkerFemale$pct.1 - sexMarkerFemale$pct.2
    sexMarkerFemale <- sexMarkerFemale %>% arrange(desc(avg_log2FC))
    VlnPlot(sample_sex_labeled, features = rownames(sexMarkerFemale)[2])
    
    geneListFemale <- unlist(split(sexMarkerFemale$p_val_adj, rownames(sexMarkerFemale)))
    geneListFemale <- geneListFemale[geneListFemale < 0.01]
    geneListFemale <- names(geneListFemale)
    write.csv(x = geneListFemale, file = file.path(base_path, "femaleGenes.csv"), 
              quote = FALSE, row.names = FALSE)
  } else {
    cat("Skipping sex-specific marker analysis in sample", s, "because one or both sex identities are missing.\n")
  }
  
  # Comparative marker analysis by time point within each sex group.
  Idents(sample_sex_labeled) <- sample_sex_labeled$sample
  for (sexGroup in c("Male", "Female")) {
    cat("Performing comparative marker analysis for sex:", sexGroup, "\n")
    group_subset <- subset(sample_sex_labeled, subset = sex == sexGroup)
    
    # Check if the subset returned any cells.
    if (ncol(group_subset) == 0) {
      cat("No cells found for sex group:", sexGroup, "in sample", s, "\n")
      next
    }
    
    available_time_levels <- levels(Idents(group_subset))
    if(length(unique(Idents(group_subset))) < 2) {
      cat("Only one time identity present for", sexGroup, "in sample", s, "; skipping comparative marker analysis.\n")
    } else {
      for (time_point in available_time_levels) {
        cell_count <- table(Idents(group_subset))[time_point]
        if(!is.na(cell_count) && cell_count > 0){
          cat("Comparing time point:", time_point, "for", sexGroup, "\n")
          markers <- tryCatch(
            {
              FindMarkers(group_subset, ident.1 = time_point, only.pos = TRUE)
            },
            error = function(e){
              cat("Error in FindMarkers for", sexGroup, "at", time_point, ":", conditionMessage(e), "\n")
              return(NULL)
            }
          )
          if(!is.null(markers) && nrow(markers) > 0){
            markers$pctDiff <- markers$pct.1 - markers$pct.2
            markers <- markers %>% arrange(desc(avg_log2FC))
            cat("Top markers for", sexGroup, "at time", time_point, ":\n")
            print(head(markers))
            outFile <- file.path(base_path, paste0(sexGroup, "_", time_point, "_markers.csv"))
            write.csv(markers, file = outFile, quote = FALSE)
          } else {
            cat("No markers found for", sexGroup, "at time", time_point, "\n")
          }
        } else {
          cat("Time point", time_point, "has zero cells for", sexGroup, "\n")
        }
      }
    }
  }
}




# First, check if your Seurat object's metadata already contains a "sex" column.
print(colnames(seuObj.integrated_Hou@meta.data))

# If not, compute cell sex. For example, using Xist and Ddx3y expression:
values <- FetchData(seuObj.integrated_Hou, vars = c("Xist", "Ddx3y"))
head(values)

# Function to determine cell sex based on gene expression,
# using a threshold to account for low-level expression.
determineSex <- function(xist, ddx3y, threshold = 1) {
  if(xist >= threshold & ddx3y < threshold) {
    "Female"
  } else if(ddx3y >= threshold & xist < threshold) {
    "Male"
  } else if(xist >= threshold & ddx3y >= threshold) {
    "Both"
  } else {
    "Undetermined"
  }
}

# Fetch expression data for "Xist" and "rna_Ddx3y" (correct column for Ddx3y)
expression_values <- FetchData(seuObj.integrated_Hou, vars = c("Xist", "rna_Ddx3y"))

# Calculate cell sex for each cell using the new function.
# Adjust the threshold below if needed.
cell_sex <- sapply(1:nrow(expression_values), function(i) {
  xist_val <- as.numeric(expression_values[i, "Xist"])
  ddx_val  <- as.numeric(expression_values[i, "rna_Ddx3y"])
  determineSex(xist_val, ddx_val, threshold = 1)
})

# Add the computed sex information to the Seurat object's metadata.
seuObj.integrated_Hou <- AddMetaData(seuObj.integrated_Hou, metadata = cell_sex, col.name = "sex")

# Now, re-create the table of sex counts by sample.
sample_sex_counts <- table(seuObj.integrated_Hou$sample, seuObj.integrated_Hou$sex)
print(sample_sex_counts)

################################################################################
#                                   Neutrotime
################################################################################
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rstatix)

# Join layers and set default assay.
Neutrophils_Hou[["joined"]] <- JoinLayers(Neutrophils_Hou[["RNA"]])
DefaultAssay(Neutrophils_Hou) <- "joined"

# Add module scores for the bm gene set.
Neutrophils_Hou <- AddModuleScore(object = Neutrophils_Hou, features = bmGenes, search = TRUE, assay = "joined", name = "neutrotimebm", ctrl = 20, nbin = 10)


# Add module scores for early and late gene lists.
Neutrophils_Hou <- AddModuleScore(object = Neutrophils_Hou, features = earlyGenes_list, search = TRUE, assay = "joined", name = "neutrotimeEarly", ctrl = 20, nbin = 10)
Neutrophils_Hou <- AddModuleScore(object = Neutrophils_Hou, features = lateGenes_list, search = TRUE, assay = "joined", name = "neutrotimeLate", ctrl = 20, nbin = 10)


# Check if UMAP reduction exists; if not, compute UMAP using PCA reduction.
if (!"umap" %in% names(Neutrophils_Hou@reductions)) {
  if (!"pca" %in% names(Neutrophils_Hou@reductions)) {
    Neutrophils_Hou <- ScaleData(Neutrophils_Hou, verbose = FALSE)
    Neutrophils_Hou <- RunPCA(Neutrophils_Hou, verbose = FALSE)
  }
  Neutrophils_Hou <- RunUMAP(Neutrophils_Hou, reduction = "pca", dims = 1:10)
}

# UMAP plot for Early Neutrotime module score.
umap_early <- FeaturePlot(
  Neutrophils_Hou,
  features = "neutrotimeEarly1",
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5
) + ggtitle("UMAP: Early Neutrotime Module Score by Sample")
print(umap_early)

# UMAP plot for Late Neutrotime module score.
umap_late <- FeaturePlot(
  Neutrophils_Hou,
  features = "neutrotimeLate1",
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5
) + ggtitle("UMAP: Late Neutrotime Module Score")
print(umap_late)

# Combined UMAP plot blending early and late scores.
combined_umap <- FeaturePlot(
  Neutrophils_Hou,
  features = c("neutrotimeEarly1", "neutrotimeLate1"),
  reduction = "umap",
  blend = TRUE,
  label = FALSE,
  blend.threshold = 0.2,  # Adjust threshold as needed.
  pt.size = 0.5
) + ggtitle("Combined UMAP: Early (Red) and Late (Green) Neutrotime Module Scores")
print(combined_umap)


############################violin


# Define gene lists for module scores.

# Module genes for BM score
bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", "Retnlg", "Lcn2", "Olfm4", "Chil3", 
                  "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", "Lyz1"))

# Module genes for Early Neutrotime score
earlyGenes <- list(c("Lyz2", "Ngp", "Wfdc21", "Camp", "Ifitm6", "Lcn2", "Ltf", 
                     "Arhgdib", "Mmp8", "Anxa1", "Prdx5", "Pglyrp1", "Cybb", 
                     "Ly6c2", "Chil3", "Dstn", "Cd177", "Pfn1", "Lgals3", "Ly6g", 
                     "Retnlg", "Adpgk", "Mgst1", "mt-Co1", "Serpinb1a", "Tkt", 
                     "Aldh2", "Capg", "Lamtor4", "Hmgn2"))

# Module genes for Late Neutrotime score
lateGenes <- list(c("Rpl41", "Cebpb", "Fgl2", "Lst1", "Junb", "Ifitm1", "Hbb-bt", 
                    "Jund", "Rps9", "Ccl6", "Dusp1", "Csf3r", "H2-D1", "Hba-a2", 
                    "Ftl1", "Wfdc17", "Fau", "Fth1", "Hba-a1", "Tyrobp", "Srgn", 
                    "Msrb1", "Btg1", "Ifitm2", "Malat1", "Rps27", "Fxyd5", "Hbb-bs"))

# Read in Blood, BM, and Spleen objects.
Blood  <- readRDS("Blood.rds")
BM     <- readRDS("BM.rds")
Spleen <- readRDS("Spleen.rds")

# Define "sample" metadata if not already set.
Blood$sample   <- "Blood"
BM$sample      <- "BoneMarrow"
Spleen$sample  <- "Spleen"

# Rename cells to ensure unique cell names in each dataset
Blood   <- RenameCells(Blood, new.names = paste0("Blood_", colnames(Blood)))
BM      <- RenameCells(BM, new.names = paste0("BM_", colnames(BM)))
Spleen  <- RenameCells(Spleen, new.names = paste0("Spleen_", colnames(Spleen)))

# Merge the additional objects with your existing 'Neutrophils_Hou' object.
Neutrophils_Hou <- merge(Neutrophils_Hou, y = c(Blood, BM, Spleen))

# Define the desired order of samples.
sample_order <- c("Uninjured", "0.5d", "1d", "3d", 
                  "7d", "14d", "60d", "90d", 
                  "BoneMarrow", "Spleen", "Blood")
Neutrophils_Hou$sample <- factor(Neutrophils_Hou$sample, levels = sample_order)

# OPTIONAL: Re-run normalization and scaling after merging
Neutrophils_Hou <- NormalizeData(Neutrophils_Hou)
Neutrophils_Hou <- FindVariableFeatures(Neutrophils_Hou)
Neutrophils_Hou <- ScaleData(Neutrophils_Hou)
Neutrophils_Hou <- RunPCA(Neutrophils_Hou, npcs = 30)
Neutrophils_Hou <- RunUMAP(Neutrophils_Hou, dims = 1:20)

# Recalculate module scores to ensure all cells get a value.
Neutrophils_Hou <- AddModuleScore(Neutrophils_Hou, features = bmGenes, name = "neutrotimebm")
Neutrophils_Hou <- AddModuleScore(Neutrophils_Hou, features = earlyGenes, name = "neutrotimeEarly")
Neutrophils_Hou <- AddModuleScore(Neutrophils_Hou, features = lateGenes, name = "neutrotimeLate")

# Violin plots for bm, early, and late module scores.
p_bm <- VlnPlot(Neutrophils_Hou, features = "neutrotimebm1", group.by = "sample", pt.size = 0.1) + 
  ggtitle("bm Scores") + xlab("Time (days)") + ylim(c(-0.5, 4))
p_early <- VlnPlot(Neutrophils_Hou, features = "neutrotimeEarly1", group.by = "sample", pt.size = 0.1) + 
  ggtitle("Early Neutrotime Scores") + xlab("Time (days)") + ylim(c(-0.5, 3))
p_late <- VlnPlot(Neutrophils_Hou, features = "neutrotimeLate1", group.by = "sample", pt.size = 0.08) + 
  ggtitle("Late Neutrotime Scores") + xlab("Time (days)") + ylim(c(0.0, 3.0))

# Print the violin plots.
print(p_bm)
print(p_early)
print(p_late)



# Optional: Save UMAP and violin plots.
ggsave("Neutrotime_Early_Violin.png", plot = p_early, width = 8, height = 6, dpi = 300)
ggsave("Neutrotime_Late_Violin.png", plot = p_late, width = 8, height = 6, dpi = 300)
ggsave("UMAP_Early_Neutrotime.png", plot = umap_early, width = 8, height = 6, dpi = 300)
ggsave("UMAP_Late_Neutrotime.png", plot = umap_late, width = 8, height = 6, dpi = 300)

# Exclude BoneMarrow, Spleen, and Blood for UMAP plots.
cells_to_keep <- colnames(Neutrophils_Hou)[!(Neutrophils_Hou$sample %in% c("BoneMarrow", "Blood", "Spleen"))]
Idents(Neutrophils_Hou) <- "sample"
################################################################################
Neutrophils_Hou[["joined"]] <- JoinLayers(Neutrophils_Hou[["RNA"]])
DefaultAssay(Neutrophils_Hou) <- "joined"
Neu.Timemarkers <- FindAllMarkers(Neutrophils_Hou, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
Neu.Timemarkers$pct.diff <- Neu.Timemarkers$pct.1 - Neu.Timemarkers$pct.2

# Bone marrow score
Neutrophils_Hou$sample <- factor(Neutrophils_Hou$sample, levels = c(0,1,3), labels = c(0,1,3))
Idents(Neutrophils_Hou) <- "sample"
bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", "Retnlg", "Lcn2", "Olfm4", "Chil3", "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", "Lyz1"))
Neutrophils_Hou <- AddModuleScore(object = Neutrophils_Hou, features = bmGenes, search = TRUE, assay = "joined", name = "neutrotimebm", ctrl = 20, nbin = 10)
VlnPlot(Neutrophils_Hou, features = "BoneMarrow1") + ggtitle("Bone Marrow Hou")