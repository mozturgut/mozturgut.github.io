###############################################################################
#                    Demo CellCall Arg1 Analysis 
#              (For demonstration with simulated data)
###############################################################################

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

# Source the main cellcall analysis script
source("cellcall_arg1_analysis.R")

###############################################################################
#                        Create Simulated Dataset
###############################################################################

create_simulated_seurat_object <- function() {
  
  message("Creating simulated dataset for CellCall analysis demonstration...")
  
  # Set random seed for reproducibility
  set.seed(123)
  
  # Define parameters
  n_genes <- 2000
  n_cells_per_sample <- c(500, 800, 600, 400)  # Naive, 3dpi, 7dpi, 14dpi
  sample_names <- c("young_naive", "young_3d", "young_7d", "young_14d")
  time_points <- c("Naive", "3 dpi", "7 dpi", "14 dpi")
  
  # Create gene names
  gene_names <- c("Arg1", paste0("Gene_", 1:(n_genes-1)))
  
  # Initialize combined data
  all_counts <- NULL
  all_metadata <- NULL
  
  for (i in 1:length(sample_names)) {
    n_cells <- n_cells_per_sample[i]
    sample_name <- sample_names[i]
    time_point <- time_points[i]
    
    # Create base expression matrix
    counts_matrix <- matrix(
      rpois(n_genes * n_cells, lambda = 2), 
      nrow = n_genes, 
      ncol = n_cells
    )
    rownames(counts_matrix) <- gene_names
    colnames(counts_matrix) <- paste0(sample_name, "_cell_", 1:n_cells)
    
    # Simulate Arg1 expression with time-specific patterns
    if (time_point == "Naive") {
      # Low Arg1 expression in naive
      arg1_prob <- 0.1
      arg1_high_cells <- sample(1:n_cells, size = round(n_cells * arg1_prob))
      counts_matrix["Arg1", arg1_high_cells] <- rpois(length(arg1_high_cells), lambda = 5)
      
    } else if (time_point == "3 dpi") {
      # High Arg1 expression at 3dpi (the pattern mentioned)
      arg1_prob <- 0.4  # 40% of cells express Arg1
      arg1_high_cells <- sample(1:n_cells, size = round(n_cells * arg1_prob))
      counts_matrix["Arg1", arg1_high_cells] <- rpois(length(arg1_high_cells), lambda = 8)
      
    } else if (time_point == "7 dpi") {
      # Moderate Arg1 expression
      arg1_prob <- 0.25
      arg1_high_cells <- sample(1:n_cells, size = round(n_cells * arg1_prob))
      counts_matrix["Arg1", arg1_high_cells] <- rpois(length(arg1_high_cells), lambda = 6)
      
    } else {  # 14 dpi
      # Low Arg1 expression returning to baseline
      arg1_prob <- 0.15
      arg1_high_cells <- sample(1:n_cells, size = round(n_cells * arg1_prob))
      counts_matrix["Arg1", arg1_high_cells] <- rpois(length(arg1_high_cells), lambda = 4)
    }
    
    # Create metadata
    metadata <- data.frame(
      sample = sample_name,
      time = time_point,
      nCount_RNA = colSums(counts_matrix),
      nFeature_RNA = colSums(counts_matrix > 0),
      row.names = colnames(counts_matrix)
    )
    
    # Add some variation in sequencing depth to simulate technical effects
    depth_variation <- rnorm(n_cells, mean = 1, sd = 0.3)
    depth_variation[depth_variation < 0.3] <- 0.3  # Minimum depth
    metadata$nCount_RNA <- round(metadata$nCount_RNA * depth_variation)
    metadata$nFeature_RNA <- round(metadata$nFeature_RNA * depth_variation)
    
    # Combine data
    if (is.null(all_counts)) {
      all_counts <- counts_matrix
      all_metadata <- metadata
    } else {
      all_counts <- cbind(all_counts, counts_matrix)
      all_metadata <- rbind(all_metadata, metadata)
    }
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = as(all_counts, "dgCMatrix"),
    meta.data = all_metadata,
    project = "CellCall_Demo"
  )
  
  # Add mitochondrial percentage (simulated)
  seurat_obj[["percent.mt"]] <- runif(ncol(seurat_obj), min = 1, max = 15)
  
  # Basic processing
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
  
  message("Simulated Seurat object created successfully!")
  message(paste("Total cells:", ncol(seurat_obj)))
  message(paste("Total genes:", nrow(seurat_obj)))
  message(paste("Samples:", paste(unique(seurat_obj$sample), collapse = ", ")))
  
  return(seurat_obj)
}

###############################################################################
#                        Run Demo Analysis
###############################################################################

run_demo_cellcall_analysis <- function() {
  
  message("CELLCALL ARG1 ANALYSIS DEMONSTRATION")
  message("====================================")
  
  # Create simulated data
  demo_seurat <- create_simulated_seurat_object()
  
  # Save the demo object
  saveRDS(demo_seurat, "demo_seurat_object.rds")
  
  # Run the cellcall analysis
  demo_results <- run_cellcall_arg1_analysis(seurat_object = demo_seurat)
  
  # Print summary results
  message("\nDEMO RESULTS SUMMARY:")
  message("====================")
  print(demo_results$summary_data)
  
  # Create additional demo-specific visualizations
  create_demo_visualizations(demo_seurat, demo_results)
  
  return(demo_results)
}

# Function to create demo-specific visualizations
create_demo_visualizations <- function(demo_seurat, demo_results) {
  
  demo_dir <- "demo_cellcall_plots"
  if (!dir.exists(demo_dir)) {
    dir.create(demo_dir, recursive = TRUE)
  }
  
  # UMAP with Arg1 expression
  p_umap <- FeaturePlot(demo_seurat, features = "Arg1", pt.size = 1) +
    ggtitle("Demo: Arg1 Expression on UMAP") +
    theme_minimal()
  
  ggsave(file.path(demo_dir, "demo_umap_arg1.png"), p_umap, width = 8, height = 6)
  
  # UMAP colored by time point
  p_umap_time <- DimPlot(demo_seurat, group.by = "time", pt.size = 1) +
    ggtitle("Demo: Cells by Time Point") +
    theme_minimal()
  
  ggsave(file.path(demo_dir, "demo_umap_timepoints.png"), p_umap_time, width = 8, height = 6)
  
  # Violin plot of Arg1 by time
  p_violin <- VlnPlot(demo_seurat, features = "Arg1", group.by = "time", pt.size = 0.5) +
    ggtitle("Demo: Arg1 Expression by Time Point") +
    theme_minimal()
  
  ggsave(file.path(demo_dir, "demo_violin_arg1.png"), p_violin, width = 8, height = 6)
  
  message(paste("Demo visualizations saved to:", demo_dir))
}

###############################################################################
#                        Example Usage Instructions
###############################################################################

print_usage_instructions <- function() {
  
  cat("\n")
  cat("CELLCALL ARG1 ANALYSIS - USAGE INSTRUCTIONS\n")
  cat("============================================\n\n")
  
  cat("This implementation provides comprehensive analysis of Arg1 expression patterns\n")
  cat("across different time points with special focus on the 3dpi pattern.\n\n")
  
  cat("FOR DEMONSTRATION (with simulated data):\n")
  cat("----------------------------------------\n")
  cat("demo_results <- run_demo_cellcall_analysis()\n\n")
  
  cat("FOR REAL DATA ANALYSIS:\n")
  cat("-----------------------\n")
  cat("# Option 1: Load from file\n")
  cat("results <- run_cellcall_arg1_analysis(seurat_object_path = 'your_seurat_object.rds')\n\n")
  
  cat("# Option 2: Use object in memory\n")
  cat("results <- run_cellcall_arg1_analysis(seurat_object = your_seurat_object)\n\n")
  
  cat("# Option 3: Integrated analysis with existing pipeline\n")
  cat("source('integrate_cellcall_analysis.R')\n")
  cat("integrated_results <- run_integrated_cellcall_analysis()\n\n")
  
  cat("KEY FEATURES:\n")
  cat("-------------\n")
  cat("✓ Analyzes Arg1 expression patterns across time points\n")
  cat("✓ Includes Arg1=0 cells to account for sequencing depth\n")
  cat("✓ Highlights 3dpi pattern as mentioned in requirements\n")
  cat("✓ Generates comprehensive visualizations and reports\n")
  cat("✓ Performs differential expression analysis\n")
  cat("✓ Correlates with sequencing quality metrics\n\n")
  
  cat("OUTPUT FILES:\n")
  cat("-------------\n")
  cat("• cellcall_plots/: Main visualization files\n")
  cat("• cellcall_arg1_report.txt: Basic analysis report\n")
  cat("• enhanced_cellcall_plots/: Additional enhanced plots\n")
  cat("• comprehensive_cellcall_arg1_report.txt: Detailed report\n")
  cat("• arg1_positive_vs_negative_markers.csv: Differential expression\n\n")
}

###############################################################################
#                        Script Initialization
###############################################################################

message("Demo CellCall Analysis Script Loaded!")
print_usage_instructions()

# Uncomment the next line to run the demo automatically:
# demo_results <- run_demo_cellcall_analysis()