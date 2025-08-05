###############################################################################
#            Integration Script for CellCall Arg1 Analysis
#        with Existing Single-Cell RNA-seq Analysis Pipeline
###############################################################################

# Source the cellcall analysis script
source("cellcall_arg1_analysis.R")

# Load additional required libraries for integration
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

###############################################################################
#                     Integration with Existing Pipeline
###############################################################################

# Function to integrate cellcall analysis with existing neutrophil analysis
integrate_cellcall_with_neutrophil_analysis <- function() {
  
  message("Integrating CellCall Arg1 Analysis with Existing Pipeline")
  message("========================================================")
  
  # Try to load existing integrated objects
  integrated_files <- c("seuObj_s_integrated.rds", "combined_seu.rds", "Neutrophils_Hou.rds")
  seurat_obj <- NULL
  
  for (file in integrated_files) {
    if (file.exists(file)) {
      message(paste("Loading", file))
      seurat_obj <- readRDS(file)
      break
    }
  }
  
  if (is.null(seurat_obj)) {
    message("No existing integrated object found. Please run the main analysis first.")
    return(NULL)
  }
  
  # Check if the object has the required metadata
  required_columns <- c("sample")
  missing_columns <- setdiff(required_columns, colnames(seurat_obj@meta.data))
  
  if (length(missing_columns) > 0) {
    message(paste("Missing required metadata columns:", paste(missing_columns, collapse = ", ")))
    return(NULL)
  }
  
  # Run cellcall analysis
  cellcall_results <- run_cellcall_arg1_analysis(seurat_object = seurat_obj)
  
  return(cellcall_results)
}

# Function to create enhanced visualizations combining cellcall with existing analysis
create_enhanced_arg1_visualizations <- function(seurat_obj, cellcall_results) {
  
  message("Creating enhanced visualizations...")
  
  # Create enhanced plots directory
  enhanced_dir <- "enhanced_cellcall_plots"
  if (!dir.exists(enhanced_dir)) {
    dir.create(enhanced_dir, recursive = TRUE)
  }
  
  # 1. UMAP plot with Arg1 expression overlay
  if ("umap" %in% names(seurat_obj@reductions)) {
    p_umap_arg1 <- FeaturePlot(seurat_obj, features = "Arg1", 
                               reduction = "umap", pt.size = 0.5) +
      scale_color_viridis_c(name = "Arg1\nExpression") +
      ggtitle("Arg1 Expression on UMAP") +
      theme_minimal()
    
    ggsave(file.path(enhanced_dir, "umap_arg1_expression.png"), 
           p_umap_arg1, width = 8, height = 6)
  }
  
  # 2. UMAP split by time points focusing on 3dpi
  if ("time" %in% colnames(seurat_obj@meta.data) && "umap" %in% names(seurat_obj@reductions)) {
    p_umap_time <- DimPlot(seurat_obj, reduction = "umap", 
                           group.by = "time", split.by = "time", ncol = 2) +
      ggtitle("UMAP by Time Points - Focus on 3dpi Pattern")
    
    ggsave(file.path(enhanced_dir, "umap_by_timepoints.png"), 
           p_umap_time, width = 12, height = 8)
  }
  
  # 3. Violin plot comparing Arg1 expression across cell types (if cell_type exists)
  if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
    p_violin_celltype <- VlnPlot(seurat_obj, features = "Arg1", 
                                 group.by = "cell_type", pt.size = 0.1) +
      ggtitle("Arg1 Expression Across Cell Types") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(enhanced_dir, "arg1_expression_by_celltype.png"), 
           p_violin_celltype, width = 10, height = 6)
  }
  
  # 4. Heatmap of Arg1 expression by sample and time
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    # Calculate average Arg1 expression by sample
    avg_expression <- seurat_obj@meta.data %>%
      mutate(Arg1_expr = FetchData(seurat_obj, vars = "Arg1")$Arg1) %>%
      group_by(sample) %>%
      summarise(
        avg_arg1 = mean(Arg1_expr, na.rm = TRUE),
        median_arg1 = median(Arg1_expr, na.rm = TRUE),
        percent_positive = sum(Arg1_expr > 0) / n() * 100,
        .groups = 'drop'
      )
    
    # Create heatmap-style plot
    p_heatmap <- avg_expression %>%
      pivot_longer(cols = c("avg_arg1", "median_arg1", "percent_positive"),
                   names_to = "metric", values_to = "value") %>%
      mutate(metric = case_when(
        metric == "avg_arg1" ~ "Average Expression",
        metric == "median_arg1" ~ "Median Expression",
        metric == "percent_positive" ~ "% Positive Cells"
      )) %>%
      ggplot(aes(x = sample, y = metric, fill = value)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Value") +
      labs(title = "Arg1 Expression Metrics by Sample",
           subtitle = "Highlighting 3dpi pattern",
           x = "Sample", y = "Metric") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(enhanced_dir, "arg1_heatmap_by_sample.png"), 
           p_heatmap, width = 10, height = 6)
  }
  
  # 5. Correlation with sequencing depth
  qc_data <- seurat_obj@meta.data %>%
    mutate(Arg1_expr = FetchData(seurat_obj, vars = "Arg1")$Arg1)
  
  p_depth_correlation <- ggplot(qc_data, aes(x = nCount_RNA, y = Arg1_expr)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = "Arg1 Expression vs Sequencing Depth",
         subtitle = "Assessing technical vs biological zeros",
         x = "Total UMI Count", y = "Arg1 Expression") +
    theme_minimal()
  
  p_feature_correlation <- ggplot(qc_data, aes(x = nFeature_RNA, y = Arg1_expr)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = "Arg1 Expression vs Feature Count",
         subtitle = "Evaluating detection sensitivity",
         x = "Number of Features", y = "Arg1 Expression") +
    theme_minimal()
  
  combined_correlation <- p_depth_correlation + p_feature_correlation
  ggsave(file.path(enhanced_dir, "arg1_sequencing_depth_correlation.png"), 
         combined_correlation, width = 12, height = 5)
  
  message(paste("Enhanced visualizations saved to:", enhanced_dir))
}

# Function to perform differential expression analysis for Arg1+ vs Arg1- cells
perform_arg1_differential_analysis <- function(seurat_obj) {
  
  message("Performing differential expression analysis for Arg1+ vs Arg1- cells...")
  
  # Add Arg1 status to metadata
  arg1_expression <- FetchData(seurat_obj, vars = "Arg1")$Arg1
  seurat_obj$Arg1_status <- ifelse(arg1_expression > 0, "Arg1_positive", "Arg1_negative")
  
  # Set identity to Arg1 status
  Idents(seurat_obj) <- "Arg1_status"
  
  # Find markers for Arg1+ cells
  tryCatch({
    arg1_markers <- FindMarkers(seurat_obj, 
                                ident.1 = "Arg1_positive", 
                                ident.2 = "Arg1_negative",
                                min.pct = 0.1,
                                logfc.threshold = 0.25)
    
    # Add gene names as a column
    arg1_markers$gene <- rownames(arg1_markers)
    
    # Save results
    write.csv(arg1_markers, "arg1_positive_vs_negative_markers.csv", row.names = FALSE)
    
    # Create volcano plot
    p_volcano <- ggplot(arg1_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(alpha = 0.6, size = 1) +
      geom_point(data = subset(arg1_markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
                 color = "red", size = 1.5) +
      labs(title = "Differential Expression: Arg1+ vs Arg1- Cells",
           subtitle = "Red points: significant genes (padj < 0.05, |logFC| > 0.5)",
           x = "Average Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
      theme_minimal()
    
    ggsave("enhanced_cellcall_plots/arg1_differential_volcano.png", 
           p_volcano, width = 8, height = 6)
    
    message("Differential expression analysis completed. Results saved to arg1_positive_vs_negative_markers.csv")
    return(arg1_markers)
    
  }, error = function(e) {
    message(paste("Error in differential expression analysis:", e$message))
    return(NULL)
  })
}

# Function to create a comprehensive analysis report
create_comprehensive_report <- function(seurat_obj, cellcall_results, de_results = NULL) {
  
  report_file <- "comprehensive_cellcall_arg1_report.txt"
  
  sink(report_file)
  
  cat("COMPREHENSIVE CELLCALL ARG1 ANALYSIS REPORT\n")
  cat("==========================================\n\n")
  
  cat("Dataset Overview:\n")
  cat("-----------------\n")
  cat(paste("Total cells in dataset:", ncol(seurat_obj), "\n"))
  cat(paste("Total genes:", nrow(seurat_obj), "\n"))
  cat(paste("Samples:", paste(unique(seurat_obj$sample), collapse = ", "), "\n\n"))
  
  if ("time" %in% colnames(seurat_obj@meta.data)) {
    cat(paste("Time points:", paste(unique(seurat_obj$time), collapse = ", "), "\n\n"))
  }
  
  cat("Arg1 Expression Analysis:\n")
  cat("-------------------------\n")
  arg1_expr <- FetchData(seurat_obj, vars = "Arg1")$Arg1
  cat(paste("Cells with detectable Arg1 expression:", sum(arg1_expr > 0), 
            sprintf("(%.1f%%)", sum(arg1_expr > 0)/length(arg1_expr)*100), "\n"))
  cat(paste("Cells with Arg1 = 0:", sum(arg1_expr == 0),
            sprintf("(%.1f%%)", sum(arg1_expr == 0)/length(arg1_expr)*100), "\n"))
  cat(paste("Mean Arg1 expression:", sprintf("%.3f", mean(arg1_expr)), "\n"))
  cat(paste("Median Arg1 expression:", sprintf("%.3f", median(arg1_expr)), "\n"))
  cat(paste("Max Arg1 expression:", sprintf("%.3f", max(arg1_expr)), "\n\n"))
  
  cat("3DPI Pattern Analysis:\n")
  cat("----------------------\n")
  if ("3 dpi" %in% names(cellcall_results$results)) {
    result_3dpi <- cellcall_results$results[["3 dpi"]]
    cat("The 3dpi time point shows distinct Arg1 expression patterns:\n")
    cat(paste("- Arg1+ cells:", result_3dpi$arg1_positive_cells, 
              sprintf("(%.1f%%)", result_3dpi$arg1_positive_percentage), "\n"))
    cat(paste("- Mean expression:", sprintf("%.3f", result_3dpi$mean_arg1_expression), "\n"))
    cat("This pattern may indicate specific biological processes at this time point.\n\n")
  }
  
  cat("Sequencing Depth Impact:\n")
  cat("------------------------\n")
  cat("Analysis of Arg1=0 cells suggests potential technical factors:\n")
  for (tp in names(cellcall_results$results)) {
    result <- cellcall_results$results[[tp]]
    if (!is.na(result$avg_features_arg1_zero)) {
      cat(paste("- ", tp, ": Avg features in Arg1=0 cells =", 
                sprintf("%.1f", result$avg_features_arg1_zero), "\n"))
    }
  }
  cat("\nCells with lower feature counts may have Arg1=0 due to technical limitations.\n\n")
  
  if (!is.null(de_results)) {
    cat("Differential Expression Summary:\n")
    cat("-------------------------------\n")
    sig_genes <- subset(de_results, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
    cat(paste("Significantly different genes between Arg1+ and Arg1- cells:", nrow(sig_genes), "\n"))
    if (nrow(sig_genes) > 0) {
      cat("Top upregulated genes in Arg1+ cells:\n")
      top_up <- head(sig_genes[order(-sig_genes$avg_log2FC), ], 5)
      for (i in 1:nrow(top_up)) {
        cat(paste("  ", top_up$gene[i], " (logFC:", sprintf("%.2f", top_up$avg_log2FC[i]), ")\n"))
      }
    }
  }
  
  cat("\nRecommendations:\n")
  cat("----------------\n")
  cat("1. Focus on 3dpi time point for further mechanistic studies\n")
  cat("2. Consider sequencing depth when interpreting Arg1=0 cells\n")
  cat("3. Validate Arg1 expression patterns using alternative methods\n")
  cat("4. Investigate biological pathways enriched in Arg1+ cells\n\n")
  
  cat("Output Files Generated:\n")
  cat("-----------------------\n")
  cat("- cellcall_plots/: Visualization files\n")
  cat("- enhanced_cellcall_plots/: Enhanced visualization files\n")
  cat("- cellcall_arg1_report.txt: Basic cellcall report\n")
  cat("- comprehensive_cellcall_arg1_report.txt: This comprehensive report\n")
  cat("- arg1_positive_vs_negative_markers.csv: Differential expression results\n")
  
  sink()
  
  message(paste("Comprehensive report saved to:", report_file))
}

###############################################################################
#                            Main Integration Function
###############################################################################

# Main function to run the complete integrated analysis
run_integrated_cellcall_analysis <- function(load_existing = TRUE) {
  
  message("Starting Integrated CellCall Arg1 Analysis")
  message("==========================================")
  
  # Step 1: Load or find existing data
  if (load_existing) {
    cellcall_results <- integrate_cellcall_with_neutrophil_analysis()
    if (is.null(cellcall_results)) {
      message("Could not load existing data. Please check file paths.")
      return(NULL)
    }
  } else {
    message("Please provide a Seurat object to analyze")
    return(NULL)
  }
  
  # Step 2: Load the Seurat object for additional analysis
  seurat_obj <- NULL
  integrated_files <- c("seuObj_s_integrated.rds", "combined_seu.rds", "Neutrophils_Hou.rds")
  
  for (file in integrated_files) {
    if (file.exists(file)) {
      seurat_obj <- readRDS(file)
      break
    }
  }
  
  if (is.null(seurat_obj)) {
    message("Could not load Seurat object for enhanced analysis")
    return(cellcall_results)
  }
  
  # Step 3: Create enhanced visualizations
  create_enhanced_arg1_visualizations(seurat_obj, cellcall_results)
  
  # Step 4: Perform differential expression analysis
  de_results <- perform_arg1_differential_analysis(seurat_obj)
  
  # Step 5: Create comprehensive report
  create_comprehensive_report(seurat_obj, cellcall_results, de_results)
  
  message("Integrated CellCall analysis completed successfully!")
  message("All results and visualizations have been generated.")
  
  return(list(
    cellcall_results = cellcall_results,
    differential_expression = de_results,
    seurat_object = seurat_obj
  ))
}

###############################################################################
#                            Script Execution
###############################################################################

message("Integration script loaded successfully!")
message("Use run_integrated_cellcall_analysis() to start the complete analysis")

# Uncomment the following line to run the analysis automatically:
# integrated_results <- run_integrated_cellcall_analysis()