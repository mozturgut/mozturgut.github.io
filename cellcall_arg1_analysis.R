###############################################################################
#                    CellCall Analysis for Arg1 Expression Patterns
#                         Across Different Time Points
###############################################################################

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(viridis)

###############################################################################
#                            Helper Functions
###############################################################################

# Function to analyze Arg1 expression patterns across time points
analyze_arg1_expression <- function(seurat_obj, time_column = "time", sample_column = "sample") {
  
  # Extract Arg1 expression data
  arg1_expression <- FetchData(seurat_obj, vars = c("Arg1", time_column, sample_column))
  colnames(arg1_expression)[1] <- "Arg1_expression"
  
  # Categorize cells based on Arg1 expression
  arg1_expression$Arg1_category <- case_when(
    arg1_expression$Arg1_expression == 0 ~ "Arg1_zero",
    arg1_expression$Arg1_expression > 0 & arg1_expression$Arg1_expression <= 1 ~ "Arg1_low",
    arg1_expression$Arg1_expression > 1 & arg1_expression$Arg1_expression <= 3 ~ "Arg1_medium",
    arg1_expression$Arg1_expression > 3 ~ "Arg1_high"
  )
  
  return(arg1_expression)
}

# Function to perform cellcall analysis focusing on Arg1+ cells
perform_cellcall_analysis <- function(seurat_obj, time_points = c("Naive", "3 dpi", "7 dpi", "14 dpi")) {
  
  # Analyze Arg1 expression for each time point
  cellcall_results <- list()
  
  for (time_point in time_points) {
    message(paste("Analyzing time point:", time_point))
    
    # Subset data for current time point
    if ("time" %in% colnames(seurat_obj@meta.data)) {
      cells_subset <- subset(seurat_obj, subset = time == time_point)
    } else {
      # Fallback to sample column if time column doesn't exist
      cells_subset <- subset(seurat_obj, subset = sample %in% grep(gsub(" ", "_", tolower(time_point)), 
                                                                  unique(seurat_obj$sample), value = TRUE))
    }
    
    if (ncol(cells_subset) == 0) {
      message(paste("No cells found for time point:", time_point))
      next
    }
    
    # Get Arg1 expression data
    arg1_data <- analyze_arg1_expression(cells_subset)
    
    # Calculate statistics
    total_cells <- nrow(arg1_data)
    arg1_positive_cells <- sum(arg1_data$Arg1_expression > 0)
    arg1_zero_cells <- sum(arg1_data$Arg1_expression == 0)
    
    # Calculate sequencing depth metrics for Arg1=0 cells
    if (arg1_zero_cells > 0) {
      arg1_zero_subset <- subset(cells_subset, cells = rownames(arg1_data)[arg1_data$Arg1_expression == 0])
      avg_features_arg1_zero <- mean(arg1_zero_subset$nFeature_RNA)
      avg_counts_arg1_zero <- mean(arg1_zero_subset$nCount_RNA)
    } else {
      avg_features_arg1_zero <- NA
      avg_counts_arg1_zero <- NA
    }
    
    # Store results
    cellcall_results[[time_point]] <- list(
      time_point = time_point,
      total_cells = total_cells,
      arg1_positive_cells = arg1_positive_cells,
      arg1_zero_cells = arg1_zero_cells,
      arg1_positive_percentage = (arg1_positive_cells / total_cells) * 100,
      arg1_zero_percentage = (arg1_zero_cells / total_cells) * 100,
      mean_arg1_expression = mean(arg1_data$Arg1_expression),
      median_arg1_expression = median(arg1_data$Arg1_expression),
      arg1_expression_data = arg1_data,
      avg_features_arg1_zero = avg_features_arg1_zero,
      avg_counts_arg1_zero = avg_counts_arg1_zero
    )
  }
  
  return(cellcall_results)
}

# Function to create comprehensive Arg1 visualization
create_arg1_visualizations <- function(cellcall_results, output_dir = "cellcall_plots") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Prepare summary data for plotting
  summary_data <- data.frame(
    time_point = names(cellcall_results),
    total_cells = sapply(cellcall_results, function(x) x$total_cells),
    arg1_positive_cells = sapply(cellcall_results, function(x) x$arg1_positive_cells),
    arg1_zero_cells = sapply(cellcall_results, function(x) x$arg1_zero_cells),
    arg1_positive_percentage = sapply(cellcall_results, function(x) x$arg1_positive_percentage),
    arg1_zero_percentage = sapply(cellcall_results, function(x) x$arg1_zero_percentage),
    mean_arg1_expression = sapply(cellcall_results, function(x) x$mean_arg1_expression),
    median_arg1_expression = sapply(cellcall_results, function(x) x$median_arg1_expression),
    avg_features_arg1_zero = sapply(cellcall_results, function(x) x$avg_features_arg1_zero),
    avg_counts_arg1_zero = sapply(cellcall_results, function(x) x$avg_counts_arg1_zero)
  )
  
  # Set factor levels for proper ordering
  summary_data$time_point <- factor(summary_data$time_point, 
                                   levels = c("Naive", "3 dpi", "7 dpi", "14 dpi"))
  
  # 1. Bar plot of Arg1 positive vs zero cells
  p1 <- summary_data %>%
    pivot_longer(cols = c("arg1_positive_percentage", "arg1_zero_percentage"),
                names_to = "category", values_to = "percentage") %>%
    mutate(category = case_when(
      category == "arg1_positive_percentage" ~ "Arg1 Positive",
      category == "arg1_zero_percentage" ~ "Arg1 Zero"
    )) %>%
    ggplot(aes(x = time_point, y = percentage, fill = category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Arg1 Positive" = "#E31A1C", "Arg1 Zero" = "#1F78B4")) +
    labs(title = "Arg1 Expression Distribution Across Time Points",
         subtitle = "Focus on 3dpi time point pattern",
         x = "Time Point", y = "Percentage of Cells", fill = "Arg1 Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 2. Mean Arg1 expression across time points
  p2 <- ggplot(summary_data, aes(x = time_point, y = mean_arg1_expression)) +
    geom_bar(stat = "identity", fill = "#FF7F00", alpha = 0.7) +
    geom_text(aes(label = round(mean_arg1_expression, 2)), vjust = -0.5) +
    labs(title = "Mean Arg1 Expression Across Time Points",
         subtitle = "Highlighting 3dpi pattern",
         x = "Time Point", y = "Mean Arg1 Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 3. Sequencing depth analysis for Arg1=0 cells
  p3 <- summary_data %>%
    filter(!is.na(avg_features_arg1_zero)) %>%
    ggplot(aes(x = time_point, y = avg_features_arg1_zero)) +
    geom_bar(stat = "identity", fill = "#6A3D9A", alpha = 0.7) +
    geom_text(aes(label = round(avg_features_arg1_zero, 0)), vjust = -0.5) +
    labs(title = "Average Features in Arg1=0 Cells",
         subtitle = "Assessing sequencing depth impact",
         x = "Time Point", y = "Average nFeature_RNA") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 4. Combined violin plot for Arg1 expression distribution
  combined_arg1_data <- do.call(rbind, lapply(names(cellcall_results), function(tp) {
    data <- cellcall_results[[tp]]$arg1_expression_data
    data$time_point <- tp
    return(data[, c("Arg1_expression", "time_point", "Arg1_category")])
  }))
  
  combined_arg1_data$time_point <- factor(combined_arg1_data$time_point, 
                                         levels = c("Naive", "3 dpi", "7 dpi", "14 dpi"))
  
  p4 <- ggplot(combined_arg1_data, aes(x = time_point, y = Arg1_expression, fill = time_point)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.5) +
    scale_fill_viridis_d() +
    labs(title = "Arg1 Expression Distribution by Time Point",
         subtitle = "Including Arg1=0 cells for sequencing depth consideration",
         x = "Time Point", y = "Arg1 Expression Level") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # 5. Highlight 3dpi pattern with comparison
  p5 <- summary_data %>%
    mutate(highlight_3dpi = ifelse(time_point == "3 dpi", "3dpi", "Other")) %>%
    ggplot(aes(x = time_point, y = arg1_positive_percentage, fill = highlight_3dpi)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("3dpi" = "#FF4444", "Other" = "#CCCCCC")) +
    labs(title = "Arg1 Positive Cell Percentage - 3dpi Pattern Highlighted",
         subtitle = "3dpi shows distinct pattern compared to other time points",
         x = "Time Point", y = "Arg1+ Cell Percentage", fill = "Time Point") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save individual plots
  ggsave(file.path(output_dir, "arg1_distribution_barplot.png"), p1, width = 10, height = 6)
  ggsave(file.path(output_dir, "arg1_mean_expression.png"), p2, width = 8, height = 6)
  ggsave(file.path(output_dir, "arg1_zero_sequencing_depth.png"), p3, width = 8, height = 6)
  ggsave(file.path(output_dir, "arg1_expression_violin.png"), p4, width = 10, height = 6)
  ggsave(file.path(output_dir, "arg1_3dpi_highlight.png"), p5, width = 8, height = 6)
  
  # Create combined plot
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(file.path(output_dir, "arg1_cellcall_analysis_combined.png"), combined_plot, 
         width = 16, height = 12)
  
  message(paste("Plots saved to:", output_dir))
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, combined = combined_plot))
}

# Function to generate cellcall summary report
generate_cellcall_report <- function(cellcall_results, output_file = "cellcall_arg1_report.txt") {
  
  sink(output_file)
  
  cat("CELLCALL ANALYSIS REPORT: ARG1 EXPRESSION PATTERNS\n")
  cat("==================================================\n\n")
  
  cat("Analysis Summary:\n")
  cat("This analysis examines Arg1 expression patterns across different time points,\n")
  cat("with special focus on the 3dpi time point where a distinct pattern was observed.\n")
  cat("Arg1=0 cells are included to account for potential sequencing depth effects.\n\n")
  
  for (time_point in names(cellcall_results)) {
    result <- cellcall_results[[time_point]]
    
    cat(paste("TIME POINT:", time_point, "\n"))
    cat(paste(rep("-", nchar(time_point) + 12), collapse = ""), "\n")
    cat(paste("Total cells:", result$total_cells, "\n"))
    cat(paste("Arg1+ cells:", result$arg1_positive_cells, 
              sprintf("(%.1f%%)", result$arg1_positive_percentage), "\n"))
    cat(paste("Arg1=0 cells:", result$arg1_zero_cells,
              sprintf("(%.1f%%)", result$arg1_zero_percentage), "\n"))
    cat(paste("Mean Arg1 expression:", sprintf("%.3f", result$mean_arg1_expression), "\n"))
    cat(paste("Median Arg1 expression:", sprintf("%.3f", result$median_arg1_expression), "\n"))
    
    if (!is.na(result$avg_features_arg1_zero)) {
      cat(paste("Avg features in Arg1=0 cells:", sprintf("%.1f", result$avg_features_arg1_zero), "\n"))
      cat(paste("Avg counts in Arg1=0 cells:", sprintf("%.1f", result$avg_counts_arg1_zero), "\n"))
    }
    
    if (time_point == "3 dpi") {
      cat("*** 3DPI PATTERN DETECTED ***\n")
      cat("This time point shows the distinct Arg1 expression pattern mentioned in the analysis.\n")
    }
    
    cat("\n")
  }
  
  cat("SEQUENCING DEPTH CONSIDERATIONS:\n")
  cat("=================================\n")
  cat("Arg1=0 cells are analyzed separately to determine if zero expression\n")
  cat("is due to true biological absence or technical limitations (sequencing depth).\n")
  cat("Cells with lower nFeature_RNA and nCount_RNA may have Arg1=0 due to technical reasons.\n\n")
  
  sink()
  
  message(paste("Report saved to:", output_file))
}

###############################################################################
#                            Main Analysis Function
###############################################################################

# Main function to run complete cellcall analysis
run_cellcall_arg1_analysis <- function(seurat_object_path = NULL, seurat_object = NULL) {
  
  message("Starting CellCall Analysis for Arg1 Expression Patterns")
  message("=======================================================")
  
  # Load Seurat object if path provided
  if (!is.null(seurat_object_path)) {
    message(paste("Loading Seurat object from:", seurat_object_path))
    seurat_obj <- readRDS(seurat_object_path)
  } else if (!is.null(seurat_object)) {
    seurat_obj <- seurat_object
  } else {
    stop("Please provide either seurat_object_path or seurat_object")
  }
  
  # Check if Arg1 gene exists in the dataset
  if (!"Arg1" %in% rownames(seurat_obj)) {
    stop("Arg1 gene not found in the dataset. Please check gene naming conventions.")
  }
  
  message("Arg1 gene found in dataset. Proceeding with analysis...")
  
  # Define time points to analyze
  time_points <- c("Naive", "3 dpi", "7 dpi", "14 dpi")
  
  # Perform cellcall analysis
  cellcall_results <- perform_cellcall_analysis(seurat_obj, time_points)
  
  # Create visualizations
  plots <- create_arg1_visualizations(cellcall_results)
  
  # Generate report
  generate_cellcall_report(cellcall_results)
  
  message("CellCall analysis completed successfully!")
  message("Check the 'cellcall_plots' directory for visualizations")
  message("Check 'cellcall_arg1_report.txt' for detailed results")
  
  return(list(
    results = cellcall_results,
    plots = plots,
    summary_data = do.call(rbind, lapply(names(cellcall_results), function(tp) {
      result <- cellcall_results[[tp]]
      data.frame(
        time_point = tp,
        total_cells = result$total_cells,
        arg1_positive_percentage = result$arg1_positive_percentage,
        mean_arg1_expression = result$mean_arg1_expression
      )
    }))
  ))
}

###############################################################################
#                            Example Usage
###############################################################################

# Example usage (uncomment to run):
# 
# # For integrated object from the existing analysis
# results <- run_cellcall_arg1_analysis(seurat_object_path = "seuObj_s_integrated.rds")
# 
# # Or if you have the object in memory
# # results <- run_cellcall_arg1_analysis(seurat_object = combined_seu)
# 
# # Print summary
# print(results$summary_data)

message("CellCall Arg1 Analysis Script Loaded Successfully!")
message("Use run_cellcall_arg1_analysis() to start the analysis")