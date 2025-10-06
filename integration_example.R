###############################################################################
#                    Integration Example for Existing Pipeline
#              Adding CellCall Arg1 Analysis to 3-17-25.Rmd
###############################################################################

# This file shows how to integrate the CellCall Arg1 analysis 
# into the existing single-cell analysis pipeline from 3-17-25.Rmd

###############################################################################
#        Example: Add to the end of your existing 3-17-25.Rmd
###############################################################################

# Add this section after the neutrophil analysis in your existing pipeline

# ```{r cellcall-arg1-analysis, error = TRUE}
# 
# ## CellCall Analysis for Arg1 Expression Patterns
# 
# # Source the cellcall analysis script
# source("cellcall_arg1_analysis.R")
# 
# # Load the integrated object (should already be available)
# if (!exists("combined_seu")) {
#   combined_seu <- readRDS("seuObj_s_integrated.rds")
# }
# 
# # Check if Arg1 is present in the dataset
# if ("Arg1" %in% rownames(combined_seu)) {
#   message("Arg1 gene found. Running CellCall analysis...")
#   
#   # Run the cellcall analysis
#   cellcall_results <- run_cellcall_arg1_analysis(seurat_object = combined_seu)
#   
#   # Display summary results
#   cat("CellCall Arg1 Analysis Summary:\n")
#   print(cellcall_results$summary_data)
#   
#   # Show the key plots
#   print(cellcall_results$plots$p1)  # Distribution barplot
#   print(cellcall_results$plots$p5)  # 3dpi highlight plot
#   
# } else {
#   message("Arg1 gene not found in dataset. Checking alternative names...")
#   
#   # Check for alternative gene names
#   possible_names <- c("ARG1", "Arginase1", "arginase", "arg1")
#   found_name <- NULL
#   
#   for (name in possible_names) {
#     if (name %in% rownames(combined_seu)) {
#       found_name <- name
#       break
#     }
#   }
#   
#   if (!is.null(found_name)) {
#     message(paste("Found Arg1 as:", found_name))
#     # You would need to rename the gene or modify the analysis accordingly
#   } else {
#     message("Arg1 gene not found under any alternative names.")
#   }
# }
# ```

###############################################################################
#        Example: For Neutrophil-Specific Analysis
###############################################################################

# Add this to analyze Arg1 specifically in neutrophils

# ```{r neutrophil-arg1-analysis, error = TRUE}
# 
# ## Arg1 Analysis in Neutrophils Only
# 
# # Check if neutrophils object exists
# if (exists("neutrophils")) {
#   neutrophil_obj <- neutrophils
# } else if (exists("Neutrophils_Hou")) {
#   neutrophil_obj <- Neutrophils_Hou
# } else {
#   # Subset neutrophils from combined object
#   if ("cell_type" %in% colnames(combined_seu@meta.data)) {
#     neutrophil_obj <- subset(combined_seu, subset = cell_type == "Neutrophils")
#   } else {
#     message("No neutrophil subset available for Arg1 analysis")
#     neutrophil_obj <- NULL
#   }
# }
# 
# if (!is.null(neutrophil_obj) && "Arg1" %in% rownames(neutrophil_obj)) {
#   
#   message("Running Arg1 analysis specifically in neutrophils...")
#   
#   # Run cellcall analysis on neutrophils only
#   neutrophil_cellcall_results <- run_cellcall_arg1_analysis(seurat_object = neutrophil_obj)
#   
#   # Create neutrophil-specific visualizations
#   cat("Neutrophil Arg1 Analysis Results:\n")
#   print(neutrophil_cellcall_results$summary_data)
#   
#   # Show neutrophil-specific plots
#   print(neutrophil_cellcall_results$plots$p4)  # Violin plot
#   
#   # Add UMAP overlay for neutrophils
#   if ("umap" %in% names(neutrophil_obj@reductions)) {
#     p_neutro_umap_arg1 <- FeaturePlot(neutrophil_obj, features = "Arg1", 
#                                       reduction = "umap", pt.size = 1) +
#       ggtitle("Arg1 Expression in Neutrophils (UMAP)") +
#       scale_color_viridis_c(name = "Arg1\nExpression")
#     print(p_neutro_umap_arg1)
#   }
# }
# ```

###############################################################################
#        Example: Integration with Existing Violin Plots
###############################################################################

# Add Arg1 to your existing neutrophil marker genes analysis

# ```{r enhanced-neutrophil-markers, error = TRUE}
# 
# ## Enhanced Neutrophil Marker Analysis Including Arg1
# 
# # Original neutrophil markers from your analysis
# neutrophil_marker_genes <- c("Ly6g", "S100a8", "S100a9", "Mmp9", "Il1r2", "Camp", "Ltf")
# 
# # Add Arg1 to the marker gene list
# if ("Arg1" %in% rownames(neutrophils)) {
#   enhanced_neutrophil_markers <- c(neutrophil_marker_genes, "Arg1")
# } else {
#   enhanced_neutrophil_markers <- neutrophil_marker_genes
#   message("Arg1 not found, using original marker set")
# }
# 
# # Create enhanced violin plots including Arg1
# for (gene in enhanced_neutrophil_markers) {
#   if (gene %in% rownames(neutrophils)) {
#     
#     # Fetch data for the current gene
#     expr_df <- FetchData(neutrophils, vars = c(gene, "sample"))
#     colnames(expr_df)[1] <- "expression"
#     
#     # Set sample factor levels
#     expr_df$sample <- factor(expr_df$sample, 
#                              levels = c("Young_Naive", "Young_3", "Young_7", "Young_14",
#                                         "Aged_Naive", "Aged_3", "Aged_7", "Aged_14"))
#     
#     # Create violin plot with enhanced styling for Arg1
#     if (gene == "Arg1") {
#       # Special highlighting for Arg1 3dpi pattern
#       p_gene <- ggplot(expr_df, aes(x = sample, y = expression, fill = sample)) +
#         geom_violin(trim = TRUE, scale = "width") +
#         geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
#         labs(title = paste("*** ARG1 EXPRESSION PATTERN ***"),
#              subtitle = "Note the distinct pattern at 3dpi time points",
#              x = "Age-Time Group",
#              y = paste(gene, "expression")) +
#         theme_classic() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1),
#               plot.title = element_text(hjust = 0.5, color = "red", face = "bold"),
#               plot.subtitle = element_text(hjust = 0.5, color = "blue"))
#       
#       # Highlight 3dpi samples
#       p_gene <- p_gene + 
#         annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, 
#                  alpha = 0.1, fill = "red") +
#         annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf, 
#                  alpha = 0.1, fill = "red")
#       
#     } else {
#       # Standard plot for other genes
#       p_gene <- ggplot(expr_df, aes(x = sample, y = expression, fill = sample)) +
#         geom_violin(trim = TRUE, scale = "width") +
#         geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
#         labs(title = paste("Expression of", gene, "by Age-Time Group"),
#              x = "Age-Time Group",
#              y = paste(gene, "expression")) +
#         theme_classic() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1),
#               plot.title = element_text(hjust = 0.5))
#     }
#     
#     # Print and save the plot
#     print(p_gene)
#     ggsave(filename = paste0("Enhanced_Violin_", gene, ".png"),
#            plot = p_gene, width = 10, height = 7)
#   }
# }
# ```

###############################################################################
#        Example: Summary Statistics Table
###############################################################################

# Create a summary table including Arg1 statistics

# ```{r arg1-summary-table, error = TRUE}
# 
# ## Create Enhanced Summary Table with Arg1 Statistics
# 
# if (exists("cellcall_results") && "Arg1" %in% rownames(combined_seu)) {
#   
#   # Extract summary data
#   cellcall_summary <- cellcall_results$summary_data
#   
#   # Create a comprehensive summary table
#   enhanced_summary <- cellcall_summary %>%
#     mutate(
#       `Time Point` = time_point,
#       `Total Cells` = total_cells,
#       `Arg1+ Cells (%)` = paste0(round(arg1_positive_percentage, 1), "%"),
#       `Mean Arg1 Expression` = round(mean_arg1_expression, 3),
#       `3dpi Pattern` = ifelse(time_point == "3 dpi", "*** PATTERN ***", "")
#     ) %>%
#     select(`Time Point`, `Total Cells`, `Arg1+ Cells (%)`, 
#            `Mean Arg1 Expression`, `3dpi Pattern`)
#   
#   # Display the table
#   cat("ARG1 EXPRESSION SUMMARY TABLE\n")
#   cat("=============================\n")
#   print(enhanced_summary)
#   
#   # Save as CSV
#   write.csv(enhanced_summary, "Arg1_Expression_Summary.csv", row.names = FALSE)
#   
#   # Highlight the 3dpi pattern
#   cat("\n*** 3DPI PATTERN DETECTED ***\n")
#   cat("The 3dpi time point shows distinct Arg1 expression characteristics.\n")
#   cat("This pattern warrants further investigation for biological significance.\n")
# }
# ```

###############################################################################
#                        Instructions for Integration
###############################################################################

message("Integration Example Script Loaded!")
message("=====================================")
message("This script shows examples of how to integrate CellCall Arg1 analysis")
message("into your existing single-cell analysis pipeline.")
message("")
message("To integrate:")
message("1. Copy the relevant code blocks into your existing .Rmd file")
message("2. Adjust variable names to match your pipeline")
message("3. Run the analysis sections after your main analysis")
message("")
message("Key integration points:")
message("- After loading/creating your integrated Seurat object")
message("- After neutrophil subset analysis")
message("- After existing marker gene analysis")
message("- Before final summary and conclusions")