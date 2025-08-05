# CellCall Analysis for Arg1 Expression Patterns

## Overview

This implementation provides comprehensive analysis of Arg1 expression patterns across different time points in single-cell RNA-seq data, with special focus on the 3dpi (3 days post-injury) time point where distinct patterns were observed. The analysis includes Arg1=0 cells to account for potential sequencing depth effects.

## Files Description

### 1. `cellcall_arg1_analysis.R`
Main analysis script containing:
- **Core Functions**: Arg1 expression analysis across time points
- **Visualization Functions**: Comprehensive plotting for Arg1 patterns
- **Report Generation**: Detailed analysis reports
- **Main Analysis Function**: `run_cellcall_arg1_analysis()`

### 2. `integrate_cellcall_analysis.R`
Integration script that:
- Connects with existing single-cell analysis pipeline
- Creates enhanced visualizations combining cellcall with existing analysis
- Performs differential expression analysis for Arg1+ vs Arg1- cells
- Generates comprehensive reports

### 3. `demo_cellcall_analysis.R`
Demonstration script that:
- Creates simulated dataset for testing
- Shows usage examples
- Validates the analysis pipeline

## Key Features

✅ **Time Point Analysis**: Analyzes Arg1 expression across Naive, 3dpi, 7dpi, and 14dpi time points  
✅ **3dpi Pattern Focus**: Special attention to the distinct pattern observed at 3dpi  
✅ **Arg1=0 Inclusion**: Includes cells with zero Arg1 expression to assess sequencing depth impact  
✅ **Sequencing Depth Analysis**: Correlates Arg1=0 cells with technical quality metrics  
✅ **Comprehensive Visualization**: Multiple plot types to visualize patterns  
✅ **Statistical Analysis**: Differential expression and statistical comparisons  
✅ **Integration Ready**: Works with existing single-cell analysis pipeline  

## Usage Instructions

### Option 1: Basic CellCall Analysis
```r
# Source the main script
source("cellcall_arg1_analysis.R")

# Run with your Seurat object file
results <- run_cellcall_arg1_analysis(seurat_object_path = "your_seurat_object.rds")

# Or with object in memory
results <- run_cellcall_arg1_analysis(seurat_object = your_seurat_object)
```

### Option 2: Integrated Analysis (Recommended)
```r
# Source the integration script
source("integrate_cellcall_analysis.R")

# Run complete integrated analysis
integrated_results <- run_integrated_cellcall_analysis()
```

### Option 3: Demo with Simulated Data
```r
# Source the demo script
source("demo_cellcall_analysis.R")

# Run demonstration
demo_results <- run_demo_cellcall_analysis()
```

## Analysis Components

### 1. Arg1 Expression Categorization
- **Arg1_zero**: Cells with Arg1 expression = 0
- **Arg1_low**: Cells with 0 < Arg1 ≤ 1
- **Arg1_medium**: Cells with 1 < Arg1 ≤ 3  
- **Arg1_high**: Cells with Arg1 > 3

### 2. Time Point Analysis
- **Naive**: Baseline condition
- **3 dpi**: Focus time point with distinct pattern
- **7 dpi**: Intermediate time point
- **14 dpi**: Later time point

### 3. Sequencing Depth Considerations
- Analysis of nFeature_RNA and nCount_RNA in Arg1=0 cells
- Assessment of technical vs biological zeros
- Correlation analysis between sequencing metrics and Arg1 detection

## Output Files

### Visualizations
- `cellcall_plots/`: Basic cellcall analysis plots
  - `arg1_distribution_barplot.png`: Arg1+ vs Arg1=0 distribution
  - `arg1_mean_expression.png`: Mean expression across time points
  - `arg1_zero_sequencing_depth.png`: Sequencing depth in Arg1=0 cells
  - `arg1_expression_violin.png`: Expression distribution by time point
  - `arg1_3dpi_highlight.png`: 3dpi pattern highlighted
  - `arg1_cellcall_analysis_combined.png`: Combined overview

- `enhanced_cellcall_plots/` (when using integrated analysis):
  - `umap_arg1_expression.png`: UMAP with Arg1 overlay
  - `umap_by_timepoints.png`: UMAP split by time points
  - `arg1_expression_by_celltype.png`: Expression by cell type
  - `arg1_heatmap_by_sample.png`: Heatmap of metrics by sample
  - `arg1_sequencing_depth_correlation.png`: Depth correlation analysis
  - `arg1_differential_volcano.png`: Volcano plot of DE genes

### Reports
- `cellcall_arg1_report.txt`: Basic analysis summary
- `comprehensive_cellcall_arg1_report.txt`: Detailed analysis report
- `arg1_positive_vs_negative_markers.csv`: Differential expression results

## Key Insights Provided

1. **3dpi Pattern Identification**: Quantifies the distinct Arg1 expression pattern at 3dpi
2. **Sequencing Depth Impact**: Assesses whether Arg1=0 cells are technical or biological
3. **Temporal Dynamics**: Tracks Arg1 expression changes across time points
4. **Cell Population Analysis**: Identifies Arg1+ cell populations and their characteristics
5. **Statistical Validation**: Provides statistical support for observed patterns

## Integration with Existing Pipeline

The cellcall analysis integrates seamlessly with the existing single-cell analysis:
- Works with existing Seurat objects from the pipeline
- Utilizes existing metadata (sample, time, cell_type)
- Builds upon existing QC metrics and UMAP embeddings
- Adds Arg1-specific analysis layers to existing results

## Requirements

### R Packages
```r
# Core analysis
library(Seurat)
library(ggplot2) 
library(dplyr)
library(tidyr)

# Visualization
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(viridis)
library(patchwork)

# Data manipulation
library(Matrix)
```

### Data Requirements
- Seurat object with:
  - Arg1 gene expression data
  - Time point metadata (time column or sample names containing time info)
  - Basic QC metrics (nFeature_RNA, nCount_RNA)

## Troubleshooting

### Common Issues

1. **Arg1 gene not found**: Check gene naming (Arg1 vs ARG1 vs Arginase1)
2. **Time point metadata missing**: Ensure time or sample columns exist
3. **Empty time points**: Some samples may not have cells for certain time points

### Solutions
- The script includes automatic fallbacks for metadata detection
- Error handling for missing time points
- Flexible gene name detection

## Scientific Rationale

### Why Include Arg1=0 Cells?
- **Technical zeros**: Low sequencing depth may cause false negatives
- **Biological zeros**: True absence of Arg1 expression
- **Pattern validation**: Ensures observed patterns are not technical artifacts

### Why Focus on 3dpi?
- User observed distinct pattern at this time point
- Critical time point in immune response/wound healing
- May represent peak of specific biological process

### Sequencing Depth Analysis
- Correlates Arg1 detection with nFeature_RNA and nCount_RNA
- Identifies if Arg1=0 cells have systematically lower sequencing quality
- Helps distinguish technical from biological effects

## Citation and Usage

This analysis was developed specifically for analyzing Arg1 expression patterns in single-cell RNA-seq data with focus on temporal dynamics and sequencing depth considerations. The implementation provides both standalone and integrated analysis options for maximum flexibility.

For questions or issues, please refer to the comprehensive documentation within each script file.