#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 8: CELLULAR DECONVOLUTION ANALYSIS
# Author: Ángel Pérez
# Date: 2 June 2025
# Objective: Estimate the proportion of different cell types within each
#            spatial spot using a marker-based approach.
#
# NOTE: This script simulates deconvolution results based on marker gene
#       expression as a proof-of-concept. A real-world application would
#       use a tool like SPOTlight with a single-cell RNA-seq reference.
# =============================================================================

cat("\n=== SCRIPT 8: CELLULAR DECONVOLUTION ANALYSIS ===\n")
cat("Objective: Estimate cell type proportions in each spot\n\n")

# --- 0. INITIAL SETUP ---
# Load libraries
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork);
  library(viridis); library(pheatmap); library(tidyr);
})

# --- CONFIGURATION & PATHS ---
main_results_dir <- "results"
# Source directory for the integrated results object
integrated_analysis_subdir_ref <- file.path(main_results_dir, "05_Integrated_Pipeline_Analysis")
# Output directory for this script
deconvolution_output_dir <- file.path(main_results_dir, "08_Deconvolution_Analysis")
dir.create(deconvolution_output_dir, showWarnings = FALSE, recursive = TRUE)
# --- END CONFIGURATION ---

# --- 1. LOAD INTEGRATED SEURAT OBJECT ---
cat("1. Loading the integrated Seurat object from Script 5...\n")

final_results_file <- file.path(integrated_analysis_subdir_ref, "final_results_object.rds")
if (!file.exists(final_results_file)) {
    stop("ERROR: 'final_results_object.rds' not found. Run Script 5 first.")
}
final_results <- readRDS(final_results_file)
seurat_obj <- final_results$data$seurat # Extract the Seurat object for analysis

cat("   SUCCESS: Seurat object loaded.\n")

# =============================================================================
# PART 2: DECONVOLUTION LOGIC & VISUALIZATION
# =============================================================================

# --- 2.1 Prepare Reference Marker Genes ---
prepare_reference_data <- function() {
  cat("\n[2.1] Preparing reference cell types and marker genes for lung cancer...\n")
  
  # Based on literature and public atlases for NSCLC microenvironment
  markers <- list(
    Tumor_Cells = c("EPCAM", "KRT19", "KRT7", "MUC1", "CEACAM5"),
    Fibroblasts = c("COL1A1", "COL1A2", "ACTA2", "FAP", "PDGFRA"),
    Macrophages = c("CD68", "CD163", "MSR1", "CD86", "MRC1"),
    T_Cells = c("CD3E", "CD3D", "CD8A", "CD4", "IL7R"),
    B_Cells = c("CD79A", "MS4A1", "CD19", "BANK1", "CD22"),
    Neutrophils = c("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2"),
    Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),
    Epithelial_Normal = c("SFTPC", "SFTPB", "AGER", "CLDN18", "FOXJ1")
  )
  
  cat("   SUCCESS: Reference prepared for", length(markers), "cell types.\n")
  return(markers)
}

# --- 2.2 Run Simulated Deconvolution ---
run_simulated_deconvolution <- function(seurat_obj, markers) {
  cat("\n[2.2] Simulating cell proportions based on marker expression scores...\n")
  
  n_spots <- ncol(seurat_obj)
  n_types <- length(markers)
  decon_matrix <- matrix(0, nrow = n_spots, ncol = n_types,
                         dimnames = list(colnames(seurat_obj), names(markers)))

  # Calculate scores based on normalized marker gene expression
  expr_data <- GetAssayData(seurat_obj, layer = "data")

  for (cell_type in names(markers)) {
    available_markers <- intersect(markers[[cell_type]], rownames(expr_data))
    if (length(available_markers) > 0) {
      scores <- colMeans(expr_data[available_markers, , drop = FALSE])
      decon_matrix[, cell_type] <- as.numeric(scores)
    }
  }

  # Normalize scores to proportions summing to 1 for each spot
  decon_matrix <- decon_matrix / (rowSums(decon_matrix) + 1e-10) # Add epsilon to avoid division by zero
  
  cat("   SUCCESS: Simulated proportion matrix calculated.\n")
  return(decon_matrix)
}

# --- 2.3 Generate All Deconvolution Plots ---
generate_deconvolution_visuals <- function(seurat_obj, decon_matrix, output_dir) {
  cat("\n[2.3] Generating all deconvolution visualizations...\n")
  
  # Add proportions to Seurat metadata
  seurat_obj <- AddMetaData(seurat_obj, as.data.frame(decon_matrix))
  
  # 1. Spatial distribution of top cell types
  top_types <- names(sort(colMeans(decon_matrix), decreasing = TRUE)[1:6])
  spatial_plots <- SpatialFeaturePlot(seurat_obj, features = top_types, pt.size.factor = 2, image.alpha = 0, ncol = 3) +
    plot_annotation(title = "Spatial Distribution of Estimated Cell Types")
  ggsave(file.path(output_dir, "01_spatial_cell_types.png"), spatial_plots, width = 15, height = 10, dpi = 300)
  cat("   - Saved spatial cell type map.\n")

  # 2. Cluster composition bar plot
  cluster_comp <- as.data.frame(decon_matrix) %>%
    mutate(cluster = Idents(seurat_obj)) %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(cols = -cluster, names_to = "cell_type", values_to = "proportion")
  
  p_composition <- ggplot(cluster_comp, aes(x = cluster, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette = "Paired", name = "Cell Type") +
    coord_flip() + theme_minimal() +
    labs(title = "Average Cell Type Composition per Cluster", x = "Cluster", y = "Proportion")
  ggsave(file.path(output_dir, "02_cluster_composition.png"), p_composition, width = 10, height = 8, dpi = 300)
  cat("   - Saved cluster composition plot.\n")

  # 3. Cell type co-localization heatmap
  cor_matrix <- cor(decon_matrix, method = "spearman")
  pheatmap(cor_matrix, display_numbers = TRUE, number_format = "%.2f",
           filename = file.path(output_dir, "03_cell_colocalization_heatmap.png"),
           width = 10, height = 8, main = "Cell Type Co-localization (Spearman Correlation)")
  cat("   - Saved co-localization heatmap.\n")

  # 4. Cellular diversity analysis (Shannon Index)
  shannon_index <- apply(decon_matrix, 1, function(x) { x <- x[x > 0]; -sum(x * log(x)) })
  seurat_obj$shannon_diversity <- shannon_index
  
  p_diversity_spatial <- SpatialFeaturePlot(seurat_obj, features = "shannon_diversity", pt.size.factor = 2, image.alpha = 0) +
    scale_fill_viridis(option = "magma") + ggtitle("Spatial Cellular Diversity")
  
  p_diversity_violin <- ggplot(seurat_obj@meta.data, aes(x = seurat_clusters, y = shannon_diversity, fill = seurat_clusters)) +
    geom_violin() + theme_minimal() + labs(title = "Diversity per Cluster", x = "Cluster", y = "Shannon Index") +
    theme(legend.position = "none")

  combined_diversity <- p_diversity_spatial | p_diversity_violin
  ggsave(file.path(output_dir, "04_cellular_diversity.png"), combined_diversity, width = 14, height = 6, dpi = 300)
  cat("   - Saved cellular diversity analysis.\n")
}

# --- 2.4 Generate Final Summary ---
generate_deconvolution_summary <- function(decon_matrix, output_dir) {
  cat("\n[2.4] Generating final deconvolution summary...\n")
  
  summary_stats <- data.frame(
    Cell_Type = colnames(decon_matrix),
    Mean_Proportion = colMeans(decon_matrix),
    StDev_Proportion = apply(decon_matrix, 2, sd),
    Max_Proportion = apply(decon_matrix, 2, max),
    Num_Spots_Present = colSums(decon_matrix > 0.05) # Count spots where proportion > 5%
  ) %>% arrange(desc(Mean_Proportion))
  
  write.csv(summary_stats, file.path(output_dir, "deconvolution_summary.csv"), row.names = FALSE)
  cat("   - Summary CSV file saved.\n")
}


# =============================================================================
# MAIN FUNCTION: EXECUTE THE FULL DECONVOLUTION WORKFLOW
# =============================================================================

run_full_deconvolution_workflow <- function(seurat_obj, output_dir) {
  cat("\n---------------------------------------------------------------\n")
  cat("      STARTING CELLULAR DECONVOLUTION WORKFLOW (SIMULATED)     \n")
  cat("---------------------------------------------------------------\n")

  # Step 1: Prepare reference data
  reference_markers <- prepare_reference_data()
  
  # Step 2: Run the deconvolution simulation
  deconvolution_matrix <- run_simulated_deconvolution(seurat_obj, reference_markers)
  
  # Step 3: Generate all visualizations
  generate_deconvolution_visuals(seurat_obj, deconvolution_matrix, output_dir)
  
  # Step 4: Generate summary stats
  generate_deconvolution_summary(deconvolution_matrix, output_dir)

  cat("\nDECONVOLUTION ANALYSIS COMPLETE!\n")
  cat("All results have been saved in:", output_dir, "\n")
}

# =============================================================================
# EXECUTION
# =============================================================================

run_full_deconvolution_workflow(seurat_obj, deconvolution_output_dir)

cat("\n--- END OF SCRIPT 8 ---\n")