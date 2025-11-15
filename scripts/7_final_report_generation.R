#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 7: FINAL REPORT GENERATION AND VISUALIZATION
# Author: Ángel Pérez
# Date: 2 June 2025
# Objective: Generate final comparative plots with correct orientation,
#            visualize cell signatures for each pipeline, and compile a
#            comprehensive final text report.
# =============================================================================

cat("\n=== SCRIPT 7: FINAL REPORT GENERATION ===\n")

# --- 0. INITIAL SETUP ---
# Load libraries
suppressPackageStartupMessages({
  library(Seurat); library(Giotto); library(ggplot2); library(dplyr);
  library(patchwork); library(RColorBrewer); library(viridis);
  library(data.table); library(gridExtra); library(grid);
})

# --- CONFIGURATION & PATHS ---
main_results_dir <- "results"
# Source directory for integrated results
integrated_analysis_subdir_ref <- file.path(main_results_dir, "05_Integrated_Pipeline_Analysis")
# Output directory for this script
final_report_output_dir <- file.path(main_results_dir, "07_Final_Report_and_Visuals")
dir.create(final_report_output_dir, showWarnings = FALSE, recursive = TRUE)
# --- END CONFIGURATION ---

# --- 1. LOAD REQUIRED DATA OBJECTS ---
cat("1. Loading 'final_results' and 'cancer_signatures' objects...\n")

# Define file paths
final_results_file_path <- file.path(integrated_analysis_subdir_ref, "final_results_object.rds")
cancer_signatures_file_path <- file.path(integrated_analysis_subdir_ref, "cancer_signatures_object.rds")

# Load final_results
if (!file.exists(final_results_file_path)) stop("ERROR: 'final_results_object.rds' not found. Run Script 5 first.")
final_results <- readRDS(final_results_file_path)
cat("   SUCCESS: 'final_results' object loaded.\n")

# Load cancer_signatures
if (!file.exists(cancer_signatures_file_path)) stop("ERROR: 'cancer_signatures_object.rds' not found. Run Script 5 first.")
cancer_signatures <- readRDS(cancer_signatures_file_path)
cat("   SUCCESS: 'cancer_signatures' object loaded.\n")

# Verify data integrity
if (!is.list(final_results) || is.null(final_results$data) || is.null(final_results$markers)) {
  stop("ERROR: 'final_results' object is malformed. Missing 'data' or 'markers' elements.")
}

# =============================================================================
# PART 2: CREATE FINAL SPATIAL COMPARISON (CORRECTED ORIENTATION)
# =============================================================================

create_final_spatial_comparison <- function(data, output_dir) {
  cat("\n[SPATIAL] Creating final comparison with corrected orientation...\n")

  # --- Seurat Plot ---
  p_seurat <- SpatialDimPlot(data$seurat, image.alpha = 0, pt.size.factor = 1.6, stroke = 0.1) +
    ggtitle("Seurat", subtitle = paste0(data$n_clusters["seurat"], " clusters")) +
    theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

  # --- Giotto Plot (Y-axis flipped manually) ---
  giotto_plot_data <- merge(pDataDT(data$giotto), getSpatialLocations(data$giotto, output = "data.table"), by = "cell_ID")
  giotto_plot_data$sdimy_flipped <- max(giotto_plot_data$sdimy) - giotto_plot_data$sdimy
  p_giotto <- ggplot(giotto_plot_data, aes(x = sdimx, y = sdimy_flipped, color = factor(leiden_clus))) +
    geom_point(size = 0.8, alpha = 0.8) + coord_fixed() + theme_void() +
    ggtitle("Giotto", subtitle = paste0(data$n_clusters["giotto"], " clusters")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

  # --- Optimized Plot (conditionally) ---
  if (data$opt_type == "giotto") {
    opt_plot_data <- merge(pDataDT(data$optimized), getSpatialLocations(data$optimized, output = "data.table"), by = "cell_ID")
    opt_plot_data$sdimy_flipped <- max(opt_plot_data$sdimy) - opt_plot_data$sdimy
    p_optimized <- ggplot(opt_plot_data, aes(x = sdimx, y = sdimy_flipped, color = factor(leiden_clus))) +
      geom_point(size = 0.8, alpha = 0.8) + coord_fixed() + theme_void() +
      ggtitle("Optimized", subtitle = paste0(data$n_clusters["optimized"], " clusters")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  } else { # seurat
    p_optimized <- SpatialDimPlot(data$optimized, image.alpha = 0, pt.size.factor = 1.6, stroke = 0.1) +
      ggtitle("Optimized", subtitle = paste0(data$n_clusters["optimized"], " clusters")) +
      theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  }

  # --- Combine and Save ---
  combined_plot <- p_seurat + p_giotto + p_optimized +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(
      title = "Final Spatial Comparison - Consistent Orientation",
      subtitle = "All pipelines visualized without background image for direct comparison",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
    )

  ggsave(file.path(output_dir, "01_final_spatial_comparison.png"),
         combined_plot, width = 18, height = 6, dpi = 300)

  cat("   Final spatial comparison plot saved to:", output_dir, "\n")
  return(combined_plot)
}

# =============================================================================
# PART 3: VISUALIZE SIGNATURES PER PIPELINE
# =============================================================================

visualize_signatures_per_pipeline <- function(data, signatures, output_dir) {
  cat("\n[SIGNATURES] Visualizing cancer signatures for each pipeline...\n")

  # --- Helper function for Seurat-based objects ---
  visualize_seurat_signatures <- function(seurat_obj, pipeline_name) {
    plots_list <- list()
    for (sig_name in names(signatures)) {
      available_genes <- intersect(signatures[[sig_name]], rownames(seurat_obj))
      if (length(available_genes) > 0) {
        score_col <- paste0(sig_name, "_score")
        seurat_obj <- AddModuleScore(seurat_obj, features = list(available_genes), name = score_col, ctrl = 100)
        plots_list[[sig_name]] <- SpatialFeaturePlot(seurat_obj, features = paste0(score_col, "1"), pt.size.factor = 2, image.alpha = 0) +
          scale_fill_viridis(option = "plasma") + ggtitle(sig_name)
      }
    }
    return(plots_list)
  }

  # --- Helper function for Giotto-based objects ---
  visualize_giotto_signatures <- function(giotto_obj, pipeline_name) {
    plots_list <- list()
    expr_matrix <- getExpression(giotto_obj, values = "normalized", output = "matrix")
    plot_data_base <- merge(pDataDT(giotto_obj), getSpatialLocations(giotto_obj, output = "data.table"), by = "cell_ID")
    plot_data_base$sdimy_flipped <- max(plot_data_base$sdimy) - plot_data_base$sdimy

    for (sig_name in names(signatures)) {
      available_genes <- intersect(signatures[[sig_name]], rownames(expr_matrix))
      if (length(available_genes) > 0) {
        score_values <- colMeans(expr_matrix[available_genes, , drop = FALSE])
        plot_data <- plot_data_base
        plot_data$score <- score_values[match(plot_data$cell_ID, names(score_values))]
        plots_list[[sig_name]] <- ggplot(plot_data, aes(x = sdimx, y = sdimy_flipped, color = score)) +
          geom_point(size = 1.2) + scale_color_viridis(option = "plasma") + coord_fixed() + theme_void() + ggtitle(sig_name)
      }
    }
    return(plots_list)
  }

  # --- Generate plots for all pipelines ---
  all_plots <- list(
    Seurat = visualize_seurat_signatures(data$seurat, "Seurat"),
    Giotto = visualize_giotto_signatures(data$giotto, "Giotto"),
    Optimized = if (data$opt_type == "giotto") visualize_giotto_signatures(data$optimized, "Optimized") else visualize_seurat_signatures(data$optimized, "Optimized")
  )

  # --- Save combined plots ---
  for (pipeline_name in names(all_plots)) {
    if (length(all_plots[[pipeline_name]]) > 0) {
      combined_plot <- wrap_plots(all_plots[[pipeline_name]], ncol = 2) +
        plot_annotation(title = paste(pipeline_name, "- Cell Type Signatures"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
      
      file_name <- paste0("02_signatures_", tolower(pipeline_name), ".png")
      ggsave(file.path(output_dir, file_name), combined_plot, width = 12, height = 10, dpi = 300)
      cat("   Signature plot saved:", file_name, "\n")
    }
  }

  cat("   All signature visualizations saved to:", output_dir, "\n")
}

# =============================================================================
# PART 4: GENERATE FINAL TEXT REPORT
# =============================================================================

generate_final_report <- function(data, markers_data, output_dir) {
  cat("\n[REPORT] Generating final comprehensive text report...\n")

  n_markers <- sapply(markers_data, function(df) if(is.data.frame(df)) nrow(df) else 0)

  # --- Report Content ---
  report_text <- paste0(
    "===============================================================\n",
    "             FINAL REPORT - MASTER'S THESIS PROJECT            \n",
    "===============================================================\n\n",
    "Title: Optimization of Spatial Transcriptomics Pipelines for Solid Tumor Analysis\n",
    "Author: Ángel Ignacio Pérez Santiago\n",
    "Date: ", Sys.Date(), "\n\n",
    "---------------------------------------------------------------\n",
    "1. EXECUTIVE SUMMARY\n",
    "---------------------------------------------------------------\n\n",
    "This project presents a comprehensive comparison of three analysis pipelines\n",
    "for spatial transcriptomics applied to the study of lung cancer:\n",
    "1. Seurat Pipeline (Industry Standard)\n",
    "2. Giotto Pipeline (Specialized Spatial Tool)\n",
    "3. Optimized Pipeline (Hybrid strategy developed in this project)\n\n",
    "---------------------------------------------------------------\n",
    "2. KEY RESULTS & PIPELINE COMPARISON\n",
    "---------------------------------------------------------------\n\n",
    "SEURAT:\n",
    "  - Clusters identified: ", data$n_clusters["seurat"], "\n",
    "  - Markers detected: ", n_markers["seurat"], "\n",
    "  - Strengths: Fast, efficient, well-established.\n\n",
    "GIOTTO:\n",
    "  - Clusters identified: ", data$n_clusters["giotto"], " (+", data$n_clusters["giotto"] - data$n_clusters["seurat"], " vs Seurat)\n",
    "  - Markers detected: ", n_markers["giotto"], "\n",
    "  - Strengths: High granularity, advanced spatial analytics.\n\n",
    "OPTIMIZED (Type: ", data$opt_type, "):\n",
    "  - Clusters identified: ", data$n_clusters["optimized"], " (+", data$n_clusters["optimized"] - data$n_clusters["seurat"], " vs Seurat)\n",
    "  - Markers detected: ", n_markers["optimized"], "\n",
    "  - Strengths: Automated decision-making, optimal balance of speed and detail.\n\n",
    "---------------------------------------------------------------\n",
    "3. CONCLUSIONS & RECOMMENDATIONS\n",
    "---------------------------------------------------------------\n\n",
    "The Optimized Pipeline represents a significant improvement over individual\n",
    "methods by intelligently combining the strengths of both Seurat and Giotto.\n",
    "Its adaptive strategy ensures efficient use of computational resources while\n",
    "maximizing analytical granularity, enabling a more detailed characterization\n",
    "of the tumor microenvironment.\n\n",
    "RECOMMENDED USAGE:\n",
    "- For exploratory analysis: Use the Optimized Pipeline.\n",
    "- For rapid screening: Use the Seurat Pipeline.\n",
    "- For deep spatial analysis: Use the Giotto Pipeline.\n\n",
    "===============================================================\n",
    "                         END OF REPORT                         \n",
    "===============================================================\n"
  )

  # --- Save Report ---
  writeLines(report_text, file.path(output_dir, "FINAL_REPORT.txt"))
  cat("   Final text report saved to:", output_dir, "\n")
}

# =============================================================================
# MAIN FUNCTION: EXECUTE ALL FINAL TASKS
# =============================================================================

run_final_analysis_suite <- function(final_results_obj, signatures_obj, output_dir) {
  cat("\n-----------------------------------------------------------------\n")
  cat("              RUNNING FINAL ANALYSIS & REPORTING SUITE           \n")
  cat("-----------------------------------------------------------------\n\n")

  # Phase 1: Spatial Comparison
  create_final_spatial_comparison(final_results_obj$data, output_dir)

  # Phase 2: Signatures per Pipeline
  visualize_signatures_per_pipeline(final_results_obj$data, signatures_obj, output_dir)

  # Phase 3: Final Report
  generate_final_report(final_results_obj$data, final_results_obj$markers, output_dir)

  cat("\nFINAL ANALYSIS SUITE COMPLETED SUCCESSFULLY\n")
  cat("\nAll final outputs generated in:", output_dir, "\n")
}

# =============================================================================
# EXECUTION
# =============================================================================

run_final_analysis_suite(final_results, cancer_signatures, final_report_output_dir)

cat("\n--- END OF SCRIPT 7 ---\n")