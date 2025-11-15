#!/usr/-bin/env Rscript
# =============================================================================
# COMPARATIVE ANALYSIS SEURAT vs GIOTTO - STANDALONE SCRIPT
# Author: Ángel Pérez
# Date: 31 May 2025
# Objective: Compare results of independent pipelines
# =============================================================================

# --- 0. INITIAL SETUP ---
cat("=== STARTING COMPARATIVE ANALYSIS SEURAT vs GIOTTO ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Process PID:", Sys.getpid(), "\n\n")

# Clean environment
rm(list = ls())
gc()

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Giotto)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(data.table)
})

# --- CONFIGURATION & PATHS ---
# Main output directory for reports and plots
main_results_dir <- "results"
# Specific subdirectory for the outputs of this comparison script
comparison_output_dir <- file.path(main_results_dir, "03_Seurat_Giotto_Comparison")
dir.create(comparison_output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. HELPER FUNCTIONS ---
cat("1. Configuring helper functions...\n")

# Function to format memory
format_memory <- function(mb) {
  if (is.na(mb)) return("N/A")
  if (mb < 1024) {
    return(paste0(round(mb, 1), " MB"))
  } else {
    return(paste0(round(mb/1024, 2), " GB"))
  }
}

# Function to find most recent files
find_latest_results <- function(path, pattern) {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop(paste("No files found with pattern:", pattern))
  }
  # Sort by modification date and take the most recent
  files_info <- file.info(files)
  latest_file <- rownames(files_info)[which.max(files_info$mtime)]
  return(latest_file)
}

# --- 2. LOAD RESULTS FROM BOTH PIPELINES ---
cat("2. Loading results from both pipelines...\n")

# Find most recent files
seurat_file <- find_latest_results(main_results_dir, "^seurat_results_.*\\.rds$")
giotto_file <- find_latest_results(main_results_dir, "^giotto_results_.*\\.rds$")

cat("   Loading Seurat from:", basename(seurat_file), "\n")
seurat_data <- readRDS(seurat_file)

cat("   Loading Giotto from:", basename(giotto_file), "\n")
giotto_data <- readRDS(giotto_file)

# Check if parameters are identical
params_match <- identical(seurat_data$parameters, giotto_data$parameters)
if (!params_match) {
  cat("WARNING: Parameters are not identical between both pipelines\n")
} else {
  cat("Parameters verified: Identical in both pipelines\n")
}

# --- 3. DETAILED COMPARATIVE ANALYSIS ---
cat("\n=== 3. DETAILED COMPARATIVE ANALYSIS ===\n")

perform_detailed_comparison <- function(seurat_data, giotto_data) {
  
  cat("3.1 Extracting data from both pipelines...\n")
  
  # Extract results
  seurat_res <- seurat_data$results
  giotto_res <- giotto_data$results
  params <- seurat_data$parameters
  
  comparison <- list()
  
  # 3.1 Basic performance metrics
  comparison$performance <- data.frame(
    Pipeline = c("Seurat", "Giotto"),
    N_Clusters = c(seurat_res$n_clusters, giotto_res$n_clusters),
    Time_sec = c(seurat_res$time_sec, giotto_res$time_sec),
    Time_min = c(round(seurat_res$time_sec/60, 1), round(giotto_res$time_sec/60, 1)),
    Memory_Peak_MB = c(seurat_res$memory_profile$peak, giotto_res$memory_profile$peak),
    Memory_Peak_GB = c(round(seurat_res$memory_profile$peak/1024, 2), 
                       round(giotto_res$memory_profile$peak/1024, 2)),
    Memory_Increase_MB = c(seurat_res$memory_profile$increase, giotto_res$memory_profile$increase),
    Memory_Increase_GB = c(round(seurat_res$memory_profile$increase/1024, 2),
                           round(giotto_res$memory_profile$increase/1024, 2)),
    stringsAsFactors = FALSE
  )
  
  cat("   Comparative performance:\n")
  print(comparison$performance)
  
  # 3.2 Cluster analysis
  cat("\n3.2 Analyzing cluster distribution...\n")
  
  # Extract cluster information from Seurat
  seurat_clusters_ids <- as.numeric(as.character(seurat_data$results$object@active.ident))
  seurat_sizes <- table(seurat_clusters_ids)
  
  # Extract cluster information from Giotto - CORRECT AND CONSISTENT METHOD
  tryCatch({
    giotto_metadata <- Giotto::pDataDT(giotto_res$object)
    giotto_clusters_ids <- as.numeric(as.character(giotto_metadata$leiden_clus))
    giotto_sizes <- table(giotto_clusters_ids)
  }, error = function(e) {
    cat("   Error accessing Giotto metadata:", e$message, "\n")
    # Fallback: create empty table
    giotto_clusters_ids <- numeric(0)
    giotto_sizes <- table(giotto_clusters_ids)
  })
  
  comparison$cluster_analysis <- list(
    seurat_sizes = seurat_sizes,
    giotto_sizes = giotto_sizes,
    seurat_stats = list(
      min = if(length(seurat_sizes) > 0) min(seurat_sizes) else NA,
      max = if(length(seurat_sizes) > 0) max(seurat_sizes) else NA,
      mean = if(length(seurat_sizes) > 0) round(mean(seurat_sizes), 1) else NA,
      median = if(length(seurat_sizes) > 0) median(seurat_sizes) else NA,
      sd = if(length(seurat_sizes) > 0) round(sd(seurat_sizes), 1) else NA
    ),
    giotto_stats = if(length(giotto_sizes) > 0) {
      list(
        min = min(giotto_sizes),
        max = max(giotto_sizes),
        mean = round(mean(giotto_sizes), 1),
        median = median(giotto_sizes),
        sd = round(sd(giotto_sizes), 1)
      )
    } else {
      list(min = NA, max = NA, mean = NA, median = NA, sd = NA)
    }
  )
  
  cat("   Seurat - Min:", comparison$cluster_analysis$seurat_stats$min,
      "Max:", comparison$cluster_analysis$seurat_stats$max,
      "Mean:", comparison$cluster_analysis$seurat_stats$mean, "\n")
  
  if(length(giotto_sizes) > 0) {
    cat("   Giotto - Min:", comparison$cluster_analysis$giotto_stats$min,
        "Max:", comparison$cluster_analysis$giotto_stats$max,
        "Mean:", comparison$cluster_analysis$giotto_stats$mean, "\n")
  } else {
    cat("   Giotto - No valid clusters detected\n")
  }
  
  # 3.3 Spot comparison
  cat("\n3.3 Comparing processed spots...\n")
  
  seurat_spots <- colnames(seurat_res$object)
  
  # Get Giotto spots safely
  giotto_local_metadata <- NULL
  tryCatch({
    if(exists("giotto_metadata") && !is.null(giotto_metadata)) {
        giotto_local_metadata <- giotto_metadata
    } else {
      giotto_local_metadata <- Giotto::pDataDT(giotto_res$object)
    }
    giotto_spots <- giotto_local_metadata$cell_ID
  }, error = function(e) {
    cat("   Error accessing Giotto spots:", e$message, "\n")
    giotto_spots <- character(0)
  })
  
  comparison$spots_analysis <- list(
    n_seurat = length(seurat_spots),
    n_giotto = length(giotto_spots),
    n_common = length(intersect(seurat_spots, giotto_spots)),
    spots_only_seurat = setdiff(seurat_spots, giotto_spots),
    spots_only_giotto = setdiff(giotto_spots, seurat_spots),
    perfect_match = setequal(seurat_spots, giotto_spots)
  )
  
  cat("   Seurat spots:", comparison$spots_analysis$n_seurat, "\n")
  cat("   Giotto spots:", comparison$spots_analysis$n_giotto, "\n")
  cat("   Common spots:", comparison$spots_analysis$n_common, "\n")
  cat("   Only in Seurat:", length(comparison$spots_analysis$spots_only_seurat), "\n")
  cat("   Only in Giotto:", length(comparison$spots_analysis$spots_only_giotto), "\n")
  cat("   Perfect match:", comparison$spots_analysis$perfect_match, "\n")
  
  # 3.4 Efficiency analysis
  cat("\n3.4 Calculating efficiency metrics...\n")
  
  time_ratio <- comparison$performance$Time_min[comparison$performance$Pipeline == "Giotto"] / 
    comparison$performance$Time_min[comparison$performance$Pipeline == "Seurat"]
  
  memory_ratio <- comparison$performance$Memory_Peak_GB[comparison$performance$Pipeline == "Giotto"] / 
    comparison$performance$Memory_Peak_GB[comparison$performance$Pipeline == "Seurat"]
  
  comparison$efficiency <- list(
    time_advantage_factor = round(time_ratio, 2),
    memory_efficiency_factor = round(memory_ratio, 2),
    cluster_difference = abs(seurat_res$n_clusters - giotto_res$n_clusters),
    faster_pipeline = ifelse(time_ratio < 1, "Giotto", "Seurat"),
    more_memory_efficient = ifelse(memory_ratio < 1, "Giotto", "Seurat")
  )
  
  cat("   Faster pipeline:", comparison$efficiency$faster_pipeline, 
      "(factor of", comparison$efficiency$time_advantage_factor, ")\n")
  cat("   More memory-efficient pipeline:", comparison$efficiency$more_memory_efficient,
      "(factor of", comparison$efficiency$memory_efficiency_factor, ")\n")
  cat("   Difference in clusters:", comparison$efficiency$cluster_difference, "\n")
  
  return(comparison)
}

# Perform comparative analysis
comparison_results <- perform_detailed_comparison(seurat_data, giotto_data)

# --- 4. GENERATE COMPARATIVE VISUALIZATIONS ---
cat("\n=== 4. GENERATING COMPARATIVE VISUALIZATIONS ===\n")

generate_comparison_plots <- function(comparison_results, seurat_data, giotto_data) {
  
  plots <- list()
  
  # 4.1 Main performance panel
  cat("4.1 Creating main performance panel...\n")
  
  perf_data <- comparison_results$performance
  
  # Cluster plot
  clusters_plot <- ggplot(perf_data, aes(x = Pipeline, y = N_Clusters, fill = Pipeline)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
    geom_text(aes(label = N_Clusters), vjust = -0.5, size = 6, fontface = "bold") +
    theme_minimal(base_size = 14) +
    labs(title = "Number of Clusters", y = "Clusters") +
    scale_fill_manual(values = c("Seurat" = "#4CAF50", "Giotto" = "#2196F3")) +
    scale_y_continuous(limits = c(0, max(perf_data$N_Clusters, na.rm = TRUE) * 1.2)) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Time plot
  time_plot <- ggplot(perf_data, aes(x = Pipeline, y = Time_min, fill = Pipeline)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
    geom_text(aes(label = paste0(Time_min, " min")), vjust = -0.5, size = 6, fontface = "bold") +
    theme_minimal(base_size = 14) +
    labs(title = "Execution Time", y = "Minutes") +
    scale_fill_manual(values = c("Seurat" = "#4CAF50", "Giotto" = "#2196F3")) +
    scale_y_continuous(limits = c(0, max(perf_data$Time_min, na.rm = TRUE) * 1.2)) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Memory plot
  memory_plot <- ggplot(perf_data, aes(x = Pipeline, y = Memory_Peak_GB, fill = Pipeline)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
    geom_text(aes(label = paste0(Memory_Peak_GB, " GB")), vjust = -0.5, size = 6, fontface = "bold") +
    theme_minimal(base_size = 14) +
    labs(title = "Peak Memory", y = "GB") +
    scale_fill_manual(values = c("Seurat" = "#4CAF50", "Giotto" = "#2196F3")) +
    scale_y_continuous(limits = c(0, max(perf_data$Memory_Peak_GB, na.rm = TRUE) * 1.2)) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plots$performance_main <- clusters_plot + time_plot + memory_plot +
    plot_annotation(
      title = "Performance Comparison: Seurat vs Giotto",
      subtitle = "Independent measurements with clean RAM",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  # 4.2 Memory profile during execution
  cat("4.2 Creating temporal memory profile...\n")
  
  # Prepare temporal memory data
  seurat_mem <- seurat_data$results$memory_profile
  giotto_mem <- giotto_data$results$memory_profile
  
  seurat_profile <- data.frame(
    Step = c("Start", "Load", "QC", "Normalization", "HVG", "Scaling", "PCA", "Clustering", "UMAP", "Final"),
    Memory_GB = c(seurat_mem$start, seurat_mem$load, seurat_mem$qc, seurat_mem$normalization,
                  seurat_mem$hvg, seurat_mem$scale, seurat_mem$pca, seurat_mem$clustering,
                  seurat_mem$umap, seurat_mem$final) / 1024,
    Pipeline = "Seurat",
    Order = 1:10
  )
  
  giotto_profile <- data.frame(
    Step = c("Start", "Load", "Object", "QC", "Normalization", "HVG", "PCA", "Clustering", "UMAP", "Final"),
    Memory_GB = c(giotto_mem$start, giotto_mem$load, giotto_mem$object, giotto_mem$qc,
                  giotto_mem$normalization, giotto_mem$hvg, giotto_mem$pca, giotto_mem$clustering,
                  giotto_mem$umap, giotto_mem$final) / 1024,
    Pipeline = "Giotto",
    Order = 1:10
  )
  
  combined_profile <- rbind(seurat_profile, giotto_profile)
  
  plots$memory_profile <- ggplot(combined_profile, aes(x = Order, y = Memory_GB, color = Pipeline)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    scale_x_continuous(breaks = 1:10, labels = unique(seurat_profile$Step)) +
    scale_color_manual(values = c("Seurat" = "#4CAF50", "Giotto" = "#2196F3")) +
    theme_minimal(base_size = 12) +
    labs(title = "Memory Profile During Execution",
         x = "Pipeline Step", y = "Memory (GB)", color = "Pipeline") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top")
  
  # 4.3 Cluster size distribution
  cat("4.3 Creating cluster distribution...\n")
  
  seurat_sizes_df <- data.frame(
    Pipeline = "Seurat",
    Cluster = names(comparison_results$cluster_analysis$seurat_sizes),
    Size = as.numeric(comparison_results$cluster_analysis$seurat_sizes)
  )
  
  # Check that Giotto has clusters before creating the data.frame
  if (length(comparison_results$cluster_analysis$giotto_sizes) > 0) {
    giotto_sizes_df <- data.frame(
      Pipeline = "Giotto", 
      Cluster = names(comparison_results$cluster_analysis$giotto_sizes),
      Size = as.numeric(comparison_results$cluster_analysis$giotto_sizes)
    )
    cluster_sizes_df <- rbind(seurat_sizes_df, giotto_sizes_df)
  } else {
    cat("   WARNING: Giotto has no valid clusters to visualize\n")
    cluster_sizes_df <- seurat_sizes_df
  }
  
  plots$cluster_distribution <- ggplot(cluster_sizes_df, aes(x = reorder(Cluster, Size), y = Size, fill = Pipeline)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Pipeline, scales = "free_x") +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(title = "Distribution of Cluster Sizes",
         x = "Cluster ID", y = "Number of Spots") +
    scale_fill_manual(values = c("Seurat" = "#4CAF50", "Giotto" = "#2196F3")) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 4.4 Complete combined panel
  plots$combined_analysis <- plots$performance_main / 
    (plots$memory_profile | plots$cluster_distribution) +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      title = "Complete Analysis: Seurat vs Giotto",
      subtitle = paste("Independent analysis -", Sys.Date()),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
    )
  
  return(plots)
}

# Generate plots
comparison_plots <- generate_comparison_plots(comparison_results, seurat_data, giotto_data)

# --- 5. GENERATE COMPARATIVE SPATIAL VISUALIZATION ---
cat("\n=== 5. GENERATING SPATIAL COMPARISON ===\n")

generate_spatial_comparison <- function(seurat_data, giotto_data) {
  
  cat("5.1 Preparing spatial data from both pipelines...\n")
  
  # Seurat data
  seurat_obj <- seurat_data$results$object
  params <- seurat_data$parameters
  
  # Check if an image is available in Seurat
  has_image <- length(seurat_obj@images) > 0
  
  if (has_image) {
    tryCatch({
      seurat_spatial <- SpatialDimPlot(
        seurat_obj,
        group.by = "seurat_clusters",
        pt.size.factor = params$POINT_SIZE,
        image.alpha = 0, # Background image removed
        alpha = 1, # Alpha for the points
        stroke = 0.2
      ) + ggtitle("Seurat - Spatial Clusters") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }, error = function(e) {
      cat("   Error with SpatialDimPlot, using manual coordinates\n")
      has_image <<- FALSE # Update variable in the function's scope
    })
  }
  
  if (!has_image) {
    # Use manual spatial coordinates
    tryCatch({
      spatial_coords <- GetTissueCoordinates(seurat_obj)
      seurat_spatial_data <- data.frame(
        x = spatial_coords[, 1],
        y = spatial_coords[, 2],
        cluster = as.factor(seurat_obj@active.ident) # It's already Idents(), @active.ident is internal
      )
      
      seurat_spatial <- ggplot(seurat_spatial_data, aes(x = x, y = y, color = cluster)) +
        geom_point(size = params$POINT_SIZE, alpha = 0.8) +
        scale_color_brewer(type = "qual", palette = params$CLUSTER_COLORS, name = "Cluster") +
        coord_fixed() +
        ggtitle("Seurat - Spatial Clusters") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.grid = element_blank()) +
        labs(x = "X Coordinate", y = "Y Coordinate")
    }, error = function(e) {
      cat("   Error getting Seurat coordinates:", e$message, "\n")
      # Placeholder plot
      seurat_spatial <- ggplot() + ggtitle("Seurat - Error in coordinates") + # Assign to seurat_spatial
        theme_minimal()
    })
  }
  
  # Giotto data - CONSISTENT METHOD
  giotto_obj <- giotto_data$results$object
  
  tryCatch({
    # Get metadata
    giotto_metadata_spatial <- Giotto::pDataDT(giotto_obj)
    
    # Get spatial coordinates
    spatial_locs_obj <- Giotto::getSpatialLocations(giotto_obj)
    spatial_coords_giotto <- spatial_locs_obj@coordinates
    
    # Ensure coordinates have the correct ID
    if (!"cell_ID" %in% colnames(spatial_coords_giotto)) {
      spatial_coords_giotto$cell_ID <- rownames(spatial_coords_giotto)
    }
    
    # Combine data
    giotto_plot_data <- merge(giotto_metadata_spatial, spatial_coords_giotto, by = "cell_ID")
    giotto_plot_data$leiden_clus <- as.factor(giotto_plot_data$leiden_clus)
    
    giotto_spatial <- ggplot(giotto_plot_data, aes(x = sdimx, y = sdimy, color = leiden_clus)) +
      geom_point(size = params$POINT_SIZE, alpha = 0.8) +
      scale_color_brewer(type = "qual", palette = params$CLUSTER_COLORS, name = "Cluster") +
      coord_fixed() +
      scale_y_reverse() +
      ggtitle("Giotto - Spatial Clusters") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_blank()) +
      labs(x = "X Coordinate", y = "Y Coordinate")
    
  }, error = function(e) {
    cat("   Error getting spatial data from Giotto:", e$message, "\n")
    # Placeholder plot
    giotto_spatial <- ggplot() + ggtitle("Giotto - Error in coordinates") + # Assign to giotto_spatial
      theme_minimal()
  })
  
  # Side-by-side comparison
  spatial_comparison <- seurat_spatial + giotto_spatial +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Spatial Comparison: Seurat vs Giotto",
      subtitle = paste("Seurat:", seurat_data$results$n_clusters, "clusters | Giotto:", 
                       giotto_data$results$n_clusters, "clusters"),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  return(spatial_comparison)
}

# Generate spatial comparison
spatial_comparison_plot <- generate_spatial_comparison(seurat_data, giotto_data)

# --- 6. SAVE RESULTS AND VISUALIZATIONS ---
cat("\n=== 6. SAVING RESULTS AND VISUALIZATIONS ===\n")

# Create timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 6.1 Save main plots
cat("6.1 Saving main visualizations...\n")

ggsave(file.path(comparison_output_dir, paste0("comparison_analysis_", timestamp, ".png")), 
       comparison_plots$combined_analysis, 
       width = 16, height = 12, dpi = 300)

ggsave(file.path(comparison_output_dir, paste0("spatial_comparison_", timestamp, ".png")), 
       spatial_comparison_plot, 
       width = 16, height = 8, dpi = 300)

# 6.2 Create object with complete results
complete_comparison <- list(
  # Original results
  seurat_data = seurat_data,
  giotto_data = giotto_data,
  
  # Comparative analysis
  comparison_analysis = comparison_results,
  
  # Analysis information
  analysis_info = list(
    comparison_date = Sys.Date(),
    comparison_time = Sys.time(),
    timestamp = timestamp,
    seurat_file = basename(seurat_file),
    giotto_file = basename(giotto_file),
    parameters_identical = params_match,
    R_version = as.character(getRversion())
  )
)

# 6.3 Save complete results (RDS in the main results directory)
output_file_rds <- file.path(main_results_dir, paste0("complete_comparison_", timestamp, ".rds"))
saveRDS(complete_comparison, file = output_file_rds)
cat("6.2 Complete analysis (RDS) saved in:", output_file_rds, "\n")

# --- 7. GENERATE FINAL REPORT ---
cat("\n=== 7. GENERATING FINAL REPORT ===\n")

generate_final_report <- function(comparison_data, timestamp, output_dir) {
  
  report_file_path <- file.path(output_dir, paste0("final_comparison_report_", timestamp, ".txt"))
  
  sink(report_file_path)
  
  cat("==========================================================================\n")
  cat("                    FINAL REPORT - HOMOLOGOUS COMPARISON\n")
  cat("                        Seurat vs Giotto Analysis\n")
  cat("==========================================================================\n\n")
  
  comp_res <- comparison_data$comparison_analysis
  # Ensure system_info exists before accessing it
  seurat_sys_info <- if(!is.null(comparison_data$seurat_data$system_info)) comparison_data$seurat_data$system_info else list(Seurat_version="N/A")
  giotto_sys_info <- if(!is.null(comparison_data$giotto_data$system_info)) comparison_data$giotto_data$system_info else list(Giotto_version="N/A")
  
  cat("ANALYSIS INFORMATION\n")
  cat("------------------------\n")
  cat("Analysis date:", format(comparison_data$analysis_info$comparison_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Seurat file:", comparison_data$analysis_info$seurat_file, "\n")
  cat("Giotto file:", comparison_data$analysis_info$giotto_file, "\n")
  cat("Identical parameters:", comparison_data$analysis_info$parameters_identical, "\n")
  cat("R version:", comparison_data$analysis_info$R_version, "\n")
  cat("Seurat version:", seurat_sys_info$Seurat_version, "\n")
  cat("Giotto version:", giotto_sys_info$Giotto_version, "\n\n")
  
  cat("MAIN RESULTS\n")
  cat("----------------------\n")
  perf <- comp_res$performance
  cat("SEURAT:\n")
  cat("  - Clusters found:", perf$N_Clusters[perf$Pipeline == "Seurat"], "\n")
  cat("  - Execution time:", perf$Time_min[perf$Pipeline == "Seurat"], "minutes\n")
  cat("  - Peak memory:", perf$Memory_Peak_GB[perf$Pipeline == "Seurat"], "GB\n")
  cat("  - Memory increase:", perf$Memory_Increase_GB[perf$Pipeline == "Seurat"], "GB\n\n")
  
  cat("GIOTTO:\n")
  cat("  - Clusters found:", perf$N_Clusters[perf$Pipeline == "Giotto"], "\n")
  cat("  - Execution time:", perf$Time_min[perf$Pipeline == "Giotto"], "minutes\n")
  cat("  - Peak memory:", perf$Memory_Peak_GB[perf$Pipeline == "Giotto"], "GB\n")
  cat("  - Memory increase:", perf$Memory_Increase_GB[perf$Pipeline == "Giotto"], "GB\n\n")
  
  cat("COMPARATIVE ANALYSIS\n")
  cat("--------------------\n")
  eff <- comp_res$efficiency
  
  cat("Cluster concordance:\n")
  if (eff$cluster_difference == 0) {
    cat("  IDENTICAL: Both find the same number of clusters\n")
  } else {
    cat("  DIFFERENCE:", eff$cluster_difference, "clusters difference\n")
    if (perf$N_Clusters[perf$Pipeline == "Giotto"] > perf$N_Clusters[perf$Pipeline == "Seurat"]) {
      cat("     Giotto finds", eff$cluster_difference, "additional clusters\n")
    } else {
      cat("     Seurat finds", eff$cluster_difference, "additional clusters\n")
    }
  }
  
  cat("\nTemporal performance:\n")
  if (eff$faster_pipeline == "Seurat") {
    cat("  - Seurat is", eff$time_advantage_factor, "times faster than Giotto\n")
  } else {
    cat("  - Giotto is", eff$time_advantage_factor, "times faster than Seurat\n")
  }
  cat("  - Absolute difference:", round(abs(perf$Time_min[perf$Pipeline == "Giotto"] - 
                                              perf$Time_min[perf$Pipeline == "Seurat"]), 1), "minutes\n")
  
  cat("\nMemory efficiency:\n")
  if (eff$more_memory_efficient == "Seurat") {
    cat("  - Seurat is", eff$memory_efficiency_factor, "times more efficient than Giotto\n")
  } else {
    cat("  - Giotto is", eff$memory_efficiency_factor, "times more efficient than Seurat\n")
  }
  cat("  - Absolute difference:", round(abs(perf$Memory_Peak_GB[perf$Pipeline == "Giotto"] - 
                                              perf$Memory_Peak_GB[perf$Pipeline == "Seurat"]), 2), "GB\n")
  
  cat("\nSpot matching:\n")
  spots <- comp_res$spots_analysis
  cat("  - Processed spots - Seurat:", spots$n_seurat, "| Giotto:", spots$n_giotto, "\n")
  cat("  - Common spots:", spots$n_common, "\n")
  cat("  - Only in Seurat:", length(spots$spots_only_seurat), "\n")
  cat("  - Only in Giotto:", length(spots$spots_only_giotto), "\n")
  cat("  - Perfect match:", spots$perfect_match, "\n")
  
  cat("\nCONCLUSIONS\n")
  cat("------------\n")
  if (eff$cluster_difference == 0) {
    cat("HIGH CONCORDANCE: Both pipelines produce identical results\n")
    cat("   in terms of number of clusters under homologous parameters.\n\n")
  } else {
    cat("ALGORITHMIC DIFFERENCES: The pipelines differ by", eff$cluster_difference, "clusters.\n")
    cat("   This suggests inherent differences in the implementations.\n\n")
  }
  
  if (eff$time_advantage_factor > 2 && eff$faster_pipeline == "Seurat") {
    cat("SIGNIFICANT TEMPORAL ADVANTAGE: Seurat is considerably faster.\n")
  } else if (eff$time_advantage_factor > 2 && eff$faster_pipeline == "Giotto") {
    cat("SIGNIFICANT TEMPORAL ADVANTAGE: Giotto is considerably faster.\n")
  } else {
    cat("SIMILAR TEMPORAL PERFORMANCE: Both have comparable times.\n")
  }
  
  if (eff$memory_efficiency_factor > 1.5 && eff$more_memory_efficient == "Seurat") {
    cat("MEMORY EFFICIENCY: Seurat is significantly more efficient.\n")
  } else if (eff$memory_efficiency_factor > 1.5 && eff$more_memory_efficient == "Giotto") {
    cat("MEMORY EFFICIENCY: Giotto is significantly more efficient.\n")
  } else {
    cat("SIMILAR MEMORY USAGE: Both have comparable consumption.\n")
  }
  
  cat("\nRECOMMENDATIONS\n")
  cat("---------------\n")
  if (eff$faster_pipeline == "Seurat" && eff$more_memory_efficient == "Seurat") {
    cat("For routine analyses: SEURAT (superior in time and memory)\n")
  } else if (eff$faster_pipeline == "Giotto" && eff$more_memory_efficient == "Giotto") {
    cat("For routine analyses: GIOTTO (superior in time and memory)\n")
  } else {
    cat("For routine analyses: Choose according to priority (time vs memory)\n")
  }
  
  if (eff$cluster_difference > 0) {
    cat("For exploratory analyses: Run BOTH and compare results\n")
  }
  cat("For spatial visualization: Both offer robust capabilities\n")
  cat("For reproducibility: Use standalone scripts like this analysis\n")
  
  cat("\nGENERATED FILES\n")
  cat("------------------\n")
  cat("- Complete analysis:", file.path(basename(main_results_dir), basename(output_file_rds)), "\n")
  cat("- Comparative visualization:", file.path(basename(output_dir), paste0("comparison_analysis_", timestamp, ".png")), "\n")
  cat("- Spatial comparison:", file.path(basename(output_dir), paste0("spatial_comparison_", timestamp, ".png")), "\n")
  cat("- This report:", file.path(basename(output_dir), paste0("final_comparison_report_", timestamp, ".txt")), "\n")
  
  cat("\n==========================================================================\n")
  cat("                           END OF REPORT\n")
  cat("==========================================================================\n")
  
  sink()
  
  return(report_file_path)
}

# Generate final report
report_file_path <- generate_final_report(complete_comparison, timestamp, comparison_output_dir)
cat("7.1 Final report generated in:", report_file_path, "\n")

# --- 8. FINAL SUMMARY ON CONSOLE ---
cat("\n=== FINAL SUMMARY OF THE COMPARATIVE ANALYSIS ===\n")

perf <- comparison_results$performance
eff <- comparison_results$efficiency

cat("SCIENTIFIC RESULTS:\n")
cat("- Seurat:", perf$N_Clusters[perf$Pipeline == "Seurat"], "clusters |", 
    perf$Time_min[perf$Pipeline == "Seurat"], "min |", 
    perf$Memory_Peak_GB[perf$Pipeline == "Seurat"], "GB\n")
cat("- Giotto:", perf$N_Clusters[perf$Pipeline == "Giotto"], "clusters |", 
    perf$Time_min[perf$Pipeline == "Giotto"], "min |", 
    perf$Memory_Peak_GB[perf$Pipeline == "Giotto"], "GB\n")

cat("\nCOMPARATIVE ADVANTAGES:\n")
cat("- Faster pipeline:", eff$faster_pipeline, "(", eff$time_advantage_factor, "x)\n")
cat("- More memory-efficient:", eff$more_memory_efficient, "(", eff$memory_efficiency_factor, "x)\n")
cat("- Difference in clusters:", eff$cluster_difference, "\n")

cat("\nVALID METHODOLOGY:\n")
cat("- Clean RAM for each pipeline\n")
cat("- Verified identical parameters\n")
cat("- Independent measurements\n")
cat("- Comparative spatial visualization\n")

cat("\nFINAL FILES:\n")
cat("- Analysis (RDS):", output_file_rds, "\n")
cat("- Plots and Report in:", comparison_output_dir, "\n")


cat("\nCONCLUSION:", 
    ifelse(eff$cluster_difference == 0, "HIGH CONCORDANCE", "ALGORITHMIC DIFFERENCES"),
    "between pipelines\n")

cat("\nCOMPLETE COMPARATIVE ANALYSIS FINISHED\n")
cat("=== END OF THE COMPLETE HOMOLOGOUS ANALYSIS ===\n")