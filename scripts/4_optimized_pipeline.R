#!/usr/bin/env Rscript
# =============================================================================
# OPTIMIZED HYBRID PIPELINE: SEURAT-GIOTTO
# Author: Ángel Pérez
# Date: 1 June 2025
# Objective: Combine advantages of Seurat (efficiency) and Giotto (granularity)
# =============================================================================

# --- 0. PHILOSOPHY OF THE OPTIMIZED PIPELINE ---
cat("=== OPTIMIZED SEURAT-GIOTTO PIPELINE ===\n")
cat("STRATEGY: Seurat's Efficiency + Giotto's Granularity\n")
cat("DECISION: Quick initial analysis - Selective in-depth analysis\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Giotto)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(future)
  library(pryr)
  library(arrow)
  library(data.table)
  library(scales)
  library(RColorBrewer)
})

# Configuration
set.seed(12345)
plan("sequential")

# --- CONFIGURATION & PATHS ---
# Main output directory
main_results_dir <- "results"
# Specific subdirectory for the outputs of this script (reports and plots)
optimized_output_subdir <- file.path(main_results_dir, "04_Optimized_Pipeline_Outputs")
dir.create(optimized_output_subdir, showWarnings = FALSE, recursive = TRUE)
# Input data directory
data_dir <- "data/visium_dataset"
# --- END CONFIGURATION ---


# --- 1. MONITORING FUNCTIONS ---
cat("1. Configuring monitoring system...\n")

get_memory_usage <- function() {
  if (Sys.info()["sysname"] == "Linux") {
    status_file <- "/proc/self/status"
    if (file.exists(status_file)) {
      status_lines <- readLines(status_file)
      vmrss_line <- grep("^VmRSS:", status_lines, value = TRUE)
      if (length(vmrss_line) > 0) {
        mem_kb <- as.numeric(gsub(".*VmRSS:\\s*(\\d+)\\s*kB.*", "\\1", vmrss_line))
        return(mem_kb / 1024)
      }
    }
  }
  tryCatch({
    mem_info <- pryr::mem_used()
    return(as.numeric(mem_info) / 1024^2)
  }, error = function(e) return(NA))
}

format_memory <- function(mb) {
  if (is.na(mb)) return("N/A")
  if (mb < 1024) return(paste0(round(mb, 1), " MB"))
  else return(paste0(round(mb/1024, 2), " GB"))
}

memory_monitor <- function(label = "Checkpoint") {
  current_mem <- get_memory_usage()
  cat("   [MEMORY]", label, ":", format_memory(current_mem), "\n")
  return(current_mem)
}

# --- 2. OPTIMIZED PIPELINE PARAMETERS ---
cat("2. Defining optimized pipeline parameters...\n")

optimized_params <- list(
  # sample_name
  sample_name = "Visium_HD_Human_Lung_Cancer_16um",
  
  # Common QC
  MIN_FEATURES = 500,
  MIN_COUNTS = 1000,
  MAX_MITO_PERCENT = 20,
  MIN_GENES_PER_SPOT = 3,
  
  # Normalization
  NORMALIZATION_METHOD = "LogNormalize",
  SCALE_FACTOR = 10000,
  N_HVGS = 3000,
  
  # PCA
  N_PCS = 50,
  N_PCS_USE = 20,
  
  # KEY DECISION: Clustering parameters
  SEURAT_RESOLUTION = 0.5,      # Conservative for initial analysis
  GIOTTO_RESOLUTION = 0.8,      # More granular for detailed analysis
  K_NEIGHBORS = 20,
  
  # DECISION CRITERIA
  MAX_MEMORY_GB = 15,           # Memory limit to decide pipeline
  GRANULARITY_THRESHOLD = 12,   # If Seurat < 12 clusters, use Giotto
  
  # Visualization
  POINT_SIZE = 0.8,
  IMAGE_ALPHA = 0.7, 
  CLUSTER_COLORS = "Set3"
)

cat("Optimized parameters defined\n")
print(optimized_params[c("SEURAT_RESOLUTION", "GIOTTO_RESOLUTION", 
                         "MAX_MEMORY_GB", "GRANULARITY_THRESHOLD")])

# --- FUNCTION FOR SPATIAL PLOT OF THE OPTIMIZED RESULT (Simplified Version) ---
generate_optimized_spatial_plot <- function(final_obj, obj_type, params, title_prefix = "Optimized Pipeline") {
  cat("Generating spatial plot for the result of the", title_prefix, "...\n")
  plot_title <- paste(title_prefix, "-", obj_type)
  
  if (obj_type == "seurat") {
    n_clusters_seurat <- length(unique(Idents(final_obj))) # Assumes final_obj is Seurat and valid
    plot_title <- paste(plot_title, "-", n_clusters_seurat, "clusters")
    
    spatial_plot <- SpatialDimPlot(
      final_obj,
      label = FALSE, 
      pt.size.factor = params$POINT_SIZE * 2, 
      image.alpha = 0, 
      alpha = 1, 
      stroke = 0.1
    ) + 
      ggtitle(plot_title) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right" 
      )
  } else if (obj_type == "giotto") {
    giotto_metadata <- Giotto::pDataDT(final_obj) # Assumes final_obj is Giotto and valid
    spatial_locs_obj <- Giotto::getSpatialLocations(final_obj)
    spatial_coords <- spatial_locs_obj@coordinates
    if (!"cell_ID" %in% colnames(spatial_coords)) {
      spatial_coords$cell_ID <- rownames(spatial_coords)
    }
    plot_data <- merge(giotto_metadata, spatial_coords, by = "cell_ID")
    plot_data$leiden_clus <- as.factor(plot_data$leiden_clus) # Assumes 'leiden_clus' exists
    
    n_clusters_giotto <- length(unique(plot_data$leiden_clus))
    plot_title <- paste(plot_title, "-", n_clusters_giotto, "clusters")

    spatial_plot <- ggplot(plot_data, aes(x = sdimx, y = sdimy, color = leiden_clus)) +
      geom_point(size = params$POINT_SIZE, alpha = 0.8) +
      scale_color_brewer(type = "qual", palette = params$CLUSTER_COLORS, name = "Cluster") +
      coord_fixed() +
      scale_y_reverse() + 
      ggtitle(plot_title) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank()
      ) +
      labs(x = "X Coordinate", y = "Y Coordinate")
  } else { # Shouldn't get here if the selection logic is correct
    cat("Unknown object type for plotting:", obj_type, "\n")
    spatial_plot <- ggplot() + ggtitle(paste(title_prefix,"- Error: Unknown object type")) + theme_minimal()
  }
  return(spatial_plot)
}

# --- 3. MAIN OPTIMIZED PIPELINE ---
cat("\n=== 3. EXECUTING OPTIMIZED PIPELINE ===\n")

run_optimized_pipeline <- function(params, data_path) {
  
  cat("STARTING OPTIMIZED HYBRID PIPELINE\n")
  start_time <- Sys.time()
  initial_memory <- memory_monitor("Start optimized pipeline")
  
  results <- list()
  
  # ========================================================================
  # PHASE 1: INITIAL ANALYSIS WITH SEURAT (FAST AND EFFICIENT)
  # ========================================================================
  cat("\nPHASE 1: INITIAL ANALYSIS WITH SEURAT\n")
  cat("Objective: Quick initial clustering and quality control\n")
  
  seurat_obj <- Load10X_Spatial(
    data.dir = data_path,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = params$sample_name,
    filter.matrix = TRUE
  )
  memory_monitor("Data loaded")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  n_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = 
                         nFeature_Spatial > params$MIN_FEATURES & 
                         nCount_Spatial > params$MIN_COUNTS & 
                         percent.mt < params$MAX_MITO_PERCENT)
  n_after <- ncol(seurat_obj)
  cat("   Spots filtered:", n_before - n_after, "removed,", n_after, "remaining\n")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = params$NORMALIZATION_METHOD,
                              scale.factor = params$SCALE_FACTOR, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                     nfeatures = params$N_HVGS, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj),
                       npcs = params$N_PCS, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:params$N_PCS_USE,
                              k.param = params$K_NEIGHBORS, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, algorithm = 4, 
                             resolution = params$SEURAT_RESOLUTION, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:params$N_PCS_USE, verbose = FALSE)
  seurat_clusters <- length(unique(Idents(seurat_obj)))
  seurat_memory <- memory_monitor("Seurat completed")
  cat("   SEURAT: ", seurat_clusters, "clusters identified\n")
  results$seurat <- list(
    object = seurat_obj,
    n_clusters = seurat_clusters,
    memory_used = seurat_memory - initial_memory
  )
  
  # ========================================================================
  # PHASE 2: EVALUATION AND INTELLIGENT DECISION
  # ========================================================================
  cat("\nPHASE 2: EVALUATION AND DECISION\n")
  current_memory_gb <- seurat_memory / 1024
  use_giotto <- FALSE
  decision_reason <- ""
  if (current_memory_gb > params$MAX_MEMORY_GB) {
    decision_reason <- paste0("Insufficient memory (", round(current_memory_gb, 1), "GB > ", 
                              params$MAX_MEMORY_GB, "GB limit)")
    cat("   DECISION: Seurat only -", decision_reason, "\n")
  } else if (seurat_clusters < params$GRANULARITY_THRESHOLD) {
    use_giotto <- TRUE
    decision_reason <- paste0("Low granularity detected (", seurat_clusters, " < ", 
                              params$GRANULARITY_THRESHOLD, " clusters)")
    cat("   DECISION: Add Giotto -", decision_reason, "\n")
  } else {
    decision_reason <- paste0("Seurat analysis sufficient (", seurat_clusters, " clusters, ",
                              round(current_memory_gb, 1), "GB memory)")
    cat("   DECISION: Seurat only -", decision_reason, "\n")
  }
  results$decision <- list(
    use_giotto = use_giotto,
    reason = decision_reason,
    seurat_clusters = seurat_clusters,
    memory_at_decision = current_memory_gb
  )
  
  # ========================================================================
  # PHASE 3: DETAILED ANALYSIS WITH GIOTTO (CONDITIONAL)
  # ========================================================================
  if (use_giotto) {
    cat("\nPHASE 3: DETAILED ANALYSIS WITH GIOTTO\n")
    cat("Objective: Higher granularity and advanced spatial analysis\n")
    tryCatch({
      h5_file_path <- file.path(data_path, "filtered_feature_bc_matrix.h5")
      spatial_info_dir <- file.path(data_path, "spatial")
      parquet_file_path <- file.path(spatial_info_dir, "tissue_positions.parquet")
      temp_expression_data <- Giotto::get10Xmatrix_h5(path_to_data = h5_file_path)
      if (is.list(temp_expression_data) && "Gene Expression" %in% names(temp_expression_data)) {
        expression_matrix <- temp_expression_data[["Gene Expression"]]
      } else if (inherits(temp_expression_data, "Matrix") || is.matrix(temp_expression_data)) {
        expression_matrix <- temp_expression_data
      } else if (is.list(temp_expression_data) && length(temp_expression_data) == 1) {
        expression_matrix <- temp_expression_data[[1]]
      }
      spatial_locations_df <- as.data.frame(arrow::read_parquet(parquet_file_path))
      seurat_spots_ids <- colnames(seurat_obj)
      common_barcodes <- intersect(colnames(expression_matrix), seurat_spots_ids)
      expression_matrix <- expression_matrix[, common_barcodes]
      spatial_locations_df <- spatial_locations_df[spatial_locations_df$barcode %in% common_barcodes, ]
      spatial_locations_df <- spatial_locations_df[match(common_barcodes, spatial_locations_df$barcode), ]
      rownames(spatial_locations_df) <- spatial_locations_df$barcode
      colnames(spatial_locations_df)[colnames(spatial_locations_df) == "pxl_col_in_fullres"] <- "sdimx"
      colnames(spatial_locations_df)[colnames(spatial_locations_df) == "pxl_row_in_fullres"] <- "sdimy"
      spatial_locations_df <- spatial_locations_df[, c("barcode", "sdimx", "sdimy")]
      giotto_obj <- Giotto::createGiottoObject(
        expression = expression_matrix,
        spatial_locs = spatial_locations_df
      )
      giotto_obj <- Giotto::addStatistics(giotto_obj, expression_values = 'raw',
                                          stats = c("cell"), verbose = FALSE)
      spots_to_keep <- common_barcodes
      giotto_obj <- subsetGiotto(giotto_obj, cell_ids = spots_to_keep)
      giotto_obj <- Giotto::filterGiotto(giotto_obj, expression_values = 'raw',
                                         expression_threshold = 1,
                                         feat_det_in_min_cells = params$MIN_GENES_PER_SPOT,
                                         min_det_feats_per_cell = NULL, verbose = FALSE)
      giotto_obj <- Giotto::normalizeGiotto(giotto_obj, expression_values = 'raw',
                                            norm_methods = 'standard',
                                            scalefactor = params$SCALE_FACTOR,
                                            log_norm = TRUE, scale_feats = FALSE,
                                            scale_cells = FALSE, verbose = FALSE)
      giotto_obj <- Giotto::calculateHVF(giotto_obj, expression_values = 'normalized',
                                         method = 'var_p_resid', nr_expression_groups = 20,
                                         var_number = params$N_HVGS, return_gobject = TRUE,
                                         verbose = FALSE)
      giotto_obj <- Giotto::runPCA(giotto_obj, expression_values = 'normalized',
                                   genes_to_use = 'hvf', scale_unit = TRUE,
                                   center = TRUE, ncp = params$N_PCS, verbose = FALSE)
      giotto_obj <- Giotto::createNearestNetwork(giotto_obj, dimensions_to_use = 1:params$N_PCS_USE,
                                                 k = params$K_NEIGHBORS, name = "snn_network",
                                                 verbose = FALSE)
      giotto_obj <- Giotto::doLeidenCluster(giotto_obj, network_name = "snn_network",
                                            resolution = params$GIOTTO_RESOLUTION,
                                            name = "leiden_clus")
      giotto_obj <- Giotto::runUMAP(giotto_obj, dimensions_to_use = 1:params$N_PCS_USE,
                                    name = "umap", verbose = FALSE)
      giotto_clusters <- length(unique(Giotto::pDataDT(giotto_obj)$leiden_clus))
      giotto_memory <- memory_monitor("Giotto completed")
      cat("   GIOTTO: ", giotto_clusters, "clusters identified\n")
      results$giotto <- list(
        object = giotto_obj,
        n_clusters = giotto_clusters,
        memory_used = giotto_memory - seurat_memory
      )
    }, error = function(e) {
      cat("   ERROR in Giotto:", e$message, "\n")
      cat("   Continuing with Seurat results only\n")
      results$giotto <- NULL
    })
  } else {
    cat("\nPHASE 3: SKIPPED (Seurat sufficient)\n")
    results$giotto <- NULL
  }
  
  # ========================================================================
  # PHASE 4: INTEGRATION AND FINAL RESULTS
  # ========================================================================
  cat("\nPHASE 4: RESULTS INTEGRATION\n")
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "secs")
  final_memory <- memory_monitor("Optimized pipeline completed")
  results$summary <- list(
    total_time_sec = as.numeric(total_time),
    total_time_min = round(as.numeric(total_time)/60, 1),
    initial_memory_mb = initial_memory,
    final_memory_mb = final_memory,
    memory_increase_mb = final_memory - initial_memory,
    seurat_clusters = results$seurat$n_clusters,
    giotto_clusters = if(!is.null(results$giotto) && !is.null(results$giotto$n_clusters)) results$giotto$n_clusters else NA,
    strategy_used = if(!is.null(results$giotto)) "Hybrid (Seurat + Giotto)" else "Seurat Only",
    decision_reason = results$decision$reason,
    efficiency_gain = if(!is.null(results$giotto) && !is.null(results$giotto$n_clusters)) {
      paste0("Seurat: ", results$seurat$n_clusters, " clusters in ", 
             round(seurat_memory/1024, 1), "GB - Giotto: +", 
             results$giotto$n_clusters - results$seurat$n_clusters, " clusters")
    } else {
      paste0("Seurat Only: ", results$seurat$n_clusters, " clusters in ", 
             round(seurat_memory/1024, 1), "GB (sufficient)")
    }
  )
  return(results)
}

# --- 4. EXECUTION OF THE PIPELINE ---
cat("\n=== 4. EXECUTING OPTIMIZED PIPELINE ===\n")
optimized_results <- run_optimized_pipeline(optimized_params, data_dir)

# --- 5. ANALYSIS OF RESULTS ---
cat("\n=== 5. ANALYSIS OF OPTIMIZED PIPELINE RESULTS ===\n")
analyze_optimized_results <- function(results) {
  summary_data_local <- results$summary
  cat("OPTIMIZED PIPELINE SUMMARY:\n")
  cat("Strategy:", summary_data_local$strategy_used, "\n")
  cat("Decision:", summary_data_local$decision_reason, "\n")
  cat("Total time:", summary_data_local$total_time_min, "minutes\n")
  cat("Total memory:", format_memory(summary_data_local$memory_increase_mb), "\n")
  cat("Efficiency:", summary_data_local$efficiency_gain, "\n")
  cat("\nDETAILED ANALYSIS:\n")
  
  comparison_type_local <- if(summary_data_local$strategy_used == "Hybrid (Seurat + Giotto)") "hybrid" else "seurat_only"

  if (comparison_type_local == "hybrid") {
    cat("HYBRID MODE EXECUTED:\n")
    cat("- Seurat (initial):", results$seurat$n_clusters, "clusters\n")
    cat("- Giotto (detailed):", results$giotto$n_clusters, "clusters\n")
    cat("- Granularity gain:", results$giotto$n_clusters - results$seurat$n_clusters, "additional clusters\n")
  } else {
    cat("SEURAT ONLY MODE:\n")
    cat("- Clusters identified:", results$seurat$n_clusters, "\n")
    cat("- Reason:", results$decision$reason, "\n")
  }
  cat("\nADVANTAGES OF THE OPTIMIZED PIPELINE:\n")
  cat("Intelligent decision based on resources\n")
  cat("Maximum computational efficiency\n")
  cat("Adaptive granularity as needed\n")
  cat("Robust fallback for errors\n")
  return(list(
    type = comparison_type_local,
    summary = summary_data_local, 
    recommendation = if(comparison_type_local == "hybrid") {
      "Successful hybrid pipeline - combines efficiency and granularity"
    } else {
      "Efficient pipeline - Seurat sufficient for this dataset"
    }
  ))
}
analysis <- analyze_optimized_results(optimized_results)

# --- 6. SAVE RESULTS ---
cat("\n=== 6. SAVING OPTIMIZED PIPELINE RESULTS ===\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
optimized_complete <- list(
  parameters = optimized_params,
  results = optimized_results,
  analysis = analysis,
  execution_info = list(
    timestamp = timestamp,
    date = Sys.Date(),
    R_version = as.character(getRversion()),
    Seurat_version = as.character(packageVersion("Seurat")),
    Giotto_version = as.character(packageVersion("Giotto"))
  )
)
output_file_rds <- file.path(main_results_dir, paste0("optimized_pipeline_results_", timestamp, ".rds"))
saveRDS(optimized_complete, file = output_file_rds)
cat("Results saved in:", output_file_rds, "\n")

# --- 7. FINAL SUMMARY (CONSOLE, PLOT AND FILE) ---
cat("\n=== 7. FINAL SUMMARY OF THE OPTIMIZED PIPELINE ===\n")

# Determine the final object and type for the plot, based on the STRATEGY USED
final_plot_object <- NULL
final_object_type <- ""
plot_file_path <- "Plot not generated (final object not determined)."

if (optimized_results$summary$strategy_used == "Hybrid (Seurat + Giotto)") {
    final_plot_object <- optimized_results$results$giotto$object 
    final_object_type <- "giotto"
    cat("The Giotto object will be used for the final plot (hybrid strategy).\n")
} else { # "Seurat Only"
    final_plot_object <- optimized_results$results$seurat$object
    final_object_type <- "seurat"
    cat("The Seurat object will be used for the final plot (Seurat only strategy).\n")
}

# Generate the spatial plot of the final result
optimizado_spatial_plot <- NULL
if (!is.null(final_plot_object)) {
    optimizado_spatial_plot <- generate_optimized_spatial_plot(
      final_obj = final_plot_object, 
      obj_type = final_object_type, 
      params = optimized_params
    )
}

# Save the plot
if (!is.null(optimizado_spatial_plot) && inherits(optimizado_spatial_plot, "ggplot")) {
    plot_file_name <- paste0("plot_optimized_pipeline_", timestamp, ".png")
    plot_file_path <- file.path(optimized_output_subdir, plot_file_name)
    ggsave(plot_file_path, optimizado_spatial_plot, width = 10, height = 8, dpi = 300)
    cat("Optimized pipeline plot saved in:", plot_file_path, "\n")
} else {
    cat("Could not generate or save the optimized pipeline spatial plot.\n")
    if(is.null(final_plot_object)) {
        plot_file_path <- "Plot not generated (final object is null)."
    } else if (is.null(optimizado_spatial_plot) || !inherits(optimizado_spatial_plot, "ggplot")) {
         plot_file_path <- "Plot not generated (error in plotting function)."
    }
}

# Prepare report content
summary_data_final <- optimized_results$summary 
report_content <- paste0(
  "=== FINAL SUMMARY OF THE OPTIMIZED PIPELINE ===\n\n",
  "STRATEGY USED: ", summary_data_final$strategy_used, "\n",
  "TOTAL TIME: ", summary_data_final$total_time_min, " minutes\n",
  "MEMORY USED (increase): ", format_memory(summary_data_final$memory_increase_mb), "\n",
  "FINAL CLUSTERS (Seurat initial): ", summary_data_final$seurat_clusters, "\n",
  "FINAL CLUSTERS (Giotto detailed, if applicable): ", if(!is.na(summary_data_final$giotto_clusters)) summary_data_final$giotto_clusters else "N/A", "\n",
  "FINAL CLUSTERS (Reported by the pipeline): ", if(summary_data_final$strategy_used == "Hybrid (Seurat + Giotto)" && !is.na(summary_data_final$giotto_clusters)) summary_data_final$giotto_clusters else summary_data_final$seurat_clusters, "\n",
  "DECISION: ", summary_data_final$decision_reason, "\n",
  "RECOMMENDATION: ", analysis$recommendation, "\n\n",
  "OPTIMIZED PIPELINE COMPLETED SUCCESSFULLY\n",
  "RDS results in: ", output_file_rds, "\n",
  "Spatial plot in: ", plot_file_path, "\n"
)

# Print to console
cat(report_content)

# Save report to text file
report_file_txt_path <- file.path(optimized_output_subdir, paste0("report_optimized_pipeline_", timestamp, ".txt"))
writeLines(report_content, con = report_file_txt_path)
cat("Report saved in:", report_file_txt_path, "\n")

cat("\n=== END OF THE OPTIMIZED PIPELINE ===\n")