#!/usr/bin/env Rscript
# =============================================================================
# GIOTTO ANALYSIS WITH MEMORY MONITORING - STANDALONE SCRIPT
# Author: Ángel Pérez
# Date: 31 May 2025
# Objective: Homologous Giotto pipeline with precise resource measurement
# =============================================================================

# --- 0. INITIAL SETUP ---
cat("=== STARTING STANDALONE GIOTTO ANALYSIS ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Process PID:", Sys.getpid(), "\n\n")

# Clean environment completely at the beginning
rm(list = ls())
gc()

# Load libraries
suppressPackageStartupMessages({
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
set.seed(12345)  # Reproducibility
plan("sequential")  # No parallelization for fair comparison

# --- CONFIGURATION & PATHS ---
# This script assumes it is executed from the project's root directory.
data_dir <- "data/visium_dataset"
main_results_dir <- "results"
giotto_plots_dir <- file.path(main_results_dir, "02_Individual_Giotto_Results")
dir.create(giotto_plots_dir, showWarnings = FALSE, recursive = TRUE)
# --- END CONFIGURATION ---

# --- 1. MEMORY MONITORING FUNCTIONS ---
cat("1. Configuring memory monitoring system...\n")

# Function to get current memory in MB
get_memory_usage <- function() {
  if (Sys.info()["sysname"] == "Linux") {
    # On Linux, use /proc/self/status
    status_file <- "/proc/self/status"
    if (file.exists(status_file)) {
      status_lines <- readLines(status_file)
      vmrss_line <- grep("^VmRSS:", status_lines, value = TRUE)
      if (length(vmrss_line) > 0) {
        # Extract value in kB and convert to MB
        mem_kb <- as.numeric(gsub(".*VmRSS:\\s*(\\d+)\\s*kB.*", "\\1", vmrss_line))
        return(mem_kb / 1024)  # Convert to MB
      }
    }
  }
  
  # Fallback using pryr (less precise but universal)
  tryCatch({
    mem_info <- pryr::mem_used()
    return(as.numeric(mem_info) / 1024^2)  # Convert to MB
  }, error = function(e) {
    return(NA)
  })
}

# Function to format memory
format_memory <- function(mb) {
  if (is.na(mb)) return("N/A")
  if (mb < 1024) {
    return(paste0(round(mb, 1), " MB"))
  } else {
    return(paste0(round(mb/1024, 2), " GB"))
  }
}

# Memory monitor during execution
memory_monitor <- function(label = "Checkpoint") {
  current_mem <- get_memory_usage()
  cat("   [MEMORY]", label, ":", format_memory(current_mem), "\n")
  return(current_mem)
}

# Initial memory
initial_memory <- memory_monitor("Clean script start")

# --- 2. DEFINE PARAMETERS ---
cat("2. Defining parameters for Giotto...\n")

params <- list(
  # sample_name is used for metadata but data_dir is the key path
  sample_name = "Visium_HD_Human_Lung_Cancer_16um",
  
  # QC Parameters
  MIN_FEATURES = 500,
  MIN_COUNTS = 1000,
  MAX_MITO_PERCENT = 20,
  MIN_GENES_PER_SPOT = 3,  # For gene filtering in Giotto
  
  # Normalization
  NORMALIZATION_METHOD = "LogNormalize", # Although Giotto uses its 'standard', we maintain nominal consistency
  SCALE_FACTOR = 10000,
  
  # Highly variable genes
  N_HVGS = 3000,
  
  # PCA
  N_PCS = 50,  # Calculate 50
  N_PCS_USE = 20,  # Use 20 for clustering
  
  # Clustering
  CLUSTERING_METHOD = "leiden",
  RESOLUTION = 0.5,
  K_NEIGHBORS = 20,
  
  # Spatial visualization
  POINT_SIZE = 0.8,
  IMAGE_ALPHA = 0.7, # Does not apply directly to Giotto::spatPlot without image, but kept for consistency
  CLUSTER_COLORS = "Set3"
)

# Print parameters
cat("\nDefined parameters:\n")
print(params[c("MIN_FEATURES", "MIN_COUNTS", "MAX_MITO_PERCENT", 
               "N_HVGS", "N_PCS_USE", "CLUSTERING_METHOD", "RESOLUTION")])

# --- 3. GIOTTO PIPELINE WITH MONITORING ---
cat("\n=== 3. EXECUTING GIOTTO PIPELINE ===\n")

run_giotto_analysis <- function(params, data_path) {
  
  cat("3.1 Loading data with Giotto...\n")
  start_time <- Sys.time()
  mem_start <- memory_monitor("Start Giotto pipeline")
  
  # Define paths
  h5_file_path <- file.path(data_path, "filtered_feature_bc_matrix.h5")
  spatial_info_dir <- file.path(data_path, "spatial")
  parquet_file_path <- file.path(spatial_info_dir, "tissue_positions.parquet")
  
  # 3.1 Load expression matrix
  temp_expression_data <- Giotto::get10Xmatrix_h5(path_to_data = h5_file_path)
  
  # Handle output as a list
  if (is.list(temp_expression_data) && "Gene Expression" %in% names(temp_expression_data)) {
    expression_matrix <- temp_expression_data[["Gene Expression"]]
  } else if (inherits(temp_expression_data, "Matrix") || is.matrix(temp_expression_data)) {
    expression_matrix <- temp_expression_data
  } else if (is.list(temp_expression_data) && length(temp_expression_data) == 1) {
    expression_matrix <- temp_expression_data[[1]]
  } else {
    stop("Unrecognized matrix format")
  }
  
  mem_load <- memory_monitor("Matrix loaded")
  cat("   Matrix loaded:", ncol(expression_matrix), "spots,", nrow(expression_matrix), "genes\n")
  
  # 3.2 Load spatial positions
  spatial_locations_df <- as.data.frame(arrow::read_parquet(parquet_file_path))
  
  # 3.3 Synchronize barcodes
  cat("3.2 Synchronizing barcodes...\n")
  barcodes_expr <- colnames(expression_matrix)
  barcodes_spatial <- spatial_locations_df$barcode
  common_barcodes <- intersect(barcodes_expr, barcodes_spatial)
  
  # Filter both
  expression_matrix <- expression_matrix[, common_barcodes]
  spatial_locations_df <- spatial_locations_df[spatial_locations_df$barcode %in% common_barcodes, ]
  spatial_locations_df <- spatial_locations_df[match(colnames(expression_matrix), spatial_locations_df$barcode), ]
  
  # Prepare positions
  rownames(spatial_locations_df) <- spatial_locations_df$barcode
  colnames(spatial_locations_df)[colnames(spatial_locations_df) == "pxl_col_in_fullres"] <- "sdimx"
  colnames(spatial_locations_df)[colnames(spatial_locations_df) == "pxl_row_in_fullres"] <- "sdimy"
  
  # Keep only necessary columns
  spatial_locations_df <- spatial_locations_df[, c("barcode", "sdimx", "sdimy")]
  
  # 3.4 Create Giotto object
  cat("3.3 Creating Giotto object...\n")
  giotto_obj <- Giotto::createGiottoObject(
    expression = expression_matrix,
    spatial_locs = spatial_locations_df
  )
  mem_object <- memory_monitor("Giotto object created")
  
  # 3.5 Calculate statistics
  cat("3.4 Calculating QC statistics...\n")
  giotto_obj <- Giotto::addStatistics(
    gobject = giotto_obj,
    expression_values = 'raw',
    stats = c("cell"),
    verbose = FALSE
  )
  
  # 3.6 Calculate mitochondrial % manually
  raw_matrix <- Giotto::getExpression(giotto_obj, values = 'raw', output = 'matrix')
  mito_genes <- grep(pattern = "^MT-", x = rownames(raw_matrix), value = TRUE)
  total_counts <- Matrix::colSums(raw_matrix)
  
  if(length(mito_genes) > 0) {
    mito_counts <- Matrix::colSums(raw_matrix[mito_genes, , drop = FALSE])
    mito_percentage <- (mito_counts / total_counts) * 100
    mito_percentage[is.na(mito_percentage)] <- 0
    mito_percentage[is.infinite(mito_percentage)] <- 0
  } else {
    mito_percentage <- rep(0, ncol(raw_matrix))
    names(mito_percentage) <- colnames(raw_matrix)
  }
  
  # Add to metadata
  mito_dt <- data.table::data.table(
    cell_ID = names(mito_percentage), 
    percent_mt = mito_percentage
  )
  
  giotto_obj <- Giotto::addCellMetadata(
    gobject = giotto_obj,
    new_metadata = mito_dt,
    by_column = TRUE,
    column_cell_ID = 'cell_ID'
  )
  
  # 3.7 Filter spots
  cat("3.5 Filtering spots by QC...\n")
  metadata <- Giotto::pDataDT(giotto_obj)
  n_before <- nrow(metadata)
  
  # Apply homologous filters
  spots_to_keep <- metadata[
    nr_feats > params$MIN_FEATURES & 
      total_expr > params$MIN_COUNTS & 
      percent_mt < params$MAX_MITO_PERCENT
  ]$cell_ID
  
  giotto_obj <- subsetGiotto(
    gobject = giotto_obj,
    cell_ids = spots_to_keep
  )
  
  n_after <- length(spots_to_keep)
  mem_qc <- memory_monitor("QC applied")
  cat("   Spots filtered:", n_before - n_after, "removed,", n_after, "remaining\n")
  
  # 3.8 Filter genes
  cat("3.6 Filtering genes...\n")
  giotto_obj <- Giotto::filterGiotto(
    gobject = giotto_obj,
    expression_values = 'raw',
    expression_threshold = 1,
    feat_det_in_min_cells = params$MIN_GENES_PER_SPOT,
    min_det_feats_per_cell = NULL,
    verbose = FALSE
  )
  
  # 3.9 Normalization
  cat("3.7 Normalizing with standard method (LogNormalize)...\n")
  giotto_obj <- Giotto::normalizeGiotto(
    gobject = giotto_obj,
    expression_values = 'raw',
    norm_methods = 'standard',
    scalefactor = params$SCALE_FACTOR,
    log_norm = TRUE,
    scale_feats = FALSE,  # Correct for v4.2.1
    scale_cells = FALSE,
    verbose = FALSE
  )
  mem_norm <- memory_monitor("Normalization")
  
  # 3.10 Variable genes - CORRECT SYNTAX FOR v4.2.1
  cat("3.8 Identifying", params$N_HVGS, "variable genes...\n")
  giotto_obj <- Giotto::calculateHVF(
    gobject = giotto_obj,
    expression_values = 'normalized',
    method = 'var_p_resid',
    nr_expression_groups = 20,
    var_number = params$N_HVGS,  # Correct argument for v4.2.1
    return_gobject = TRUE,
    verbose = FALSE
  )
  mem_hvg <- memory_monitor("Variable genes")
  
  # 3.11 PCA
  cat("3.9 Running PCA...\n")
  giotto_obj <- Giotto::runPCA(
    gobject = giotto_obj,
    expression_values = 'normalized',
    genes_to_use = 'hvf',
    scale_unit = TRUE,
    center = TRUE,
    ncp = params$N_PCS,
    verbose = FALSE
  )
  mem_pca <- memory_monitor("PCA")
  
  # 3.12 kNN network
  cat("3.10 Creating nearest neighbor network (k =", params$K_NEIGHBORS, ")...\n")
  giotto_obj <- Giotto::createNearestNetwork(
    gobject = giotto_obj,
    dimensions_to_use = 1:params$N_PCS_USE,
    k = params$K_NEIGHBORS,
    name = "snn_network",
    verbose = FALSE
  )
  
  # 3.13 Clustering with Leiden (without verbose)
  cat("3.11 Clustering with Leiden (resolution", params$RESOLUTION, ")...\n")
  giotto_obj <- Giotto::doLeidenCluster(
    gobject = giotto_obj,
    network_name = "snn_network",
    resolution = params$RESOLUTION,
    name = "leiden_clus"
  )
  mem_cluster <- memory_monitor("Clustering")
  
  # 3.14 UMAP
  cat("3.12 Calculating UMAP...\n")
  giotto_obj <- Giotto::runUMAP(
    gobject = giotto_obj,
    dimensions_to_use = 1:params$N_PCS_USE,
    name = "umap",
    verbose = FALSE
  )
  mem_umap <- memory_monitor("UMAP")
  
  # Summary
  metadata_final <- Giotto::pDataDT(giotto_obj)
  n_clusters <- length(unique(metadata_final$leiden_clus))
  cat("\n   GIOTTO COMPLETE:", n_clusters, "clusters identified\n")
  
  end_time <- Sys.time()
  time_taken <- difftime(end_time, start_time, units = "secs")
  mem_final <- memory_monitor("Final Giotto")
  
  # Collect memory measurements
  memory_profile <- list(
    start = mem_start,
    load = mem_load,
    object = mem_object,
    qc = mem_qc,
    normalization = mem_norm,
    hvg = mem_hvg,
    pca = mem_pca,
    clustering = mem_cluster,
    umap = mem_umap,
    final = mem_final,
    peak = max(c(mem_start, mem_load, mem_object, mem_qc, mem_norm, mem_hvg, 
                 mem_pca, mem_cluster, mem_umap, mem_final), na.rm = TRUE),
    increase = mem_final - mem_start
  )
  
  return(list(
    object = giotto_obj,
    n_clusters = n_clusters,
    time_sec = as.numeric(time_taken),
    memory_profile = memory_profile
  ))
}

# Run Giotto analysis
giotto_results <- run_giotto_analysis(params, data_dir)

# --- 4. GENERATE SPATIAL VISUALIZATION ---
cat("\n=== 4. GENERATING SPATIAL VISUALIZATION ===\n")

generate_giotto_spatial_plot <- function(giotto_obj, params) {
  
  cat("4.1 Generating Giotto spatial plot...\n")
  
  # Get spatial data from Giotto - CORRECT METHOD
  giotto_metadata <- Giotto::pDataDT(giotto_obj)
  
  # Extract spatial coordinates from the Giotto object
  spatial_locs_obj <- Giotto::getSpatialLocations(giotto_obj)
  
  # Convert the spatial object to data.table
  spatial_coords <- spatial_locs_obj@coordinates
  
  # Ensure coordinates have the correct ID
  if (!"cell_ID" %in% colnames(spatial_coords)) {
    spatial_coords$cell_ID <- rownames(spatial_coords)
  }
  
  # Combine metadata with spatial locations
  plot_data <- merge(giotto_metadata, spatial_coords, by = "cell_ID")
  
  # CONVERT leiden_clus to factor for discrete scale
  plot_data$leiden_clus <- as.factor(plot_data$leiden_clus)
  
  spatial_plot <- ggplot(plot_data, aes(x = sdimx, y = sdimy, color = leiden_clus)) +
    geom_point(size = params$POINT_SIZE, alpha = 0.8) +
    scale_color_brewer(type = "qual", palette = params$CLUSTER_COLORS, name = "Cluster") +
    coord_fixed() +
    scale_y_reverse() +  # Invert Y to match typical image orientation
    ggtitle("Giotto - Spatial Clusters") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid = element_blank()
    ) +
    labs(x = "X Coordinate", y = "Y Coordinate")
  
  return(spatial_plot)
}

# Generate spatial plot
giotto_spatial_plot <- generate_giotto_spatial_plot(giotto_results$object, params)

# --- 5. SAVE RESULTS ---
cat("\n=== 5. SAVING RESULTS ===\n")

# Create timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 5.1 Create object with complete results
giotto_complete_results <- list(
  # Parameters used
  parameters = params,
  
  # Pipeline results
  results = giotto_results,
  
  # System information
  system_info = list(
    pipeline = "Giotto",
    date = Sys.Date(),
    time = Sys.time(),
    timestamp = timestamp,
    pid = Sys.getpid(),
    R_version = as.character(getRversion()),
    Giotto_version = as.character(packageVersion("Giotto")),
    initial_memory_mb = initial_memory,
    session_info = sessionInfo()
  )
)

# 5.2 Save RDS results (in the main results directory)
output_file_rds <- file.path(main_results_dir, paste0("giotto_results_", timestamp, ".rds"))
saveRDS(giotto_complete_results, file = output_file_rds)
cat("5.1 Giotto results (RDS) saved in:", output_file_rds, "\n")

# 5.3 Save spatial plot (in the designated subfolder)
output_file_plot <- file.path(giotto_plots_dir, paste0("giotto_spatial_", timestamp, ".png"))
ggsave(output_file_plot, giotto_spatial_plot, 
       width = 10, height = 8, dpi = 300)
cat("5.2 Spatial plot saved in:", output_file_plot, "\n")

# 5.4 Generate quick summary
cat("\n=== GIOTTO SUMMARY ===\n")
cat("Clusters found:", giotto_results$n_clusters, "\n")
cat("Execution time:", round(giotto_results$time_sec/60, 1), "minutes\n")
cat("Peak memory:", format_memory(giotto_results$memory_profile$peak), "\n")
cat("Memory increase:", format_memory(giotto_results$memory_profile$increase), "\n")
cat("Final spots:", length(Giotto::pDataDT(giotto_results$object)$cell_ID), "\n")

# 5.5 Information for next script
cat("\nRDS file for comparison:", output_file_rds, "\n")
cat("Giotto pipeline completed successfully\n")
cat("Execute now: 3_Compare_pipelines.R\n")
cat("\n=== END OF GIOTTO ANALYSIS ===\n")