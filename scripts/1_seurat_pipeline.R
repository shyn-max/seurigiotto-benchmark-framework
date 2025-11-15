#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 1: SEURAT ANALYSIS WITH MEMORY MONITORING - STANDALONE SCRIPT
# Author: Ángel Pérez
# Date: 31 May 2025
# Objective: Homologous Seurat pipeline with precise resource measurement
# =============================================================================

# --- 0. INITIAL SETUP ---
cat("=== SCRIPT 1: STARTING STANDALONE SEURAT ANALYSIS ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Process PID:", Sys.getpid(), "\n\n")

# Clean environment completely at the beginning
rm(list = ls())
gc()

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(future)
  library(pryr)
  library(scales)
  library(RColorBrewer)
})

# --- CONFIGURATION & PATHS ---
# This script assumes it is executed from the project's root directory.
# All paths are relative to that root.

# Input Data Directory
# IMPORTANT: Place your Visium data inside the folder specified here.
data_dir <- "data/visium_dataset"

# Output Directories
main_results_dir <- "results"
seurat_plots_dir <- file.path(main_results_dir, "01_Individual_Seurat_Results")
dir.create(seurat_plots_dir, showWarnings = FALSE, recursive = TRUE)
# --- END CONFIGURATION ---

# General settings
set.seed(12345)  # Reproducibility
plan("sequential")  # No parallelization for fair comparison


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
cat("2. Defining parameters for Seurat...\n")

params <- list(
  # sample_name is used for the Seurat object slice name
  sample_name = "Visium_HD_Human_Lung_Cancer_16um",
  
  # QC Parameters
  MIN_FEATURES = 500,
  MIN_COUNTS = 1000,
  MAX_MITO_PERCENT = 20,
  
  # Normalization
  NORMALIZATION_METHOD = "LogNormalize",
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
  IMAGE_ALPHA = 0.7,
  CLUSTER_COLORS = "Set3"
)

# Print parameters
cat("\nDefined parameters:\n")
print(params[c("MIN_FEATURES", "MIN_COUNTS", "MAX_MITO_PERCENT", 
               "N_HVGS", "N_PCS_USE", "CLUSTERING_METHOD", "RESOLUTION")])

# --- 3. SEURAT PIPELINE WITH MONITORING ---
cat("\n=== 3. EXECUTING SEURAT PIPELINE ===\n")

run_seurat_analysis <- function(params, data_path) {
  
  cat("3.1 Loading data with Seurat...\n")
  start_time <- Sys.time()
  mem_start <- memory_monitor("Start Seurat pipeline")
  
  # Load data
  seurat_obj <- Load10X_Spatial(
    data.dir = data_path,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = params$sample_name,
    filter.matrix = TRUE
  )
  
  mem_load <- memory_monitor("Data loaded")
  cat("   Data loaded:", ncol(seurat_obj), "spots,", nrow(seurat_obj), "genes\n")
  
  # 3.2 Quality control
  cat("3.2 Applying quality control...\n")
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter
  n_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = 
                         nFeature_Spatial > params$MIN_FEATURES & 
                         nCount_Spatial > params$MIN_COUNTS & 
                         percent.mt < params$MAX_MITO_PERCENT
  )
  n_after <- ncol(seurat_obj)
  mem_qc <- memory_monitor("QC applied")
  
  cat("   Spots filtered:", n_before - n_after, "removed,", n_after, "remaining\n")
  
  # 3.3 Normalization
  cat("3.3 Normalizing with LogNormalize...\n")
  seurat_obj <- NormalizeData(
    seurat_obj, 
    normalization.method = params$NORMALIZATION_METHOD,
    scale.factor = params$SCALE_FACTOR,
    verbose = FALSE
  )
  mem_norm <- memory_monitor("Normalization")
  
  # 3.4 Variable genes
  cat("3.4 Identifying", params$N_HVGS, "variable genes...\n")
  seurat_obj <- FindVariableFeatures(
    seurat_obj, 
    selection.method = "vst",
    nfeatures = params$N_HVGS,
    verbose = FALSE
  )
  mem_hvg <- memory_monitor("Variable genes")
  
  # 3.5 Scaling
  cat("3.5 Scaling data...\n")
  seurat_obj <- ScaleData(
    seurat_obj, 
    features = VariableFeatures(seurat_obj),
    verbose = FALSE
  )
  mem_scale <- memory_monitor("Scaling")
  
  # 3.6 PCA
  cat("3.6 Running PCA...\n")
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = params$N_PCS,
    verbose = FALSE
  )
  mem_pca <- memory_monitor("PCA")
  
  # 3.7 Clustering with LEIDEN
  cat("3.7 Clustering with Leiden (resolution", params$RESOLUTION, ")...\n")
  seurat_obj <- FindNeighbors(
    seurat_obj, 
    dims = 1:params$N_PCS_USE,
    k.param = params$K_NEIGHBORS,
    verbose = FALSE
  )
  
  seurat_obj <- FindClusters(
    seurat_obj,
    algorithm = 4,  # 4 = Leiden
    resolution = params$RESOLUTION,
    verbose = FALSE
  )
  mem_cluster <- memory_monitor("Clustering")
  
  # 3.8 UMAP
  cat("3.8 Calculating UMAP...\n")
  seurat_obj <- RunUMAP(
    seurat_obj,
    dims = 1:params$N_PCS_USE,
    verbose = FALSE
  )
  mem_umap <- memory_monitor("UMAP")
  
  # Summary
  n_clusters <- length(unique(Idents(seurat_obj)))
  cat("\n   SEURAT COMPLETE:", n_clusters, "clusters identified\n")
  
  end_time <- Sys.time()
  time_taken <- difftime(end_time, start_time, units = "secs")
  mem_final <- memory_monitor("Final Seurat")
  
  # Collect memory measurements
  memory_profile <- list(
    start = mem_start,
    load = mem_load,
    qc = mem_qc,
    normalization = mem_norm,
    hvg = mem_hvg,
    scale = mem_scale,
    pca = mem_pca,
    clustering = mem_cluster,
    umap = mem_umap,
    final = mem_final,
    peak = max(c(mem_start, mem_load, mem_qc, mem_norm, mem_hvg, 
                 mem_scale, mem_pca, mem_cluster, mem_umap, mem_final), na.rm = TRUE),
    increase = mem_final - mem_start
  )
  
  return(list(
    object = seurat_obj,
    n_clusters = n_clusters,
    time_sec = as.numeric(time_taken),
    memory_profile = memory_profile
  ))
}

# Run Seurat analysis
seurat_results <- run_seurat_analysis(params, data_dir)

# --- 4. GENERATE SPATIAL VISUALIZATION ---
cat("\n=== 4. GENERATING SPATIAL VISUALIZATION ===\n")

generate_seurat_spatial_plot <- function(seurat_obj, params) {
  
  cat("4.1 Generating Seurat spatial plot...\n")
  
  # Check if an image is available
  has_image <- length(seurat_obj@images) > 0
  
  if (has_image) {
    # With tissue image
    spatial_plot <- SpatialDimPlot(
      seurat_obj,
      group.by = "seurat_clusters",
      pt.size.factor = params$POINT_SIZE,
      alpha = c(params$IMAGE_ALPHA, 1),
      stroke = 0.2
    ) +
      ggtitle("Seurat - Spatial Clusters") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
  } else {
    # Without image, use spatial coordinates
    spatial_coords <- GetTissueCoordinates(seurat_obj)
    cluster_data <- data.frame(
      x = spatial_coords[, 1],
      y = spatial_coords[, 2],
      cluster = Idents(seurat_obj)
    )
    
    spatial_plot <- ggplot(cluster_data, aes(x = x, y = y, color = cluster)) +
      geom_point(size = params$POINT_SIZE, alpha = 0.8) +
      scale_color_brewer(type = "qual", palette = params$CLUSTER_COLORS, name = "Cluster") +
      coord_fixed() +
      ggtitle("Seurat - Spatial Clusters") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank()
      ) +
      labs(x = "X Coordinate", y = "Y Coordinate")
  }
  
  return(spatial_plot)
}

# Generate spatial plot
seurat_spatial_plot <- generate_seurat_spatial_plot(seurat_results$object, params)

# --- 5. SAVE RESULTS ---
cat("\n=== 5. SAVING RESULTS ===\n")

# Create timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 5.1 Create object with complete results
seurat_complete_results <- list(
  # Parameters used
  parameters = params,
  
  # Pipeline results
  results = seurat_results,
  
  # System information
  system_info = list(
    pipeline = "Seurat",
    date = Sys.Date(),
    time = Sys.time(),
    timestamp = timestamp,
    pid = Sys.getpid(),
    R_version = as.character(getRversion()),
    Seurat_version = as.character(packageVersion("Seurat")),
    initial_memory_mb = initial_memory,
    session_info = sessionInfo()
  )
)

# 5.2 Save RDS results (in the main results directory)
output_file_rds <- file.path(main_results_dir, paste0("seurat_results_", timestamp, ".rds"))
saveRDS(seurat_complete_results, file = output_file_rds)
cat("5.1 Seurat results (RDS) saved in:", output_file_rds, "\n")

# 5.3 Save spatial plot (in the designated subfolder)
output_file_plot <- file.path(seurat_plots_dir, paste0("seurat_spatial_", timestamp, ".png"))
ggsave(output_file_plot, seurat_spatial_plot, 
       width = 10, height = 8, dpi = 300)
cat("5.2 Spatial plot saved in:", output_file_plot, "\n")

# 5.4 Generate quick summary
cat("\n=== SEURAT SUMMARY ===\n")
cat("Clusters found:", seurat_results$n_clusters, "\n")
cat("Execution time:", round(seurat_results$time_sec/60, 1), "minutes\n")
cat("Peak memory:", format_memory(seurat_results$memory_profile$peak), "\n")
cat("Memory increase:", format_memory(seurat_results$memory_profile$increase), "\n")
cat("Final spots:", ncol(seurat_results$object), "\n")

# 5.5 Information for next script
cat("\nRDS file for comparison:", output_file_rds, "\n")
cat("Seurat pipeline completed successfully\n")
cat("Execute now: 2_Giotto_pipeline.R\n")
cat("\n=== END OF SEURAT ANALYSIS ===\n")