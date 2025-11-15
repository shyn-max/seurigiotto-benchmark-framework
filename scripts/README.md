# Scripts Workflow Documentation

This document provides a detailed breakdown of each R script in the `seurigiotto-benchmark` project. The scripts are designed to be run sequentially, especially for the full integrated analysis.

---

### Workflow Overview

The project is structured as a series of modular scripts:
1.  **Benchmarking (Scripts 1-3):** Independent execution of Seurat and Giotto with harmonized parameters, followed by a direct performance comparison.
2.  **Advanced Pipelines (Script 4):** A novel hybrid pipeline that makes an intelligent choice between Seurat and Giotto based on real-time analysis criteria.
3.  **Integrated Analysis (Scripts 5-9):** A suite of scripts that load the results from the initial pipelines to perform deep comparative analysis, biological interpretation, and final report generation.

---

### Script Details

#### `01_seurat_pipeline.R`
-   **Purpose:** Runs a complete, standalone Spatial Transcriptomics analysis using the **Seurat** package.
-   **Key Features:** Includes detailed memory and execution time monitoring for benchmarking purposes.
-   **Inputs:** Raw spatial data from `data/visium_dataset`.
-   **Outputs:** An `.rds` file (`seurat_results_[timestamp].rds`) in the `results/` directory containing the Seurat object and performance metrics. A spatial plot is saved in `results/01_Individual_Seurat_Results/`.

#### `02_giotto_pipeline.R`
-   **Purpose:** Runs a complete, standalone analysis using the **Giotto** package, with parameters homologous to the Seurat script.
-   **Key Features:** Includes identical memory and time monitoring for a fair benchmark.
-   **Inputs:** Raw spatial data from `data/visium_dataset`.
-   **Outputs:** An `.rds` file (`giotto_results_[timestamp].rds`) and a spatial plot in their respective output directories.

#### `03_compare_pipelines.R`
-   **Purpose:** Loads the `.rds` outputs from scripts 1 and 2 to perform a direct **performance comparison**.
-   **Key Features:** Generates comparative plots for cluster count, execution time, and memory usage. Produces a final automated text report summarizing the benchmark.
-   **Inputs:** The latest `seurat_results_...rds` and `giotto_results_...rds` files from the `results/` directory.
-   **Outputs:** Comparative plots and a final report in `results/03_Seurat_Giotto_Comparison/`.

#### `4_optimized_pipeline.R`
-   **Purpose:** Executes the novel **Optimized Hybrid Pipeline**.
-   **Key Features:** First runs a fast Seurat analysis, then intelligently decides whether to proceed with a more granular Giotto analysis based on predefined memory and cluster count thresholds.
-   **Inputs:** Raw spatial data from `data/visium_dataset`.
-   **Outputs:** A detailed `.rds` object (`optimized_pipeline_results_...rds`) containing the results of the chosen path (Seurat-only or Seurat+Giotto) and the decision logic.

#### `5_integrated_analysis.R`
-   **Purpose:** Acts as the **main integration script**. It loads the results from the Seurat, Giotto, and Optimized pipelines to perform a comprehensive scientific comparison.
-   **Key Features:** Runs marker gene analysis and cell signature scoring across all three results.
-   **Inputs:** The latest `.rds` files from the Seurat, Giotto, and Optimized pipelines.
-   **Outputs:** A master `.rds` object (`final_results_object.rds`) containing all data, plus comparative plots (heatmaps, bar charts) in `results/05_Integrated_Pipeline_Analysis/`.

#### `6_native_visualizations.R`
-   **Purpose:** Compares the **native plotting functions** (`SpatialDimPlot` vs. `spatPlot`) of each pipeline.
-   **Inputs:** The `final_results_object.rds` from script 5.
-   **Outputs:** Side-by-side plots showing how each tool visualizes its own data, saved in `results/06_Native_Pipeline_Visualizations/`.

#### `7_final_report_generation.R`
-   **Purpose:** Generates the **final summary figures and a comprehensive text report** that synthesizes all findings.
-   **Key Features:** Creates plots with corrected spatial orientation and visualizes cell signatures for each pipeline.
-   **Inputs:** The `final_results_object.rds` and `cancer_signatures_object.rds` from script 5.
-   **Outputs:** Final comparative plots and a `FINAL_REPORT.txt` in `results/07_Final_Report_and_Visuals/`.

#### `8_deconvolution_analysis.R`
-   **Purpose:** Performs a **simulated deconvolution** to estimate cell type proportions in each spot based on marker gene expression.
-   **Key Features:** Demonstrates a proof-of-concept for deconvolution when single-cell reference data is not available. Generates spatial maps of cell types, cluster composition plots, and diversity analysis.
-   **Inputs:** The `final_results_object.rds` from script 5.
-   **Outputs:** A suite of deconvolution plots in `results/08_Deconvolution_Analysis/`.

#### `9_balanced_deconvolution.R`
-   **Purpose:** A **refined deconvolution analysis** using a smaller, balanced set of markers to analyze the Tumor Microenvironment (TME) without the overwhelming signal from fibroblasts.
-   **Key Features:** Compares full tissue composition vs. TME-only, and calculates Immune-to-Tumor ratios.
-   **Inputs:** The `final_results_object.rds` from script 5.
-   **Outputs:** Advanced comparative plots focusing on the TME in `results/09_Balanced_Deconvolution/`.