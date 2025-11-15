\# Data Directory



This directory contains the necessary input data for the `seurigiotto-benchmark` project.



\## Included Files



\-   `Tumor.csv`

\-   `Fibroblasts.csv`

\-   `Macrophages.csv`

\-   `Neutrophils.csv`



These files contain lists of marker genes used for cell type signature scoring and simulated deconvolution analyses.



\## Required Data (Not Included)



For the analysis pipelines to run, you must download the main spatial transcriptomics dataset from 10x Genomics.



1\.  \*\*Dataset:\*\* Visium HD CytAssist Gene Expression of Human Lung Cancer (Fixed Frozen)

2\.  \*\*Download Link:\*\* \[10x Genomics Datasets](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-lung-cancer-fixed-frozen)

3\.  \*\*Action:\*\* Download the "Gene Expression" data package. After unzipping, you will have a folder containing the `filtered\\\_feature\\\_bc\\\_matrix.h5` file and a `spatial` subdirectory.

4\.  \*\*Placement:\*\* Rename this main folder to `visium\\\_dataset` and place it \*\*inside this `data/` directory\*\*.



The final structure should look like this:



seurigiotto-benchmark/

└── data/

├── visium\_dataset/

│ ├── filtered\_feature\_bc\_matrix.h5

│ └── spatial/

│ └── ... (spatial data files)

├── Tumor.csv

└── ... (other signature files)

