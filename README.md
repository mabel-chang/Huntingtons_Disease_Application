# R Shiny Huntington’s Disease Application BF591 (F23)
## Overview:
The application is an interactive R Shiny dashboard designed to facilitate the exploration and analysis of Huntington’s Disease (HD) transcriptomic data from Labadorf et al. (GEO accession: GSE64810). It offers a user-friendly interface for examining sample metadata, normalized gene expression counts, differential expression (DE) results, and gene set enrichment analysis (GSEA).

The study utilized mRNA sequencing of post-mortem BA5 human brain tissue samples, including:
* 20 samples from individuals diagnosed with Huntington’s disease
* 45 age- and sex-matched neurologically healthy controls

## Features
### 1. Samples Tab
* Upload sample metadata (CSV format)
* Summary: Displays column names, data types, and numeric column means
* Table: Interactive, searchable, and sortable data table
* Plots: Histogram of user-selected numeric metadata variable

### 2. Counts Tab
* Upload normalized counts matrix (CSV format from GEO)
* Apply filters with sliders for:
  * Gene variance percentile
  * Minimum number of non-zero samples
* Summary: Table of passing vs. non-passing genes after filtering
* Diagnostics: Scatterplots of variance and zero-count thresholds
* Heatmap: Log-transformed gene counts after filtering
* PCA: Principal Component Analysis (interactive axis selection)

### 3. Differential Expression (DE) Tab
* Upload DE results (CSV from GEO supplementary files)
* Input: Interactive table of raw DE results
* Result:
  * Volcano Plot with customizable axes, highlight colors, and significance threshold
  * Filtered Table based on adjusted p-value cutoff

### 4. Gene Set Enrichment Analysis (GSEA) Tab
* Upload preprocessed GSEA results (CSV file)
* Barplot: Top positive and negative enriched pathways
* Table:
  * Filter by adjusted p-value threshold
  * Filter by positive/negative enrichment scores (NES)
  * Download filtered table as CSV
* Scatterplot: NES vs. -log10(adjusted p-value) with threshold highlighting

## Data Sources
* Sample metadata: SRA Run Selector metadata for Huntington’s Disease dataset
* Normalized counts matrix: Downloaded from GEO accession (supplementary files)
* Differential expression results: GEO supplementary data
* GSEA analysis: Performed with MSigDB C2 curated gene sets

## Challenges
During development, the major challenge was integrating standalone code into the Shiny environment.
For example, while tableOutput worked in testing, it failed inside the Shiny framework for DE results.
This issue was resolved by switching to DT (DataTables) for rendering interactive tables.

## Requirements
The application was built in R using the following packages:
* shiny
* DT
* tidyverse
* colourpicker
* ggplot2
* reshape2
* pheatmap
* fgsea

## Running the Application
Clone this repository:

git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>

Install required R packages (if not already installed):
install.packages(c("shiny", "DT", "tidyverse", "colourpicker", "ggplot2", "reshape2", "pheatmap"))
BiocManager::install("fgsea")

Launch the Shiny app:
shiny::runApp("app.R")

## Usage
* Download the GSE64810 dataset from NCBI GEO
* Prepare your metadata, normalized counts, DE results, and GSEA results as CSV files.
* Upload the relevant files into each corresponding tab.
* Explore summaries, visualizations, and statistical analyses interactively.

## Acknowledgements:
* Data from
    Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC, Hadzi TC, Bregu J, et al. (2015) RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. PLoS ONE 10(12): e0143563. https://doi.org/10.1371/journal.pone.0143563
* Huntington’s Disease mRNA-seq study (GSE64810)
* MSigDB C2 curated gene sets for enrichment analysis

## Author:
Mabel Chang
Developed as the final project for BF591: Bioinformatics Programming & Applications at Boston University
