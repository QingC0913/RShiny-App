## RShiny-App
Name: Qing Cheng

Date: Dec 15, 2024

BF591 R for Biological Sciences 

Boston University 

## Data:
Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls

This dataset profiled gene expression with RNASeq in post-mortem human dorsolateral prefrontal cortex from patients who died from
Huntington’s Disease and age- and sex-matched neurologically healthy controls.

[Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)

## Preprocessing Steps: 

1. Get relevant samples information from series matrix using python

  **Input**: `GSE64810_series_matrix.csv`\
  **Processing**: `series.py`\
  **Output**: `processed_series_matrix.csv`

  **Input**: `processed_series_matrix.csv.csv`\
  **Processing**: `series.py`\
  **Output**: `processed_samples.csv`
  
2. formatting in R

  **Input**: `processed_samples.csv`\
  **Processing**: `preprocessing.Rmd`\
  **Output**: `samples_matrix.csv`

## Content:
### Sample Information Exploration tab

**Input**: `samples_matrix.csv`\
**Functionality**: 

* Data Information 
* Sample Summary
* Data Table
  * sortable columns 
* Sample Exploration Plots
  * User can choose which column to plot
  
  
### Counts Matrix Exploration tab

**Input**: `norm_counts.csv` or `norm_DESeq2.csv`\
**Functionality**: 

* Data Summary
* Diagnostic Plots
  * filter by percentile variance 
* Counts Heatmap
  * Clustered heatmap of filtered genes
* Principal Component Analysis 
  * User can choose which two PCs to plot 

### Differential Expression Analysis tab
**Input**: `DESeq2_diffexp.csv`\
**Functionality**: 

* Differential Expression Results 
  * Filter results table by log2FoldChange and p*adj
* Visualizations
  * Filter significant results on volcano plot by log2FC and p*adj
  
### Correlation Network Analysis tab
**Input**: same as COUNTS tab\
**Functionality**: 

* Selected Genes Heatmap
  * User can upload a list of genes of interest
  * Clustered heatmap of correlation among genes of interest 
  * Clustered heatmap of gene counts 
* Correlation Network 
  * Network graph of genes 
  * User can select graph layout for better visualization 
  * User can adjust minimum correlation threshold for edges
  * Shortest path between any two vertices
* Network Metrics 
  
  
