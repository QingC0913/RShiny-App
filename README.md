# RShiny-App

## proprocessing steps: 

1. get relevant samples information from series matrix using python

  **input**: `GSE64810_series_matrix.csv`\
  **processing**: `series.py`\
  **output**: `processed_samples.csv`
  
2. formatting in R

  **input**: `processed_samples.csv`\
  **processing**: `preprocessing.Rmd`\
  **output**: `samples.csv`

## content
SAMPLES tab

  **input**: `samples.csv`
  
COUNTS tab

  **input**: `norm_counts.csv` or `GSE64810_mlhd_DESeq2_norm_counts_adjust.csv`
  
DE tab

  **input**: `DESeq2_diffexp.csv`
  
DE tab

  same as COUNTS tab
  
  