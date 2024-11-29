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

  **input**: `norm_small.csv` (genes 1-1000 only bc original file too large)

DE tab

  **input**: `DESeq2_diffexp.csv`
  
  