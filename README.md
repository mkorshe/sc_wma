
# Optimizing integration of single-cell eQTL summary statistic #

This project focuses on optimizing the integration of single-cell eQTL summary statistics by exploring alternative parameters for meta-analysis besides the number of donors, such as the number of cells per donor or the number of reads per cohort. The main goal of this project is to improve the accuracy and reliability of single-cell eQTL meta-analysis results.

## Project Structure
The project is organized into three main folders:

### 1. Preprocessing 
This folder contains scripts for processing and filtering single-cell eQTL summary statistics.

### 2. Weighted Metaanalysis
This folder contains scripts for conducting weighted meta-analysis using different sample-size-like parameters.

### 3. Postprocessing of WMA results
This folder contains scripts for post-processing and visualizing the results of the weighted meta-analysis.

## Weighted Metaanalysis Concept
In single-cell studies, the number of donors (sample size) is a common parameter used for conducting meta-analysis. However, there are several other parameters that could be considered as alternatives to the number of donors. For instance, the number of cells per donor or the number of reads per cohort could be tested as alternative parameters for conducting meta-analysis of single-cell summary statistics.

Different studies vary in terms of the number of cells collected and reads (or counts) per cell. Therefore, this project suggests testing several sample-size-like parameters such as the average number of cells per donor, average counts per cell, average counts per donor, total number of cells per cohort, total number of counts per cohort as an alternative to the number of donors in a single-cell eQTL meta-analysis (see fig.1 in the manuscript).

## Manuscript
For more information about the project, please refer to the manuscript available on bioRxiv.
