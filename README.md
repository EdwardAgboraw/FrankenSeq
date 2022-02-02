# FrankenSeq

## Installation

FrankenSeq can be installed from Github using the following R commands:
```
# install the package
library("devtools")
install_github("EdwardAgboraw/FrankenSeq")

# launch the application
library(FrankenSeq)
launch_FrankenSeq()
```

Note: FrankenSeq was built in R version 4.1.1.

## About

FrankenSeq is a comprehensive and modular analysis platform for the clustering of scRNA-seq data.

FrankenSeq accepts either SingleCellExperiment/Seurat object RDS files, or Gene x Cell matrix CSV files as input.

FrankenSeq is composed of eight tabs - ‘Data Processing’, ‘Feature Selection’, ‘Dimension Reduction’, ‘Cluster Validation’, ‘Cluster Analysis’, ‘Results’, 'DEG Analysis', and 'Deep Learning'

These tabs are intended to be used in sequence, with each accepting as input the output of the tab before it. The output of the first five tabs is a Seurat object containing the combined results of the analysis so far. The output of ‘Results’ is a CSV file containing either a summary of the final cluster assignments, the full cluster assignments, or the specific genes selected by the chosen feature selection algorithm.

The 'Feature Selection', 'Dimension Reduction', and 'Cluster Analysis' tabs each provide multiple different options for conduction the titular analysis.

The FrankenSeq pipeline is reactive - any changes made to earlier steps in the analysis cause alterations in the results of any steps below them. For example, changing the Feature Selection option results in the recalculation of both the dimension reduction and cluster analysis results.

There are two exceptions to the above. The "Cluster Validation" Tab is fully optional - "Cluster Analysis" can be done with or without it having been run first. And the "Deep Learning" tab is an alternative to "Cluster Analysis", and as such takes its input from the "Dimension Reduction" tab (instead of the tabs directly preceeding it).

The basic structure of the FrankenSeq Analysis Pipeline is shown below:

