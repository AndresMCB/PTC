# PTC
PTC package contains functions to determine causal miRNAs-mRNA relationships


## Introduction
Inspired by the pseudo-time concept we develop a novel approach, called the pseudo-time causality (PTC) based approach, to elucidate the miRNA-mRNA interactions during biological processes, using gene expression data with the expression profiles of matched miRNAs and mRNAs in the same cells or samples. Given a biological process, PTC firstly transforms the matched miRNA and mRNA single cell gene expression data to pseudo-time data using the marker genes of the biological process. PTC relies on the causal invariance property [(Peters et al., 2015](https://doi.org/10.1111/rssb.12167); [Pfister et al., 2018)](https://doi.org/10.1080/01621459.2018.1491403) to find the causal relationships
between miRNAs andmRNAs from the pseudo-time data. 

We have applied PTC to the single cell dataset from [Wang et al., 2019](https://doi.org/10.1038/s41467-018-07981-6) and the bulk data from [Pham et al., 2019](https://doi.org/10.1186/s12859-019-2668-x). In both cases, VIM, an EMT marker, is used to define the pseudo-time.
The results have shown that PTC significantly outperforms the benchmark methods in identifying experimentally confirmed miRNA-mRNA interactions using either single cell or bulk data. The results suggest that the temporal information during a biological process is useful for revealing the miRNA-mRNA interactions characterising the biological process.

## Installation 
PTC runs in the R statistical computing environment. R version 3.6.1 or higher and Bioconductor version 3.11 or higher are required.
1. Please install Bioconductor, you can use the following code in R

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.11")
```
2. Install Bioconductor dependencies required by PTC 
```R
BiocManager::install(c('miRBaseConverter', 'CancerSubtypes'))
```
3. Install PTC package from github repository 
```R
devtools::install_github('AndresMCB/PTC')
```
## Documentation 
Detailed information about the functions implemented in PTC can be found in the [user manual](PTC_1.0.0.pdf)

Please find the datasets employed in our paper in the folder [data](data/)

