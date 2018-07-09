# sc_batch_remove

## About
sc_batch_remove is wrapper for scRNA-seq analysis batch removing software.

## Requirements
- R (>= 3.5.0)
  - Seurat
  - sva
  
## Example
### datasets
- pbmc_4k_filtered
- pbmc_8k_filtered

datasets are available from [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets)


### Command
```
$ Rscript seurat_combat.R pbmc_4k_filtered pbmc_8k_filtered
```

### Standard Output
```
$ Rscript seurat_combat.R pbmc_4k_filtered pbmc_8k_filtered
Loading required package: ggplot2
Loading required package: cowplot

Attaching package: ‘cowplot’

The following object is masked from ‘package:ggplot2’:

    ggsave

Loading required package: Matrix

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Loading required package: mgcv
Loading required package: nlme

Attaching package: ‘nlme’

The following object is masked from ‘package:dplyr’:

    collapse

This is mgcv 1.8-23. For overview type 'help("mgcv-package")'.
Loading required package: genefilter
Loading required package: BiocParallel
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|

========== batch distribution before filtering ===========
pbmc_4k_filtered pbmc_8k_filtered
            4340             8381
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene means
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variance to mean ratios
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
[1] 436
Scaling data matrix
  |======================================================================| 100%
pdf
  2
  |======================================================================| 100%
Time Elapsed:  37.0334763526917 secs

An object of class seurat in project pbmc_4k_filtered
 33694 genes across 5224 samples.
Warning message:
Removed 6881 rows containing missing values (geom_point).
pdf
  2
Parameters used in latest FindClusters calculation run on: 2018-07-09 14:50:55
=============================================================================
Resolution: 0.6
-----------------------------------------------------------------------------
Modularity Function    Algorithm         n.start         n.iter
     1                   1                 100             10
-----------------------------------------------------------------------------
Reduction used          k.param          prune.SNN
     pca                 30                0.0667
-----------------------------------------------------------------------------
Dims used in calculation
=============================================================================
1 2 3 4 5 6 7 8 9 10

pdf
  2
pdf
  2

========== batch distribution after filtering ===========
   1    2
3500 1724

pdf
  2
Regressing out: batchid
  |======================================================================| 100%
Time Elapsed:  24.1224076747894 secs
Scaling data matrix
  |======================================================================| 100%
pdf
  2
Found2batches
Adjusting for0covariate(s) or covariate level(s)
Standardizing Data across genes
Fitting L/S model and finding priors
Finding parametric adjustments
Adjusting the Data

Scaling data matrix
  |======================================================================| 100%
pdf
  2
```

### Outputs

