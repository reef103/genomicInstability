# genomicInstability R package

**Version:** 1.0.0  
**Dependencies:** mixtools  

---
This document describes the basic capabilities of the *genomicInstability* package for R.
The package is devoted to the estimation of Copy Number Variation (CNV) from RNA-Seq data, by estimating the association between gene expression co-regulation and coding-genes chromosomal location.
This information can then be integrated across all loci and chromosomes to infer an index of genomic instability--or Genomic Instability Score (GIS).

The package includes functions to infer CNVs from gene expression data, infer the GIS, and plot the results.

## Use-case examples

### Loading the package and example data
The following code loads the package and 2 datasets:

1. scRNA-Seq expression profile for a syngeneic model of melanoma (B16F10) implanted s.c. in C57BL6 mice.

2. scRNA-Seq data for a compendium of cells from 8 mouse organs which is intended to be used as reference (null model) for the semi-supervised approach (see below).

```r
# Load the package
library(genomicInstability)

# Loading the expression data
dset <- readRDS("melanoma-scrnaseq.rds")

# Loading the expression data for the NULL model
dnull <- readRDS("Null_model_mat.rds")
```

### Unsupervised CNV inference

**CNVs can be inferred** using a completely unsupervised approach, where only the expression matrix to be analyzed and the species, either mouse or human, have to be provided to the `inferCNV()` function:

```r
# Infering the CNVs without use of a NULL model
results <- inferCNV(dset, species="mouse")
```

This function creates an *inferCNV*-class object, which is a list with a slot named `nes`, containing the matrix of inferred CNV results, with samples/cells in columns and chromosome fragments in rows.


The **genomic instability score (GIS)** can then be estimated with the function `genomicInstabilityScore()`, and the results visualized with the function `giDensityPlot`.

```r
# Computing the genomic instability score
results <- genomicInstabilityScore(results)

# Plotting results
giDensityPlot(results)
```

The function `genomicInstabilityScore()` updates the *inferCNV*-class object including slotS called `gis` and `gi_likelihood`.

`genomicInstabilityScore()` also computes the genomic instability likelihood, as the relative likelihood between the first and the last distributions of a mixture Gaussian model fit to the genomic instability score results (Fig. 1).

![Unsupervised GIS density](./density-unsupervised.png)  
Figure 1: Density distribution plot for the genomic instability score and mixture of 3 Gaussian distributions fitted to the data.


The results for the inferredCNVs can can also be plotted with the `plot()` function.
An example plot is shown in Fig. 2.

```r
# Plotting results
plot(results, output="cnv-unsupervised.png", gamma=1.2, resolution=48)
```
![Unsupervised GIS density](./cnv-unsupervised.png)  
Figure 2. Plot of the inferred CNVs.

### Semi-supervised inference of CNVs

**CNVs can be inferred** using a semi-supervised approach, where an additional expression matrix, representing normal tissues, can be used as reference.
Below is an example using a compendium of normal cells from 8 mouse organs, including the inference of the CNVs and estimation of the genomic instability score.

```r
# Infering the CNVs
results <- inferCNV(dset, dnull, species="mouse")

# Computing the genomic instability score
results <- genomicInstabilityScore(results)
```

The genomic instability results visualized with the function `giDensityPlot`.

```r
# Plotting results
giDensityPlot(results)
```

When a null model (reference) is provided, the genomic instability likelihood is estimated as the relative likelihood between the null model (shown in black, Fig. 3) and the last distribution of a mixture Gaussian model fit to the genomic instability score results.

![Supervised GIS density](./density-supervised.png)  
Figure 3: Density distribution plot for the genomic instability score and mixture of 3 Gaussian distributions fitted to the data. The Distribution of the reference dataset is shown in black.


The results for the inferredCNVs can can also be plotted with the `plot()` function (Fig. 4):

```r
# Plotting results
plot(results, output="cnv-supervised.png", resolution=48)
```

![Supervised GIS density](./cnv-supervised.png)  
Figure 4. Plot of the inferred CNVs.

## Advanced functionality
The functions provided by package *genomicInstability* are organized in three functional categories:

1. CNV inference, with the function `inferCNV()`
2. Genomic Instability Score estimation, with the function `genomicInstabilityScore()`
3. Genomic unstable (tumor) and stable (normal) samples or cells inference, with the function `giLikelihood()`
4. Graphical display of the results, with the functions `plot()` and `giDensityPlot()`

### 1. CNV inference
The function **inferCNV()** predicts loci with Copy Number Variation (CNV) based on the co-expression of transcripts coded by genes located closely in the chromatin.

#### Arguments

- **expmat**: An expression matrix, of genes (rows) by samples (columns), where the samples can also be single cells.
It is recommended to correct classical RNASeq data for gene-length, but since the analysis is performed for each sample/cell independently, no sequencing depth normalization is required.

- **nullmat**: An optional expression matrix to be used as normal reference.
This can be a set of tissue lineage-matching normal samples or cells, or a pool of normal cells from different tissue lineages.
The matrix should be normalized by gene-length in the case of RNASeq and represent genes (rows) by samples/cells (columns).

- **species**: This is a string indicating the genome species, and can be either "human" or "mouse".

- **k**: Integer indicating the number of genes per chromosomal window to analyze. The default is 100.

- **skip**: Integer indicating the distance between consecutive genome windows expressed as number of genes. The default is 25.

- **min_geneset**: Integer >1 indicating the minimum size allowed for the genomic windows.

#### Return

Object of class *inferCNV*, which is a list containing the following slots:

- **nes**: Matrix of Normalized Enrichment Scores (NES) for chromosomal windows (rows) by samples/cells (columns) corresponding to the expression matrix "expmat".
- **null**: Matrix of NES for the chromosomal windows (rows) by samples/cells (columns) corresponding to the normal samples/cells used as reference ("nullmat").
- **param**: List of input parameters.

### 2. Genomic Instability Score (GIS) estimation

The function `genomicInstabilityScore()` takes an object of class *genomicInstability* with a valid *nes* slot and computes a score of genomic instability as the log2(variance) for each sample/cell across all chromosome windows (loci).
The results are stored in two additional slots:

- **gis**: Vector of genomic instability scores for the samples represented in "expmat".

- **gisnull**: Vector of genomic instability scores for the samples represented in "nullmat".

### 3. Genomic Instability Likelihood

The genomic instability likelihood is computed by fitting a mixture of *k* gaussian models to the GIS data.
As default, the function will try to fit between 1 and 3 models, but the behavior can be fine-tuned.

#### Arguments

- **inferCNV**: Object of class inferCNV with a *nes* slot.

- **recompute**: Logical whether the stored gaussian-fit information, if present, should be ignored and re-computed based on the data.

- **distros**: Vector of 2 numbers indicating the minimum and maximum number of gaussian models to fit. Default is c(1, 3).
It is recommended to check the models fitted to the data with the `giDensityPlot()` function (see below).
This will assist on the election of the number of models to fit.

- **tumor**: Optional vector of integers indicating the positions for the models assumed to correspond to the tumors.
Such positions are sorted based on the mean of each model.
It is recommended to check the models fitted to the data with the `giDensityPlot()` function (see below).
This will assist on the election of the models most likely representing genomically unstable (tumor) populations of samples/cells.

- **normal**: Optional vector of integers indicating the positions for the normal models.
The distribution corresponding to the data-set provided as reference by the "nullmat" argument of `inferNV()` function will be used as default if provided.
If not provided, the model with the lowest mean will be used.
Same as before, it is recommended to use `giDensityPlot() to assist on the selection of the models best representing genomically stable (normal) populations of samples/cells.

#### Details

This function fit a mixture of gaussian models to the distribution of GIS.
By default, the function fits up to 3 models, and considers the right-most (higher mean) as representing the tumors and the left-most (lower mean) as representing the normal samples or cells.
This behavior can be modified by:

1. Selecting the number of gaussian models to fit
2. Selecting the distributions representing the tumor and normal samples

### 4. Display of results

The `genomicInstability` package provides the function `plot()`--plot.inferCNV()--for plotting the inferred CNV on the different chromosomes, and `giDensityPlot()` to plot the distribution of the GIS and the fitted Gaussian models.

#### Density distribution plot for GIS and models

This can be obtained by the `giDensityPlot()` function.
It is recommended to generate a diagnostic plot with `giDensityPlot()`, so the results of `giLikelihood()` can be better interpreted and decisions can be make, based on the models fitted to the distribution of the data, regarding the models considered as representing genomic stable (normal) and genomic unstable (tumor) samples/cells.

##### Arguments

- **inferCNV**: An object of class *genomicInstability* containing the results generated by the function `giLikelihood()`.

- **legend**: String indicating where to locate the legend, either "top", "topleft", or "topright", It also accepts "none" as indication that the legend should not be included.

- **Additional optional arguments passed to the plot() function**, including for example, defining the `ylim`.

#### Chromosomal display of the inferred CNVs

This plot, generated by the `plot()` function when using a *genomic Instability*-class object as argument, display the areas of the genome (loci) showing enrichment on co-activated transcription (potentially the result of genetic amplifications), or co-repressed transcription (potentially genetic deletions).
Potential amplifications and deletions are shown in red and blue color, respectively, in a heatmap representing samples or cells in rows, and genomic locations (loci) in columns.
The color intensity (gamma) as well as the automatic identification of "stable", "intermediate", and genomically "unstable" samples/cells parameters can be specified.

##### Arguments

- **inferCNV**: Object of class *genomicInstability*.

- **output**: Optional string indicating the name of the output file for the plot.
If provided, it should include the file extension, which could be ".pdf", ".png", or ".jpg", and which indicate the format of the output file.
If omitted, the plot will be generated on the default output device.

- **threshold**: Numeric between 0 and 1 indicating the likelihood threshold to define the "stable", "intermediate" and genomically "unstable" populations.

- **gamma**: Positive number indicating the gamma transformation for the colors in the plot.

- **resolution**: Positive integer indicating the resolution, in ppi, for the png and jpg file formats.