#' @details
#' The basic functionality of this package can be performed by inferCNV(),
#' to infer the enrichment of loci-blocks on gene expresion;
#' genomicInstabilityScore(), to estimate the genomic instability
#' for each of the cells in the scRNASeq dataset; giLikelihood(), to estimate
#' the relative likelihood for each cell to be normal (low genomic instability)
#' or tumor (high genomic instability); plot() and giDensityPlot() to plot the
#' scores per loci-block and the distribution of the genomic instability score,
#' respectively.
#' @seealso [inferCNV()] for estimating loci-block enrichment,
#' [genomicInstabilityScore()] for estimating the genomic instability of each
#' cell in the dataset, [giLikelihood()] for estimating the relative
#' likelihood for the cells to be normal or neoplastic, [plot.inferCNV()] and
#' [giDensityPlot()] to plot the results.
#' @keywords internal
"_PACKAGE"
