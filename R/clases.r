# Initializationa and dependencies

#' @importFrom mixtools normalmixEM
#' @importFrom SummarizedExperiment assays
#' @importFrom checkmate assertInt
#' @importFrom checkmate assertMatrix
#' @importFrom checkmate assertLogical
#' @importFrom checkmate assertIntegerish
#' @importFrom checkmate assertClass
#' @importFrom checkmate assertNumeric
NULL

#' Chromosomal coordinate of human and mouse known genes
#' 
#' A dataset containing the chromosomal coordinate for known human
#' and mouse genes
#' 
#' @format data.frame with 2 columns: Chromosome and Coordinate.
#' To access this data use:
#' \describe{
#'   \item{data(hg38)}{Human}
#'   \item{data(mm10)}{Mouse}
#' }
"genePosition"

#' Average length of human and mouse known genes
#' 
#' A dataset containing the average length for known mouse
#' and human genes
#' 
#' @format Vector of integers indicating the average length in bp
#'   for each gene, indicated with EntrezIDs as name argument.
#' To access this data use:
#' \describe{
#'   \item{data(hg38)}{Human}
#'   \item{data(mm10)}{Mouse}
#' }
"geneLength"
