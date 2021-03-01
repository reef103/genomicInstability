# Functions required for the analysis and inference of genetic intability from
# RNA-Seq

#' Topological gene sets
#'
#' This function generates a list of sets of k genes encoded by neighbor loci
#'
#' @param species Character string indicating the species, either human or mouse
#' @param k Integer indicating the number of genes per set
#' @param skip Interger indicating the displacement of the window for selecting
#' the k genes
#'
#' @return List of topoligically-close gene sets
#'
#' @examples
#' chrom_set <- generateChromosomeGeneSet('human')
#' length(chrom_set)
#' chrom_set[1:2]
#'
#' @export
generateChromosomeGeneSet <- function(species = c("human", "mouse"), k = 100,
skip = 25) {
    # validate options for species
    species <- match.arg(species)
    # Validate input conditions
    checkmate::assertInt(k, lower = 10, upper = 1000)
    checkmate::assertInt(skip, lower = 1, upper = k)
    # Get the gene coordinates
    switch(species, human = data(hg38, package = "genomicInstability",
        envir = environment()),
        mouse = data(mm10, package = "genomicInstability",
        envir = environment()))
    # Generate vector of positions per gene
    gene_position <- as.numeric(as.vector(genePosition[, 2]))
    names(gene_position) <- rownames(genePosition)
    # Sorted list of gene position per chromosome
    genes_per_chromosome <- tapply(gene_position, genePosition[, 1], sort)
    # Sorted list of genes per chromosome
    genes_per_chromosome <- lapply(genes_per_chromosome, names)
    # Remove alternative and control chromosomes
    pos <- unique(c(grep("_alt", names(genes_per_chromosome)),
    grep("_random", names(genes_per_chromosome))))
    # Chromosomes with less than 10 genes
    pos <- unique(c(pos, which(sapply(genes_per_chromosome, length) < k)))
    genes_per_chromosome <- genes_per_chromosome[-pos]
    # Genesets per chromosome
    geneset_chromosome <- lapply(genes_per_chromosome, splitVectorWindow, k = k,
        skip = skip)
    # Flatten the chromosomal structure for the genesets
    geneset <- unlist(geneset_chromosome, recursive = FALSE, use.names = FALSE)
    # Name the genesets by chromosome-set
    names(geneset) <- paste(rep(names(geneset_chromosome),
        sapply(geneset_chromosome, length)), unlist(lapply(geneset_chromosome,
        names), use.names = FALSE), sep = "-")
    return(geneset)
}

#' Inference of CNV from expression data
#'
#' This function estimates the CNV score based on expression data
#'
#' @param expmat Matrix of gene expression profiles or signatures with genes
#' '(entrezID) in rows and samples in columns
#' @param nullmat Optional matrix with same number of rows as \code{expmat} to
#' be used as null model
#' @param species Character string indicating the species, either human or mouse
#' @param k Integer indicating the number of genes per set
#' @param skip Interger indicating the displacement of the window for selecting
#' the k genes
#' @param min_geneset Integer indicating the minimum size for the genesets
#' @param verbose Logical, whether progress should be reported
#'
#' @return Object of class inferCNV, which is a list containing matrix of nes,
#' and parameters (param), including species, window (k) and skip
#'
#' @examples
#' data(GSE103322, package='HNSCgenomicInstability')
#' cnv <- inferCNV(tpm_matrix)
#' class(cnv)
#' names(cnv)
#' cnv$nes[1:5, 1:3]
#'
#' @export
inferCNV <- function(expmat, nullmat = NULL, species = c("human", "mouse"),
    k = 100, skip = 25, min_geneset = 10, verbose = TRUE) {
    # Check values for species
    species <- match.arg(species)
    # Validate input
    checkmate::assertMatrix(expmat, mode = "numeric", all.missing = FALSE,
        min.rows = 1000, min.cols = 1, row.names = "named", col.names = "named")
    checkmate::assertMatrix(nullmat, mode = "numeric", all.missing = FALSE,
        min.rows = 1000, min.cols = 1, row.names = "named", null.ok = TRUE)
    checkmate::assertInt(k, lower = 10, upper = 1000)
    checkmate::assertInt(skip, lower = 1, upper = k)
    checkmate::assertInt(min_geneset, lower = 2, upper = k)
    checkmate::assertLogical(verbose, len = 1)
    # Compatibilize null model
    if (!is.null(nullmat)) {
        genes <- intersect(rownames(expmat), rownames(nullmat))
        if (length(genes) < 100)
            stop("Genes in expmat and nullmat do not match")
        expmat <- expmat[match(genes, rownames(expmat)), , drop = FALSE]
        nullmat <- nullmat[match(genes, rownames(nullmat)), , drop = FALSE]
    }
    # Generate genesets
    if (verbose)
        message("Generating the genesets from the genome information")
    geneset <- generateChromosomeGeneSet(species, k = k, skip = skip)
    # filter represented genes
    geneset <- lapply(geneset, keepVectorElements, elements = rownames(expmat))
    # keep genesets with at least min_geneset genes
    geneset <- geneset[sapply(geneset, length) >= min_geneset]
    if (length(geneset) == 0)
        stop(paste0("No geneset with at least ", min_geneset, " genes"))
    # Enrichment for the expmat
    if (verbose)
        message("Computing the enrichment for the genesets in the expression
            matrix")
    expmat_nes <- sREA(expmat, geneset)
    # Enrichment of the null model
    nullnes <- NULL  # Initialize the null nes variable
    if (!is.null(nullmat)) {
        if (verbose)
            message("Computing null model")
        expmat_null <- sREA(nullmat, geneset)
        # Estimating NES
        if (verbose)
            message("Estimating the normalized enrichment scores")
        expmat_nes <- t(sapply(seq_len(nrow(expmat_nes)),
            computeNesForMatrixRow, nesmat = expmat_nes, nullmat = expmat_null))
        nullnes <- t(sapply(seq_len(nrow(expmat_null)), computeNesForMatrixRow,
            nesmat = expmat_null, nullmat = expmat_null))
        rownames(expmat_nes) <- rownames(nullnes) <- names(geneset)
    }
    # Returning the results
    res <- list(nes = expmat_nes, null = nullnes,
        param = list(species = species, k = k, skip = skip))
    class(res) <- "inferCNV"
    return(res)
}


#' Genomic Instability Analysis
#'
#' This function computes the genomic instability for an object of class
#' inferCNV
#'
#' @param cnv Object of class inferCNV generated by inferCNV() function
#' @param likelihood Logical, whether the genomic instability likelihood should
#' be estimated
#'
#' @return Object of class inferCNV with updated slots for gis and gisnull
#'
#' @examples
#'
#' data(GSE103322, package='HNSCgenomicInstability')
#' cnv <- inferCNV(tpm_matrix)
#' cnv <- genomicInstabilityScore(cnv)
#' plot(density(cnv$gis))
#'
#' @seealso [inferCNV()] to infer the enrichment of loci-blocks in the gene
#' expression data.
#'
#' @export
genomicInstabilityScore <- function(cnv, likelihood = FALSE) {
    # Validating inputs
    checkmate::assertLogical(likelihood, len = 1)
    validateInferCNV(cnv, "nes")
    # Compute genomic instability score
    gis <- log2(colVars(cnv[["nes"]]))
    gisnull <- NULL  # Initialize gisnull
    if (!is.null(cnv[["null"]])) {
        gisnull <- log2(colVars(cnv[["null"]]))
    }
    cnv[["gis"]] <- gis
    cnv[["gisnull"]] <- gisnull
    if (likelihood) {
        cnv <- giLikelihood(cnv)
    }
    return(cnv)
}

#' Genomic instability likelihood
#'
#' This function computes the genomic instability likelihood
#'
#' @param inferCNV InferCNV-class object
#' @param recompute Logical, whether the model fits should be re-computed
#' @param distros Vector of 2 integers indicating the minimum and maximum number
#' of Gaussian models to fit
#' @param tumor Optional vector of integers indicating the Gaussians considered
#' as tumors
#' @param normal Optional vector of integers indicating the Gaussians considered
#' as normal. This is only useful when no null model has been provided for the
#' analysis
#'
#' @return Updated inferCNV-class object with gi_likelihood slot
#'
#' @examples
#'
#' data(GSE103322, package='HNSCgenomicInstability')
#' cnv <- inferCNV(tpm_matrix)
#' cnv <- genomicInstabilityScore(cnv)
#' cnv <- giLikelihood(cnv, distros=c(3, 3), tumor=2:3)
#' print(cnv$gi_fit)
#' plot(density(cnv$gi_likelihood, from=0, to=1))
#'
#' @seealso [genomicInstabilityScore()] to estimate the genomic instability
#' score for each cell in the dataset, and [inferCNV()] to infer the enrichment
#' of loci-blocks in the gene expression data.
#' @export
giLikelihood <- function(inferCNV, recompute = TRUE, distros = c(1, 3),
    tumor = NULL, normal = NULL) {
    # Validate inputs
    validateInferCNV(inferCNV, "nes")
    checkmate::assertLogical(recompute, len = 1)
    checkmate::assertIntegerish(distros, lower = 1, upper = 10, len = 2,
        any.missing = FALSE)
    checkmate::assertIntegerish(tumor, lower = 1, upper = max(distros),
        any.missing = FALSE, null.ok = TRUE)
    checkmate::assertIntegerish(normal, lower = 1, upper = max(distros),
        any.missing = FALSE, null.ok = TRUE)
    # sort distros indexes
    distros <- sort(distros)
    # Compute GIS if not done previously
    if (is.null(inferCNV[["gis"]]))
        inferCNV <- genomicInstabilityScore(inferCNV)
    # Fit mixture gaussians to the results
    if (is.null(inferCNV[["gi_fit"]]) | recompute) {
        results_fit <- mixGaussianFit(inferCNV[["gis"]], min = distros[1],
            max = distros[2])
        inferCNV[["gi_fit"]] <- results_fit
    }
    results_fit <- inferCNV[["gi_fit"]]
    # If no tumor models are selected, asign the last gaussian
    if (is.null(tumor))
        tumor <- 2:length(results_fit[["mu"]])
    if (any(tumor > length(results_fit[["mu"]])))
        stop("Tumor selection is larger than the number of models",
            call. = FALSE)
    # If there is a null model
    if (!is.null(inferCNV[["gisnull"]]) & is.null(normal)) {
        # Keep only the null model and the ones selected as tumor
        results_fit[["mu"]] <- c(mean(inferCNV[["gisnull"]], na.rm = TRUE),
            results_fit[["mu"]][tumor])
        results_fit[["sigma"]] <- c(sd(inferCNV[["gisnull"]], na.rm = TRUE),
            results_fit[["sigma"]][tumor])
        # Adjust the tumor positions to account for the removed models and the
        # null model
        tumor <- (seq_len(length(tumor))) + 1
    } else {
        if (is.null(normal))
            normal <- 1
        if (any(normal > length(results_fit[["mu"]])))
            stop("Normal selection is larger than the number of models",
                call. = FALSE)
        # Keep only selected models
        results_fit[["mu"]] <- results_fit[["mu"]][c(normal, tumor)]
        results_fit[["sigma"]] <- results_fit[["sigma"]][c(normal, tumor)]
        # Adjust normal and tumor indexes
        normal <- seq_len(length(normal))
        tumor <- (seq_len(length(tumor))) + length(normal)
    }
    # Compute relative likelihood
    inferCNV[["gi_likelihood"]] <- rowSums(predict(results_fit,
        inferCNV[["gis"]], tumor))
    return(inferCNV)
}
