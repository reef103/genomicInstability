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
#' chrom_set[seq_len(2)]
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
    pos <- unique(c(pos, which(vapply(genes_per_chromosome, length,
                                      numeric(1)) < k)))
    genes_per_chromosome <- genes_per_chromosome[-pos]
    # Genesets per chromosome
    geneset_chromosome <- lapply(genes_per_chromosome, splitVectorWindow, k = k,
        skip = skip)
    # Flatten the chromosomal structure for the genesets
    geneset <- unlist(geneset_chromosome, recursive = FALSE, use.names = FALSE)
    # Name the genesets by chromosome-set
    names(geneset) <- paste(rep(names(geneset_chromosome),
        vapply(geneset_chromosome, length, numeric(1))),
        unlist(lapply(geneset_chromosome, names), use.names = FALSE), sep = "-")
    return(geneset)
}

#' Inference of CNV from expression data
#'
#' This function estimates the CNV score based on expression data
#'
#' @param expmat Matrix or sparse matrix of gene expression profiles or
#'  signatures with genes '(entrezID) in rows and samples in columns
#' @param nullmat Optional matrix or sparse matrix with same number of rows as \code{expmat} to
#' be used as null model or vector of index positions for the cells in
#' expmat that can be considered normal diploid
#' @param species Character string indicating the species, either human or mouse
#' @param k Integer indicating the number of genes per set
#' @param skip Interger indicating the displacement of the window for selecting
#' the k genes
#' @param min_geneset Integer indicating the minimum size for the genesets
#' @param batch_size Number indicating the maximum ize of the batch of cells
#' to be analyzed by aREA per iteration. This parameter is used to limit
#' memory load when transforming the sparse matrix normalized expression into
#' full matrix ranked signature
#' @param verbose Logical, whether progress should be reported
#'
#' @return Object of class inferCNV, which is a list containing matrix of nes,
#' and parameters (param), including species, window (k) and skip
#'
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' dset <- eh[["EH5419"]]
#' tpm_matrix <- SummarizedExperiment::assays(dset)$TPM
#' set.seed(1)
#' tpm_matrix <- tpm_matrix[, sample(ncol(tpm_matrix), 500)]
#' cnv <- inferCNV(tpm_matrix)
#' class(cnv)
#' names(cnv)
#' cnv$nes[1:5, 1:3]
#'
#' @export
inferCNV <- function(expmat, nullmat = NULL, species = c("human", "mouse"),
    k = 100, skip = 25, min_geneset = 10, batch_size = 1000, verbose = TRUE) {
    # Check values for species
    species <- match.arg(species)
    # Validate input
    checkmate::assert(checkmate::testMatrix(expmat, mode = "numeric",
        all.missing = FALSE, min.rows = 1000, min.cols = 1, row.names = "named",
        col.names = "named"),
        checkmate::testClass(expmat, "Matrix"))
    checkmate::assert(checkmate::testMatrix(nullmat, mode = "numeric",
        all.missing = FALSE, min.rows = 1000, min.cols = 1, row.names = "named",
        null.ok = TRUE),
        checkmate::testNumeric(nullmat, lower = 1, upper = ncol(expmat)),
        checkmate::testClass(nullmat, "Matrix"))
    checkmate::assertInt(k, lower = 10, upper = 1000)
    checkmate::assertInt(skip, lower = 1, upper = k)
    checkmate::assertInt(min_geneset, lower = 2, upper = k)
    checkmate::assertNumber(batch_size, lower = 1, upper = Inf)
    checkmate::assertLogical(verbose, len = 1)
    batch_size <- round(batch_size)
    # Compatibilize null model
    if (length(nullmat) > 0 && !is.null(nrow(nullmat))) {
        genes <- intersect(rownames(expmat), rownames(nullmat))
        if (length(genes) < 100)
            stop("Genes in expmat and nullmat do not match")
        expmat <- expmat[match(genes, rownames(expmat)), , drop = FALSE]
        nullmat <- nullmat[match(genes, rownames(nullmat)), , drop = FALSE]
    }
    # Generate genesets
    if (verbose) {
        message("Generating the genesets from the genome information")
    }
    geneset <- generateChromosomeGeneSet(species, k = k, skip = skip)
    # filter represented genes
    geneset <- lapply(geneset, keepVectorElements, elements = rownames(expmat))
    # keep genesets with at least min_geneset genes
    geneset <- geneset[vapply(geneset, length, numeric(1)) >= min_geneset]
    if (length(geneset) == 0) {
        stop(paste0("No geneset with at least ", min_geneset, " genes"))
    }
    genes <- unique(unlist(geneset, use.names = FALSE))
    # Enrichment for the expmat
    if (verbose) {
        message("Computing the enrichment for the genesets in the expression matrix")
    }
    expmat <- expmat[rownames(expmat) %in% genes, , drop = FALSE]
    expmat_nes <- sREAsm(expmat, geneset, batch_size = batch_size)
    # Enrichment of the null model
    nullnes <- NULL  # Initialize the null nes variable
    if (length(nullmat) > 0) {
        if (verbose) {
            message("Computing null model")
        }
        if (!is.null(nrow(nullmat))) {
            nullmat <- nullmat[rownames(nullmat) %in% genes, , drop = FALSE]
            expmat_null <- sREAsm(nullmat, geneset, batch_size = batch_size)
        } else {
            expmat_null <- expmat_nes[, round(nullmat), drop = FALSE]
        }
        # Estimating NES
        if (verbose) {
            message("Estimating the normalized enrichment scores")
        }
        expmat_nes <- t(vapply(seq_len(nrow(expmat_nes)),
            computeNesForMatrixRow, numeric(ncol(expmat_nes)),
            nesmat = expmat_nes, nullmat = expmat_null))
        nullnes <- t(vapply(seq_len(nrow(expmat_null)), computeNesForMatrixRow,
            numeric(ncol(expmat_null)), nesmat = expmat_null,
            nullmat = expmat_null))
        rownames(expmat_nes) <- rownames(nullnes) <- names(geneset)
    }
    # Returning the results
    res <- list(nes = expmat_nes, null = nullnes,
        param = list(species = species, k = k, skip = skip))
    class(res) <- "inferCNV"
    return(res)
}

#' Inference of CNV from single-cell sparse raw-count matices
#'
#' This function estimates the CNV score based on single-cell expression data
#'
#' @param expmat Sparse matrix of raw-counts
#' @param nullmat Optional sparse matrix with same number of rows as
#' \code{expmat} to be used as null model or index for the cells in expmat
#' considered a diploid reference
#' @param species Character string indicating the species, either human or mouse
#' @param k Integer indicating the number of genes per set
#' @param skip Interger indicating the displacement of the window for selecting
#' the k genes
#' @param min_geneset Integer indicating the minimum size for the genesets
#' @param center Logical, whether to subtract the mean expression of each gene
#' @param scale Logical, whether to divide by each gene standard deviation
#' @param batch_size Number indicating the maximum ize of the batch of cells
#' to be analyzed by aREA per iteration. This parameter is used to limit
#' memory load when transforming the sparse matrix normalized expression into
#' full matrix ranked signature
#' @param verbose Logical, whether progress should be reported
#'
#' @return Object of class inferCNV, which is a list containing matrix of nes,
#' and parameters (param), including species, window (k) and skip
#'
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' dset <- eh[["EH5419"]]
#' tpm_matrix <- SummarizedExperiment::assays(dset)$TPM
#' set.seed(1)
#' tpm_matrix <- tpm_matrix[, sample(ncol(tpm_matrix), 500)]
#' cnv <- inferCNV(tpm_matrix)
#' class(cnv)
#' names(cnv)
#' cnv$nes[1:5, 1:3]
#'
#' @export
inferCNVsc <- function(expmat, nullmat = NULL, species = c("human", "mouse"),
    k = 100, skip = 25, min_geneset = 10, center = TRUE, scale = FALSE,
    batch_size = 1000, verbose = TRUE) {
    # Check values for species
    species <- match.arg(species)
    # Validate input
    checkmate::assertClass(expmat, "Matrix")
    checkmate::assert(checkmate::testClass(nullmat, "Matrix", null.ok = TRUE),
        checkmate::testNumeric(nullmat, lower = 1, upper = ncol(expmat)))
    checkmate::assertInt(k, lower = 10, upper = 1000)
    checkmate::assertInt(skip, lower = 1, upper = k)
    checkmate::assertInt(min_geneset, lower = 2, upper = k)
    checkmate::assertFlag(center)
    checkmate::assertFlag(scale)
    checkmate::assertNumber(batch_size, lower = 1, upper = Inf)
    checkmate::assertLogical(verbose, len = 1)
    batch_size <- round(batch_size)
    # Compatibilize null model
    if (length(nullmat) > 0 && !is.null(nrow(nullmat))) {
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
    geneset <- geneset[vapply(geneset, length, numeric(1)) >= min_geneset]
    if (length(geneset) == 0) {
        stop(paste0("No geneset with at least ", min_geneset, " genes"))
    }
    genes <- unique(unlist(geneset, use.names = FALSE))
    # Normalizing the gene expression
    expmat <- cpm(expmat)
    expmat@x <- log2(expmat@x + 1)
    # Enrichment for the expmat
    if (verbose)
        message("Computing the enrichment for the genesets in the expression matrix")
    expmat <- expmat[rownames(expmat) %in% genes, , drop = FALSE]
    ref <- expmat
    if (length(nullmat) > 0) {
        if (!is.null(nrow(nullmat))) {
            nullmat <- cpm(nullmat)
            nullmat@x <- log2(nullmat@x + 1)
            genes <- intersect(rownames(expmat), rownames(nullmat))
            expmat <- expmat[match(genes, rownames(expmat)), , drop = FALSE]
            nullmat <- nullmat[match(genes, rownames(nullmat)), , drop = FALSE]
            ref <- nullmat
        } else {
            ref <- expmat[, round(nullmat), drop = FALSE]
        }
    }
    m1 <- 0
    sd1 <- 1
    if (center) {
        m1 <- Matrix::rowMeans(ref, na.rm = TRUE)
    }
    if (scale) {
        sd1 <- sqrt(rowVarsSM(ref))
    }
    expmat_nes <- sREAsm(expmat, geneset, m1 = m1, sd1 = sd1,
        batch_size = batch_size)
    # Enrichment of the null model
    nullnes <- NULL  # Initialize the null nes variable
    if (length(nullmat) > 0) {
        if (verbose) {
            message("Computing null model")
        }
        if (!is.null(nrow(nullmat))) {
            expmat_null <- sREAsm(nullmat, geneset, m1 = m1, sd1 = sd1,
                batch_size = batch_size)
        } else {
            expmat_null <- expmat_nes[, round(nullmat), drop = FALSE]
        }
        # Estimating NES
        if (verbose) {
            message("Estimating the normalized enrichment scores")
        }
        expmat_nes <- t(vapply(seq_len(nrow(expmat_nes)),
            computeNesForMatrixRow, numeric(ncol(expmat_nes)),
            nesmat = expmat_nes, nullmat = expmat_null))
        nullnes <- t(vapply(seq_len(nrow(expmat_null)), computeNesForMatrixRow,
            numeric(ncol(expmat_null)), nesmat = expmat_null,
            nullmat = expmat_null))
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
#' eh <- ExperimentHub::ExperimentHub()
#' dset <- eh[["EH5419"]]
#' tpm_matrix <- SummarizedExperiment::assays(dset)$TPM
#' set.seed(1)
#' tpm_matrix <- tpm_matrix[, sample(ncol(tpm_matrix), 500)]
#' cnv <- inferCNV(tpm_matrix)
#' cnv <- genomicInstabilityScore(cnv)
#' plot(density(cnv$gis))
#'
#' @seealso [inferCNV()] to infer the enrichment of loci-blocks in the gene
#' expression data.
#'
#' @export
genomicInstabilityScore <- function(cnv, method = c("var", "meansq"),
    likelihood = FALSE) {
    # Validating inputs
    checkmate::assertLogical(likelihood, len = 1)
    method <- match.arg(method)
    validateInferCNV(cnv, "nes")
    # Compute genomic instability score
    gisnull <- NULL  # Initialize gisnull
    switch(method,
    var = {
        gis <- log2(colVars(cnv[["nes"]]))
        if (!is.null(cnv[["null"]])) {
            gisnull <- log2(colVars(cnv[["null"]]))
        }
    },
    meansq = {
        gis <- log2(colMeans(cnv[["nes"]] ^ 2, na.rm = TRUE))
        if (!is.null(cnv[["null"]])) {
            gisnull <- log2(colMeans(cnv[["null"]] ^ 2, na.rm = TRUE))
        }
    })
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
#' as tumors, accepts first and last
#' @param normal Optional vector of integers indicating the Gaussians considered
#' as normal. This is only useful when no null model has been provided for the
#' analysis, accepts firt and last
#' @param adjust Method for multiple hypohesis test adjustment to be used when
#' null model is provided but there is no clear separation between distributions
#' @param pval Number indicating the p-value threshold
#'
#' @return Updated inferCNV-class object with gi_likelihood slot
#'
#' @examples
#'
#' eh <- ExperimentHub::ExperimentHub()
#' dset <- eh[["EH5419"]]
#' tpm_matrix <- SummarizedExperiment::assays(dset)$TPM
#' set.seed(1)
#' tpm_matrix <- tpm_matrix[, sample(ncol(tpm_matrix), 500)]
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
giLikelihood <- function(inferCNV, recompute = TRUE, distros = c(1, 10),
    tumor = NULL, normal = NULL, adjust = c("fdr", "none", "holm", "hochberg",
    "hommel", "bonferroni", "BH", "BY"), pval = .05) {
    # Validate inputs
    validateInferCNV(inferCNV, "nes")
    checkmate::assertLogical(recompute, len = 1)
    checkmate::assertIntegerish(distros, lower = 1, upper = 10, len = 2,
        any.missing = FALSE)
    checkmate::assert(
        checkmate::testIntegerish(tumor, lower = 1, upper = max(distros),
        any.missing = FALSE, null.ok = TRUE),
        checkmate::testString(tumor, pattern = "first|last"))
    checkmate::assert(
        checkmate::testIntegerish(normal, lower = 1, upper = max(distros),
        any.missing = FALSE, null.ok = TRUE),
        checkmate::testString(normal, pattern = "first|last"))
    checkmate::assertNumber(pval, lower = 0, upper = 1)
    adjust <- match.arg(adjust)
    # sort distros indexes
    distros <- sort(distros)
    # Compute GIS if not done previously
    if (is.null(inferCNV[["gis"]])) {
        inferCNV <- genomicInstabilityScore(inferCNV)
    }
    # Fit mixture gaussians to the results
    if (is.null(inferCNV[["gi_fit"]]) | recompute) {
        results_fit <- mixGaussianFit(inferCNV[["gis"]], min = distros[1],
            max = distros[2])
        inferCNV[["gi_fit"]] <- results_fit
    }
    results_fit <- inferCNV[["gi_fit"]]

    # If gisnull and no tumor and normal then return likelihood based on p-value
    if (length(inferCNV[["gisnull"]]) > 0 & length(tumor) == 0 & length(normal) == 0) {
        null_fit <- mixGaussianFit(inferCNV[["gisnull"]], min = distros[1],
            max = distros[2])
        if (max(null_fit$mu) >= max(results_fit$mu)) {
            cdf <- n1platform:::aecdf(inferCNV$gisnull)
            p <- cdf(inferCNV$gis, alternative = "greater")$p.value %>%
                p.adjust(method = adjust)
            inferCNV[["gi_likelihood"]] <- n1platform:::sigT(-log10(p),
                slope = 2, inflection = -log10(pval))
            return(inferCNV)
        }
        pos <- which(results_fit$mu > max(null_fit$mu))
        results_fit$mu <- c(null_fit$mu, results_fit$mu[pos])
        results_fit$sigma <- c(null_fit$sigma, results_fit$sigma[pos])
        results_fit$lambda <- c(null_fit$lambda, results_fit$lambda[pos])
        results_fit$lambda <- results_fit$lambda / sum(results_fit$lambda)
        normal <- seq_len(length(null_fit$mu))
        tumor <- seq_len(length(results_fit$mu))[-normal]
    }

    # Adjust tumor and normal
    if (length(tumor) > 0) {
        tumor[tummor == "fist"] <- 1
        tumor[tumor == "last"] <- length(results_fit[["mu"]])
        tumor <- round(as.numeric(tumor))
    }
    if (length(normal) > 0) {
        normal[normal == "first"] <- 1
        normal[normal == "last"] <- length(results_fit[["mu"]])
        normal <- round(as.numeric(normal))
    }

    # If no tumor and normal models are selected, asign the last gaussian
    if (length(tumor) == 0 & length(normal) == 0) {
        normal <- 1
        if (length(results_fit[["mu"]]) > 1) {
            normal <- seq_len(length(results_fit[["mu"]]) - 1)
        }
        tumor <- length(results_fit[["mu"]])
    }
    if (length(normal) == 0) {
        pos <- seq_len(length(results_fit[["mu"]]))
        normal <- pos[pos != tumor]
    }
    if (length(tumor) == 0) {
        pos <- seq_len(length(results_fit[["mu"]]))
        tumor <- pos[pos != normal]
    }
    if (any(tumor > length(results_fit[["mu"]]))) {
        stop("Tumor selection is larger than the number of models",
            call. = FALSE)
    }
    if (any(normal > length(results_fit[["mu"]])))
        stop("Normal selection is larger than the number of models",
            call. = FALSE)

    # Keep only selected models
    results_fit[["mu"]] <- results_fit[["mu"]][c(normal, tumor)]
    results_fit[["sigma"]] <- results_fit[["sigma"]][c(normal, tumor)]
    # Adjust normal and tumor indexes
    normal <- seq_len(length(normal))
    tumor <- (seq_len(length(tumor))) + length(normal)
    results_fit$normal <- normal
    results_fit$tumor <- tumor

    # Compute relative likelihood
    inferCNV[["gi_fit"]] <- results_fit
    inferCNV[["gi_likelihood"]] <- rowSums(predict(results_fit,
        inferCNV[["gis"]], tumor))
    return(inferCNV)
}
