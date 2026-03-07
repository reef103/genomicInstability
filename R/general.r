# General functions The functions below are all internal for the package

#' Genomic instability functions for debugging
#' 
#' @param path String indicating the path to the genomicInstability package
#' 
#' @export
sourceGenomicInstability <- function(path = "~/code/genomicInstability") {
    checkmate::assertString(path)
    source(file.path(path, "R/analysis.r"))
    source(file.path(path, "R/enrichment.r"))
    source(file.path(path, "R/general.r"))
    source(file.path(path, "R/mixGaussianFit.r"))
    source(file.path(path, "R/plots.r"))
}


# Split a vector based on a window with and displacement This function split a
# vector into a list of vectors with a given window width and displacement
# @param x Vector to split @param k Integer indicating the window width
# @param skip Integer indicating the displacement @return List of vectors with
# width = k
splitVectorWindow <- function(x, k = 200, skip = 50) {
    # Set initial positions
    start_pos <- seq(1, length(x) - k + 1, by = skip)
    res <- lapply(start_pos, subsetVectorRange, x = x, k = k)
    names(res) <- seq_len(length(res))
    return(res)
}

# Range of elements from vector This function returns a subset of k elements of
# a vector starting at position start_pos @param start_pos Integer indicating
# the starting position for the subset @param x Vector @param k Integer
# indicating the number of elements to include in the subset @return Vector
subsetVectorRange <- function(start_pos, x, k = 200) {
    end <- start_pos + k - 1
    if (start_pos > length(x))
        stop("start greater than length(x)")
    x[start_pos:(min(end, length(x)))]
}

# Variance of rows for arrays with NA values This function computes the variance
# by rows ignoring NA values @param x Numeric matrix @return Vector with the
# variance by row results
rowVars <- function(x) {
    ave <- rowMeans(x, na.rm = TRUE)
    pos <- which(is.na(x))
    largo <- rowSums(!is.na(x))
    x[pos] <- rep(ave, ncol(x))[pos]
    res <- (x - ave)^2 %*% rep(1, ncol(x))/(largo - 1)
    return(res[, 1])
}

# Variance of rows for sparse matrices with NA values This function computes the variance
# by rows ignoring NA values @param x Numeric matrix @return Vector with the
# variance by row results
rowVarsSM <- function(x) {
    ave <- Matrix::rowMeans(x, na.rm = TRUE)
    pos <- which(is.na(x))
    largo <- Matrix::rowSums(!is.na(x))
    x[pos] <- rep(ave, ncol(x))[pos]
    res <- (x - ave)^2 %*% rep(1, ncol(x))/(largo - 1)
    return(res[, 1])
}


# Variance of columns for arrays with NA values This function computes the
# variance by columns ignoring NA values @param x Numeric matrix @return Vector
# with the variance by column results
colVars <- function(x) rowVars(t(x))


# Keep vector elements This function filter a vector keeping only selected
# elements @param x Vector @param elements vector of elements to keep @return
# Filtered vector
keepVectorElements <- function(x, elements) x[x %in% elements]

# Get elements from string This function select the selected positions from a
# delimited string and return them as a matrix @param x Vector of character
# strings @param sep Character string indicating the separation character/s
# @param pos Vector of integers indicating the positions of interest @return
# Vector or Matrix with same number of rows as elements in the input x vector
# and columns as positions selected
getElementsFromString <- function(x, sep = "-", pos = 1) {
    vapply(strsplit(x, sep), subsetVector, character(length(pos)), pos = pos)
}

# Get positions from vector This function returns a slected subset of a vector
# @param x Vector @param pos Integer indicating the positios to return @return
# Verctor
subsetVector <- function(x, pos = 1) x[pos]

# Subset Matrix by columns @param i Vector indicating the columns to select
# @param x Matrix @return Matrix with selected columns i
subsetMatrixByColumns <- function(i, x) x[, i, drop = FALSE]

# Integration with trapezoid method This function integrate over a numerical
# range using the trapezoid method @param x Numeric vector of x values @param y
# Numeric vector of y values @return Number
integrateTZ <- function(x, y) {
    pos <- order(x)
    x <- x[pos]
    y <- y[pos]
    idx = 2:length(x)
    return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}

# Validate genomicInstability objects This function asserts the validity of
# inferCNV-class objects @param x Instance of class inferCNV to check @param
# slots String vector with slots to check @return x, error is trigered if the
# test is not suscessful
validateInferCNV <- function(x, slots = "nes") {
    checkmate::assertClass(x, "inferCNV")
    if ("nes" %in% slots)
        checkmate::assertMatrix(x[["nes"]], mode = "numeric",
            all.missing = FALSE,  min.rows = 1, min.cols = 1,
            row.names = "named", col.names = "named")
}

#' Counts per million implementation for sparse matrices
#' 
#' @param x Numeric sparse matrix
#' 
#' @return sparse matrix of CPM
#' @export
cpm <- function(x) {
    checkmate::assert(testMatrix(x), testClass(x, "Matrix"))
    col_sum <- Matrix::colSums(x, na.rm = TRUE)
    f <- Matrix::sparseMatrix(i = seq_len(length(col_sum)),
        j = seq_len(length(col_sum)), x = 1e+06 / col_sum)
    res <- x %*% f
    colnames(res) <- colnames(x)
    return(res)
}

#' Maximum overlap between reference and test gaussians
#' 
#' Overlap between two normal distributions
#'
#' Compute the overlap integral \eqn{\int \min(f_1(x), f_2(x))\,dx} where
#' \eqn{f_1,f_2} are normal densities with given means and sds.
#'
#' @param mu1 Numeric mean for the first normal
#' @param mu2 Numeric mean for the second normal
#' @param sd1 Numeric standard deviation for the first normal
#' @param sd2 Numeric standard deviation for the second normal
#' @param tol Numeric tolerance for numeric checks
#' @return Numeric overlap in [0,1]
normalOverlap <- function(mu1, mu2, sd1 = 1, sd2 = 1, tol = 1e-12) {
    checkmate::assertNumeric(c(mu1, mu2), len = 2)
    checkmate::assertNumeric(c(sd1, sd2), lower = 0, any.missing = FALSE)
    if (sd1 <= 0 || sd2 <= 0)
        stop("sd must be positive")

    s1 <- sd1; s2 <- sd2
    A <- 1/(2 * s2^2) - 1/(2 * s1^2)
    B <- -mu2/(s2^2) + mu1/(s1^2)
    C <- mu2^2/(2 * s2^2) - mu1^2/(2 * s1^2) - log(s2/s1)
    roots <- numeric(0)

    if (abs(A) < tol) {
        if (abs(B) < tol) {
            roots <- numeric(0)
        } else {
            roots <- -C / B
        }
    } else {
        D <- B^2 - 4 * A * C
        if (D < -tol) {
            roots <- numeric(0)
        } else {
            D <- max(0, D)
            r1 <- (-B - sqrt(D)) / (2 * A)
            r2 <- (-B + sqrt(D)) / (2 * A)
            roots <- sort(c(r1, r2))
            if (length(roots) == 2 && abs(roots[1] - roots[2]) < tol)
                roots <- roots[1]
        }
    }

    if (length(roots) == 0) {
        f <- function(x) pmin(dnorm(x, mu1, s1), dnorm(x, mu2, s2))
        val <- stats::integrate(f, -Inf, Inf, rel.tol = 1e-8)$value
        return(as.double(val))
    }

    pts <- c(-Inf, as.numeric(roots), Inf)
    overlap <- 0
    for (i in seq_len(length(pts) - 1)) {
        a <- pts[i]; b <- pts[i + 1]
        mid <- if (is.finite(a) && is.finite(b)) (a + b)/2 else if (is.finite(a)) a + 1 else if (is.finite(b)) b - 1 else 0
        if (dnorm(mid, mu1, s1) <= dnorm(mid, mu2, s2)) {
            overlap <- overlap + (pnorm(b, mu1, s1) - pnorm(a, mu1, s1))
        } else {
            overlap <- overlap + (pnorm(b, mu2, s2) - pnorm(a, mu2, s2))
        }
    }
    as.double(overlap)
}

#' Maximum overlap between reference and test gaussians
#'
#' Accepts either a numeric vector `c(mu, sd)` or a list with elements `mu`/`mean` and `sd`.
#' @param null_fit First gaussian (numeric or list)
#' @param test_fit Second gaussian (numeric or list)
#' @return Numeric overlap in [0,1]
maxGaussianOverlap <- function(null_fit, test_fit) {
    # Expect lists with numeric vectors `mu` and `sigma`
    if (!is.list(null_fit) || !is.list(test_fit))
        stop("null_fit and test_fit must be lists with elements 'mu' and 'sigma'")

    # try common name variants
    get_vec <- function(f, name1, name2 = NULL) {
        if (!is.null(f[[name1]])) return(as.numeric(f[[name1]]))
        if (!is.null(name2) && !is.null(f[[name2]])) return(as.numeric(f[[name2]]))
        NULL
    }

    mu_null <- get_vec(null_fit, "mu", "mean")
    sd_null <- get_vec(null_fit, "sigma", "sd")
    mu_test <- get_vec(test_fit, "mu", "mean")
    sd_test <- get_vec(test_fit, "sigma", "sd")

    if (is.null(mu_null) || is.null(sd_null))
        stop("null_fit must contain 'mu' and 'sigma' numeric vectors")
    if (is.null(mu_test) || is.null(sd_test))
        stop("test_fit must contain 'mu' and 'sigma' numeric vectors")

    mu_null <- as.numeric(mu_null)
    sd_null <- as.numeric(sd_null)
    mu_test <- as.numeric(mu_test)
    sd_test <- as.numeric(sd_test)

    if (length(mu_null) != length(sd_null))
        stop("null_fit: 'mu' and 'sigma' must have the same length")
    if (length(mu_test) != length(sd_test))
        stop("test_fit: 'mu' and 'sigma' must have the same length")

    m <- length(mu_null)
    k <- length(mu_test)

    res <- numeric(k)
    for (j in seq_len(k)) {
        overlaps <- vapply(seq_len(m), function(i) {
            normalOverlap(mu_null[i], mu_test[j], sd_null[i], sd_test[j])
        }, numeric(1))
        res[j] <- max(overlaps)
    }
    res
}