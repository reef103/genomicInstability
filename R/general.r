# General functions
# The functions below are all internal for the package

# Split a vector based on a window with and displacement
# 
# This function split a vector into a list of vectors with a given window width and displacement
# 
# @param x Vector to split
# @param k Integer indicating the window width
# @param skip Integer indicating the displacement
# @return List of vectors with width = k
splitVectorWindow <- function(x, k=200, skip=50) {
    # Set initial positions
    start_pos <- seq(1, length(x)-k+1, by=skip)
    res <- lapply(start_pos, subsetVectorRange, x=x, k=k)
    names(res) <- 1:length(res)
    return(res)
}

# Range of elements from vector
# 
# This function returns a subset of k elements of a vector starting at position start_pos
# 
# @param start_pos Integer indicating the starting position for the subset
# @param x Vector
# @param k Integer indicating the number of elements to include in the subset
# @return Vector
subsetVectorRange <- function(start_pos, x, k=200) {
    end <- start_pos+k-1
    if (start_pos > length(x)) stop("start greater than length(x)")
    x[start_pos:(min(end, length(x)))]
}

# Variance of rows for arrays with NA values
#
# This function computes the variance by rows ignoring NA values
#
# @param x Numeric matrix
# @return Vector with the variance by row results
rowVars <- function(x) {
    ave <- rowMeans(x, na.rm=TRUE)
    pos <- which(is.na(x))
    largo <- rowSums(!is.na(x))
    x[pos] <- rep(ave, ncol(x))[pos]
    res <- (x-ave)^2 %*% rep(1,ncol(x))/(largo-1)
    return(res[, 1])
}

# Variance of columns for arrays with NA values
#
# This function computes the variance by columns ignoring NA values
#
# @param x Numeric matrix
# @return Vector with the variance by column results
colVars <- function(x) rowVars(t(x))


# Keep vector elements
# 
# This function filter a vector keeping only selected elements
# 
# @param x Vector
# @param elements vector of elements to keep
# @return Filtered vector
keepVectorElements <- function(x, elements) x[x %in% elements]

# Get elements from string
# 
# This function select the selected positions from a delimited string and return them as a matrix
# 
# @param x Vector of character strings
# @param sep Character string indicating the separation character/s
# @param pos Vector of integers indicating the positions of interest
# @return Vector or Matrix with same number of rows as elements in the input x vector and columns as positions selected
getElementsFromString <- function(x, sep="-", pos=1) {
    sapply(strsplit(x, sep), subsetVector, pos=pos)
}

# Get positions from vector
# 
# This function returns a slected subset of a vector
# 
# @param x Vector
# @param pos Integer indicating the positios to return
# @return Verctor
subsetVector <- function(x, pos=1) x[pos]

# Subset Matrix by columns
# 
# @param i Vector indicating the columns to select
# @param x Matrix
# @return Matrix with selected columns i
subsetMatrixByColumns <- function(i, x) x[, i, drop=FALSE]

# Integration with trapezoid method
# 
# This function integrate over a numerical range using the trapezoid method
# 
# @param x Numeric vector of x values
# @param y Numeric vector of y values
# @return Number
integrateTZ <- function(x, y) {
    pos <- order(x)
    x <- x[pos]
    y <- y[pos]
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

# Validate genomicInstability objects
# 
# This function asserts the validity of inferCNV-class objects
# 
# @param x Instance of class inferCNV to check
# @param slots String vector with slots to check
# @return x, error is trigered if the test is not suscessful
validateInferCNV <- function(x, slots="nes") {
    checkmate::assertClass(x, "inferCNV")
    if ("nes" %in% slots)
        checkmate::assertMatrix(x[["nes"]], mode="numeric", all.missing=FALSE, min.rows=1, min.cols=1, row.names="named", col.names="named")
}
