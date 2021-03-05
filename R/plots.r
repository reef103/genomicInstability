# Plot functions

#' Genomic instability plot
#'
#' This function plot the genomic instability distribution, gaussian fits and
#' null distribution if available
#'
#' @param inferCNV Object of class inferCNV
#' @param legend Character string indicating the location of the legend. none
#' to not include it
#' @param ... Additional parameters for plot()
#'
#' @return None, a figure is created in the default output device
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
#' giDensityPlot(cnv)
#'
#' @seealso [giLikelihood()] to estimate the relative likelihood,
#' [genomicInstabilityScore()] to estimate the genomic instability score for
#' each cell in the dataset, and [inferCNV()] to infer the enrichment
#' of loci-blocks in the gene expression data.
#' @export
giDensityPlot <- function(inferCNV, legend = c("topleft", "top", "topright",
    "none"),
    ...) {
    # Validate inputs
    validateInferCNV(inferCNV, "nes")
    legend <- match.arg(legend)
    # Compute what is missing
    if (is.null(inferCNV[["gi_fit"]]))
        inferCNV <- giLikelihood(inferCNV)
    results_fit <- inferCNV[["gi_fit"]]
    plot(results_fit, inferCNV[["gis"]], xlab = "Genomic instability score",
        ...)
    if (!is.null(inferCNV[["gisnull"]])) {
        dnull <- density(inferCNV[["gisnull"]], na.rm = TRUE, adj = 2)
        dnull$y <- dnull$y * max(density(inferCNV[["gis"]], na.rm = TRUE)$y) /
            max(dnull$y)
        lines(dnull, lwd = 2)
        if (legend != "none")
            legend(legend, c("Samples", "Null model"),
                fill = c("grey", "black"), bty = "n")
    }
}

#' Plot chromosome map
#'
#' This function generates a chromosomes map plot for the inferred CNVs
#'
#' @param x Object of class inferCNV
#' @param output Optional output PDF file name (with extension)
#' @param threshold Likelihood threshold for identifying genomically inestable
#' cells/samples, 0 disables this filter
#' @param gamma Number indicating the gamma transformation for the colors
#' @param ... Additional parameters for plot
#' @param resolution Integer indicating the ppi for the png and jpg output files
#'
#' @return Nothing, a plot is generated in the default output devise
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
#' \donttest{plot(cnv, output='test.png')}
#'
#' @seealso [giLikelihood()] to estimate the relative likelihood,
#' [genomicInstabilityScore()] to estimate the genomic instability score for
#' each cell in the dataset, and [inferCNV()] to infer the enrichment
#' of loci-blocks in the gene expression data.
#' @method plot inferCNV
#' @export
plot.inferCNV <- function(x, output = NULL, threshold = 0.2, gamma = 1.5,
    resolution = 150, ...) {
    inferCNV <- x
    # Validate inputs
    validateInferCNV(inferCNV, "nes")
    checkmate::assertCharacter(output, min.chars = 1, any.missing = FALSE,
        len = 1, null.ok = TRUE)
    checkmate::assertNumeric(threshold, lower = 0, upper = 1,
        any.missing = FALSE, len = 1)
    checkmate::assertNumeric(gamma, lower = 0, upper = 10, any.missing = FALSE,
        len = 1)
    checkmate::assertInt(resolution, lower = 10, upper = 10000)
    # If threshold is between 0 and 1 and gi_likelihood is missing, compute
    # gi_likelihood
    if (is.null(inferCNV[["gi_likelihood"]]) & !(threshold > 0 & threshold < 1))
        inferCNV <- giLikelihood(inferCNV)
    # define nes object with a single nes matrix
    nes <- list(` ` = t(inferCNV[["nes"]]))
    # If threshold is between 0 and 1, split the nes matrix
    if (threshold > 0 & threshold < 1) {
        threshold <- sort(c(threshold, 1 - threshold))
        nes <- list()
        nes[["Unstable"]] <- t(inferCNV[["nes"]])[inferCNV[["gi_likelihood"]] >=
            threshold[2], , drop = FALSE]
        nes[["Intermediate"]] <- t(inferCNV[["nes"]])[inferCNV[["gi_likelihood"]] > threshold[1]
            & inferCNV[["gi_likelihood"]] < threshold[2], , drop = FALSE]
        if (threshold[1] == threshold[2]) {
            nes[["Stable"]] <- t(inferCNV[["nes"]])[inferCNV[["gi_likelihood"]]
                < threshold[1], , drop = FALSE]
        } else {
            nes[["Stable"]] <- t(inferCNV[["nes"]])[inferCNV[["gi_likelihood"]]
                <= threshold[1], , drop = FALSE]
        }
        # Remove empty matrices
        nes <- nes[sapply(nes, nrow) > 0]
    }
    # Sort matrices by hierarchical cluster
    nes <- lapply(nes, sortRowsByHclust)
    # Split matrices by chromosome
    nes <- lapply(nes, splitByChromosome)
    # Prepare for the plot
    nesplot <- transformNesForPlot(nes, method = "scale")
    # Size of the chromosome blocks
    chrom_size <- sapply(nes[[1]], ncol)
    # Size of the groups
    group_size <- sapply(nes, function(x) nrow(x[[1]]))
    # Figure panels widths and heights
    widths <- chrom_size/100 + 0.02
    widths[1] <- widths[1] + 0.4 - 0.01
    widths[length(widths)] <- widths[length(widths)] + 0.1 - 0.01
    heights <- group_size/200 + 0.02
    heights[1] <- heights[1] + 0.4 - 0.01
    heights[length(heights)] <- heights[length(heights)] + 0.1 - 0.01
    # Generating the plot
    output_flag <- FALSE  # Flag to whether the graphic device should be closed
    if (!is.null(output)) {
        if (length(grep("\\.pdf", output)) > 0) {
            pdf(output, w = sum(widths), h = sum(heights), pointsize = 10,
                useD = FALSE)
            output_flag <- TRUE  # Flag to close the device
        } else if (length(grep("\\.png", output)) > 0) {
            png(output, w = sum(widths), h = sum(heights), units = "in",
                res = resolution)
            output_flag <- TRUE  # Flag to close the device
        } else if (length(grep("\\.jpg", output)) > 0) {
            jpeg(output, w = sum(widths), h = sum(heights), units = "in",
                res = resolution)
            output_flag <- TRUE  # Flag to close the device
        }
    }
    layout(matrix(seq_len(length(chrom_size) * length(group_size)),
        length(group_size), length(chrom_size), byrow = TRUE), widths = widths,
        heights = heights)
    for (group in seq_len(length(group_size))) {
        for (chrom in seq_len(length(chrom_size))) {
            # Defining the margins
            if (group == 1 & chrom == 1) {
                par(mai = c(0.01, 0.4, 0.4, 0.01))
            } else if (group == length(group_size) & chrom == 1) {
                par(mai = c(0.1, 0.4, 0.01, 0.01))
            } else if (chrom == 1) {
                par(mai = c(0.01, 0.4, 0.01, 0.01))
            } else if (group == 1 & chrom == length(chrom_size)) {
                par(mai = c(0.01, 0.01, 0.4, 0.1))
            } else if (group == length(group_size)
                & chrom == length(chrom_size)) {
                par(mai = c(0.1, 0.01, 0.01, 0.1))
            } else if (chrom == length(chrom_size)) {
                par(mai = c(0.01, 0.01, 0.01, 0.1))
            } else if (group == 1) {
                par(mai = c(0.01, 0.01, 0.4, 0.01))
            } else if (group == length(group_size)) {
                par(mai = c(0.1, 0.01, 0.01, 0.01))
            } else {
                par(mai = c(0.01, 0.01, 0.01, 0.01))
            }
            # Drawing the heatmaps
            plothm(nesplot[[group]][[chrom]], grid = FALSE, scmax = 1,
                gama = gamma, ...)
            # Ading the axis labels
            if (group == 1)
                axis(3, ncol(nes[[group]][[chrom]])/2, names(nes[[group]])
                    [chrom], tick = FALSE, las = 1, line = -0.5)
            if (chrom == 1)
                axis(2, nrow(nes[[group]][[chrom]])/2, names(nes)[group],
                    tick = FALSE, las = 3, line = -0.5)
        }
    }
    if (output_flag)
        dev <- dev.off()
}

# Sort the rows of a matrix based on hierarchical cluster analysis @param x
# Numeric matrix @param metric Character string indicating the distance metric
# @param method Character string indicating the method for the hierarchical
# cluster analysis @return Numerical matrix
sortRowsByHclust <- function(x, metric = c("euclidean", "maximum", "manhattan",
    "canberra", "binary", "minkowski"), method = c("complete", "ward.D",
    "ward.D2", "single", "complete", "average", "mcquitty", "median",
    "centroid")) {
    # Validate inputs
    checkmate::assertMatrix(x, mode = "numeric", any.missing = FALSE,
        min.rows = 2, min.cols = 2)
    metric <- match.arg(metric)
    method <- match.arg(method)
    # Distance
    dd <- dist(x, method = metric)
    pos <- hclust(dd, method = method)$order
    x[pos, , drop = FALSE]
}

# split matrix by chromosome @param x matrix of NES with cells or samples in
# rows and chromosome fragments in columns @return List of matrices
splitByChromosome <- function(x) {
    # Get the chromosomes from the names of the columns
    chr <- getElementsFromString(colnames(x), sep = "-", pos = 1)
    # Split matrix by chromosomes
    res <- tapply(seq_len(ncol(x)), chr, subsetMatrixByColumns, x = x)
    # Sort the chromosomes
    chr_name <- gsub("chr", "", names(res))
    # replace non-numeric names by consecutive numbers
    suppressWarnings(chr_numeric <- !is.na(as.numeric(chr_name)))
    pos <- which(!chr_numeric)
    chr_name[pos] <- seq_len(length(pos)) + length(which(chr_numeric))
    res <- res[order(as.numeric(chr_name))]
    res
}

# colorScale This function generates a color scale
#
# @param x Vector or matrix of numeric values
# @param color Vector of character strings indicating the colors for the scale.
# Up to three colors can be defined. While is used for the missing color
# @param gama Number indicating the gama transformation
# @param alpha Number between 0 and 1 indicating the transparency of the color
# (1 for absolute color)
# @param scmax Number indicating the maximum value for the scale
# @param nacol Character string indicating the color for missing values
# @return Vector of colors
colorScale <- function(x, color = c("royalblue", "firebrick2"), gama = 1,
    alpha = 1, scmax = 0, nacol = "grey80") {
    if (length(color) == 1)
        color <- c(color, "white", color)
    if (length(color) == 2)
        color <- c(color[1], "white", color[2])
    if (scmax == 0)
        scmax <- max(abs(x), na.rm = TRUE)
    pos <- which(abs(x) > scmax)
    if (length(pos) > 0)
        x[pos] <- scmax * sign(x[pos])
    x <- abs(x/scmax)^gama * sign(x)
    color <- t(col2rgb(color))
    col <- sapply(x, function(x, color) {
        colSums(color * c(abs(x) * (x < 0), 1 - abs(x), x * (x > 0)))
    }, color = color/255)
    pos <- which(colSums(is.na(col)) > 0)
    col[is.na(col)] <- 0
    col <- apply(col, 2, function(x, alpha) rgb(x[1], x[2], x[3],
        alpha = alpha), alpha = alpha)
    col[pos] <- nacol
    return(col)
}

# Plot heatmap This function produce a heatmap plot from a numerical matrix
#
# @param x Numerical matrix
# @param color Two character strings vector describing the colors for the
# heatmap
# @param gama Number, indicating the exponential transformation for the color
# scale
# @param cex Number indicating the magnification factor for the labels
# @param grid Logical, whether a grid should be ploted
# @param scale Number between 0 and .9 indicating the proportion of vertical
# space used to draw the color scale
# @param scmax Optional number indicating the maximum value to be allowed for
# the heatmap
# @param box Logical, whether to draw a box around the plot
# @param ... Additional parameters to pass to the plot function
#
# @return Nothing, a heatmap is produced in the default output device
plothm <- function(x, color = c("royalblue", "firebrick2"), gama = 1, cex = 1,
    grid = TRUE, scale = FALSE, scmax = 0, box = TRUE, ...) {
    if (scale > 0) {
        if (scale == 1)
            ff <- 6/(nrow(x) + 5) else ff <- scale
        pari <- par("mai")
        layout(matrix(seq_len(2), 2, 1), h = c(1 - ff, ff))
        if (round(sum(pari - c(1.02, 0.82, 0.82, 0.42)), 2) == 0)
            pari <- c(0.2, 0.2, 1.2, 1.2)
        par(mai = pari)
        plothm.(x, color = color, gama = gama, scmax = scmax, box = box, ...)
        axis(4, nrow(x):1, rownames(x), tick = FALSE, line = 0, las = 2,
            adj = 0, cex.axis = cex)
        axis(3, seq_len(ncol(x)), colnames(x), tick = FALSE, line = 0, las = 2,
            adj = 0, cex.axis = cex)
        ra <- seq(-1, 1, length = 100)
        coli <- colorScale(x = ra, color = color, gama = gama, scmax = scmax)
        par(mai = pari * c(0, 1, 0, 3) + c(0.5, 0, 0.1, 0))
        image(seq_len(length(ra)), 1, matrix(seq_len(length(ra)),
            length(ra), 1), col = coli, ylab = "", xlab = "", axes = FALSE)
        if (scmax == 0)
            scmax <- max(abs(x), na.rm = TRUE)
        axis(1, seq(1, length(ra), length = 5), round(seq(-scmax, scmax,
            length = 5), 1), cex.axis = cex)
    } else plothm.(x = x, color = color, gama = gama, grid = grid,
        scmax = scmax, box = box, ...)
}

# Plotting functions
plothm. <- function(x, color = c("royalblue", "firebrick2"), gama = 1,
    grid = TRUE, scmax = 0, box = TRUE, ...) {
    coli <- colorScale(x = x[nrow(x):1, , drop = FALSE], color = color,
        gama = gama, scmax = scmax)
    image(seq_len(ncol(x)), seq_len(nrow(x)), t(matrix(seq_len(ncol(x) *
        nrow(x)), nrow(x), ncol(x))), col = coli, ylab = "", xlab = "",
        axes = FALSE, ...)
    if (box)
        box()
    if (grid)
        grid(ncol(x), nrow(x), col = "black", lty = 1)
}

# Transform data for the heatmap plot @param x either matrix or list of matrices
# @param method Character string indicating the transformation method @return
# Transformed matrix or list of matrices
transformNesForPlot <- function(x, method = c("scale", "equalize")) {
    method <- match.arg(method)
    # Transformation function
    transformation <- nesTransformationFunction(as.vector(unlist(x,
        use.names = FALSE)), method = method)
    x <- applyNesTransformation(x, transformation)
    return(x)
}

# Transformation function for NES
#
# @param x Numeric vector
# @param method Character string indicating the transformation method
#
# @return Function that takes a vector or matrix of values to be transformed
# and returns a vector or matrix of transformed values
nesTransformationFunction <- function(x, method = c("scale", "equalize")) {
    # Check method values
    method <- match.arg(method)
    # Check x is a vector
    if (!is.vector(x))
        stop("x must be a vector")
    # switch between methods
    switch(method, scale = {
        xmax <- max(abs(x), na.rm = TRUE)
        results <- function(x) x/xmax
        return(results)
    }, equalize = {
        results <- function(x) x
        return(results)
    })
}

# Apply NES transformation @param x MAtrix or list of matrices @param
# transformation Function generated by nesTRansformationFunction() function
# @return Transformed matrix or list of matrices
applyNesTransformation <- function(x, transformation) {
    if (is.list(x)) {
        x <- lapply(x, applyNesTransformation, transformation)
    } else {
        x <- transformation(x)
    }
    return(x)
}
