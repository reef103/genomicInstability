# Functions for enrichment analysis All the functions below are internal for the
# package

sREA <- function(signatures, groups) {
    # Check signatures is a matrix
    if (is.null(nrow(signatures)))
        signatures <- matrix(signatures, length(signatures), 1,
            dimnames = list(names(signatures), "sample1"))
    # Copula transformation of the signature matrix by colulmns
    sig <- qnorm(apply(signatures, 2, rank)/(nrow(signatures) + 1))
    gr <- vapply(groups, function(x, samp) {
        samp %in% x
    }, logical(nrow(sig)), samp = rownames(sig))
    gr <- t(gr)
    nn <- rowSums(gr)
    gr <- gr/nn
    es <- gr %*% sig
    return(es * sqrt(nn))
}

aecdf1 <- function(dnull, symmetric = FALSE, x, alternative = c("two.sided",
    "greater", "less")) {
    dnull <- dnull[is.finite(dnull)]
    if (symmetric) {
        iqr <- quantile(abs(dnull), c(0.5, 1 - 5/length(dnull)))
        pd <- ecdf(abs(dnull))
        a <- list(x = abs(dnull), y = pd(abs(dnull)))
        a1 <- list(x = a$x[length(a$x) - (tl2:iq2) + 1] - iqr[3],
            y = log(1 - epd(iqr[3])) - log(1 - a$y[length(a$x) -
                (tl2:iq2) + 1]))
        a1 <- lapply(a1, function(x, pos) x[pos], pos = which(is.finite(a1$y)))
        if (length(a1$x) < 3)
            stop("Not enough permutations to compute NULL distribution",
                call. = FALSE)
        a1 <- list(x = a$x[length(a$x) - (15:4)] - iqr[2],
            y = log(1 - pd(iqr[2])) - log(1 - a$y[length(a$x) - (15:4)]))
        a1 <- lapply(a1, function(x, pos) x[pos], pos = which(is.finite(a1$y)))
        if (length(a1$x) < 3)
            stop("Not enough permutations to compute NULL distribution",
                call. = FALSE)
        fit <- lm(y ~ 0 + x, data = a1)
        alternative <- match.arg(alternative)
        x1 <- abs(x)
        p <- exp(log(1 - pd(iqr[2])) - predict(fit, list(x = x1 - iqr[2])))
        p <- p * (x1 > iqr[2]) + (1 - pd(x1)) * (x1 <= iqr[2])
        nes <- qnorm(p/2, lower.tail = FALSE) * sign(x)
        switch(alternative, two.sided = {
            p <- p
        }, greater = {
            p <- p/2
            p[x < 0] <- 1 - p[x < 0]
        }, less = {
            p <- p/2
            p[x > 0] <- 1 - p[x > 0]
        })
        names(nes) <- names(p) <- names(x)
        return(list(nes = nes, p.value = p))
    }
    iqr <- quantile(dnull, c(5/length(dnull), 0.5, 1 - 5/length(dnull)))
    pd <- ecdf(dnull)
    a <- list(x = dnull, y = pd(dnull))
    a1 <- list(x = a$x[5:14] - iqr[1], y = log(pd(iqr[1])) - log(a$y[5:14]))
    a1 <- lapply(a1, function(x, pos) x[pos], pos = which(is.finite(a1$y)))
    if (length(a1$x) < 3)
        stop("Not enough permutations to compute NULL distribution",
            call. = FALSE)
    fit1 <- lm(y ~ 0 + x, data = a1)
    a1 <- list(x = a$x[length(a$x) - (15:4)] - iqr[3],
               y = log(1 - pd(iqr[3])) - log(1 - a$y[length(a$x) - (15:4)]))
    a1 <- lapply(a1, function(x, pos) x[pos], pos = which(is.finite(a1$y)))
    if (length(a1$x) < 3)
        stop("Not enough permutations to compute NULL distribution",
            call. = FALSE)
    fit2 <- lm(y ~ 0 + x, data = a1)
    alternative <- match.arg(alternative)
    p1 <- exp(log(pd(iqr[1])) - predict(fit1, list(x = x - iqr[1])))
    p2 <- exp(log(1 - pd(iqr[3])) - predict(fit2, list(x = x - iqr[3])))
    p <- p1 * (x < iqr[1]) + p2 * (x > iqr[3]) + pd(x) * (x >= iqr[1] &
        x < iqr[2]) + (1 - pd(x)) * (x >= iqr[2] & x <= iqr[3])
    nes <- qnorm(p, lower.tail = FALSE) * sign(x - iqr[2])
    switch(alternative, two.sided = {
        p <- p * 2
    }, greater = {
        p[x < iqr[2]] <- 1 - p[x < iqr[2]]
    }, less = {
        p[x >= iqr[2]] <- 1 - p[x >= iqr[2]]
    })
    names(nes) <- names(p) <- names(x)
    return(list(nes = nes, p.value = p))
}

computeNesForMatrixRow <- function(i, nesmat, nullmat) {
    aecdf1(nullmat[i, ], x = nesmat[i, ])$nes
}
