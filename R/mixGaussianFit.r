# Functionsto fit mixtures of Gaussians

#' Find the modes for a multimodal distribution
#' 
#' This function returns the modes of a multimodal distribution
#' 
#' @param x Numeric vector
#' @param adj Number indicating the adjust parameter for bandwidth of the density estimation
#' @param thr Threshold for lambda, distrivutions with lambda below this threshold will be discarded
#' @return list of three elements: mean, sd and lambda
getPeaks3 <- function(x, adj=1.2, thr=1e-2) {
    den <- density(x, adj=adj, na.rm=TRUE, n=512*50)
    sp <- smooth.spline(den)
    sp2 <- predict(sp, den$x, deriv=2)$y
    pos <- which(sp2[-length(sp2)]*sp2[-1]<0)
    x2 <- (den$x[pos]+den$x[pos+1])/2
    sp3 <- predict(sp, x2, deriv=3)$y
    posi <- sapply(which(sp3<0), function(i, x2, posi) {
        if (length(x2)<(i+1)) return(NULL)
        pos <- which(posi>x2[i] & posi<x2[i+1])
        if (length(pos)==0) return((x2[i+1]+x2[i])/2)
        return(posi[pos[1]])
    }, x2=x2, posi=getPeaks(x, adj=adj))
    posi <- unlist(posi[sapply(posi, length)>0], use.names=FALSE)
    posi1 <- sapply(which(sp3>0), function(i, x2) {
        if (length(x2)<(i+1)) return(NULL)
        (x2[i]+x2[i+1])/2
    }, x2=x2)
    posi1 <- unlist(posi1[sapply(posi1, length)>0], use.names=FALSE)
    posi <- sapply(posi, function(posi, x) which(x>posi)[1], x=den$x)
    posi1 <- sapply(posi1, function(posi, x) which(x>posi)[1], x=den$x)
    m <- sort(den$x[posi])
    if (length(posi1)==0) b <- c(min(den$x)-diff(range(den$x)),  max(den$x)+diff(range(den$x)))
    else b <- c(min(den$x)-diff(range(den$x)), sort(den$x[posi1]), max(den$x)+diff(range(den$x)))
    sigma <- sapply(m, function(x, den, b) {
        ymax <- approx(den, xout=x)$y
        x1 <- den$x[den$y<(ymax/2)]
        dd <- x1-x
        dd1 <- b-x
        x2 <- c(x1[dd<0][which.max(dd[dd<0])], x1[dd>0][which.min(dd[dd>0])])
        x3 <- c(b[dd1<0][which.max(dd1[dd1<0])], b[dd1>0][which.min(dd1[dd1>0])])
        dd2 <- x2-x
        dd3 <- x3-x
        opt <- c(dd2[1]-dd3[1], dd3[2]-dd2[2])
        pos <- which(opt>0)
        if (length(pos)>0) tmp <- dd2[pos][which.min(abs(dd2[pos]))]
        else tmp <- dd3[which.min((dd2-dd3)^2)]
        res <- -3.934e-6+.8493*abs(tmp)
        if (length(res)==0) res <- NA
        res
    }, den=den, b=b)
    lambda <- sapply(1:length(m), function(i, m, sigma, den) {
        if (is.na(sigma[i])) return(NA)
        x2 <- seq(m[i]-sigma[i], m[i]+sigma[i], length=100)
        integrateTZ(x2, approx(den, xout=x2)$y)
    }, m=m, sigma=sigma, den=den)
    lambda[is.na(lambda)] <- 0
    lambda <- lambda/sum(lambda)
    pos <- which(lambda>thr)
    list(m=m[pos], sigma=sigma[pos], lambda=lambda[pos])
}

#' Find the modes for a multimodal distribution
#' 
#' This function returns the modes of a multimodal distribution
#' 
#' @param x Numeric vector
#' @param adj Number indicating the adjust parameter for bandwidth of the density estimation
#' @return Numeric vector of modes
getPeaks <- function(x, adj=1.5) {
    den <- density(x, adj=adj, na.rm=TRUE, n=512*50)
    posi <- which(c(FALSE, diff(diff(den$y)>0)<0, FALSE))
    posi <- posi[order(den$y[posi], decreasing=TRUE)]
    sort(den$x[posi])
}

#' Fit a mixture of gaussian curves
#' 
#' This function fits a mixture of gaussians to the distribution of any population
#' 
#' @param x Numeric vectors
#' @param thr Minimum lambda (proportion of the distribution) to include in the analysis
#' @param min Optional minimum number of distributions
#' @param max Optional maximum number of distributions
#' @param adj Optional vector of 2 components defining the range of values for the density adj parameter
#' @return Object of class mgfit. List of fitted parameters
mixGaussianFit <- function(x, thr=1e-2, min=1, max=1e6, adj=c(.8, 2)) {
    x <- x[is.finite(x)]
    fit <- lapply(seq(min(adj), max(adj), length=10), function(adj, x, thr, den, min, max) {
        param <- getPeaks3(x, adj=adj, thr=thr)
        param <- lapply(param, function(x) c(x, 0))
        if ((length(param$m)-1)<min | (length(param$m)-1)>max) return(NULL)
        while(any(param$lambda<thr)) {
            pos <- which(param$lambda<thr)
            param <- lapply(param, function(x, pos) x[-pos], pos=pos)
            suppressMessages(fit <- normalmixEM(x, lambda=param$lambda, mu=param$m, sd=param$sigma, epsilon=1e-50))
            ye <- sapply(1:length(fit$lambda), function(i, x, fit) {
                dnorm(x, fit$mu[i], fit$sigma[i])*fit$lambda[i]
            }, x=den$x, fit=fit)
            if (is.null(dim(ye))) ye <- matrix(ye, length(ye), 1)
            mse <- mean((rowSums(ye)-den$y)^2)
            param <- list(mu=fit$mu, sigma=fit$sigma, lambda=fit$lambda)
        }
        list(mu=fit$mu, sigma=fit$sigma, lambda=fit$lambda, loglik=fit$loglik, mse=mse)
    }, x=x, thr=thr, den=density(x, n=512*20), min=min, max=max)
    fit <- fit[sapply(fit, length)>0]
    if (length(fit)==0) stop("No sucessful filt with provided parameters", call.=FALSE)
    mse <- sapply(fit, function(x) x$mse)
    mse[!is.finite(mse)] <- 1e6
    res <- fit[[which.min(mse)]]
    class(res) <- "mgfit"
    return(res)
}

#' Predict relative likelihood for mixGaussianFit
#' 
#' This function computes the relative likelihood based on parameters fitted with mixGaussianFit function
#' 
#' @param fit Object generated by mixGaussianFit function
#' @param x Numerical vector of observations
#' @param k Optional vector of integers indicating the distributions to use
#' @return Matrix of relative likelihood, observations in rowa and distributions in columns
predict.mgfit <- function(fit, x, k=NULL) {
    # fit should have at least 2 gaussians
    if (length(fit$mu)<2) stop("Mixture of at least 2 gaussians is required", call.=FALSE)
    # If k is not specified then fit for all gaussians
    if (is.null(k)) k <- 1:length(fit$mu)
    # Ensure 1 <= k <= gaussians
    k <- k[k %in% (1:length(fit$mu))]
    # Check there are still usable k
    if (length(k)==0) stop("Selected gaussians not present in the fit object", call.=FALSE)
    prob_density <- sapply(1:length(fit$mu), function(i, x, fit) dnorm(x, fit$mu[i], fit$sigma[i]), x=x, fit=fit)
    prob_density <- prob_density/rowSums(prob_density)
    colnames(prob_density) <- 1:ncol(prob_density)
    return(prob_density[, k, drop=FALSE])
}

#' Plot mixGaissuinFit objects
#' 
#' Plot a mixture of distributions fitted with mixGaussianFit
#' 
#' @param fit Object generated by the function mixGaussianFit
#' @param x Vector of values to plot the distribution
#' @param col Color for the background distribution
#' @param lwd Line weight for the fitted distributions
#' @param fitCol Optional vector of colors for the fitted distributions
#' @param main Character string for the main title
#' @param ... Additional parameters to pass to the plot function
#' @return Nothing, a plot is generated
plot.mgfit <- function(fit, x, col="grey75", lwd=2, fitCol=NULL, main="", ylim=NULL, ...) {
    if (is.null(fitCol)) fitCol <- rainbow(length(fit$mu), .8, .8, end=.8)
    den <- density(x, from=min(x, na.rm=TRUE), to=max(x, na.rm=TRUE), na.rm=TRUE)
    denmax <- 0
    for (i in 1:length(fit$mu)) denmax <- max(denmax, dnorm(den$x, fit$mu[i], fit$sigma[i])*fit$lambda[i])
    denmax <- max(denmax, den$y)
    if (is.null(ylim))
        ylim <- c(0, denmax)
    plot(den, type="n", axes=FALSE, main=main, ylim=ylim, ...)
    axis(1)
    axis(2)
    polygon(c(min(den$x), den$x, max(den$x)), c(0, den$y, 0), col=col, border=NA)
    for (i in 1:length(fit$mu)) lines(den$x, dnorm(den$x, fit$mu[i], fit$sigma[i])*fit$lambda[i], lwd=lwd, col=fitCol[i])
}
