#' R / C++ programs for Abbring, Campbell, Tilly Yang (2017) Very Simple Markov
#' Perfect Industry Dynamics: Empirics
#'
#' @author Jan Tilly
#' @seealso  A Matlab sandbox is
#'   available at \url{https://github.com/jtilly/acty}. An extensive
#'   documentation of the Matlab sandbox is available at
#'   \url{https://jtilly.github.io/acty}.
#' @docType package
#' @name actyR-package
NULL

#' Population data for muSAs from 2000 to 2009
#'
#' This data set contains population data for muSAs from 2000 to 2009 in
#' panel data format
#'
#' @docType data
#' @name population
#' @format A data frame with 573 rows and 12 variables
#' \describe{
#'   \item{id}{CBSA id}
#'   \item{name}{CBSA name}
#'   \item{pop2000}{population in 2000}
#'   \item{pop2001}{population in 2001}
#'   \item{pop2002}{population in 2002}
#'   \item{pop2003}{population in 2003}
#'   \item{pop2004}{population in 2004}
#'   \item{pop2005}{population in 2005}
#'   \item{pop2006}{population in 2006}
#'   \item{pop2007}{population in 2007}
#'   \item{pop2008}{population in 2008}
#'   \item{pop2009}{population in 2009}
#' }
#' @source \url{https://www.census.gov/popest/data/metro/totals/2009/tables/CBSA-EST2009-01.csv}
NULL

#' Firm count data for muSAs from 2000 to 2009
#'
#' This data set contains firm counts for Movie theaters for muSAs from 2000 to
#' 2009 in panel data format
#'
#' @docType data
#' @name firms
#' @format A data frame with 573 rows and 12 variables
#' \describe{
#'   \item{id}{CBSA id}
#'   \item{name}{CBSA name}
#'   \item{firms2000}{firms in 2000}
#'   \item{firms2001}{firms in 2001}
#'   \item{firms2002}{firms in 2002}
#'   \item{firms2003}{firms in 2003}
#'   \item{firms2004}{firms in 2004}
#'   \item{firms2005}{firms in 2005}
#'   \item{firms2006}{firms in 2006}
#'   \item{firms2007}{firms in 2007}
#'   \item{firms2008}{firms in 2008}
#'   \item{firms2009}{firms in 2009}
#' }
#' @source \url{http://www.census.gov/econ/cbp/}
NULL

#' Covariates for muSAs in the year 2000
#'
#' This data set contains covariates for the muSAs in the year 2000
#'
#' @docType data
#' @name cov
#' @format A data frame with 573 rows and 12 variables
#' \describe{
#'   \item{id}{CBSA id}
#'   \item{name}{CBSA name}
#'   \item{inc}{median income in 2000}
#'   \item{dummy1}{Middle Atlantic Census Region Dummy}
#'   \item{dummy2}{East North Central Census Region Dummy}
#'   \item{dummy3}{West North Central Census Region Dummy}
#'   \item{dummy4}{South Atlantic Census Region Dummy}
#'   \item{dummy5}{East South Central Census Region Dummy}
#'   \item{dummy6}{West South Central Census Region Dummy}
#'   \item{dummy7}{Mountain Census Region Dummy}
#'   \item{dummy8}{Pacific Census Region Dummy}
#'   \item{dummy9}{Low Diversity Dummy}
#' }
#' @source \url{http://factfinder.census.gov/faces/nav/jsf/pages/index.xhtml}
NULL


#' First Step Likelihood Function
#'
#' This function computes the likelihood used in the first step of
#' the estimation procedure.
#'
#' @param mu_sigma is a vector whose elements are \code{mu} and \code{sigma}
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state
#' @param Settings is an R list with program Settings, most notably cCheck
#' @param Pi is the transition probabillity matrix if \code{Pi} is not
#' provided the transition probability matrix is computed based on \code{x}
#' times length(myGrid).
#' @return Returns the scalar valued negative log-likelihood
nfxpLikelihoodStep1 = function(mu_sigma, Data, Settings, Pi = NULL) {

    if (is.null(Pi)) {
        Pi = tauchen(Settings$logGrid, mu_sigma[1], mu_sigma[2])
    }

    # use linear indexing to pick the right values out of Pi
    from = as.vector(Data$C[, 1:ncol(Data$C) - 1])
    to = as.vector(Data$C[,2:ncol(Data$C)])

    if (any(is.na(Pi[Settings$cCheck * (to - 1) + from])) | any(Pi[Settings$cCheck * (to - 1) + from] <= 0, na.rm = TRUE)) {
        return(Inf)
    } else {
        return(-sum(log(Pi[Settings$cCheck * (to - 1) + from])))
    }
}

#' First Step Gradient Function
#'
#' This function computes the gradient used in the first step of
#' the estimation procedure.
#'
#' @param mu_sigma is a vector whose elements are \code{mu} and \code{sigma}
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state
#' @param Settings is an R list with program Settings, most notably cCheck
#' @param Pi is the transition probabillity matrix if \code{Pi} is not
#' provided the transition probability matrix is computed based on \code{mu_sigma}
#' times length(myGrid).
#' @return Returns the vector valued gradient of the negative log-likelihood
nfxpGradientStep1 = function(mu_sigma, Data, Settings, Pi = NULL) {

    if (is.null(Pi)) {
        Pi = tauchen(Settings$logGrid, mu_sigma[1], mu_sigma[2])
    }
    dPi.dMu = dPidMu(Pi, Settings$logGrid, mu_sigma[1], mu_sigma[2])
    dPi.dSigma  = dPidSigma(Pi, Settings$logGrid, mu_sigma[1], mu_sigma[2])
    from = as.vector(Data$C[, 1:ncol(Data$C) - 1])
    to = as.vector(Data$C[,2:ncol(Data$C)])
    likelihoodContributions = Pi[Settings$cCheck * (to - 1) + from]
    gradientContributions = cbind(dPi.dMu[Settings$cCheck * (to - 1) + from],
                                  dPi.dSigma[Settings$cCheck * (to - 1) + from]) / 100.0

    return(-colSums(gradientContributions/rep(likelihoodContributions, 2)))
}

#' Second Step Likelihood Function
#'
#' This function computes the second step likelihood function
#'
#' @param x is a vector of parameters
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state and a matrix \code{N} with number of active firms
#' @param Param is an R list with parameters
#' @param Settings is an R list with settings
#' @return Returns a vector with likelihood contributions
nfxpLikelihoodStep2 = function(x, Data, Settings, Param) {

    sX = 2 * Settings$nCheck + 1
    lX = length(x)

    x = matrix(x, nrow = lX, ncol = 1)

    if (length((Param$specificationMatrix[1:sX, 1:lX] %*% x)) != (2 * Settings$nCheck + 1)) {
        stop("Dimensions of specificationMatrix and x don't match!")
    }

    Param$k = (Param$specificationMatrix[1:sX, 1:lX] %*% x)[seq(from = 1, to = Settings$nCheck)]
    Param$phi = (Param$specificationMatrix[1:sX, 1:lX] %*% x)[seq(from=Settings$nCheck + 1, to = 2 * Settings$nCheck)]
    Param$omega = (Param$specificationMatrix[1:sX, 1:lX] %*% x)[2 * Settings$nCheck + 1]
    ll = likelihood(Settings, Param, Data)
    if (any(is.na(ll)) | any(ll <= 0, na.rm = TRUE)) {
        return(Inf)
    } else {
        return(-sum(log(ll)))
    }
}



#' Second Step Gradient Function
#'
#' This function computes the gradient of the second step likelihood function
#'
#' @param x is a vector of parameters
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state and a matrix \code{N} with number of active firms
#' @param Settings is an R list with settings
#' @param Param is an R list with parameters
#' @return Returns a vector of gradients
nfxpGradientStep2 = function(x, Data, Settings, Param) {

    sX = 2 * Settings$nCheck + 1
    lX = length(x)

    x = matrix(x, nrow = lX, ncol = 1)

    if (length((Param$specificationMatrix[1:sX, 1:lX] %*% x)) != (2 * Settings$nCheck + 1)) {
        stop("Dimensions of specificationMatrix and x don't match!")
    }

    Param$k = (Param$specificationMatrix[1:sX, 1:lX] %*% x)[seq(from = 1, to = Settings$nCheck)]
    Param$phi = (Param$specificationMatrix[1:sX, 1:lX] %*% x)[seq(from=Settings$nCheck + 1, to = 2 * Settings$nCheck)]
    Param$omega = (Param$specificationMatrix[1:sX, 1:lX] %*% x)[2 * Settings$nCheck + 1]

    GradientList = likelihoodGradient(Settings, Param, Data)
    return(- colSums(GradientList$gradientContributions[, 1:sX] %*% Param$specificationMatrix[1:sX, 1:lX] / rep(GradientList$likelihoodContributions, lX)))

}

#' Third Step Likelihood Function
#'
#' This function computes the third step likelihood function
#'
#' @param x is a vector of parameters
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state and a matrix \code{N} with number of active firms
#' @param Settings is an R list with settings
#' @param Param is an R list with parameters
#' @param Pi is the transition probabillity matrix if \code{Pi} is not
#' provided the transition probability matrix is computed based on \code{x}
#' times length(myGrid).
#' @return Returns a vector with likelihood contributions
nfxpLikelihoodStep3 = function(x, Data, Settings, Param, Pi = NULL) {

    Param$mu = x[length(x) - 1]
    Param$sigma = x[length(x)]

    if (is.null(Pi)) {
        Param$Pi = tauchen(Settings$logGrid, Param$mu, Param$sigma)
    }  else {
        Param$Pi = Pi
    }

    return(nfxpLikelihoodStep1(c(Param$mu, Param$sigma), Data, Settings, Param$Pi) +
            nfxpLikelihoodStep2(x[-c(length(x) - 1, length(x))], Data, Settings, Param))
}

#' Second Step Gradient Function
#'
#' This function computes the gradient of the second step likelihood function
#'
#' @param x is a vector of parameters
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state and a matrix \code{N} with number of active firms
#' @param Settings is an R list with settings
#' @param Param is an R list with parameters
#' @param Pi is the transition probabillity matrix if \code{Pi} is not
#' provided the transition probability matrix is computed based on \code{x}
#' times length(myGrid).
#' @return Returns a vector of gradients
nfxpGradientStep3 = function(x, Data, Settings, Param, Pi = NULL) {

    lX = length(x)
    x = matrix(x, nrow = lX, ncol = 1)

    if (length((Param$specificationMatrix %*% x)) != (2 * Settings$nCheck + 3)) {
        stop("Dimensions of specificationMatrix and x don't match!")
    }

    Param$k =     (Param$specificationMatrix %*% x)[seq(from = 1, to = Settings$nCheck)]
    Param$phi =   (Param$specificationMatrix %*% x)[seq(from=Settings$nCheck + 1, to = 2 * Settings$nCheck)]
    Param$omega = (Param$specificationMatrix %*% x)[2 * Settings$nCheck + 1]
    Param$mu =    (Param$specificationMatrix %*% x)[2 * Settings$nCheck + 2]
    Param$sigma = (Param$specificationMatrix %*% x)[2 * Settings$nCheck + 3]

    if (is.null(Pi)) {
        Param$Pi = tauchen(Settings$logGrid, Param$mu, Param$sigma)
    }  else {
        Param$Pi = Pi
    }

    # second step
    GradientList = likelihoodGradient(Settings, Param, Data)
    gradient = - colSums(GradientList$gradientContributions %*% Param$specificationMatrix /
                             rep(GradientList$likelihoodContributions, lX))


    gradient[c(lX - 1, lX)] = gradient[c(lX - 1, lX)]/100 +
        nfxpGradientStep1(c(Param$mu, Param$sigma) , Data, Settings, Param$Pi)

    return(gradient)

}

#' Third Step Variance-Covariance Matrix
#'
#' This function computes the variance-covariance matrix that corresponds
#' to the third-step likelihood. The function uses the outer-product-of-the-gradient
#' estimator of the hessian.
#'
#' @param x is a vector of parameters
#' @param Data is an R list that contains a matrix \code{C} with indices
#' of the demand state and a matrix \code{N} with number of active firms
#' @param Settings is an R list with settings
#' @param Param is an R list with parameters
#' @return Returns a vector with likelihood contributions
nfxpCovarianceStep3 = function(x, Data, Settings, Param) {

    lX = length(x)
    x = matrix(x, nrow = lX, ncol = 1)

    if (length((Param$specificationMatrix %*% x)) != (2 * Settings$nCheck + 3)) {
        stop("Dimensions of specificationMatrix and x don't match!")
    }

    Param$k = (Param$specificationMatrix %*% x)[seq(from = 1, to = Settings$nCheck)]
    Param$phi = (Param$specificationMatrix %*% x)[seq(from=Settings$nCheck + 1, to = 2 * Settings$nCheck)]
    Param$omega = (Param$specificationMatrix %*% x)[2 * Settings$nCheck + 1]
    Param$mu = (Param$specificationMatrix %*% x)[2 * Settings$nCheck + 2]
    Param$sigma = (Param$specificationMatrix %*% x)[2 * Settings$nCheck + 3]

    # likelhood and gradient contributions step 1
    Param$Pi = tauchen(Settings$logGrid, Param$mu, Param$sigma)
    dPi.dMu = dPidMu(Param$Pi, Settings$logGrid, Param$mu, Param$sigma)
    dPi.dSigma  = dPidSigma(Param$Pi, Settings$logGrid, Param$mu, Param$sigma)
    from = as.vector(Data$C[, 1:ncol(Data$C) - 1])
    to = as.vector(Data$C[,2:ncol(Data$C)])
    Step1.likelihoodContributions = Param$Pi[Settings$cCheck * (to - 1) + from]
    Step1.gradientContributions = cbind(dPi.dMu[Settings$cCheck * (to - 1) + from],
                                        dPi.dSigma[Settings$cCheck * (to - 1) + from])

    # gradient contributions step 2
    GradientList = likelihoodGradient(Settings, Param, Data)
    Step2.likelihoodContributions = GradientList$likelihoodContributions
    Step2.gradientContributions = GradientList$gradientContributions %*% Param$specificationMatrix
    Step2.gradientContributions[c(lX - 1, lX)] = Step2.gradientContributions[c(lX - 1, lX)]

    # add step 1 and step 2
    Step3.scores = Step2.gradientContributions / rep(Step2.likelihoodContributions, lX)
    Step3.scores[,c(lX - 1, lX)] = (Step3.scores[,c(lX - 1, lX)] + Step1.gradientContributions / rep(Step1.likelihoodContributions, 2))/100

    # create hessian using the outer product of the gradient estimator
    opg = t(Step3.scores) %*% Step3.scores
    return(solve(opg))
}
