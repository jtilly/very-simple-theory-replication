#' Compute mixing probabilities
#'
#' The function computes the mixing probability for a given state (N,C,W).
#'
#' @param N is an integer with the number of firms post-entry
#' @param C is an integer with the index of the demand grid
#' @param W is the a scalar with the current realization of the cost shock
#' @param vS is the post-survival value function
#' @return a scalar with the mixing probability
#'
mixingProbabilities = function(N, C, W, vS) {

  nE = N
  matCoef = matrix(0, nrow = N, ncol = N)
  for (jX in seq(from = nE - 1, to = 0)) {
    signCoef = (-1)^(jX - seq(from = 0, to = jX))
    nCk = factorial(nE - 1) / factorial(nE - 1 - jX) / (factorial(seq(from = 0, to = jX)) * factorial(jX - seq(from = 0, to = jX)))
    continuationValue = (-exp(W) + vS[seq(from = 1, to = jX + 1),C])
    matCoef[nE-jX,seq(from = 1, to = jX + 1)] = signCoef * nCk * continuationValue
  }

  vecCoef=rowSums(matCoef)
  # note that R's polyroot and Matlab's root take the opposite ordering of inputs
  mixprobcand = polyroot(vecCoef[nE:1])
  mixprobcand = Re(mixprobcand[abs(Im(mixprobcand))< 1e-10])
  aS = mixprobcand[mixprobcand >= 0 & mixprobcand <= 1]

  if (length(aS) != 1) {
    stop('Survival strategy had no or multiple solutions.')
  }

  return(aS)
}



#' Data Generation
#'
#' This function generates a data set and returns a list with matrices \code{N} and \code{C}
#'
#' @param Settings is an R list of settings
#' @param Param is an R list of of params
#' @param rCheck is an integer with the number of markets for which to generate data
#' @param tCheck is an integer with the number of time periods for which to generate data
#' @param tBurn burn in period
#' @param initC initial C
#' @param initN initial N
#' @return an R list with the matrices \code{N} and \code{C}
simulateData = function(Settings, Param, rCheck, tCheck, tBurn = 25, initC = NULL, initN = NULL) {

  W = matrix(rnorm(rCheck*(tBurn + tCheck), -0.5 * Param$omega ^ 2, Param$omega), nrow = rCheck, ncol = (tBurn + tCheck))

  Data = list(
    C = matrix(0, nrow = rCheck, ncol = (tBurn + tCheck)),
    N = matrix(0, nrow = rCheck, ncol = (tBurn + tCheck)),
    W = W)

  # compute the Eq at these parameters
  Eq  = valuefunction(Settings, Param)

  if (is.null(initC)) {
      Data$C[,1] = sample(1:Settings$cCheck, size = rCheck, replace = TRUE)
  } else {
      Data$C[,1] = initC
  }
  if (is.null(initN)) {
      Data$N[,1] = sample(1:Settings$nCheck, size = rCheck, prob = 1.0 / Settings$nCheck + vector(mode = "numeric", length = Settings$nCheck), replace = TRUE)
  } else {
      Data$N[,1] = initN
  }

  # loop over markets
  for(rX in seq(from = 1, to = rCheck)) {

    # loop over time
    for(tX in seq(from = 2, to = tBurn + tCheck)) {

      Data$C[rX, tX] = sample(1:Settings$cCheck, size = 1, prob = as.numeric(Param$Pi[Data$C[rX, tX - 1],]))


      if (Data$N[rX, tX - 1]>0 && (Eq$vS[1,Data$C[rX, tX]] <= exp(W[rX, tX]))) {
        # all firms leave
        Data$N[rX, tX] = 0
      } else if (Data$N[rX, tX - 1]>1 && (Eq$vS[pmax(1,Data$N[rX, tX - 1]),Data$C[rX, tX]] <= exp(W[rX, tX]))) {
        # some firms leave
        aS = mixingProbabilities(Data$N[rX, tX - 1], Data$C[rX, tX], W[rX, tX], Eq$vS)
        Data$N[rX, tX] = rbinom(1, Data$N[rX, tX - 1], aS)
      } else if (Data$N[rX, tX - 1]<Settings$nCheck && (Eq$vS[pmax(1,Data$N[rX, tX - 1]),Data$C[rX, tX]] > exp(W[rX, tX]))) {
        # all incumbents stay, there may be entry
        Data$N[rX, tX] = Data$N[rX, tX - 1] +
            sum((Eq$vS[seq(from = Data$N[rX, tX - 1] + 1, to = Settings$nCheck), Data$C[rX, tX]]
                - (1 + Param$phi[seq(from = Data$N[rX, tX - 1] + 1, to = Settings$nCheck)]) * exp(W[rX, tX])) > 0)
      } else {
        # all incumbents stay, no entry
        Data$N[rX, tX] = Data$N[rX, tX - 1]
      }


    }

  }

  if (tBurn > 0) {
    Data$C = Data$C[,-tBurn:-1]
    Data$N = Data$N[,-tBurn:-1]
    Data$W = Data$W[,-tBurn:-1]
  }

  return(Data)
}
