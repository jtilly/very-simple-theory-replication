#' Tauchen Discretization of an AR(1)
#'
#' Tauchen computes the transition probability matrix for a discretized
#' AR(1) with normally distributed innovations
#'
#' @param myGrid A vector with the equidistant support of the stochastic process
#' @param mu The mean of the underlying normal distribution
#' @param sigma The standard deviation of the underlying normal distribution
#' @return a square transition probability matrix with dimension length(myGrid)
#' times length(myGrid).
#' @examples
#' tauchen( c(1,2,3), 0, 100)
tauchen = function(myGrid, mu, sigma) {

    mu = mu / 100.0
    sigma = sigma / 100.0

    Pi = matrix(1, length(myGrid), length(myGrid))

    if (length(myGrid) == 1) {
        return(Pi)
    }

    d = myGrid[2] - myGrid[1]

    # interior transition <from, to>
    interiorTransition = function(to) {
        return(pnorm((to - myGrid + d / 2.0 - mu) / sigma) - pnorm((to - myGrid - d / 2.0 - mu) / sigma))
    }

    # apply all the interior transitions
    if (length(myGrid) > 2) {
        Pi[, -c(1, length(myGrid))] = sapply(myGrid[-c(1, length(myGrid))], interiorTransition)
    }

    # take care of the boundaries
    Pi[, 1] = pnorm((myGrid[1] - myGrid + d / 2.0 - mu) / sigma)
    Pi[,length(myGrid)] = 1 - pnorm((myGrid[length(myGrid)] - myGrid - d / 2.0 - mu) / sigma)

    if (any(abs(rowSums(Pi) - 1) > 1e-9)) {
        stop("row sums of the transition probability matrix don't sum to 1")
    }

    return(Pi)
}

#' Compute the ergodic distribution of a Markov Chain
#'
#' Computes the ergodic distribution of a discrete Markov chain
#'
#' @param transMat the transition probability matrix of the Markov chain
#' @return a vector with the ergodic distribution
ergodic = function(transMat) {
    eigenvector = Re(eigen(t(transMat))$vectors[, 1])
    return(eigenvector / sum(eigenvector))
}
