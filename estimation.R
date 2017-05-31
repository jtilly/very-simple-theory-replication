# load package.
library(actyR)

# make sure we are in the correct working directory
if (!all(c("graphs", "results", "tables") %in% list.files())) {
    stop("Please navigate to the correct working directory
         using `setwd`. Your working directory
         should include the folders graphs, results, and tables.")
}

# how many cores?
cl = makePSOCKcluster(parallel::detectCores(), outfile = file.path(tempdir(), "cluster_estimation.log"))
registerDoParallel(cl)

# Load data from package
data(cov)
data(firms)
data(population)

# Get rid of the first two columns
Data = list("C" = as.matrix(population[, -c(1, 2)]),
            "N" = as.matrix(firms[, -c(1, 2)]),
            "X" = as.matrix(cov[, -c(1, 2)]))

# demean income
Data$X[,1] = log(as.numeric(Data$X[, 1])) - log(mean(as.numeric(Data$X[, 1])))

# start the clock
time_begin = Sys.time()

# Make the Gauss-Legendre integration nodes and weights
GaussLegendre = legendre.quadrature.rules(32)[32][[1]]

# Set Settings- some of these are data dependent
Settings = list(tolInner = 1e-10,
                tolOuter = 1e-8,
                maxIterInner = 1000,
                maxIterOuter = 1000,
                nCheck = 9,
                cCheck = 200,
                lowBndC = log(min(Data$C) / 10000 / 1.25),
                uppBndC = log(max(Data$C) / 10000 * 1.25),
                logGrid = NA,
                d = NA,
                integrationNodes = (GaussLegendre$x + 1) / 2,
                integrationWeights = GaussLegendre$w / 2)

# Discretize the population data
Settings$logGrid = matrix(seq(from = Settings$lowBndC,
                              to = Settings$uppBndC,
                              length.out = Settings$cCheck),
                          nrow = 1,
                          ncol = Settings$cCheck)
Settings$d = Settings$logGrid[2] - Settings$logGrid[1]

# replace Data$C by its index
for (iX in seq(1, nrow(Data$C))) {
    for (tX in seq(1, ncol(Data$C))) {
        Data$C[iX, tX] = which.min(abs(Data$C[iX, tX] / 10000 - exp(Settings$logGrid)))
    }
}

# Initialize Param
Param = list(rho = 1 / 1.05, k = NA, phi = NA, omega = NA, Pi = NA, mu = NA, sigma = NA, kappa = 1)

## First step estimation
# Get initial values from log innovations
innovations = as.matrix(log(population[, 4:ncol(population)]) - log(population[, 3:(ncol(population) - 1)])) * 100
Startvalues.step1 = c(mean(innovations), 2 * sd(innovations))

# Do the optimization
Step1 = optim(Startvalues.step1,
              nfxpLikelihoodStep1,
              gr = nfxpGradientStep1,
              method = "BFGS",
              control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
              Settings = Settings,
              Data = Data)

# get the new transition probability matrix
Param$mu = Step1$par[1]
Param$sigma = Step1$par[2]
Param$Pi = tauchen(Settings$logGrid, Param$mu, Param$sigma)

# create specification matrix
specificationMatrix = matrix(0 , nrow = 21, ncol = 8)
specificationMatrix[1:4,1:4] = diag(rep(1, 4))
specificationMatrix[5:9,4] = rep(1, 5)
specificationMatrix[10:18,5] = rep(1, 9)
specificationMatrix[19,6] = 1
specificationMatrix[20,7] = 1
specificationMatrix[21,8] = 1

if ("specificationMatrix" %in% names(Param)) {
    Param$specificationMatrix = specificationMatrix
} else {
    Param = c(Param, list(specificationMatrix=specificationMatrix))
}

# Set starting values to 1 for all parameters
Startvalues.Step2 = vector(mode="numeric", length = 6) + 2
nfxpLikelihoodStep2(Startvalues.Step2, Data, Settings, Param)
nfxpGradientStep2(Startvalues.Step2, Data, Settings, Param)

# Do the estimation
Step2 = optim(Startvalues.Step2,
              nfxpLikelihoodStep2,
              gr = nfxpGradientStep2,
              method = "BFGS",
              control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
              Data = Data,
              Settings = Settings,
              Param = Param)

## Third step estimation without covariates
# Take the estimates from step 1 and step 2 as start values
Startvalues.step3 = c(Step2$par, Step1$par)

# Do the estimation
Step3 = optim(Startvalues.step3,
              nfxpLikelihoodStep3,
              gr = nfxpGradientStep3,
              method = "BFGS",
              control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
              Data = Data,
              Settings = Settings,
              Param = Param)

## Second step estimation with covariates but no diversity
likelihood.Step2.Covariates = function(x, Data, Settings, Param) {

    Param$beta = x[7:15]

    if (x[6] < 0) {
        return(Inf)
    }

    getMarketLikelihood = function (mX) {
        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))
        x[1:4] = x[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))
        return(nfxpLikelihoodStep2(x[1:6], DataRow, Settings, Param))
    }

    ll = foreach(mX = 1:NROW(Data$X), .combine = "+", .packages = "actyR", .export = "x") %dopar%  {
        getMarketLikelihood(mX)
    }

    return(ll)
}

## Second step estimation with covariates but no diversity
gradient.Step2.Covariates = function(x, Data, Settings, Param) {

    Param$beta = x[7:15]

    if (x[6] < 0) {
        return(rep(Inf, length = length(x)))
    }

    getMarketGradient = function (mX) {

        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        x[1:4] = x[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))

        gradVector = nfxpGradientStep2(x[1:6], DataRow, Settings, Param)
        gradVector[1:4] = gradVector[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))

        return(c(gradVector[1:6], sum(x[1:4] * gradVector[1:4]) * DataRow$X[1:9]))
    }

    gradientsList = foreach(sX = 1:NROW(Data$X), .packages = "actyR", .export = "x") %dopar%  {
        getMarketGradient(sX)
    }

    grad = vector(mode = "numeric", length = length(x))
    for (jX in 1:NROW(Data$X)) {
        grad = grad + gradientsList[[jX]]
    }

    return(grad)
}


# Take the startvalues from the estimates without covariates
Startvalues.Step2 = c(Step2$par, rep(0.0,9))

# Do the optimization
tic()
Step2.Covariates = optim(Startvalues.Step2,
                         likelihood.Step2.Covariates,
                         gradient.Step2.Covariates,
                         method = "BFGS",
                         control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
                         Data = Data, Settings = Settings, Param = Param)
toc()

## Third step estimation with covariates but no diversity
# Define the likelihood function
likelihood.Step3.Covariates = function(x, Data, Settings, Param) {
    Param$Pi = tauchen(Settings$logGrid, x[length(x) - 1], x[length(x)])
    ll = nfxpLikelihoodStep1(c(x[length(x) - 1], x[length(x)]), Data, Settings, Param$Pi) +
        likelihood.Step2.Covariates(x[1:15], Data, Settings, Param)
    return(ll)
}

# Third step gradient
gradient.Step3.Covariates = function(x, Data, Settings, Param) {

    Param$beta = x[7:15]
    Param$mu = x[16]
    Param$sigma = x[17]

    if (x[6] < 0) {
        return(rep(Inf, length = length(x)))
    }

    getMarketGradient = function (mX) {

        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        k = x[1:4]
        x[1:4] = x[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))

        gradVector = nfxpGradientStep3(c(x[1:6],Param$mu, Param$sigma), DataRow, Settings, Param)
        gradVector[1:4] = gradVector[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))

        return(c(gradVector[1:6],
                 sum(k * gradVector[1:4]) * DataRow$X[1:9],
                 gradVector[7:8]))
    }

    gradientsList = foreach(sX = 1:NROW(Data$X), .packages = "actyR", .export = "x") %dopar%  {
        getMarketGradient(sX)
    }

    grad = vector(mode = "numeric", length = length(x))

    for (jX in 1:NROW(Data$X)) {
        grad = grad + gradientsList[[jX]]
    }

    return(grad)
}

# Take the start values from steps 1 and 2
Startvalues.Step3 = c(Step2.Covariates$par, Step1$par)

# Do the optimization
tic()
Step3.Covariates = optim(Startvalues.Step3,
                         likelihood.Step3.Covariates,
                         gr = gradient.Step3.Covariates,
                         method = "BFGS",
                         control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
                         Settings = Settings,
                         Data = Data,
                         Param = Param)
toc()

# Get the covariance from the outer-product-of the gradient estimator
covariance.Step3.Covariates = function(x) {

    lX = length(x)
    Param$phi = rep(x[5], 9)
    Param$omega = x[6]
    Param$beta = x[7:15]
    Param$mu = x[16]
    Param$sigma = x[17]

    Param$Pi = tauchen(Settings$logGrid, x[16], x[17])

    # Step 1
    dPi.dMu = dPidMu(Param$Pi, Settings$logGrid, Param$mu, Param$sigma)
    dPi.dSigma  = dPidSigma(Param$Pi, Settings$logGrid, Param$mu, Param$sigma)
    from = as.vector(Data$C[,1:ncol(Data$C) - 1])
    to = as.vector(Data$C[,2:ncol(Data$C)])
    Step1.likelihoodContributions = Param$Pi[Settings$cCheck * (to - 1) + from]
    Step1.gradientContributions = cbind(dPi.dMu[Settings$cCheck * (to - 1) + from],
                                        dPi.dSigma[Settings$cCheck * (to - 1) + from])

    # Step2
    getMarketGradientContribution = function (mX) {
        k = c(x[1:4], x[4], x[4], x[4], x[4], x[4])
        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        lX = length(x)

        Param$k = k * exp(sum(DataRow$X[1:9] * Param$beta))

        GradientList = likelihoodGradient(Settings, Param, DataRow)
        Step2.likelihoodContributions = GradientList$likelihoodContributions
        Step2.gradientContributions = GradientList$gradientContributions %*% Param$specificationMatrix

        grad = cbind(
            Step2.gradientContributions[,1:4] * exp(sum(DataRow$X[1:9] * Param$beta)),
            matrix(Step2.gradientContributions[,5], nrow = 9),
            matrix(Step2.gradientContributions[,6], nrow = 9),
            matrix(rowSums(Step2.gradientContributions[,1:4] * matrix(Param$k[1:4], nrow = 9, ncol = 4, byrow = TRUE)) %x% DataRow$X[1:9], nrow = 9, ncol = 9, byrow = TRUE),
            Step2.gradientContributions[,7:8] / 100
        )

        return(list("grad" = grad, "ll" = GradientList$likelihoodContributions))
    }

    gradientsList = foreach(sX = 1:NROW(Data$X),
                            .packages = "actyR",
                            .export = c("Data", "Settings", "Param")) %dopar%  {
                                getMarketGradientContribution(sX)
                            }

    Step2.gradientContributions = NULL
    Step2.likelihoodContributions = NULL
    for (jX in 1:NROW(Data$X)) {
        Step2.gradientContributions = rbind(Step2.gradientContributions, gradientsList[[jX]]$grad)
        Step2.likelihoodContributions = rbind(Step2.likelihoodContributions, gradientsList[[jX]]$ll)
    }

    # add step 1 and step 2
    Step3.scores = Step2.gradientContributions / rep(Step2.likelihoodContributions, lX)
    Step3.scores[,c(lX - 1, lX)] = Step3.scores[,c(lX - 1, lX)] + Step1.gradientContributions / rep(Step1.likelihoodContributions, 2)

    # create hessian using the outer product of the gradient estimator
    opg = t(Step3.scores) %*% Step3.scores

    covariance = solve(opg)

    return(covariance)

}


Step3.Covariates = c(Step3.Covariates, list(covariance = covariance.Step3.Covariates(Step3.Covariates$par), se = NULL))
Step3.Covariates$se = sqrt(diag(Step3.Covariates$covariance))

## Second step estimation with covariates and diversity
likelihood.Step2.Covariates.Diversity = function(x, Data, Settings, Param) {

    Param$beta = x[11:19]

    getMarketLikelihood = function(mX) {
        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        if (Data$X[mX,10] == 0) {
            k = x[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))
        } else {
            k = x[5:8] * exp(sum(DataRow$X[1:9] * Param$beta))
        }

        return(nfxpLikelihoodStep2(c(k, x[9] , x[10]), DataRow, Settings, Param))
    }

    ll = foreach(sX = 1:NROW(Data$X), .combine = "+", .packages = "actyR") %dopar%  {
        getMarketLikelihood(sX)
    }

    return(min(1e25, ll))
}


## Second step estimation with covariates but no diversity
gradient.Step2.Covariates.Diversity = function(x, Data, Settings, Param) {

    Param$beta = x[11:19]

    getMarketGradient = function (mX) {

        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        if (Data$X[mX,10] == 0) {
            kTilde = x[1:4]
        } else {
            kTilde = x[5:8]
        }
        k = kTilde * exp(sum(DataRow$X[1:9] * Param$beta))

        gradVector = nfxpGradientStep2(c(k, x[9] , x[10]), DataRow, Settings, Param)
        gradVector[1:4] = gradVector[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))

        if (Data$X[mX,10] == 0) {
            return(c(gradVector[1:4],
                     rep(0,4),
                     gradVector[5:6],
                     sum(kTilde * gradVector[1:4]) * DataRow$X[1:9]
            ))
        } else {
            return(c(rep(0,4),
                     gradVector[1:6],
                     sum(kTilde * gradVector[1:4]) * DataRow$X[1:9]
            ))
        }

    }

    gradientsList = foreach(sX = 1:NROW(Data$X), .packages = "actyR") %dopar%  {
        getMarketGradient(sX)
    }

    grad = vector(mode="numeric", length = length(x))
    for(jX in 1:NROW(Data$X)) {
        grad = grad + gradientsList[[jX]]
    }

    return(grad)
}


# Take the startvalues from the estimates without covariates
Startvalues.Step2 = c(Step2.Covariates$par[1:4], Step2.Covariates$par[1:4], Step2.Covariates$par[5:15])

# Do the optimization
tic()
Step2.Covariates.Diversity = optim(Startvalues.Step2,
                                   likelihood.Step2.Covariates.Diversity,
                                   gradient.Step2.Covariates.Diversity,
                                   method = "BFGS",
                                   control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
                                   Data = Data, Settings = Settings, Param = Param)
toc()


## Third step estimation with covariates but no diversity
# Define the likelihood function
likelihood.Step3.Covariates.Diversity = function(x, Data, Settings, Param) {
    Param$Pi = tauchen(Settings$logGrid, x[length(x) - 1], x[length(x)])
    ll = nfxpLikelihoodStep1(c(x[length(x) - 1], x[length(x)]), Data, Settings, Param$Pi) +
        likelihood.Step2.Covariates.Diversity(x[1:19], Data, Settings, Param)
    return(ll)
}

# Third step gradient
gradient.Step3.Covariates.Diversity = function(x, Data, Settings, Param) {

    Param$beta = x[11:19]
    Param$mu = x[20]
    Param$sigma = x[21]

    getMarketGradient = function (mX) {

        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        if (Data$X[mX,10] == 0) {
            kTilde = x[1:4]
        } else {
            kTilde = x[5:8]
        }
        k = kTilde * exp(sum(DataRow$X[1:9] * Param$beta))

        gradVector = nfxpGradientStep3(c(k, x[9:10], Param$mu, Param$sigma), DataRow, Settings, Param)
        gradVector[1:4] = gradVector[1:4] * exp(sum(DataRow$X[1:9] * Param$beta))

        if (Data$X[mX,10] == 0) {

            return(c(gradVector[1:4],
                     rep(0,4),
                     gradVector[5:6],
                     sum(kTilde * gradVector[1:4]) * DataRow$X[1:9],
                     gradVector[7:8]))

        } else {

            return(c(rep(0,4),
                     gradVector[1:6],
                     sum(kTilde * gradVector[1:4]) * DataRow$X[1:9],
                     gradVector[7:8]))

        }
    }

    gradientsList = foreach(sX = 1:NROW(Data$X), .packages = "actyR") %dopar%  {
        getMarketGradient(sX)
    }

    grad = vector(mode = "numeric", length = length(x))

    for (jX in 1:NROW(Data$X)) {
        grad = grad + gradientsList[[jX]]
    }

    return(grad)
}

# Take the start values from steps 1 and 2
Startvalues.Step3 = c(Step2.Covariates.Diversity$par, Step1$par)

# Do the optimization
tic()
Step3.Covariates.Diversity = optim(Startvalues.Step3,
                                   likelihood.Step3.Covariates.Diversity,
                                   gr = gradient.Step3.Covariates.Diversity,
                                   method = "BFGS",
                                   control = list(REPORT = 1, abstol = Settings$tolOuter, reltol = Settings$tolOuter, maxit = Settings$maxIterOuter, trace = 1),
                                   Settings = Settings,
                                   Data = Data,
                                   Param = Param)
toc()

# Get the covariance from the outer-product-of the gradient estimator
covariance.Step3.Covariates.Diversity = function(x) {

    lX = length(x)
    Param$phi = rep(x[9], 9)
    Param$omega = x[10]
    Param$beta = x[11:19]
    Param$mu = x[20]
    Param$sigma = x[21]

    Param$Pi = tauchen(Settings$logGrid, x[20], x[21])

    # Step 1
    dPi.dMu = dPidMu(Param$Pi, Settings$logGrid, Param$mu, Param$sigma)
    dPi.dSigma = dPidSigma(Param$Pi, Settings$logGrid, Param$mu, Param$sigma)
    from = as.vector(Data$C[,1:ncol(Data$C) - 1])
    to = as.vector(Data$C[,2:ncol(Data$C)])
    Step1.likelihoodContributions = Param$Pi[Settings$cCheck * (to - 1) + from]
    Step1.gradientContributions = cbind(dPi.dMu[Settings$cCheck * (to - 1) + from],
                                        dPi.dSigma[Settings$cCheck * (to - 1) + from])

    # Step2
    getMarketGradientContribution = function (mX) {

        DataRow = list(N = t(as.matrix(Data$N[mX,])),
                       C = t(as.matrix(Data$C[mX,])),
                       X = t(as.matrix(Data$X[mX,])))

        if (Data$X[mX,10] == 0) {
            k = c(x[1:4], x[4], x[4], x[4], x[4], x[4])
        } else {
            k = c(x[5:8], x[8], x[8], x[8], x[8], x[8])
        }

        lX = length(x)

        Param$k = k * exp(sum(DataRow$X[1:9] * Param$beta))

        GradientList = likelihoodGradient(Settings, Param, DataRow)
        Step2.likelihoodContributions = GradientList$likelihoodContributions
        Step2.gradientContributions = GradientList$gradientContributions %*% Param$specificationMatrix

        if (Data$X[mX,10] == 0) {
            grad = cbind(
                Step2.gradientContributions[,1:4] * exp(sum(DataRow$X[1:9] * Param$beta)),
                matrix(0,ncol = 4, nrow = 9),
                matrix(Step2.gradientContributions[,5], nrow = 9),
                matrix(Step2.gradientContributions[,6], nrow = 9),
                matrix(rowSums(Step2.gradientContributions[,1:4] * matrix(Param$k[1:4], nrow = 9, ncol = 4, byrow = TRUE)) %x% DataRow$X[1:9], nrow = 9, ncol = 9, byrow = TRUE),
                Step2.gradientContributions[,7:8] / 100
            )
        } else {
            grad = cbind(
                matrix(0,ncol = 4, nrow = 9),
                Step2.gradientContributions[,1:4] * exp(sum(DataRow$X[1:9] * Param$beta)),
                matrix(Step2.gradientContributions[,5], nrow = 9),
                matrix(Step2.gradientContributions[,6], nrow = 9),
                matrix(rowSums(Step2.gradientContributions[,1:4] * matrix(Param$k[1:4], nrow = 9, ncol = 4, byrow = TRUE)) %x% DataRow$X[1:9], nrow = 9, ncol = 9, byrow = TRUE),
                Step2.gradientContributions[,7:8] / 100
            )
        }

        return(list("grad" = grad, "ll" = GradientList$likelihoodContributions))
    }

    gradientsList = foreach(sX = 1:NROW(Data$X),
                            .packages = "actyR",
                            .export = c("Data", "Settings", "Param")) %dopar%  {
                                getMarketGradientContribution(sX)
                            }

    Step2.gradientContributions = NULL
    Step2.likelihoodContributions = NULL

    for (jX in 1:NROW(Data$X)) {
        Step2.gradientContributions = rbind(Step2.gradientContributions, gradientsList[[jX]]$grad)
        Step2.likelihoodContributions = rbind(Step2.likelihoodContributions, gradientsList[[jX]]$ll)
    }

    # add step 1 and step 2
    Step3.scores = Step2.gradientContributions / rep(Step2.likelihoodContributions, lX)
    Step3.scores[,c(lX - 1,lX)] = Step3.scores[,c(lX - 1,lX)] + Step1.gradientContributions / rep(Step1.likelihoodContributions, 2)

    # create hessian using the outer product of the gradient estimator
    opg = t(Step3.scores) %*% Step3.scores

    covariance = solve(opg)

    return(covariance)

}

Step3.Covariates.Diversity = c(Step3.Covariates.Diversity, list(covariance = covariance.Step3.Covariates.Diversity(Step3.Covariates.Diversity$par), se = NULL))
Step3.Covariates.Diversity$se = sqrt(diag(Step3.Covariates.Diversity$covariance))

# Make Latex Tables with Final Estimates
empirical.results.table =
    c("\\begin{tabular}{l c c c c}  ",
      "& \\multicolumn{4}{c}{ Geographic Preference Diversity } \\\\",
      "& All $\\mu$SAs  & Diversity $>$ 13.4 miles & & Diversity $\\leq$ 13.4 miles \\\\",
      sprintf("$k(1) \\times 10^5$ & %.2f & %.2f & & %.2f \\\\", Step3.Covariates$par[1], Step3.Covariates.Diversity$par[1], Step3.Covariates.Diversity$par[5]),
      sprintf("                    & (%.2f) & (%.2f) & & (%.2f) \\\\", Step3.Covariates$se[1],  Step3.Covariates.Diversity$se[1],  Step3.Covariates.Diversity$se[5]),
      sprintf("$k(2) \\times 10^5$ & %.2f & %.2f & & %.2f \\\\", Step3.Covariates$par[2], Step3.Covariates.Diversity$par[2], Step3.Covariates.Diversity$par[6]),
      sprintf("                    & (%.2f) & (%.2f) & & (%.2f) \\\\", Step3.Covariates$se[2],  Step3.Covariates.Diversity$se[2],  Step3.Covariates.Diversity$se[6]),
      sprintf("$k(3) \\times 10^5$ & %.2f & %.2f & & %.2f \\\\", Step3.Covariates$par[3], Step3.Covariates.Diversity$par[3], Step3.Covariates.Diversity$par[7]),
      sprintf("                    & (%.2f) & (%.2f) & & (%.2f) \\\\", Step3.Covariates$se[3],  Step3.Covariates.Diversity$se[3],  Step3.Covariates.Diversity$se[7]),
      sprintf("$k(4) \\times 10^5$ & %.2f & %.2f & & %.2f \\\\", Step3.Covariates$par[4], Step3.Covariates.Diversity$par[4], Step3.Covariates.Diversity$par[8]),
      sprintf("                    & (%.2f) & (%.2f) & & (%.2f) \\\\", Step3.Covariates$se[4],  Step3.Covariates.Diversity$se[4],  Step3.Covariates.Diversity$se[8]),
      sprintf("$\\varphi$ & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[5], Step3.Covariates.Diversity$par[9]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[5], Step3.Covariates.Diversity$se[9]),
      sprintf("Median Income & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[6], Step3.Covariates.Diversity$par[10]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[6], Step3.Covariates.Diversity$se[10]),
      sprintf("Mid Atlantic & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[7], Step3.Covariates.Diversity$par[11]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[7], Step3.Covariates.Diversity$se[11]),
      sprintf("East North Central & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[8], Step3.Covariates.Diversity$par[12]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[8], Step3.Covariates.Diversity$se[12]),
      sprintf("West North Central & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[9], Step3.Covariates.Diversity$par[13]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[9], Step3.Covariates.Diversity$se[13]),
      sprintf("South Atlantic & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[10], Step3.Covariates.Diversity$par[14]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[10], Step3.Covariates.Diversity$se[14]),
      sprintf("East South Central & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[11], Step3.Covariates.Diversity$par[15]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[11], Step3.Covariates.Diversity$se[15]),
      sprintf("West South & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[12], Step3.Covariates.Diversity$par[16]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[12], Step3.Covariates.Diversity$se[16]),
      sprintf("Mountain & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[13], Step3.Covariates.Diversity$par[17]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[13], Step3.Covariates.Diversity$se[17]),
      sprintf("Pacific  & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[14], Step3.Covariates.Diversity$par[18]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[14], Step3.Covariates.Diversity$se[18]),
      sprintf("$\\omega$ & 1.75 & & 1.74 &  \\\\ ", Step3.Covariates$par[15], Step3.Covariates.Diversity$par[19]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[15], Step3.Covariates.Diversity$se[19]),
      sprintf("$\\mu_C\\times 10^{2}$ & %.2f & & %.2f &  \\\\ ", Step3.Covariates$par[16], Step3.Covariates.Diversity$par[20]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[16], Step3.Covariates.Diversity$se[20]),
      sprintf("$\\sigma_C\\times 10^{2}$ & 1.21 & & 1.21 &  \\\\ ", Step3.Covariates$par[17], Step3.Covariates.Diversity$par[21]),
      sprintf(" & (%.2f) & & (%.2f) & \\\\ ", Step3.Covariates$se[17], Step3.Covariates.Diversity$se[21]),
      sprintf("$-\\log \\mathcal L$ & %.2f & & %.2f & \\\\", Step3.Covariates$value, Step3.Covariates.Diversity$value),
      sprintf("Number of Markets & %.0f & %.0f & & %.0f \\\\", NROW(Data$X), sum(Data$X[,10]==1), sum(Data$X[,10]==0)),
      "\\end{tabular}"
    )

con = file("tables/results_fiml_estimates_divisions.tex")
writeLines(empirical.results.table, con)
close.connection(con)

# make the table with the ratios
Step3.Covariates$ratios = c(1/Step3.Covariates$par[1]*10, Step3.Covariates$par[2:4]/Step3.Covariates$par[1:3])
Step3.Covariates.Diversity$ratiosH =
    c(1/Step3.Covariates.Diversity$par[1]*10, Step3.Covariates.Diversity$par[2:4]/Step3.Covariates.Diversity$par[1:3])
Step3.Covariates.Diversity$ratiosL =
    c(1/Step3.Covariates.Diversity$par[5]*10, Step3.Covariates.Diversity$par[6:8]/Step3.Covariates.Diversity$par[5:7])

# compute standard errors using the delta method
Step3.Covariates$ratios.se = vector(mode = "numeric", length = length(Step3.Covariates$ratios))
Step3.Covariates$ratios.se[1] = Step3.Covariates$se[1] * Step3.Covariates$ratios[1]^2 / 10

for(jX in 2:4) {
    Step3.Covariates$ratios.se[jX] = sqrt(c(1/Step3.Covariates$par[jX], -Step3.Covariates$par[jX-1]/Step3.Covariates$par[jX]^2) %*%
                                              Step3.Covariates$covariance[(jX-1):jX,(jX-1):jX] %*%
                                              c(1/Step3.Covariates$par[jX], -Step3.Covariates$par[jX-1]/Step3.Covariates$par[jX]^2))
}

Step3.Covariates.Diversity$ratiosH.se = vector(mode = "numeric", length = length(Step3.Covariates.Diversity$ratiosH))
Step3.Covariates.Diversity$ratiosH.se[1] = Step3.Covariates.Diversity$se[1] * Step3.Covariates.Diversity$ratiosH[1]^2 / 10

for(jX in 2:4) {
    Step3.Covariates.Diversity$ratiosH.se[jX] = sqrt(c(1/Step3.Covariates.Diversity$par[jX], -Step3.Covariates.Diversity$par[jX-1]/Step3.Covariates.Diversity$par[jX]^2) %*%
                                                         Step3.Covariates.Diversity$covariance[(jX-1):jX,(jX-1):jX] %*%
                                                         c(1/Step3.Covariates.Diversity$par[jX], -Step3.Covariates.Diversity$par[jX-1]/Step3.Covariates.Diversity$par[jX]^2))
}

Step3.Covariates.Diversity$ratiosL.se = vector(mode = "numeric", length = length(Step3.Covariates.Diversity$ratiosL))
Step3.Covariates.Diversity$ratiosL.se[1] = Step3.Covariates.Diversity$se[5] * Step3.Covariates.Diversity$ratiosL[1]^2 / 10

for(jX in 2:4) {
    Step3.Covariates.Diversity$ratiosL.se[jX] = sqrt(c(1/Step3.Covariates.Diversity$par[4+jX], -Step3.Covariates.Diversity$par[4+jX-1]/Step3.Covariates.Diversity$par[4+jX]^2) %*%
                                                         Step3.Covariates.Diversity$covariance[4+(jX-1):jX,4+(jX-1):jX] %*%
                                                         c(1/Step3.Covariates.Diversity$par[4+jX], -Step3.Covariates.Diversity$par[4+jX-1]/Step3.Covariates.Diversity$par[4+jX]^2))
}

empirical.results.ratios.table = c(
    "\\begin{tabular}{l c c c}  ",
    "& \\multicolumn{3}{c}{ Geographic Preference Diversity } \\\\ ",
    "&  All $\\mu$SAs  & Diversity $>$ 13.4 miles & Diversity $\\leq$ 13.4 miles  \\\\ ",
    sprintf("$1/\\left(k(1)\\times 10^3\\right)$ & %.2f & %.2f & %.2f \\\\ ", Step3.Covariates$ratios[1], Step3.Covariates.Diversity$ratiosH[1], Step3.Covariates.Diversity$ratiosL[1]),
    sprintf("& (%.2f) & (%.2f) & (%.2f) \\\\ ", Step3.Covariates$ratios.se[1], Step3.Covariates.Diversity$ratiosH.se[1], Step3.Covariates.Diversity$ratiosL.se[1]),
    sprintf("$k(2)/k(1)$ & %.2f & %.2f & %.2f \\\\ ", Step3.Covariates$ratios[2], Step3.Covariates.Diversity$ratiosH[2], Step3.Covariates.Diversity$ratiosL[2]),
    sprintf("& (%.2f) & (%.2f) & (%.2f) \\\\ ",  Step3.Covariates$ratios.se[2], Step3.Covariates.Diversity$ratiosH.se[2], Step3.Covariates.Diversity$ratiosL.se[2]),
    sprintf("$k(3)/k(2)$ & %.2f & %.2f & %.2f \\\\ ", Step3.Covariates$ratios[3], Step3.Covariates.Diversity$ratiosH[3], Step3.Covariates.Diversity$ratiosL[3]),
    sprintf("& (%.2f) & (%.2f) & (%.2f) \\\\ ",  Step3.Covariates$ratios.se[3], Step3.Covariates.Diversity$ratiosH.se[3], Step3.Covariates.Diversity$ratiosL.se[3]),
    sprintf("$k(4)/k(3)$ & %.2f & %.2f & %.2f \\\\ ", Step3.Covariates$ratios[4], Step3.Covariates.Diversity$ratiosH[4], Step3.Covariates.Diversity$ratiosL[4]),
    sprintf("& (%.2f) & (%.2f) & (%.2f) \\\\ ", Step3.Covariates$ratios.se[4], Step3.Covariates.Diversity$ratiosH.se[4], Step3.Covariates.Diversity$ratiosL.se[4]),
    sprintf("Number of Markets & %.0f & %.0f & %.0f \\\\ ", NROW(Data$X), sum(Data$X[,10]==1), sum(Data$X[,10]==0)),
    "\\end{tabular} ")

con = file("tables/results_fiml_ratios_divisions.tex")
writeLines(empirical.results.ratios.table, con)
close.connection(con)

save(Step3.Covariates, Step3.Covariates.Diversity, file = "results/empirical.RData")

# end cluster
stopCluster(cl)

# stop the clock
time_end = Sys.time()
print(time_end - time_begin)
