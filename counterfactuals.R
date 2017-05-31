# load package.
library(actyR)

# make sure we are in the correct working directory
if (!all(c("graphs", "results", "tables") %in% list.files())) {
  stop("Please navigate to the correct working directory
       using `setwd`. Your working directory
       should include the folders graphs, results, and tables.")
}

# how many cores?
cl = makePSOCKcluster(parallel::detectCores(), outfile = file.path(tempdir(), "cluster_counterfactuals.log"))
registerDoParallel(cl)

# Load the data from the package
data(cov)
data(firms)
data(population)

# Get rid of the first two columns
Data = list(
    "C" = as.matrix(population[, -c(1:2)]),
    "N" = as.matrix(firms[, -c(1:2)]),
    "X" = as.matrix(cov[, -c(1:2)])
)

# demean income
Data$X[,1] = log(as.numeric(Data$X[,1])) - log(mean(as.numeric(Data$X[,1])))

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
for (iX in seq(1,nrow(Data$C))) {
    for (tX in seq(1,ncol(Data$C))) {
        Data$C[iX, tX] = which.min(abs(Data$C[iX, tX] / 10000 - exp(Settings$logGrid)))
    }
}

# Initialize Param
Param = list(rho = 1 / 1.05, k = NA, phi = NA, omega = NA, Pi = NA, mu = NA, sigma = NA, kappa = 1)

# Model Estimates
if (file.exists("results/empirical.RData")) {

    load("results/empirical.RData")

} else {

    # use defaults
    Step3.Covariates.Diversity = list(
        "par" = c(0.339919904280447,
                  0.203904920745783,
                  0.170969277093272,
                  0.135316805751601,
                  0.384977662327072,
                  0.184905095971473,
                  0.144054268197212,
                  0.0957307161137273,
                  49.0611353611962,
                  1.74401092882882,
                  0.846343564267042,
                  -0.603498627661048,
                  -0.451086899041058,
                  0.0917481983210092,
                  -0.709561212668823,
                  -0.492424728653295,
                  -0.285562904009737,
                  -0.0706475942062051,
                  -0.0843932034976282,
                  0.336201747489111,
                  1.20808383848113))

}

## Map estimates into Param structure
Param$kL = c(Step3.Covariates.Diversity$par[5:8], rep(Step3.Covariates.Diversity$par[8], 5))
Param$kH = c(Step3.Covariates.Diversity$par[1:4], rep(Step3.Covariates.Diversity$par[4], 5))
Param$k = Param$kH
Param$phi = rep(Step3.Covariates.Diversity$par[9], 9)
Param$omega = Step3.Covariates.Diversity$par[10]
Param$mu = Step3.Covariates.Diversity$par[20]
Param$sigma = Step3.Covariates.Diversity$par[21]
Param$Pi = tauchen(Settings$logGrid, Param$mu, Param$sigma)
Param$beta = Step3.Covariates.Diversity$par[11:19]

# define useful helper functions

# create extract_n_matrix
computeExtractNmatrix = function(Settings) {

    extract_n_matrix = matrix(0, nrow = ((Settings$nCheck + 1) * Settings$cCheck), ncol = Settings$nCheck + 1)
    for (n in 0:Settings$nCheck) {
        extract_n_matrix[(n * Settings$cCheck + 1):((n + 1) * Settings$cCheck), n + 1] = 1
    }
    return(extract_n_matrix)

}


# create extract_c_matrix
computeExtractCmatrix = function(Settings, Param) {

    extract_c_matrix = matrix(0, nrow = Settings$nCheck + 1, ncol = (Settings$nCheck + 1) *
                          Settings$cCheck)
    for (nX in 1:(1 + Settings$nCheck)) {
        extract_c_matrix[nX, (Settings$cCheck * (nX - 1) + 1):(Settings$cCheck * nX)] = ergodic(Param$Pi)
    }

    return(extract_c_matrix)

}


# compute the model means for t_forward periods
computeModelMeans = function(ParamPrime, Settings, t_forward, initial_c, parallel = FALSE) {

    extract_n_matrix = computeExtractNmatrix(Settings)
    trans_mat = getModelTransitionMatrix(Settings, ParamPrime)
    trans_mat[trans_mat < 1e-7] = 0
    trans_mat = trans_mat / rowSums(trans_mat)
    nX = 0

    `%op%` = ifelse(parallel, `%dopar%`, `%do%`)

    return(t(foreach(nX = 0:Settings$nCheck, .combine = cbind) %op%  {
        distribution = matrix(0, nrow = t_forward, ncol = Settings$cCheck * (Settings$nCheck + 1))
        distribution[1, nX * Settings$cCheck + 1:Settings$cCheck] = initial_c
        for (tX in 2:t_forward) {
            distribution[tX,] = distribution[tX - 1,] %*% trans_mat
            distribution[tX,] = distribution[tX,] / sum(distribution[tX,])
        }
        distributionN = distribution %*% extract_n_matrix
        means = distributionN %*% matrix(0:Settings$nCheck, ncol = 1)
        return(means)
    }))

}

# compute the model means at time 30 as a function of demand at time zero
computeModelMeansT30 = function(Settings, Param.factual, extract_n_matrix, kappa = NULL) {

    if (!is.null(kappa)) {
        Param.factual$kappa = kappa
    }
    trans_mat.factual = getModelTransitionMatrix(Settings, Param.factual)
    trans_mat.factual31 = trans_mat.factual %^% 30
    N30meansVector.factual = trans_mat.factual31 %*% extract_n_matrix %*% matrix(0:Settings$nCheck, ncol = 1)
    N30means.factual = matrix(N30meansVector.factual, nrow = 200)
    return(N30means.factual)
}


## make initial distributions for the demand state (q1, q2, q3)
initialDistr = list(q1 = NULL, q2 = NULL, q3 = NULL)

# get the first quartiles in 2009
q1 = quantile(population$pop2009, 0.25)
initialDistr$q1 = matrix(0, nrow = 1, ncol = Settings$cCheck)
initialDistr$q1[which.min(abs(q1 - exp(Settings$logGrid)*1e4))] = 1

# get the second quartiles in 2009
q2 = median(population$pop2009)
initialDistr$q2 = matrix(0, nrow = 1, ncol = Settings$cCheck)
initialDistr$q2[which.min(abs(q2 - exp(Settings$logGrid)*1e4))] = 1

# get the third quartiles in 2009
q3 = quantile(population$pop2009, 0.75)
initialDistr$q3 = matrix(0, nrow = 1, ncol = Settings$cCheck)
initialDistr$q3[which.min(abs(q3 - exp(Settings$logGrid)*1e4))] = 1

# Get extract_n_matrix
extract_n_matrix = computeExtractNmatrix(Settings)

# make different Param structures (make a function for each of counterfactual)

# Factual
factual = function(ParamPrime) {
    return(ParamPrime)
}
Param.factual = factual(Param)

# No-Sunk Cost Experiment
nosunk = function(ParamPrime) {
    Param.nosunk = ParamPrime
    Param.nosunk$phi[] = 0
    Param.nosunk$kappa = 1 + ParamPrime$phi[1] * (1 - ParamPrime$rho)
    return(Param.nosunk)
}
Param.nosunk = nosunk(Param)

# netflix
# 25% Demand Decrease
netflix = function(ParamPrime) {
    Param.netflix = Param
    Param.netflix$k = ParamPrime$k * 0.75
    return(Param.netflix)
}
Param.netflix = netflix(Param)

netflix.joa.all = function(ParamPrime) {
    # 25% Demand Decrease: with joa for all firms
    Param.netflix.joa.all = ParamPrime
    Param.netflix.joa.all$k = ParamPrime$k * 0.75
    Param.netflix.joa.all$k[] = Param.netflix.joa.all$k[1]
    return(Param.netflix.joa.all)
}
Param.netflix.joa.all = netflix.joa.all(Param)

# 25% Demand Decrease: with joa for at most two firms
netflix.joa.duo = function(ParamPrime) {
    Param.netflix.joa.duo = ParamPrime
    Param.netflix.joa.duo$k = ParamPrime$k * 0.75
    Param.netflix.joa.duo$k[2] = Param.netflix.joa.duo$k[1]
    return(Param.netflix.joa.duo)
}
Param.netflix.joa.duo = netflix.joa.duo(Param)


# number of periods to look ahead
t_forward = 31

# Part 1: factual vs no-sunk costs

# t -> E[N_t | X,C0,N0=n]
means.factual.q1 = computeModelMeans(Param.factual,
                                     Settings,
                                     t_forward,
                                     initialDistr$q1,
                                     parallel = TRUE)

means.factual.q3 = computeModelMeans(Param.factual,
                                     Settings,
                                     t_forward,
                                     initialDistr$q3,
                                     parallel = TRUE)

means.nosunk.q1 = computeModelMeans(Param.nosunk,
                                     Settings,
                                     t_forward,
                                     initialDistr$q1,
                                     parallel = TRUE)

means.nosunk.q3 = computeModelMeans(Param.nosunk,
                                     Settings,
                                     t_forward,
                                     initialDistr$q3,
                                     parallel = TRUE)

# log(C0) -> log(E[N30 | X,C0,N_0])
N30means.factual = computeModelMeansT30(Settings,
                                        Param.factual,
                                        extract_n_matrix)


N30means.nosunk = computeModelMeansT30(Settings,
                                       Param.nosunk,
                                       extract_n_matrix)


N30means.calibrated = computeModelMeansT30(Settings,
                                           Param.nosunk,
                                           extract_n_matrix,
                                           4.75)

## Part 2: Netflix shock of 25% and JOAs

# t -> E[N_t | X,C0,N0=n]
means.netflix.q1 = computeModelMeans(Param.netflix,
                                     Settings,
                                     t_forward,
                                     initialDistr$q1,
                                     parallel = TRUE)

means.netflix.q3 = computeModelMeans(Param.netflix,
                                     Settings,
                                     t_forward,
                                     initialDistr$q3,
                                     parallel = TRUE)

means.netflix.joa.all.q1 = computeModelMeans(Param.netflix.joa.all,
                                             Settings,
                                             t_forward,
                                             initialDistr$q1,
                                             parallel = TRUE)

means.netflix.joa.all.q3 = computeModelMeans(Param.netflix.joa.all,
                                             Settings,
                                             t_forward,
                                             initialDistr$q3,
                                             parallel = TRUE)

means.netflix.joa.duo.q1 = computeModelMeans(Param.netflix.joa.duo,
                                             Settings,
                                             t_forward,
                                             initialDistr$q1,
                                             parallel = TRUE)

means.netflix.joa.duo.q3 = computeModelMeans(Param.netflix.joa.duo,
                                             Settings,
                                             t_forward,
                                             initialDistr$q3,
                                             parallel = TRUE)

# log(C0) -> log(E[N30 | X,C0,N_0])
N30means.netflix = computeModelMeansT30(Settings,
                                        Param.netflix,
                                        extract_n_matrix)

N30means.netflix.joa.all = computeModelMeansT30(Settings,
                                                Param.netflix.joa.all,
                                                extract_n_matrix)

N30means.netflix.joa.duo = computeModelMeansT30(Settings,
                                                Param.netflix.joa.duo,
                                                extract_n_matrix)

# aggregate over markets in the data

# now compute means for all markets using the empirical initial values (n_0, c_0)
# this function returns a vector with t_forward elements
computeMarketMeans = function(rX, treatment) {

    ParamPrime = Param
    if (Data$X[rX, 10]==0) {
        ParamPrime$k = ParamPrime$kH * exp(sum(Data$X[rX, 1:9] * ParamPrime$beta))
    } else {
        ParamPrime$k = ParamPrime$kL * exp(sum(Data$X[rX, 1:9] * ParamPrime$beta))
    }

    ParamPrime = treatment(ParamPrime)

    init_distr = matrix(0, nrow = 1, ncol = Settings$cCheck)
    init_distr[Data$C[rX, 10]] = 1

    return(computeModelMeans(ParamPrime, Settings, t_forward, init_distr)[Data$N[rX, 10] + 1, ])
}

# t -> E[N_t | X,C0,N0=n] where we integrate over X, C0, N0 using the data
means.empirical.factual = t(foreach(rX = 1:NROW(Data$X), .combine = cbind, .packages = "actyR") %dopar%  {
    computeMarketMeans(rX, factual)
})

colMeans(means.empirical.factual)

means.empirical.netflix = t(foreach(rX = 1:NROW(Data$X), .combine = cbind, .packages = "actyR") %dopar%  {
    computeMarketMeans(rX, netflix)
})

means.empirical.netflix.joa.all = t(foreach(rX = 1:NROW(Data$X), .combine = cbind, .packages = "actyR") %dopar%  {
    computeMarketMeans(rX, netflix.joa.all)
})

means.empirical.netflix.joa.duo = t(foreach(rX = 1:NROW(Data$X), .combine = cbind, .packages = "actyR") %dopar%  {
    computeMarketMeans(rX, netflix.joa.duo)
})

# make different selections based on demand in year 2009 so that we can plot subgroups
# of the data later
c0.fivenum = fivenum(Settings$logGrid[Data$C[,10]])
select.q1 = Settings$logGrid[Data$C[,10]]<c0.fivenum[2]
select.q2 = Settings$logGrid[Data$C[,10]]>=c0.fivenum[2] & Settings$logGrid[Data$C[,1]]<c0.fivenum[3]
select.q3 = Settings$logGrid[Data$C[,10]]>=c0.fivenum[3] & Settings$logGrid[Data$C[,1]]<c0.fivenum[4]
select.q4 = Settings$logGrid[Data$C[,10]]>=c0.fivenum[4]

## Step 3: Find the Netflix Shock that would perfectly offset JOAs in expectation

# returns the transition path for a given market rX with a JOA for all firms under
# a netflix shock `shock` (shock=1 corresponds to the baseline)
computeMarketMeans.joa.all = function(rX, shock) {
    ParamPrime = Param
    if (Data$X[rX, 10]==0) {
        ParamPrime$k = ParamPrime$kH * exp(sum(Data$X[rX, 1:9] * ParamPrime$beta)) * shock
    } else {
        ParamPrime$k = ParamPrime$kL * exp(sum(Data$X[rX, 1:9] * ParamPrime$beta)) * shock
    }
    ParamPrime$k[] = ParamPrime$k[1]
    init_distr = matrix(0, nrow = 1, ncol = Settings$cCheck)
    init_distr[Data$C[rX, 10]] = 1
    return(computeModelMeans(ParamPrime, Settings, t_forward, init_distr)[Data$N[rX, 10] + 1, ])
}

shock.joa.all = 0.62 # this is the answer
means.empirical.netflix.joa.all.offset = t(foreach(rX = 1:NROW(Data$X), .combine = cbind, .packages = "actyR") %dopar%  {
    computeMarketMeans.joa.all(rX, shock.joa.all)
})
colMeans(means.empirical.netflix.joa.all.offset)

# returns the transition path for a given market rX with a JOA for two firms under
# a netflix shock `shock` (shock=1 corresponds to the baseline)
computeMarketMeans.joa.duo = function(rX, shock) {
    ParamPrime = Param
    if (Data$X[rX, 10]==0) {
        ParamPrime$k = ParamPrime$kH * exp(sum(Data$X[rX, 1:9] * ParamPrime$beta)) * shock
    } else {
        ParamPrime$k = ParamPrime$kL * exp(sum(Data$X[rX, 1:9] * ParamPrime$beta)) * shock
    }
    ParamPrime$k[2] = ParamPrime$k[1]
    init_distr = matrix(0, nrow = 1, ncol = Settings$cCheck)
    init_distr[Data$C[rX, 10]] = 1
    return(computeModelMeans(ParamPrime, Settings, t_forward, init_distr)[Data$N[rX, 10] + 1, ])
}

shock.joa.duo = 0.82 # this is the answer
means.empirical.netflix.joa.duo.offset = t(foreach(rX = 1:NROW(Data$X), .combine = cbind, .packages = "actyR") %dopar%  {
    computeMarketMeans.joa.duo(rX, shock.joa.duo)
})
colMeans(means.empirical.netflix.joa.duo.offset)

## make all the graphs
legendString = c("estimated model", "Netflix", "Netflix + JOA", "Netflix + JOA duopoly")

# Part 1
plotNumberActiveFirms = function(y, filename, title = NULL) {

    nCheck = NCOL(y) - 1
    years = NROW(y) - 1

    if (is.null(title)) {
        pdf(file = paste("graphs/", filename, ".pdf", sep = ""), width = 5, height = 5)
        par(mar = c(c(3, 3, 0, 0) + 0.1), mgp = c(2.1, 0.9, 0))
    } else {
        pdf(file = paste("graphs/", filename, ".pdf", sep = ""), width = 5, height = 5)
        par(mar = c(c(3, 3, 3, 0) + 0.1), mgp = c(2.1, 0.9, 0))
    }

    matplot(0:years, y, type = "l",
            xlim = c(-1,years), ylim = c(0,4), lty = 1, col = 1:5,  lwd = 6,
            ylab = "number of firms", xlab = "years elapsed", cex.axis = 0.75, cex.lab = 0.75)
    if (!is.null(title)) {
        title(main=list(title, cex=0.75))
    }
    if (nCheck > 0) {
        legend(0.2, 4, 0:nCheck, lty = 1, lwd = 6, col = 0:nCheck + 1, title = "initial number of firms", horiz = TRUE, bty = "n", cex = 0.75)
    }
    dev.off()

}

plotNumberActiveFirms(t(means.factual.q1)[,1:5], "means-factual-q1", title = "first demand quartile")
plotNumberActiveFirms(t(means.factual.q3)[,1:5], "means-factual-q3", title = "third demand quartile")

pdf("graphs/longrun-distribution.pdf", width = 5, height = 5)
par(mar = c(c(3, 3, 3, 0) + 0.1), mgp = c(2.1, 0.9, 0))
plot(Settings$logGrid + log(10000), log(N30means.factual[, 1]), type = "l",
     ylab = "longrun average number of firms in logs",
     xlab = "initial demand state (log population)",
     lwd = 6, col = "black", ylim = c(-2.5, 2.0), cex.axis = 0.75, cex.lab = 0.75, lty = 1, xlim = c(9, 12.5))
lines(Settings$logGrid + log(10000), log(N30means.nosunk[, 1]), type = "l", col = "orange", lwd = 6, lty = 1)
legend(10.0, -1, c("estimated model", "no sunk costs"),
       lwd = 6, col = c("black", "orange"), lty = 1, horiz = FALSE, bty = "n", cex = 0.75)
title(main=list("no sunk cost counterfactual", cex=0.75))
dev.off()

# Part 2

plotIRF = function(y, filename, ylim=c(0,0.3), legend.y = -.26, title = NULL) {

    # set element (1,1) to zero
    y[1,1] = 0

    nCheck = NCOL(y) - 1
    years = NROW(y) - 1

    pdf(file = paste("graphs/", filename, ".pdf", sep = ""), width = 5, height = 5)
    par(mar = c(c(3, 3, 3, 0) + 0.1), mgp = c(2.1, 0.9, 0))
    matplot(0:years, y,
            type = "l", xlim = c(-1,years), lty = 1, col = 0:nCheck + 1, lwd = 5,
            ylab = "relative change", xlab = "years elapsed", cex.axis = 0.75, cex.lab = 0.75, ylim = ylim)
    if (nCheck > 0) {
        legend(-1.25,
               legend.y,
               0:nCheck,
               lty = 1,
               lwd = 5,
               col = 0:nCheck + 1,
               title = "initial number of firms",
               horiz = TRUE,
               bty = "n",
               cex = 0.75)
    }
    if (!is.null(title)) {
        title(main=list(title, cex=0.75))
    }
    dev.off()

}

plotIRF(t((means.netflix.q1-means.factual.q1)/means.factual.q1)[,1:5], "irf-netflix-q1", c(-.31, 0), title = "first demand quartile")
plotIRF(t((means.netflix.q3-means.factual.q3)/means.factual.q3)[,1:5], "irf-netflix-q3", c(-.31, 0), title = "third demand quartile")

pdf("graphs/longrun-distribution-netflix.pdf", width = 5, height = 5)
par(mar = c(c(3, 3, 3, 0) + 0.1), mgp = c(2.1, 0.9, 0))
plot(Settings$logGrid + log(10000), log(N30means.factual[, 1]), type = "l",
     ylab = "longrun average number of firms in logs",
     xlab = "initial demand state (log population)",
     lwd = 5, col = 1, ylim = c(-2.5, 2.0), cex.axis = 0.75, cex.lab = 0.75, lty = 1, xlim = c(9,12.5))
lines(Settings$logGrid + log(10000), log(N30means.netflix[, 1]), type = "l", col = 2, lwd = 5, lty = 1)
lines(Settings$logGrid + log(10000), log(N30means.netflix.joa.all[, 1]), type = "l", col = 3, lwd = 5, lty = 1)
lines(Settings$logGrid + log(10000), log(N30means.netflix.joa.duo[, 1]), type = "l", col = 4, lwd = 5, lty = 1)
title(main=list("Netflix counterfactual", cex = 0.75))
legend(10.0, -1, legendString,
       lwd = 5, col = 1:4, lty = 1, horiz = FALSE, bty = "n", cex = 0.75)
dev.off()

plotNetflixAggregateTransition = function(select = rep(TRUE, NROW(Data$X)), filename, legend.y = 1, ylim = c(0.5, 1.8), title = NULL) {

    pdf(file = paste("graphs/", filename, ".pdf", sep = ""), width = 5, height = 5)
    par(mar = c(c(3, 3, 3, 0) + 0.1), mgp = c(2.1, 0.9, 0))
    plot(0:30, colMeans(means.empirical.factual[select,]), type = "l",
         ylab = "average number of active firms",
         xlab = "years elapsed",
         lwd = 5, col = 1, cex.axis = 0.75, cex.lab = 0.75, lty = 1, ylim = ylim, pch = 1)
    title(main=list(title, cex=0.75))
    lines(0:30, colMeans(means.empirical.netflix[select,]), col = 2, lwd = 5)
    lines(0:30, colMeans(means.empirical.netflix.joa.all[select,]), col = 3, lwd = 5)
    lines(0:30, colMeans(means.empirical.netflix.joa.duo[select,]), col = 4, lwd = 5)
    legend(0.2, legend.y, legendString,
           lwd = 5, col = 1:4, lty = 1, horiz = FALSE, bty = "n", cex = 0.75)
    dev.off()

}

plotNetflixAggregateTransition(filename = "netflix-aggregate-transition-q1",
                               select = select.q1,
                               legend.y = 3.0,
                               ylim = c(0.5, 3.1),
                               title = 'first demand quartile')

plotNetflixAggregateTransition(filename = "netflix-aggregate-transition-q2",
                               select = select.q2,
                               legend.y = 3.0,
                               ylim = c(0.5, 3.1),
                               title = 'second demand quartile')

plotNetflixAggregateTransition(filename = "netflix-aggregate-transition-q3",
                               select = select.q3,
                               legend.y = 3.0,
                               ylim = c(0.5, 3.1),
                               title = 'third demand quartile')

plotNetflixAggregateTransition(filename = "netflix-aggregate-transition-q4",
                               select = select.q4,
                               legend.y = 1.5,
                               ylim = c(0.5, 3.1),
                               title = 'fourth demand quartile')


# end cluster
stopCluster(cl)

# stop the clock
time_end = Sys.time()
print(time_end - time_begin)
