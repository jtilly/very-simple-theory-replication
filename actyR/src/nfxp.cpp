// nfxp.cpp

#include <cmath>
#include <math.h>
#include <actyR.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "lib.h"
#include "nfxp.h"
#include "nfxpGradient.h"

//' Compute equilibrium value functions
//'
//' This function computes the equilibrium value functions
//' and the probabilities of entry and certain survival. Value
//' functions are computed via backward recursion. For this
//' function to work, the following \code{Settings} must have been
//' set:
//' \itemize{
//'  \item{\code{nCheck} The maximum number of firms with non-zero flow profits in some state of the world.}
//'  \item{\code{logGrid} The log of the grid}
//'  \item{\code{tolInner} The convergence tolerance for the inner loop}
//'  \item{\code{maxIterInner} The maximum number of iterations for the inner loop}
//' }
//' The following elements in \code{Param} must have been set
//' \itemize{
//'  \item{\code{k} A vector of length \code{nCheck} that determines the flow surplus}
//'  \item{\code{phi} A vector of length \code{nCheck} that governs the entry costs}
//'  \item{\code{omega} A scalar that governs the variance of the cost shocks}
//'  \item{\code{rho} The discount factor}
//'  \item{\code{Pi} A transition probability matrix for the demand grid, which is computed via
//'  \code{tauchen}}
//' }
//'
//' @param Settings is an R list of settings
//' @param Param is an R list of of parameters
//' @return List with value functions, entry probabilities and
//' certain survival probabilities.
// [[Rcpp::export]]
List valuefunction(Rcpp::List Settings, Rcpp::List Param) {

    int nCheck = Settings["nCheck"];
    rowvec logGrid = Settings["logGrid"];
    int maxIterInner = Settings["maxIterInner"];
    double tolInner = Settings["tolInner"];

    vec k = Param["k"];
    vec phi = Param["phi"];
    double omega = Param["omega"];
    double rho = Param["rho"];
    double kappa = Param["kappa"];
    mat Pi = Param["Pi"];

    // find ccheck
    int cCheck = Pi.n_cols;

    // initialize the important matrices
    mat vS(nCheck + 1, cCheck), pEntry(nCheck + 1, cCheck), pEntrySet(nCheck + 1, cCheck), pStay(nCheck, cCheck);

    // init helper variable flow surplus
    rowvec flowSurplus(cCheck);

    // init helper variables
    rowvec logvSn(cCheck), vSnPrime(cCheck), partialExp(cCheck), valueSureSurvNoEntry(cCheck), valueAdditionalEntry(cCheck), vSn(cCheck);

    // get the parameters of the distribution of W
    double w_mean = -.5 * pow(omega, 2);
    double w_std = omega;

    valueAdditionalEntry.zeros(); vS.zeros(); pEntry.zeros(); pEntrySet.zeros(); pStay.zeros();

    // transpose the Pi matrix once here (so we don't have to repeatedly transpose it in a loop)
    mat Pitranspose = Pi.t();
    vSn.ones();

    // start the value function iteration beginning with nCheck
    for(int n = nCheck - 1; n >= 0; n--) {

        int iter = 0;
        double vSdiff = 1.0;

        flowSurplus = exp(logGrid) * k(n) / (n + 1);
        valueAdditionalEntry += pEntrySet.row(n + 1) % vS.row(n + 1);

        // do the value function iteration
        while (vSdiff>tolInner && iter++ < maxIterInner) {

            logvSn  = log(vSn);
            pStay.row(n)  =  mynormcdf(logvSn - log(kappa), w_mean, w_std);
            partialExp = 1.0 - mynormcdf(- logvSn + log(kappa), w_mean, w_std);
            valueSureSurvNoEntry = (pStay.row(n) - pEntry.row(n + 1)) % vSn;
            vSnPrime = rho * (flowSurplus - kappa*partialExp + valueSureSurvNoEntry + valueAdditionalEntry) * Pitranspose;
            vSdiff = as_scalar(max(abs(vSn - vSnPrime)));
            vSn = vSnPrime;

        }

        if(vSdiff>tolInner) {
            ::Rf_warning("Algorithm did not converge.");
        }

    vS.row(n) = vSn;
    pEntry.row(n) = mynormcdf(logvSn - log(kappa+phi(n)), w_mean, w_std);
    pEntrySet.row(n) = pEntry.row(n)-pEntry.row(n + 1);

    }

    return List::create(_["vS"]         = vS,
                        _["pEntry"]     = pEntry,
                        _["pEntrySet"]  = pEntrySet,
                        _["pStay"]      = pStay);

}

//' Compute the mixing density
//'
//' This function computes the mixing density and returns it as a cube.
//'
//' @param Settings is an R list of settings
//' @param Param is an R list of of parameters
//' @param vS is a matrix with the value functions
//' @return A cube with the mixing density.
// [[Rcpp::export]]
cube computeMixingDensity(Rcpp::List Settings, Rcpp::List Param, mat vS) {

    int nCheck = Settings["nCheck"];
    int cCheck = Settings["cCheck"];

    vec w = Settings["integrationWeights"];
    vec p = Settings["integrationNodes"];


    // nchoosekmatrix
    mat nchoosekMatrix(nCheck + 1, nCheck + 1);
    for(int nX = 0; nX < (int) nCheck + 1; nX++) {
        for(int jX = 0; jX < (int) nCheck + 1; jX++) {
            nchoosekMatrix(nX, jX) = nchoosek(nX, jX);
        }
    }

    cube  mixingDensity(p.n_elem, cCheck, nCheck);
    cube  aSinv(p.n_elem, cCheck, nCheck);
    cube  daSinvdP(p.n_elem, cCheck, nCheck);
    mixingDensity.zeros(); aSinv.zeros(); daSinvdP.zeros();

    mat intWeights = repmat(w, 1, cCheck);
    double w_mean = -.5 * pow((double) Param["omega"], 2);
    double w_std = Param["omega"];
    double kappa = Param["kappa"];

    // mixing weights
    for(int n = 2; n <= nCheck; n++) {

        mat expaSInv(p.n_elem, cCheck);
        expaSInv.zeros();
        mat dexpaSInvdP(p.n_elem, cCheck);
        dexpaSInvdP.zeros();

        for(int nPrime = 1; nPrime <= n; nPrime++) {

            int nCk = nchoosekMatrix(n - 1, nPrime - 1);
            mat binomialPmf     = nCk * repmat(pow(p, nPrime - 1) % pow(1-p, n-nPrime), 1, cCheck);
            mat dbinomialPmfdP  = nCk * repmat(pow(p, nPrime-2) % pow(1-p, n-nPrime - 1) % ((nPrime - 1) * (1-p) - p * (n-nPrime)), 1, cCheck);
            mat repvS =  repmat(vS.row(nPrime - 1), p.n_elem,1);
            expaSInv +=  binomialPmf % repvS ;
            dexpaSInvdP += dbinomialPmfdP % repvS;

        }

        aSinv.slice(n - 1) =  log(expaSInv) - log(kappa);
        daSinvdP.slice(n - 1) =  dexpaSInvdP / expaSInv;

        mat normaSinv = mynormpdf(aSinv.slice(n - 1), w_mean, w_std);
        mixingDensity.slice(n - 1) = daSinvdP.slice(n - 1) % normaSinv % intWeights;

    }

    return(mixingDensity);

}


//' Compute likelihood contributions
//'
//' This function computes the vector of likelihood contributions for the
//' data set in \code{Data}. The likelihood contributions are computed as follows:
//' \itemize{
//' \item{ Use the function \code{valuefunction} to compute the equilibrium value functions
//' and the entry and sure survival probabilities}
//' \item{Compute the mixing density}
//' \item{Assemble the likelihood contributions for each observations}
//' }
//'
//' @param Settings is an R list of settings
//' @param Param is an R list of parameters
//' @param Data is an R list of the data and includes the matrices \code{N} and \code{C}
//' @return A column vector of likelihood contributions for each observation
// [[Rcpp::export]]
colvec likelihood(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Data) {

    int nCheck = Settings["nCheck"];

    vec w = Settings["integrationWeights"];
    vec p = Settings["integrationNodes"];

    Rcpp::List Eq = valuefunction(Settings, Param);

    mat vS = Eq["vS"];
    mat pEntry = Eq["pEntry"];
    mat pEntrySet = Eq["pEntrySet"];
    mat pStay = Eq["pStay"];

    mat C = Data["C"];
    mat N = Data["N"];

    int M = C.n_rows;
    int T = C.n_cols;

    // nchoosekmatrix
    mat nchoosekMatrix(nCheck + 1, nCheck + 1);
    for(int nX = 0; nX < nCheck + 1; nX++) {
        for(int jX = 0;jX < nCheck + 1; jX++) {
            nchoosekMatrix(nX, jX) = nchoosek(nX, jX);
        }
    }

    // compute the mixing density
    cube mixingDensity = computeMixingDensity(Settings, Param, vS);

    // compute likelihood
    int jX = 0;
    colvec likelihoodContributions((T - 1)*M);
    likelihoodContributions.zeros();

    // loop over all markets
    for(int mX = 0; mX < M; mX++) {

        // loop over time
        for(int tX = 1; tX < T; tX++) {

            // entry
            if (N(mX, tX) > N(mX, tX - 1)) {
                likelihoodContributions(jX) = pEntrySet(N(mX, tX) - 1, C(mX, tX) - 1);
            }
            // 0 and no entry
            if (N(mX, tX) == 0 && N(mX, tX - 1) == 0) {
                likelihoodContributions(jX) = 1-pEntry(0, C(mX, tX) - 1) ;
            }
            // number of firms stays the same (yet >0)
            if (N(mX, tX) > 0 && N(mX, tX) == N(mX, tX - 1)) {
                likelihoodContributions(jX) = pStay(N(mX, tX) - 1, C(mX, tX) - 1) - pEntry(N(mX, tX), C(mX, tX) - 1);
            }
            // all firms leave
            if (N(mX, tX) == 0 && N(mX, tX - 1) > 0) {
                likelihoodContributions(jX) =  1 - pStay(0, C(mX, tX) - 1);
            }
            // mixing
            if (N(mX, tX - 1) > 1 && N(mX, tX)<=N(mX, tX - 1)) {

                likelihoodContributions(jX) += -accu((nchoosekMatrix(N(mX,tX - 1), N(mX,tX)) * pow(p, N(mX,tX)) % pow(1-p, N(mX,tX - 1)-N(mX,tX)))
                                 % mixingDensity.slice(N(mX,tX - 1) - 1).col(C(mX,tX) - 1));

            }

            jX++;

        }

    }

    return(likelihoodContributions);
}


//' Computes the model's transition matrix
//'
//' Computes the model's implied transition matrix
//' @param Settings is an R list with settings
//' @param Param is an R list with parameters
//' @return square matrix with cCheck * (nCheck + 1) rows and columns
// [[Rcpp::export]]
mat getModelTransitionMatrix(Rcpp::List Settings, Rcpp::List Param) {

    int nCheck = Settings["nCheck"];
    int cCheck = Settings["cCheck"];

    mat Pi = Param["Pi"];

    mat transMat((nCheck + 1) * cCheck, (nCheck + 1) * cCheck);
    transMat.zeros();

    vec w = Settings["integrationWeights"];
    vec p = Settings["integrationNodes"];

    Rcpp::List Eq = valuefunction(Settings, Param);

    mat vS = Eq["vS"];
    mat pEntry = Eq["pEntry"];
    mat pEntrySet = Eq["pEntrySet"];
    mat pStay = Eq["pStay"];

    // nchoosekmatrix
    mat nchoosekMatrix(nCheck + 1, nCheck + 1);
    for(int nX = 0; nX < nCheck + 1; nX++) {
        for(int jX = 0;jX < nCheck + 1; jX++) {
            nchoosekMatrix(nX, jX) = nchoosek(nX, jX);
        }
    }

    // compute the mixing density
    cube mixingDensity = computeMixingDensity(Settings, Param, vS);

    // loop over all markets
    for(int n = 0; n <= nCheck; n++) {

        for(int c = 0; c < cCheck; c++) {

            for(int nPrime = 0; nPrime <= nCheck; nPrime++) {

                for(int cPrime = 0; cPrime < cCheck; cPrime++) {

                    // entry
                    if (nPrime > n) {
                        transMat(n * cCheck + c, nPrime * cCheck + cPrime) = pEntrySet(nPrime - 1, cPrime);
                    }
                    // 0 and no entry
                    if (nPrime == 0 && n == 0) {
                        transMat(n * cCheck + c, nPrime * cCheck + cPrime) = 1.0 - pEntry(0, cPrime) ;
                    }
                    // number of firms stays the same (yet >0)
                    if (nPrime > 0 && nPrime == n) {
                        transMat(n * cCheck + c, nPrime * cCheck + cPrime) = pStay(nPrime - 1, cPrime) - pEntry(nPrime, cPrime);
                    }
                    // all firms leave
                    if (nPrime == 0 && n > 0) {
                        transMat(n * cCheck + c, nPrime * cCheck + cPrime) =  1.0 - pStay(0, cPrime);
                    }
                    // mixing
                    if (n > 1 && nPrime <= n) {
                        transMat(n * cCheck + c, nPrime * cCheck + cPrime) += -accu((nchoosekMatrix(n, nPrime) * pow(p, nPrime) % pow(1.0 - p, n - nPrime)) % mixingDensity.slice(n - 1).col(cPrime));
                    }

                    transMat(n * cCheck + c, nPrime * cCheck + cPrime) *= Pi(c, cPrime);

                }
            }

        }

    }

    return transMat;

}
