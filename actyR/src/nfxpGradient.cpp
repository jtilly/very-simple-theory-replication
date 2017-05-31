// nfxpGradient.cpp

#include <cmath>
#include <math.h>
#include <actyR.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "lib.h"
#include "nfxp.h"
#include "nfxpGradient.h"

// Compute the gradient of the value function with respect to all model parameters
//
// This function computes the gradient of the value function with respect to all model
// parameters by applying the implicit function theorem to the Bellman operator. It
// is only called from within C++ and returns a structure of gradients.
//
// @param Settings is an R list of settings
// @param Param is an R list of parameters
// @param Eq is an R list of equilibrium objects
// @return A structure of gradients
GradientStructure valuefunctionGradient(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Eq) {

    int nCheck = Settings["nCheck"];
    int cCheck = Settings["cCheck"];
    mat vS = Eq["vS"];
    mat pEntry = Eq["pEntry"];
    mat pStay = Eq["pStay"];
    rowvec logGrid = Settings["logGrid"];
    double w_mean = -.5*pow((double) Param["omega"], 2);
    double w_std = Param["omega"];
    double omega = Param["omega"];
    double rho = Param["rho"];
    vec k = Param["k"];
    vec phi = Param["phi"];
    mat Pi = Param["Pi"];

    mat Pitranspose = Pi.t();

    cube DvS(nCheck + 1, cCheck, 2*nCheck+3), DpEntry(nCheck + 1, cCheck, 2*nCheck+3),   DpStay(nCheck + 1, cCheck, 2*nCheck+3);
    DvS.zeros(); DpEntry.zeros(), DpStay.zeros();

    cube dTdVInv(cCheck, cCheck, nCheck);
    for(int n = 0; n<nCheck;n++) {
        dTdVInv.slice(n) = inv(rho * Pi % repmat(pStay.row(n) - pEntry.row(n + 1), cCheck,1)  - eye(cCheck, cCheck)).t();
    }

    // let's take care of the k's
    for(int kX=nCheck - 1; kX >= 0; kX--) {

        mat dvalueAdditionalEntrydK(1, cCheck);
        dvalueAdditionalEntrydK.zeros();

        // we know that all derivatives for kX<n are zero (because you only get there through a zero continuation value)
        for(int n = kX; n >= 0; n--) {

            // partial derivative of the operator with respect to k
            mat dTdK(1, cCheck);
            if(n == kX) {
                dTdK = rho * (exp(logGrid) / (n + 1) - DpEntry.slice(kX).row(n + 1) % vS.row(n) + dvalueAdditionalEntrydK) * Pitranspose;
            } else {
                dTdK = rho * (-DpEntry.slice(kX).row(n + 1) % vS.row(n) + dvalueAdditionalEntrydK) * Pitranspose;
            }

            // implicit function theorem
            DvS.slice(kX).row(n) = - dTdK * dTdVInv.slice(n) ;

            // compute dPentry
            DpEntry.slice(kX).row(n) = mynormpdf(log(vS.row(n)) - log((1 + phi(n))), w_mean, w_std) / vS.row(n) % DvS.slice(kX).row(n);
            DpStay.slice(kX).row(n) = mynormpdf(log(vS.row(n)), w_mean, w_std) / vS.row(n) % DvS.slice(kX).row(n);

            // piece together derivatives from additional entry
            dvalueAdditionalEntrydK = dvalueAdditionalEntrydK
            + (DpEntry.slice(kX).row(n)- DpEntry.slice(kX).row(n + 1)) % vS.row(n)
            + (pEntry.row(n)-pEntry.row(n + 1)) % DvS.slice(kX).row(n);

        }

    }

    // let's take care of the phi's
    for(int phiX=2*nCheck - 1; phiX >= nCheck; phiX--) {

        mat dvalueAdditionalEntrydPhi(1, cCheck);
        dvalueAdditionalEntrydPhi.zeros();

        // we know that all derivatives for phi(nPrime) with nPrime<n are zero (because you only get there through a zero continuation value)
        for(int n = (phiX-nCheck); n >= 0; n--) {

            // partial derivative of the operator with respect to phi
            mat dTdPhi(1, cCheck);
                dTdPhi = rho * (- DpEntry.slice(phiX).row(n + 1) % vS.row(n) + dvalueAdditionalEntrydPhi) * Pitranspose;

            // implicit function theorem
            DvS.slice(phiX).row(n) = - dTdPhi * dTdVInv.slice(n);

            // compute dPentry
            if(n == phiX-nCheck) {
                DpEntry.slice(phiX).row(n) = mynormpdf(log(vS.row(n)) - log((1 + phi(n))), w_mean, w_std)
                    % (DvS.slice(phiX).row(n)/vS.row(n) - 1 / (1 + phi(n)));
            } else {
                DpEntry.slice(phiX).row(n) = mynormpdf(log(vS.row(n)) - log((1 + phi(n))), w_mean, w_std)
                    % (DvS.slice(phiX).row(n)/vS.row(n));
            }

            DpStay.slice(phiX).row(n) = mynormpdf(log(vS.row(n)), w_mean, w_std)
                    % DvS.slice(phiX).row(n)/vS.row(n);

            // piece together derivatives from additional entry
            dvalueAdditionalEntrydPhi = dvalueAdditionalEntrydPhi
                + (DpEntry.slice(phiX).row(n)- DpEntry.slice(phiX).row(n + 1)) % vS.row(n)
                + (pEntry.row(n)-pEntry.row(n + 1)) % DvS.slice(phiX).row(n);

        }

    }

    // let's take care of the omega
    mat dvalueAdditionalEntrydOmega(1, cCheck);
    dvalueAdditionalEntrydOmega.zeros();

    int omegaX = 2*nCheck;

    for(int n = nCheck - 1; n >= 0; n--) {

        // partial derivative of the operator with respect to phi
        mat dTdOmega(1, cCheck);
        dTdOmega = rho * (
                omega * vS.row(n) %  mynormpdf(log(vS.row(n)), w_mean, w_std) % (- 1.0 / pow(omega, 2) * log(vS.row(n)) + 0.5)
                - DpEntry.slice(omegaX).row(n + 1) % vS.row(n) + dvalueAdditionalEntrydOmega
                + omega * mynormpdf(- log(vS.row(n)), w_mean, w_std) % (0.5 + 1.0 / pow(omega, 2) * log(vS.row(n)))) * Pitranspose;

        // implicit function theorem
        DvS.slice(omegaX).row(n) = - dTdOmega * dTdVInv.slice(n);

        // compute dPentry
        DpEntry.slice(omegaX).row(n) = omega * mynormpdf(log(vS.row(n)) - log((1 + phi(n))), w_mean, w_std) % (- 1.0 / pow(omega, 2) * log(vS.row(n))
            + 1.0 / omega * DvS.slice(omegaX).row(n)/vS.row(n) + log(1 + phi(n))/pow(omega, 2) + 0.5);

        // compute dPstay
        DpStay.slice(omegaX).row(n) = omega * mynormpdf(log(vS.row(n)), w_mean, w_std) % (- 1.0 / pow(omega, 2) * log(vS.row(n))
            + 1.0 / omega * DvS.slice(omegaX).row(n)/vS.row(n) + 0.5);

        // piece together derivatives from additional entry
        dvalueAdditionalEntrydOmega = dvalueAdditionalEntrydOmega
            + (DpEntry.slice(omegaX).row(n)- DpEntry.slice(omegaX).row(n + 1)) % vS.row(n)
            + (pEntry.row(n)-pEntry.row(n + 1)) % DvS.slice(omegaX).row(n);

     }


    // let's take care of the mu
    mat dvalueAdditionalEntrydMu(1, cCheck);
    dvalueAdditionalEntrydMu.zeros();

    int muX = 2*nCheck + 1;
    mat dPi = dPidMu(Pi, logGrid, Param["mu"], Param["sigma"]);

    for(int n = nCheck - 1; n >= 0; n--) {

        // partial derivative of the operator with respect to mu
        mat dTdMu(1, cCheck);
        dTdMu = rho * (-DpEntry.slice(muX).row(n + 1) % vS.row(n) + dvalueAdditionalEntrydMu) * Pitranspose +
                vS.row(n)  * inv(Pitranspose) * dPi.t();

        // implicit function theorem
        DvS.slice(muX).row(n) = - dTdMu * dTdVInv.slice(n);

        // compute dPentry
        DpEntry.slice(muX).row(n) = mynormpdf(log(vS.row(n)) - log((1 + phi(n))), w_mean, w_std) / vS.row(n) % DvS.slice(muX).row(n);
        DpStay.slice(muX).row(n) = mynormpdf(log(vS.row(n)), w_mean, w_std) / vS.row(n) % DvS.slice(muX).row(n);

        // piece together derivatives from additional entry
        dvalueAdditionalEntrydMu = dvalueAdditionalEntrydMu
          + (DpEntry.slice(muX).row(n)- DpEntry.slice(muX).row(n + 1)) % vS.row(n)
          + (pEntry.row(n)-pEntry.row(n + 1)) % DvS.slice(muX).row(n);

     }

    // let's take care of the sigma
    mat dvalueAdditionalEntrydSigma(1, cCheck);
    dvalueAdditionalEntrydSigma.zeros();

    int sigmaX = 2*nCheck+2;
    dPi = dPidSigma(Pi, logGrid, Param["mu"], Param["sigma"]);

    for(int n = nCheck - 1; n >= 0; n--) {

        // partial derivative of the operator with respect to sigma
        mat dTdSigma(1, cCheck);
        dTdSigma = rho * (-DpEntry.slice(sigmaX).row(n + 1) % vS.row(n) + dvalueAdditionalEntrydSigma) * Pitranspose +
                vS.row(n)  * inv(Pitranspose) * dPi.t();

        // implicit function theorem
        DvS.slice(sigmaX).row(n) = - dTdSigma * dTdVInv.slice(n);

        // compute dPentry
        DpEntry.slice(sigmaX).row(n) = mynormpdf(log(vS.row(n)) - log((1 + phi(n))), w_mean, w_std) / vS.row(n) % DvS.slice(sigmaX).row(n);
        DpStay.slice(sigmaX).row(n) = mynormpdf(log(vS.row(n)), w_mean, w_std) / vS.row(n) % DvS.slice(sigmaX).row(n);

        // piece together derivatives from additional entry
        dvalueAdditionalEntrydSigma = dvalueAdditionalEntrydSigma
          + (DpEntry.slice(sigmaX).row(n)- DpEntry.slice(sigmaX).row(n + 1)) % vS.row(n)
          + (pEntry.row(n)-pEntry.row(n + 1)) % DvS.slice(sigmaX).row(n);

     }

     GradientStructure Gradient;
     Gradient.DpEntry = DpEntry;
     Gradient.DpStay = DpStay;
     Gradient.DvS = DvS;
     return(Gradient);
}

//' Compute the gradient of the mixing density for a variable of interest
//'
//' This function computes the gradient of the mixing density. Note that all
//' parameters (except for omega) enter the mixing density only through the value functions. The
//' derivative of the value function is given as argument to the function.
//'
//' @param Settings is an R list of settings
//' @param Param is an R list of parameters
//' @param Eq is an R list of equilibrium objects
//' @param D_vS is a matrix of the gradient of the value function with respect to the variable of interest
//' @param OMEGA_FLAG is a flag whether the variable in question is omega or not
//' @return A cube of the mixing density \code{p^{ - 1}(a_S, c, n)}
// [[Rcpp::export]]
cube mixingDensityGradient(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Eq, mat D_vS, bool OMEGA_FLAG=false) {

    int nCheck = Settings["nCheck"];
    int cCheck = Settings["cCheck"];

    vec w = Settings["integrationWeights"];
    vec p = Settings["integrationNodes"];

    mat vS = Eq["vS"];

    // nchoosekmatrix
    mat nchoosekMatrix(nCheck + 1, nCheck + 1);
    for(int nX=0;nX<nCheck + 1;nX++) {
        for(int jX=0;jX<nCheck + 1;jX++) {
            nchoosekMatrix(nX, jX) = nchoosek(nX, jX);
        }
    }

    cube  mixingDensity(p.n_elem, cCheck, nCheck);
    cube  aSinv(p.n_elem, cCheck, nCheck);
    cube  daSinvdP(p.n_elem, cCheck, nCheck);

    mixingDensity.zeros(); aSinv.zeros(); daSinvdP.zeros();
    cube  D_mixingDensity(p.n_elem, cCheck, nCheck);
    cube  D_aSinv(p.n_elem, cCheck, nCheck);
    cube  D_daSinvdP(p.n_elem, cCheck, nCheck);
    D_mixingDensity.zeros(); D_aSinv.zeros(); D_daSinvdP.zeros();


    mat intWeights = repmat(w, 1, cCheck);
    double w_mean = -.5*pow((double) Param["omega"], 2);
    double w_std = Param["omega"];
    double omega = Param["omega"];

    // mixing weights

    for(int n = 2; n <= nCheck; n++) {

        mat expaSInv(p.n_elem, cCheck);
        expaSInv.zeros();
        mat dexpaSInvdP(p.n_elem, cCheck);
        dexpaSInvdP.zeros();

        mat D_expaSInv(p.n_elem, cCheck);
        D_expaSInv.zeros();
        mat D_dexpaSInvdP(p.n_elem, cCheck);
        D_dexpaSInvdP.zeros();

        for(int nPrime = 1; nPrime<=n; nPrime++) {

            int nCk = nchoosekMatrix(n - 1, nPrime - 1);
            mat binomialPmf     = nCk * repmat(pow(p, nPrime - 1) % pow(1-p, n-nPrime), 1, cCheck);
            mat dbinomialPmfdP  = nCk * repmat(pow(p, nPrime-2) % pow(1-p, n-nPrime - 1) % ((nPrime - 1) * (1-p) - p * (n-nPrime)), 1, cCheck);
            mat repvS =  repmat(vS.row(nPrime - 1), p.n_elem,1);
            mat D_repvS =  repmat(D_vS.row(nPrime - 1), p.n_elem,1);
            expaSInv +=  binomialPmf % repvS ;
            dexpaSInvdP += dbinomialPmfdP % repvS;
            D_expaSInv +=  binomialPmf % D_repvS ;
            D_dexpaSInvdP += dbinomialPmfdP % D_repvS;

        }

        aSinv.slice(n - 1) =  log (expaSInv);
        D_aSinv.slice(n - 1) =  D_expaSInv / expaSInv;
        daSinvdP.slice(n - 1) =  dexpaSInvdP / expaSInv;
        D_daSinvdP.slice(n - 1) =  D_dexpaSInvdP / expaSInv - dexpaSInvdP % D_expaSInv / pow(expaSInv, 2);

        mat normaSinv = mynormpdf(aSinv.slice(n - 1), w_mean, w_std);

        mat D_normaSinv(p.n_elem, cCheck);

        if(OMEGA_FLAG) {
            D_normaSinv = mynormpdf(aSinv.slice(n - 1), w_mean, w_std) % (
                         - 1 / omega
                         - (aSinv.slice(n - 1)/omega + .5 * omega) % (- aSinv.slice(n - 1)/pow(omega, 2) + D_aSinv.slice(n - 1)/omega + 0.5));
        } else {
            D_normaSinv = - mynormpdf(aSinv.slice(n - 1), w_mean, w_std) % (aSinv.slice(n - 1)/omega + .5 * omega) % D_aSinv.slice(n - 1)/omega ;
        }

        mixingDensity.slice(n - 1) = daSinvdP.slice(n - 1) % normaSinv % intWeights;

        D_mixingDensity.slice(n - 1) = D_daSinvdP.slice(n - 1) % normaSinv % intWeights
                                    + daSinvdP.slice(n - 1)  % D_normaSinv % intWeights;

    }

    return(D_mixingDensity);

}


//' Compute the gradient of the likelihood function
//'
//' This function computes the gradient of the likelihood function.
//' It also computes the equilibrium and the likelihood function and
//' then returns a list with the likelihood and gradient contributions.
//'
//' @param Settings is an R list of settings
//' @param Param is an R list of parameters
//' @param Data an R list of the data with matrices \code{C} and \code{N}
//' @return An R list with a vector of likelihood contributions and a matrix of
//' gradient contributions.
// [[Rcpp::export]]
List likelihoodGradient(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Data) {

    int nCheck = Settings["nCheck"];
    int cCheck = Settings["cCheck"];

    vec w = Settings["integrationWeights"];
    vec p = Settings["integrationNodes"];

    Rcpp::List Eq = valuefunction(Settings, Param);
    GradientStructure Gradient = valuefunctionGradient(Settings, Param, Eq);
    colvec likelihoodContributions = likelihood(Settings, Param, Data);

    cube DvS = Gradient.DvS;
    cube DpEntry = Gradient.DpEntry;
    cube DpStay = Gradient.DpStay;

    mat C = Data["C"];
    mat N = Data["N"];

    int M = C.n_rows;
    int T = C.n_cols;

    // nchoosekmatrix
    mat nchoosekMatrix(nCheck + 1, nCheck + 1);
    for(int nX=0;nX<nCheck + 1;nX++) {
        for(int jX=0;jX<nCheck + 1;jX++) {
            nchoosekMatrix(nX, jX) = nchoosek(nX, jX);
        }
    }

    cube binomialProbabilities(p.n_elem, nCheck + 1, nCheck + 1);
    binomialProbabilities.zeros();
    for(int nX=2; nX < nCheck + 1; nX++) {
        for(int nXPrime=0; nXPrime <= nX; nXPrime++) {
            for(int pX=0; pX < (int) p.n_elem; pX++) {
                binomialProbabilities(pX, nX, nXPrime) = nchoosekMatrix(nX, nXPrime) * pow(p(pX), nXPrime) * pow(1-p(pX), nX-nXPrime);
            }
        }
    }


    // gradients
    mat gradientContributions((T - 1)*M, 2*nCheck+3);
    gradientContributions.zeros();

    for(int kX=0;kX<2*nCheck+3;kX++) {

        cube DmixingDensity(p.n_elem, cCheck, nCheck);

        if(kX==2*nCheck) {
            DmixingDensity = mixingDensityGradient(Settings, Param, Eq, DvS.slice(kX), true);
        } else {
            DmixingDensity = mixingDensityGradient(Settings, Param, Eq, DvS.slice(kX), false);
        }

        int jX=0;

        // loop over all markets
        for(int mX=0;mX<M;mX++) {

            // loop over time
            for(int tX=1;tX<T;tX++) {

                // entry
                if (N(mX, tX) > N(mX, tX - 1)) {
                    gradientContributions(jX, kX) = DpEntry(N(mX, tX) - 1, C(mX, tX) - 1, kX) - DpEntry(N(mX, tX), C(mX, tX) - 1, kX);
                }

                // 0 and no entry
                if (N(mX, tX) == 0 && N(mX, tX - 1) == 0) {
                    gradientContributions(jX, kX) = -DpEntry(0, C(mX, tX) - 1, kX) ;
                }

                // number of firms stays the same (yet >0)
                if (N(mX, tX) > 0 && N(mX, tX) == N(mX, tX - 1)) {
                    gradientContributions(jX, kX) = DpStay(N(mX, tX) - 1, C(mX, tX) - 1, kX) - DpEntry(N(mX, tX), C(mX, tX) - 1, kX);
                }

                // all firms leave
                if (N(mX, tX) == 0 && N(mX, tX - 1) > 0) {
                    gradientContributions(jX, kX) =  -DpStay(0, C(mX, tX) - 1, kX);
                }

                // mixing
                if (N(mX, tX - 1) > 1 && N(mX, tX)<=N(mX, tX - 1)) {
                    for(int pX=0; pX < (int) p.n_elem; pX++) {
                        gradientContributions(jX, kX) -= binomialProbabilities(pX, N(mX,tX - 1),  N(mX,tX))
                            * DmixingDensity(pX, C(mX,tX) - 1, N(mX,tX - 1) - 1);
                    }

                }

                jX++;

            }

        }
    }

    return List::create(_["gradientContributions"]         = gradientContributions,
                        _["likelihoodContributions"]     = likelihoodContributions);
}


//' Compute the gradient of the transition probability matrix with respect to mu
//'
//' @param Pi is the transition probability matrix
//' @param logGrid is the demand grid
//' @param mu is the mean of the innovations times 100
//' @param sigma is the standard deviation of the innovations times 100
//' @return A matrix with the derivative of the transition probability matrix with respect to mu.
// [[Rcpp::export]]
mat dPidMu(mat Pi, rowvec logGrid, double mu, double sigma) {
  mu = mu/100;
  sigma = sigma/100;
  mat dPi(Pi.n_cols,Pi.n_cols);
  double d = logGrid(1) - logGrid(0);
  for(int iX = 0; iX < (int) Pi.n_cols; iX++) {
    dPi(iX,0) = -mynormpdf(logGrid(0) - logGrid(iX)  + d/2, mu, sigma);
    for(int jX = 1; jX < (int) Pi.n_cols - 1; jX++) {
      dPi(iX,jX) = -mynormpdf(logGrid(jX) - logGrid(iX)  + d/2, mu, sigma)
                   +mynormpdf(logGrid(jX) - logGrid(iX) - d/2, mu, sigma);
    }
    dPi(iX,Pi.n_cols - 1) = mynormpdf(logGrid(Pi.n_cols - 1) - logGrid(iX) - d/2, mu, sigma);
  }
  return(dPi);
}


//' Compute the gradient of the transition probability matrix with respect to sigma
//'
//' @param Pi is the transition probability matrix
//' @param logGrid is the demand grid
//' @param mu is the mean of the innovations
//' @param sigma is the standard deviation of the innovations times 100
//' @return A matrix with the derivative of the transition probability matrix with respect to sigma.
// [[Rcpp::export]]
mat dPidSigma(mat Pi, rowvec logGrid, double mu, double sigma) {
  mu = mu/100;
  sigma = sigma/100;
  mat dPi(Pi.n_cols,Pi.n_cols);
  double d = logGrid(1) - logGrid(0);
  for(int iX = 0;iX < (int) Pi.n_cols; iX++) {
    dPi(iX,0) = -mynormpdf(logGrid(0) - logGrid(iX)  + d/2, mu, sigma) * (logGrid(0) - logGrid(iX)  + d/2 - mu)/sigma;
    for(int jX = 1; jX < (int) Pi.n_cols - 1; jX++) {
      dPi(iX,jX) = -mynormpdf(logGrid(jX) - logGrid(iX)  + d/2, mu, sigma) * (logGrid(jX) - logGrid(iX)  + d/2 - mu)/sigma
                   +mynormpdf(logGrid(jX) - logGrid(iX) - d/2, mu, sigma) * (logGrid(jX) - logGrid(iX)  - d/2 - mu)/sigma;
    }
    dPi(iX,Pi.n_cols - 1) = mynormpdf(logGrid(Pi.n_cols - 1) - logGrid(iX) - d/2, mu, sigma) * (logGrid(Pi.n_cols - 1) - logGrid(iX) - d/2 - mu)/sigma;
  }
  return(dPi);
}
