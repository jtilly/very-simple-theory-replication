#ifndef nfxpGradient_h
#define nfxpGradient_h

GradientStructure valuefunctionGradient(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Eq);
cube mixingDensityGradient(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Eq, mat D_vS, bool OMEGA_FLAG);
List likelihoodGradient(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Data);
mat dPidMu(mat Pi, rowvec logGrid, double mu, double sigma  );
mat dPidSigma(mat Pi, rowvec logGrid, double mu, double sigma  );

#endif
