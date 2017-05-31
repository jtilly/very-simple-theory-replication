#ifndef nfxp_h
#define nfxp_h

List valuefunction(Rcpp::List Settings, Rcpp::List Param);
colvec likelihood(Rcpp::List Settings, Rcpp::List Param, Rcpp::List Data);
cube computeMixingDensity(Rcpp::List Settings, Rcpp::List Param, mat vS);

#endif
