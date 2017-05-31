#include <cmath>
#include <math.h> 
#include <actyR.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "lib.h"

// define a few constants
double root2pi_1 = 0.398942280401432677939946059934381868;
double root2_1 = 1.0/sqrt(2.0);

// normal cdf
rowvec mynormcdf(rowvec x, double mean, double sigma) {    
    rowvec z(x.n_elem);
    z = root2_1*(x-mean)/sigma;
    for(size_t jX=0;jX<x.n_elem;jX++) {
        z(jX) = 0.5 *(1.0 + erf(z(jX)));
    }
    return(z);
} 


double mynormpdf(double x, double mean, double sigma) {
    double z = pow((x-mean)/sigma, 2);
    return(root2pi_1 * exp(-0.5 * z)/sigma);
}


// normal pdf
mat mynormpdf(mat x, double mean, double sigma) {    
    mat z(x.n_rows, x.n_cols);
    z = pow((x-mean)/sigma, 2);
    return(root2pi_1 * exp(-0.5 * z)/sigma);
}

// nchoosek
int nchoosek(int n, int k) {
    if (k == 0)  
        return(1);
    return (n * nchoosek(n - 1, k - 1)) / k;
}
