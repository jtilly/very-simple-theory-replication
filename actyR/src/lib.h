#ifndef lib_h
#define lib_h

struct GradientStructure {
    cube DpEntry, DpStay, DvS;
};

rowvec mynormcdf(rowvec x, double mean, double sigma); 
double mynormpdf(double x, double mean, double sigma);
mat mynormpdf(mat x, double mean, double sigma);
int nchoosek(int n, int k);

#endif
