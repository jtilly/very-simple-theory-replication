# Very Simple Markov-Perfect Industry Dynamics: Theory
*Jaap Abbring, Jeff Campbell, Jan Tilly, Nan Yang*

## Introduction
In the conclusion of *Very Simple Markov-Perfect Industry Dynamics: Theory*, we refer to the computational results that we present in a companion paper, [*Very Simple Markov-Perfect Industry Dynamics: Empirics*](https://pure.uvt.nl/portal/files/16003770/2017_021.pdf). We claim that

> [o]ur companion paper also demonstrates the practicality of applying our model to such data by estimating Motion Picture Theaters' sunk costs and the toughness of competition between them. The model's maximum likelihood estimation requires calculating a separate equilibrium for each market in the data for each trial value of its parameters, but this required only about thirty minutes using two Intel Xeon E5-2699 v3 CPUs (released by Intel in 2014) on a single machine with C++ code. We were also able to conduct many policy experiments, which calculated the effects of large demand shocks and counterfactual competition policies. This experience leads us to conclude that structural investigations of oligopoly dynamics based on this paper's model can be done with few computational resources.

This package contains the code required to replicate and verify these claims. The code in this package generates all tables and figures contained in Sections 4.3 and 4.4 of *Very Simple Markov-Perfect Industry Dynamics: Empirics*. The interested reader is referred to that paper for details.

## Install and Run
The code is organized in an R package that we refer to as `actyR`. The R package contains C++ code that requires compilation. The package has the following dependencies (all available from [CRAN](https://cran.r-project.org)):
- `Rcpp`
- `RcppArmadillo`
- `tictoc`
- `gaussquad`
- `expm`
- `doParallel`

### Linux, Mac OS
`RcppArmadillo` requires a fairly recent compiler (gcc >= 4.6). Also, Mac users need to install a Fortran compiler that is available [here](http://r.research.att.com/libs/) (more [details](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)).

To install the `actyR` and all of its dependencies, open your terminal and run
```bash
make install
```
To replicate the estimation and to compute the counterfactuals, run:
```bash
make replication
```
The estimation and counterfactuals can also be replicated individually using:
```bash
make estimation
```
and
```bash
make counterfactuals
```

For maximum performance, set the compiler optimization flags to `-O3` prior to running `make install`. This can be done by setting the content of `~/.R/Makevars` to
```bash
CC = gcc
CXX = g++
FC = gfortran
CFLAGS = -O3
CXXFLAGS = -O3
```
You may need to create the file if it doesn't exist.

### Windows
Windows users need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) before they can install this package.

Install the R package and its dependencies manually (please adjust the path to the location of the folder that contains `actyR` on your system):
```r
# install dependencies
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("tictoc")
install.packages("gaussquad")
install.packages("expm")
install.packages("doParallel")

# install actyR
setwd("C:/path/to/folder_that_contains_actyR")
install.packages("actyR", repos = NULL, type = "source", INSTALL_opts = "--no-multiarch --preclean")
```

To run the estimation, open R and run
```r
source("estimation.R")
```
To run the counterfactual simulations, open R and run
```r
source("counterfactuals.R")
```

## Contents
This replication package is structured as follows:
```
- actyR: contains the R package with R and C++ source code and data sets
- actyR/data: contains the data used in the estimation
- actyR/inst/include: contains a useful header file
- actyR/man: contains the manual pages for the R package
- actyR/R: contains the R code
- actyR/src: contains the C++ code
- actyR/test: contains unit tests
- graphs: empty directory for graphs from counterfactuals
- results: empty directory for output produced by estimation
- tables: empty directory for tables with estimation results
- makefile: makefile to install package and run code
- estimation.R: contains the main file for the estimation of the model
- counterfactuals.R: contains the main file for the counterfactuals
```

## Run Time

On a Linux server with 36 cores (two Intel Xeon E5-2699 v3 CPUs), running the entire estimation procedure took 30 minutes. Computing the counterfactuals took another 30 minutes.

On a quad core Mac (Intel Core i7-3615QM), the entire estimation procedure took 2 hours 49 minutes. Computing the counterfactuals took 1 hour 53 minutes.

## R version, package versions, and compilers used

R and its packages evolve over time. Here, we provide the `sessionInfo()` from our Linux server, where we generated all of the results presented in *Very Simple Markov-Perfect Industry Dynamics: Empirics*.

```
> sessionInfo()
R version 3.3.3 (2017-03-06)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS release 6.8 (Final)

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] actyR_1.0          expm_0.999-1       Matrix_1.2-8       doParallel_1.0.10
 [5] iterators_1.0.8    foreach_1.4.3      tictoc_1.0         gaussquad_1.0-2   
 [9] orthopolynom_1.0-5 polynom_1.3-9     

loaded via a namespace (and not attached):
[1] tools_3.3.3      Rcpp_0.12.10     codetools_0.2-15 grid_3.3.3      
[5] lattice_0.20-34
```

We compiled the C++ code in our package with `gcc`:
```
$ gcc --version
gcc (GCC) 4.9.3
```

## Further Reading

In addition to this replication package, we maintain a teaching package with Matlab code at [http://verysimple.abbring.org](http://verysimple.abbring.org). The teaching package contains a simplified version of the code that is used in *Very Simple Markov-Perfect Industry Dynamics: Empirics* that is suitable for experimentation and teaching in a graduate IO course.
