all: help
install: .install
replication: estimation counterfactuals

help:
	@echo "----------------------------------------------------------------"
	@echo "Replicate Very Simple-Markov Perfect Industry Dynamics: Theory  "
	@echo "----------------------------------------------------------------"
	@echo ""
	@echo " make install          --  install the necessary R packages     "
	@echo " make replication      --  replicate all results                "
	@echo ""
	@echo " make estimation       --  estimate model                       "
	@echo " make counterfactuals  --  compute counterfactuals              "
	@echo ""

.install:
	Rscript -e 'if(!"Rcpp" %in% installed.packages()[,1]) install.packages("Rcpp", repos = "http://cran.us.r-project.org")'
	Rscript -e 'if(!"RcppArmadillo" %in% installed.packages()[,1]) install.packages("RcppArmadillo", repos = "http://cran.us.r-project.org")'
	Rscript -e 'if(!"gaussquad" %in% installed.packages()[,1]) install.packages("gaussquad", repos = "http://cran.us.r-project.org")'
	Rscript -e 'if(!"foreach" %in% installed.packages()[,1]) install.packages("foreach", repos = "http://cran.us.r-project.org")'
	Rscript -e 'if(!"tictoc" %in% installed.packages()[,1]) install.packages("tictoc", repos = "http://cran.us.r-project.org")'
	Rscript -e 'if(!"expm" %in% installed.packages()[,1]) install.packages("expm", repos = "http://cran.us.r-project.org")'
	Rscript -e 'if(!"doParallel" %in% installed.packages()[,1]) install.packages("doParallel", repos = "http://cran.us.r-project.org")'
	R CMD INSTALL actyR --preclean
	touch .install

estimation: .install
	Rscript "estimation.R"

counterfactuals: .install
	Rscript "counterfactuals.R"

clean:
	rm -f .install

.PHONY: clean  \
	install \
	replication \
	estimation \
	counterfactuals \
	help
