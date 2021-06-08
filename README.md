# Zero modified models
The R code is for the article **"A New Algorithm for Zero-Modified Models Applied to Citation Counts"**. This article proposes a new algorithm to fit zero-modified versions of **discretised log-normal**, **hooked power-law** and **Weibull** models to citation count data from 23 different Scopus categories from 2012. The new algorithm allows the standard errors of all parameter estimates to be calculated, and hence also confidence intervals and p-values. This algorithm can also estimate both negative and positive zero-modification parameters corresponding to zero-deflation and zero-inflation in the datasets.

There are three R files (ZMDLN.R, ZMHPL.R, and ZMWeibull.R) corresponding to the zero-modified models of the three mentioned distributions. The data used for running the R code is also provided (Data_new.xlsx). This excel file has 23 columns that each column is the citation counts relative to each Scopus subject. The list of the names of Scopus subjects is also mentioned in the R files. In the following, the different steps of the algorithm for one of the three distribution for example, discretised log-normal, are explained. For another two distributions, all steps will be the same. 


The algorithm for fitting the zero-modified discretised log-normal has the following steps:

1-	Reading the data and then for each subject (each colum of the excel file), the algorithm will do the following steps.

2-	Providing the zero-modified mass probability function of discretised log-normal. 

3-	Computation of the log-likelihood function of the mentioned mass function.

4-	Using “optim” function to maximise the mentioned log-likelihood function to estimate the parameters.

5-	Extracting the estimates of parameters based on the results of **“optim”** function.

6-	Computation of the standard deviation of the estimates based on the Hessian matrix returned by “optim” function, then the confidence intervals can be obtained.

7-	Computation of Akaike Information Criterion based on the log-likelihood value returned by the “optim” function.

8-	Providing Wald hypothesis tests related to the parameters, particularly for the zero-modification parameter to determine whether there is statistical evidence of zero-modification in the data.

