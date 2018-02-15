# DML-IVg-estimation
G-estimator for the linear IV model for trials with non-compliance, with treatment effect heterogeneity, using Super learner for nuisance fits

This is an R function that implements the g-estimator for a linear IV model with treatment effect heterogeneity as detailed in 
DiazOrdaz K., Daniel R., Kreif, N. (2018) Data-adaptive doubly robust instrumental variable methods for treatment effect heterogeneity. 
arXiv e-print (arXiv:1802.02821)

The user can specify whether the non-adherence is onesided or two-sided,  treatment effect modification is of interest or not
and whether parametric or data-adaptive fits (using the Super Learner SL) are to be used.

If SL is used, the user must specify separately a library of learners suitable for the type of exposure and instrument 
for the SL to use.

In the paper, we assume both exposure and instrument are binary, but if the user is happy to assume linearity on the exposure holds
for the IV model, the methods extends to continuous exposure. See Vanteelandt and Didelez 	
Robustness and efficiency of covariate adjusted linear instrumental variable estimators eprint arXiv:1510.01770
