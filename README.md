# BayesMonophyly
This script will perform Bayesian monophyl test from MrBayes or BEAST tree files.

Due to sample of trees from MCMC, so that frequency of topology equals probability of that topology, this is done rather easilly.
Script takes on the input sample of trees from MrBayes or BEAST analysis (remove burnin) and then check the number of trees
that contains monophyly of specified taxons against trees that do not. This is then corrected by prior tree frequency.
