# BayesMonophyly
This python 2.7 script will perform Bayesian monophyl test from MrBayes or BEAST tree files. It is dependent on ete2.

Run as:
```python BayesMonophyl.py -s [species to test] -i [input files] -b [burnin for all files, 20% by default]```
and specify minimum of two species. You can also specify more files, such as several runs from MrBayes analysis (by default, MrBayes is running two runs).

_____

In Bayesian analysis, MCMC is sampling parametric space, not only searching for the best solution. After sufficient number of samples are taken and "warmum" phase is removed (burnin), this sample is equal to posterior distribution.
This mean that frequency of trees with certain topology, or certain class of topologies, equals to probability (evidence for, belief that it is true) of this topology, or class of topologies.

This means that we can look at trees and count those which contain monophyly of taxons of interest and compare their frequency against those, which do not contain this monophyly.

Note however that [Suchard et al. 2005][1] consider this approach naive as it does not correctly estimate (or take into account) error rate. For now, only simple Bayes factor (equation 2 from this work) is implemented. For this reason, take results of this script as preliminary results and test your hypothesis directly by running MrBayes (or similar software) with Stepping Stone sampling and with topology restricted to monophyly and another run with topology restricted to non-monophyly. Comparing these SS-obtained likelihoods will be much more accurate.

Have also on mind that your hypotheses have to make sense with respect to your sample. You can not test monophyly of H<sup>1</sup>: group A and sequence X against H<sup>0</sup>: group A; when you have general sample where group A does not have to by monophyletic, will produce incorrect results. I think this is obvious but [Bergsten at al. 2013][2] were able to write whole paper on it.

[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2005.00352.x/pdf "Suchard at al. 2005: Models for Estimating Bayes Factors with Applications to Phylogeny and Tests of Monophyly"

[2]: http://sysbio.oxfordjournals.org/content/early/2013/06/14/sysbio.syt029.full "Bergsten at al. 2013: Bayesian Tests of Topology Hypotheses with an Example from Diving Beetles"

TODO:
* Implement more correct Bayes Factor
* Implement non-ete2 solution for testing monophyly
