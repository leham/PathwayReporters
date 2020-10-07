# PathwayReporters.jl
A Julia module for simulating noise on parameters in a stochastic model of gene expression that includes mRNA maturation, as well as protein dynamics.  The simulation code outputs the expected mean values, the noise parameters and theoretical noise strength, the sample means of the various species as well as the extrinsic noise strengths as determined by both pathway-reporter methods [1] and the dual-reporter method.  Each parameter has a fixed scaled Beta distribution with scaling to match the given mean, except for mRNA transcription, where the standard Beta distribution parameters can be chosen.  The code is easily adjusted to simulate different forms of noise (other than scaled Beta form) on the parameters.  The repository is provided in support of the pre-print “Pathway dynamics can delineate the sources of transcriptional noise in gene expression” below [1]. 

## Related Publications 
[1] L. Ham, M.Jackson & M. P. H. Stumpf, Pathway dynamics can delinate the sources of transcriptional noise in gene expression. (2020) *bioRxiv* [2020.09.30.319814](https://doi.org/10.1101/2020.09.30.319814).
