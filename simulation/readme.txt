#This repository contains the simulation framework used to evaluate MERLIN method and other standard MR methods.


1_generates phenotypes.R — Generates exposure and outcome phenotypes under predefined genetic architectures using real genotype data.

2.1_sim.R — Implements all MR and MERLIN/CHESS analyses under non-overlapping exposure and outcome GWAS samples.

2.2_sim_overlap.R — Extends the analyses to scenarios with sample overlap and incorporates overlap-induced correlation corrections.

3_Figure2plot.R — This code generates Figure 2 by summarizing simulation results, producing boxplots of causal effect estimates, evaluating type I error rates, and plotting statistical power across methods under varying genetic architectures.

4_Figure3plot.R - This code reads simulated estimation results and produces Figure 3, comparing estimator bias, type-I error, and power for MERLIN method and standard MR methods across varying genetic architectures.

5_Figureqqplot.R - This code generates QQ plots using simulated estimation results to assess inflation and calibration of p-values for MERLIN method and standard MR methods under sample-overlap scenarios across different genetic architectures.