#This repository contains the simulation framework used to evaluate MERLIN method and other standard MR methods.


1_generates phenotypes.R — Generates exposure and outcome phenotypes under predefined genetic architectures using real genotype data.

2.1_sim.R — Implements all MR and MERLIN analyses under non-overlapping exposure and outcome GWAS samples.

2.2_sim_overlap.R — Extends the analyses to scenarios with sample overlap and incorporates overlap correlation corrections.

3_Figure2plot.R — This code reads simulated estimation of beta_A results and produces Figure 2, comparing estimator bias, type-I error, and power for MERLIN method and standard MR methods.

4_Figure3plot.R - This code reads simulated estimation of beta_I results and produces Figure 3, comparing estimator bias, type-I error, and power for MERLIN method and standard MR methods.

5_Figureqqplot.R - This code generates QQ plots using simulated estimation results to assess inflation and calibration o
