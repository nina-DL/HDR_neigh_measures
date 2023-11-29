# Alternative Approaches for Computing Highest-Density Regions

This repository contains four main R scripts for reproducing the results of the paper Alternative Approaches for Computing Highest-Density Regions. These are summarized as follows:

1. "1 Simulated_scenarios.R": defines the different scenarios evaluated in the paper and simulates corresponding data.
2. "2 HDR_functions.R": defines all the different functions used in the paper including a) the generalized approach for computing HDRs, b) the different evaluated measures, c) the different performance metrics
3. "3 HDR_AllMetrics_Simulations.R": intensive simulation setup based on 1. (which generates scenarios as in "Def_scenarios.RData") and 2.
4. "4 HDR_AllMetrics_RealData.R": evaluation of the proposed metrics on real data (made available in the "magic04.data" file).
5. "Paper_figures.R": reproduces all figures of the main paper (based on "Def_scenarios.RData" and results given in "Error Evaluation" folder). 
