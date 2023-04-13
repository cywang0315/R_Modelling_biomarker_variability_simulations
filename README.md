
Files in this repository are R codes for simulation studies in our paper submitted to Biostatistics.

Article: Modelling Biomarker Variability in Joint Analysis of Longitudinal and Time-to-Event Data

Authors: Chunyu Wang, Jiaming Shen, Christiana Charalambous, and Jianxin Pan.

This folder contains R commands and functions required to run the simulations (case 1 to case 4) presented in the paper. 
Four main R scripts are listed: 'simulation_case1.R', 'simulation_case2.R', 'simulation_case3.R', and 'simulation_case4.R' 
along with their source files.

The source files are similar except for the function 'gendat' for generating the simulated data in difference cases.
Case 3 and Case 4 involve an additional function 'knot.sel' for selecting the best knot sequence among the candidate knot patterns.

Due to data protection and confidentiality, the code for analysing the MRC trial is not included. Case 2 is designed to mimic
the MRC dataset. The corresponding R file 'simulation_case2.R' and its source file illustrate the procedure for MRC analysis.
