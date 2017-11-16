# pvf-frailty-binary-regression
Example data set and R code for analysis using PVF frailty model

This repository contains a CSV file with simulated data, R code for its analysis, and results (estimated parameters) in a Word file. To perform the analysis, simply place the R file and the CSV in the same directory and run the R code. Required packages are lme4 and BHH2. 
The example data set contains an id number for each subject, 4 explanatory variables, and an outcome (y) for each of 6 timepoints. 
The covariates are age (centered and scaled by 10 years), race (three categories), education (four categories), and treatment (two categories).
