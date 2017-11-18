# pvf-frailty-binary-regression
Example data set and R code for analysis using PVF frailty model

This repository contains a CSV file with simulated data, R code for its analysis, and results (estimated parameters and standard errors) in a PDF file. 

To perform the analysis, simply place the R file and the CSV in the same directory and run the R code. Required packages are lme4 (version 1.1-13) and BHH2 (version 2016.05.31). The version of R used to perform the analyses was 3.4.1.

The example data set contains an id number for each subject, 4 explanatory variables (age appears twice: on its original scale and transformed), and an outcome (denoted y) for each of 6 timepoints (coded as an integer 1-6). 

The covariates are age (continuous; centered and scaled by 10 years), race (three categories), education (four categories), and treatment (two categories).
