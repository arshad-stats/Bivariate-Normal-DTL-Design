# Bivariate-Normal-DTL-Design
**R code for the work: Estimation Following Surrogate Based Selection in Multi-Arm Drop-the-Losers Designs**

This repository contains the R scripts required to reproduce the simulation studies and real-data applications presented in the manuscript. It provides an estimation framework for two-stage drop-the-losers (DTL) designs assuming bivariate normal responses with a common but unknown covariance matrix. Treatment selection at the first stage is based on a surrogate endpoint, while final inference is conducted on the primary endpoint.

## Dependencies

The code is written in R. To run the scripts, you will need to install the following R packages:

```R
install.packages(c("MASS", "ggplot2", "dplyr", "tidyr", "scales", "MVN", "biotools"))
