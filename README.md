# BayesPGLS: Bayesian Phylogenetic Generalized Least Squares.

BayesPGLS is an R package to carry out Bayesian phylogenetic generalized least squares, testing the level of Pagel's lambda and allowing to include weights.

## What's in BayesPGLS?

The package includes functions to run PGLS and to visualize the results. It includes:  

- Options to run PGLS vs GLS and use model choice to test the influence of the phylogeny, particularly at low levels of Pagel's lambda;
- Allows to include weights when using response variables that are calculated as summary statistics;
- Plots for trace convergence, posterior density, and diagnostics;
- Provides a list of potential outliers and influential observations;

## How to install BayesPGLS?
To install BayesPGLS from GitHub, type the following lines of code on the R console:

```R
# Install and load 'devtools':
install.packages("devtools")
library(devtools)

# Install BayesPGLS:
install_git("https://github.com/fercol/BayesPGLS", subdir = "pkg/")
```
