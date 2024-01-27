# BayesPGLS: Bayesian Phylogenetic Generalized Least Squares.

BayesPGLS is an R package to carry out Bayesian PGLS that allows including weights.

## What's in BayesPGLS?

The package includes functions to run PGLS and to visualize the results. It includes:  

- Ability to run PGLS vs GLS and use model choice to test the importance of the phylogeny;
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
