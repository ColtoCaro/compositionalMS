# Installation

In R, you can use `dev_tools():`

```
devtools::install_github("https://github.com/ColtoCaro/compositionalMS" , INSTALL_opts="--no-staged-install")
```

Or, you can pull the source and install on a command line:

```
cd ~/workspace
git clone git@github.com:ColtoCaro/compositionalMS
R CMD INSTALL compositionalMS --no-staged-install
```

# Description

A temporary guide for using the compMS package

This package fits Bayesian compostional models to mass spectrometry proteomics data utilizing isobaric tags.
A thorough description of the modelling can be found in the manuscript
https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00699

This package is still under development and both the functions and output are subject to change.
We apologize for the sorry state of the documentation and the numerous data fields that seemingly
do nothing.  We're working on it.

Quick start of using the package in R-Studio (or R terminal) can be begun via following commands:
```
> library(devtools)
> install_github('coltocaro/compMS')
> library(compMS)
```
WARNING:  rstan can be tricky to install on many systems.  See this thread,
https://discourse.mc-stan.org/t/rstan-install-problem-on-mac-nojave/6022
from the Stan Forums on some of the difficulties with Stan 2.17.4 (these problems should
be resolved with Version 2.18.1).  Often entering `Sys.setenv(USE_CXX14=1)` is enough to fix
the problem.

The installation comes with the file, sampleDat.rda.  We recommend looking at this file to see how
 your data needs to be formatted.  An explanation of the data format can be seen by entering
`?sampleDat` into the R terminal (after installing and loading the package).

The most important function in the package is the `compBayes()` function.  Calling this on your data
will fit a composititional Bayesian model to a dataset in a form similar to that in `sampleDat.rda.`
The results will return a list object where the first component contains a table with a summary of
the proportions for each specified condition.  The second component contains a table with a summary
of proportions for each specified biological replicate (assuming this has been specified).  The
third component contains the Stan model object. Starting from component 4, the list object contains
summary information on differential expression (log2 fold-changes).  A separate table is given for
each condition.  So if ten conditions were specified then `RES[[4]] - RES[[14]]` would contain the
summary data.

After instaling and loading the package, information about this
function and the parameters it takes can be seen by entering `?compBayes` into the R terminal
(without the quotation marks).

There are also two plotting functions built into the package `catterPlot()` and `precisionPlot()`.
Again, more information can be obtained with `?caterpillar` and `?precisionPlot`.

As an example, these functions can be used as follows.
```
RES <- compBayes(sampleDat)

View(RES[[1]])
View(RES[[2]])
View(RES[[4]])
View(RES[[5]])
View(RES[[6]])

catterPlot(RES, byCond = TRUE, avgCond = TRUE)

precisionPlot(RES, avgCond = TRUE)
```