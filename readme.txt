A temporary guide for using the compMS package


This package is still under development and both the functions and output are subject to change.
We apologize for the sorry state of the documentation and the numerous data fields that seemingly
do nothing.  Someday they will do something amazing! (maybe)

The installation comes with the file, sampleDat.rda.  We recommend looking at this file to see how
 your data needs to be formatted.  An explanation of the data format can be seen by entering
?ptmDat into the R terminal (after installing and loading the package).

The most important function in the package is the compBayes() function.  Calling this on your data
will fit a composititional Bayesian model to a dataset in a form similar to that in sampleDat.rda.
The results will return a list object where the first component contains a table with a summary of
the proportions for each specified condition.  The second component contains a table with a summary
of proportions for each specified biological replicate (assuming this has been specified).  The
third component contains the Stan model object. Starting from component 4, the list object contains
summary information on differential expression (log2 fold-changes).  A separate table is given for
each condition.  So if ten conditions were specified then RES[[4]] - RES[[14]] would contain the
summary data.

After instaling and loading the package, information about this
function and the parameters it takes can be seen by entering ``?compBayes'' into the R terminal
(without the quotation marks).

There are also two plotting functions built into the package caterPlot() and precisionPlot().
Again, more information can be obtained with ?caterpillar and ?precisionPlot.

As an example, these functions can be used as follows.

RES <- compBayes(sampleDat)

View(RES[[1]])
View(RES[[2]])
View(RES[[4]])
View(RES[[5]])
View(RES[[6]])

caterPlot(RES, byCond = TRUE)

precisionPlot(RES)
