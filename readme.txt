A temporary guide for using the compMS package


This package is still under development and both the functions and output are subject to change.
We apologize for the sorry state of the documentation and the numerous data fields that seemingly 
do nothing.  Someday they will do something amazing! (maybe)  

The installation comes with the file, sampleDat.rda.  We recommend looking at this file to see how
 your data needs to be formatted.  An explanation of the data format can be seen by entering
?ptmDat into the R terminal (after installing and loading the package).  

The most important function in the package is compBayes() function.  Calling this on your data
will fit a composititional Bayesian model to a dataset in a form similar to that in sampleDat.rda.
The results will return a list object where the first component contains a table with summary 
information.  The third component contains the Stan model object.  The second slot contains 
something thatwe disabled (sorry). After instaling and loading the package, information about this 
function and the parameters it takes can be seen by entering ``?compBayes'' into the R terminal 
(without the quotation marks).  

There are also two plotting functions built into the package caterpillar() and precisionPlot().
Again, more information can be obtained with ?caterpillar and ?precisionPlot.  

As an example, these functions can be used as follows.

RES <- compBayes(sampleDat) 

View(RES[[1]])

caterpillar(RES, byCond = TRUE)

precisionPlot(RES)
