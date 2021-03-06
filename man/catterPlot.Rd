% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotingFunctions.R
\name{catterPlot}
\alias{catterPlot}
\title{Caterpillar plot}
\usage{
catterPlot(RES, byCond = FALSE, plotAll = FALSE, avgCond = FALSE)
}
\arguments{
\item{RES}{The results list generated from the
function \code{\link{compBayes}}}

\item{byCond}{A boolean parameter which determines if separate plots should
be made for each condition.  The default is FALSE which results in the
creation of a single plot.}

\item{plotAll}{A boolean variable indicating whether or not all parameters
should be plotted.  The default is false, resulting in a plot for only
parameters with 95\% credible intervals that do not contain zero.}

\item{avgCond}{A boolean parameter that determines whether or not plots
will be based on heierarchical population level mean parameters or average
within sample parameters.  The default is to use the population level
parameter, however with very small sample sizes this plot may not be useful.}
}
\description{
This function takes a list which contains a Stanfit object and summary
tables.  The output is a caterpillar plot of parameters with credible
intervals that do not contain 0.  Red line segments denote a 80\% credible
interval, while the black tails show 95\%.
}
