% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ISM_ARG.R
\name{ISM_ARG}
\alias{ISM_ARG}
\title{A full ancestral recombination graph (ARG)}
\usage{
ISM_ARG(n, rho, bacteria = FALSE, delta = NULL)
}
\arguments{
\item{n}{An integer for the number of leaf lineages.}

\item{rho}{The recombination parameter.}

\item{bacteria}{logical; If TRUE, genes recombine by conversion.}

\item{delta}{numeric; If bacteria = TRUE, delta is the mean of recombinant segment length.}
}
\value{
A list containing dataframes of edges and nodes, and other information about ARG.
}
\description{
Simulate coalescent an recombination events and construct an ARG.
Assume infinite site model.
}
\examples{
ARG1 <- ISM_ARG(5L, 1)
ARG2 <- ISM_ARG(10L, 0.5, bacteria=TRUE, delta=0.1)
}
