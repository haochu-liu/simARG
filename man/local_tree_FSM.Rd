% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_tree_FSM.R
\name{local_tree_FSM}
\alias{local_tree_FSM}
\title{Select edges for local tree graph}
\usage{
local_tree_FSM(ARG, location)
}
\arguments{
\item{ARG}{A \code{FSM_ARG} object.}

\item{location}{The site location for local tree.}
}
\value{
\code{localtree} object; Local tree at a chosen site.
}
\description{
Construct a local tree from a \code{FSM_ARG} object with given site.
}
