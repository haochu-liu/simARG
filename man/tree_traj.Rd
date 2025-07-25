% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_traj.R
\name{tree_traj}
\alias{tree_traj}
\title{Local tree trajectory}
\usage{
tree_traj(child, tree)
}
\arguments{
\item{child}{The child node index.}

\item{tree}{An \code{localtree} object.}
}
\value{
A vector of edge indices for trajectory.
}
\description{
Find the trajectory from a given node to the root node.
}
