% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{plotMap}
\alias{plotMap}
\title{Plot map}
\usage{
plotMap(coords, statusMat, column = 1, bLegend = FALSE)
}
\arguments{
\item{coords}{Coordinates of units.}

\item{statusMat}{A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).}

\item{column}{Column (i.e. timepoints) to analyse. Only one column will be used.}

\item{bLegend}{Whether to add legend with statuses.}
}
\description{
Plots a map for a specified timepoint.
}
\author{
Mikhail Churakov
}
