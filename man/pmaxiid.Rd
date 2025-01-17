% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/piid.R
\name{pmaxiid}
\alias{pmaxiid}
\title{Cumulative distribution funtcion of the distribution of k-length block maxima under the assumption of temporal independence}
\usage{
pmaxiid(q, FX, k, asRP = 0, ...)
}
\arguments{
\item{q}{a scalar denoting the return levels}

\item{FX}{the cumulative distribution function of distribution of the parent process}

\item{k}{a scalar denoting the length of block (e.g., 365 for 1 year, in the case of a daily parent process)}

\item{asRP}{a boolean indicating whether or not the function will return probabilities or return periods}

\item{...}{additional arguments (parameter) that regards the distribution of the parent process}
}
\value{
vector of probabilities (or return periods, if asRP=1)
}
\description{
Cumulative distribution funtcion of the distribution of k-length block maxima under the assumption of temporal independence
}
