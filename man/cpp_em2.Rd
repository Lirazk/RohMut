% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_em2}
\alias{cpp_em2}
\title{Estimate the mutation rate using EM}
\usage{
cpp_em2(
  mu1,
  mu2,
  mutation,
  len,
  indexes,
  maxiter,
  eps,
  prob,
  estimate_intercept = FALSE,
  intercept = 0
)
}
\arguments{
\item{mu1}{Starting value for the first mutation parameter}

\item{mu2}{Starting value for the second mutation parameter}

\item{mutation}{Vector of the number of mutations in a segment}

\item{len}{Vector of the (physical) length of a segment}

\item{indexes}{List of vector, where each vector is the indexes of the segments for an individual}

\item{maxiter}{Maximum number of EM iterations}

\item{eps}{The algorithm stops when the parameters change is less than eps, |prev_mu - mu| < eps}

\item{prob}{Matrix of probabilities, where the first column is the probability of being first cousin, and the second of second cousin.}

\item{estimate_intercept}{Whether to estimate the intercept as well, intercept is used when there's genotyping error.}

\item{intercept}{Starting value for the intercept.}

\item{stochastic}{If true, an incremental version is used instead.}
}
\value{
An estimate of the male mutation rate, female mutation rate and the intercept.
}
\description{
Estimate the mutation rate using EM
}
