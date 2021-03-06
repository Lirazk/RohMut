% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{estimate_EM}
\alias{estimate_EM}
\title{Estimate mutation parameters}
\usage{
estimate_EM(
  data,
  parms,
  probs,
  maxit = 100,
  eps = 1e-15,
  estimate_intercept = F,
  male_map = NULL,
  female_map = NULL,
  male_start = NULL,
  female_start = NULL,
  male_end = NULL,
  female_end = NULL
)
}
\arguments{
\item{data}{a data.frame with (at least) the following cols: 
subject - an id of the subject from which the segment was taken from.
mutation - the number of mutation in the segment.
length - the physical length of the segment.}

\item{parms}{A vector of starting values for the parameters.
If male_map and female_map is given, the first parameter would be the female's mutation rate and the second the male's. Otherwise we can't promise which is which.
If estimate_intercept is true then a third parameter should be provided as well, and is the starting value for the genotyping error rate.
You can also provide a third parameter when the genotyping error is known and estimation isn't needed (that is, estimate_intercept is false).}

\item{probs}{Vector of probabilities, where each is the probability of an individual to be the offspring of first cousin.}

\item{maxit}{The maximum number of EM iterations.}

\item{eps}{The tolerance. The algorithm stops when the parameters change is less than eps, |prev_mu - mu| < eps.}

\item{estimate_intercept}{Whether we should estimate the genotyping error.}

\item{male_map}{A vector with the male's recombination rate (in morgans) of each segment of specific length. Can be NULL.}

\item{female_map}{A vector with the female's recombination rate (in morgans) of each segment of specific length.}

\item{male_start}{A vector with the male's recombination rate (in morgans) around the starting point of each segment.}

\item{female_start}{A vector with the female's recombination rate (in morgans) around the starting point of each segment.}

\item{male_end}{A vector with the male's recombination rate (in morgans) around the end point of each segment.}

\item{female_end}{A vector with the female's recombination rate (in morgans) around the end point of each segment.}
}
\value{
Returns the female and male mutation rate, as well as the estimated genotyping error if estimate_intercept is true.
The first parameter would be the female's mutation rate and the second the male's if the genetic maps are provided, otherwise the order depends on the starting values.
}
\description{
estimate_EM is used to estimate the male's and female's mutation rate based on ROH when the common ancestor of the parents is very close (the assumption here is offspring of either first or second cousin).
The model is based on a mixture of Poisson, where the number of mutation is Poisson with the rate influenced by the segment length, number of male and female ancestors and the sex-specific mutation rate.
}
\examples{
library(data.table)
data("fs_cousins", package = "RohMut")
# Note: a somewhat larger file so it might take a few seconds to load
data("genetic_map", package = "RohMut")
fs_cousins <- data.table(fs_cousins)
genetic_map <- data.table(genetic_map)
# Calculate the male's and female's genetic distance for each segment by interpolation
male_dist <- vector(length = nrow(fs_cousins))
female_dist <- vector(length = nrow(fs_cousins))
for (i in 1:22) {
  index <- which(fs_cousins$chr == i)
  male_dist[index] <- approx(genetic_map[chr == i, pos],
                             genetic_map[chr == i, male_cM], xout = fs_cousins[index, end])$y -
    approx(genetic_map[chr == i, pos],
           genetic_map[chr == i, male_cM], xout = fs_cousins[index, start])$y
  
  female_dist[index] <- approx(genetic_map[chr == i, pos],
                               genetic_map[chr == i, female_cM], xout = fs_cousins[index, end])$y -
    approx(genetic_map[chr == i, pos],
           genetic_map[chr == i, female_cM], xout = fs_cousins[index, start])$y
}
# Same idea but for area around the start and end of each segment
male_start <- vector(length = nrow(fs_cousins))
male_end <- vector(length = nrow(fs_cousins))
female_start <- vector(length = nrow(fs_cousins))
female_end <- vector(length = nrow(fs_cousins))
range <- 25000
for (i in 1:22) {
  index <- which(fs_cousins$chr == i)
  male_end[index] <- approx(genetic_map[chr == i, pos],
                            genetic_map[chr == i, male_cM], xout = fs_cousins[index, end + range], rule = 2)$y -
    approx(genetic_map[chr == i, pos],
           genetic_map[chr == i, male_cM], xout = fs_cousins[index, end], rule = 2)$y
  
  male_start[index] <- approx(genetic_map[chr == i, pos],
                              genetic_map[chr == i, male_cM], xout = fs_cousins[index, start], rule = 2)$y -
    approx(genetic_map[chr == i, pos],
           genetic_map[chr == i, male_cM], xout = fs_cousins[index, start - range], rule = 2)$y
  
  female_end[index] <- approx(genetic_map[chr == i, pos],
                              genetic_map[chr == i, female_cM], xout = fs_cousins[index, end + range], rule = 2)$y -
    approx(genetic_map[chr == i, pos],
           genetic_map[chr == i, female_cM], xout = fs_cousins[index, end], rule = 2)$y
  female_start[index] <- approx(genetic_map[chr == i, pos],
                                genetic_map[chr == i, female_cM], xout = fs_cousins[index, start], rule = 2)$y -
    approx(genetic_map[chr == i, pos],
           genetic_map[chr == i, female_cM], xout = fs_cousins[index, start - range], rule = 2)$y
}

# For the estimation of probability, we use a simple logistic regression which classify each 
# individual based on the number of segments he has, and trained on other simulated data. 
# A more sophisticated models can (and probably should) be used, since otherwise there's a small bias
# in the estimate.
data("prob_model", package = "RohMut")
probs <- predict(model, fs_cousins[, .N, by = subject], type = "response")
RohMut::estimate_EM(fs_cousins, c(2e-8, 1e-8, 1e-8), probs, maxit = 100, estimate_intercept = T, 
male_map = male_dist/100,  female_map = female_dist/100, male_start = male_start/100, 
female_start = female_start/100, male_end = male_end/100, female_end = female_end/100)
# It's also possible to use acceleration schemes, for example:
SQUAREM::squarem(c(2e-8, 1e-8), function(theta) RohMut::estimate_EM(fs_cousins, theta, probs, maxit = 1, 
male_map = male_dist/100, female_map = female_dist/100, male_start = male_start/100, 
female_start = female_start/100, male_end = male_end/100, female_end = female_end/100), 
control = list(tol=1e-15))
}
