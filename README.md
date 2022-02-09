# RohMut

This R packages estimates the sex-specific mutation rate based from runs of homozygosity(ROH) segments from close common ancestor (Parents should be either first or second cousins).

We assume that the number of mutations in each generation is a Poisson process, with the rate depending on the sex of the ancestor.

Looking for example at the offspring of a first cousin, there are three possible pedigrees - Either the grand grandparents had two daughters/sons and those are the grandparents of the individual, or they had one daughter and one son, and they are the grandparents of the individual. We assume exchangeability, so the probability of two sons is 0.25, one son and one daughter is 0.5, and of two daughters is 0.25. Finally, each segment came from either the grand grandmother or grand grandfather with equal probability (actually due to the different recombination rate, the probability is slightly different from 0.5, but we assume it's 0.5).

That gives us a hierarchical Poisson mixture model, which we use maximum likelihood through a variation of the Expectation maximization algorithm (EM) to estimate the mutation rates.
Note that the likelihood is symmetric with respect to the mutation rate, so after the estimation, if the recombination rate is provided, we use it to select which of the parameters is the male's mutation rate and which is the female's based on the probability of recombination.

## Installation

To install using the devtools package:

```R
# If devtools is not available, install it first
# install.packages("devtools")
devtools::install_github("Lirazk/RohMut")
```
## Quick start

The main function is estimate_EM, which returns the estimated female and male mutation rate (in that order if a genomic map is provided), and if required, also the genotyping error rate.
It has the following required parameters:

* data - Data frame, with each row in it containing a data on one ROH segment for some individual. The expected columns are:
1. subject - the id of an individual.
2. mutation - the number of mutations in the segment.
3. length - the physical length of the segment.
* parms - Vector of starting values for the parameters. Can either be a vector of 2 elements if there is no genotyping error, or 3 in the case it is.
* probs - Vector of the probability of an individual to be the offspring of a first cousin.

It also has the following optional parameters:

* maxit - The maximum number of EM iterations.
* eps - The algorithm stops when the parameters change is less than eps, |prev_mu - mu| < eps.
* estimate_intercept - Whether we should estimate the genotyping error.
* male_map - A vector with the male's recombination rate (in morgans) of each segment of specific length.
* female_map - A vector with the female's recombination rate (in morgans) of each segment of specific length.
* male_start - A vector with the male's recombination rate (in morgans) around the starting point of each segment
* female_start - A vector with the female's recombination rate (in morgans) around the starting point of each segment.
* male_end - A vector with the male's recombination rate (in morgans) around the end point of each segment.
* female_end - A vector with the female's recombination rate (in morgans) around the end point of each segment.

### A Minimal working example:

```R
library(RohMut)

data("fs_cousins", package = "RohMut")
data("prob_model", package = "RohMut")

fs_cousins <- data.table(fs_cousins)

# 'model' is a logistic regression model, trained on the number of segments to predict whether it's a first or second cousin.
probs <- predict(model, fs_cousins[, .N, by = subject], type = "response")

fs_cousins[, mutation := rpois(nrow(fs_cousins), length * (1e-8 + 2.5e-9 * (Nf - male_meiosis) + 7.5e-9 * (Nm + male_meiosis)))]
# SQUAREM is used to boost the convergence rate of the algorithm.
SQUAREM::squarem(c(1e-8, 2e-8, 1e-7), 
          function(theta) RohMut::estimate_EM(fs_cousins, theta, probs, maxit = 1, estimate_intercept = T), 
      control = list(tol=1e-15))
```

### Example using the recombination rate

```R
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

data("prob_model", package = "RohMut")
probs <- predict(model, fs_cousins[, .N, by = subject], type = "response")
SQUAREM::squarem(runif(2, 0, 1e-8), function(theta) RohMut::estimate_EM(fs_cousins, theta, probs, maxit = 1, 
                male_map = male_dist/100, female_map = female_dist/100, male_start = male_start/100, 
                female_start = female_start/100, male_end = male_end/100, female_end = female_end/100), 
            control = list(tol=1e-15))
```

Estimates from 250 runs, with known genotyping error on 500 first cousins and 478 second cousins (The mode around 5e-9 for both the mutation rate disappear as the sample size increases):

![](result.png)