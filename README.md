# RohMut

## Installation

```
devtools::install_github("Lirazk/RohMut")
```

## Use

There are two main ways to estimate the mutation rate. The first is to provide the male and female genetic distance of each segment + the genetic distance in area around the start and end of each segment, see ?estimate_EM for an example. One problem with it is that it's slightly biased (more so when the recombination model is wrong, which it is for real data).
The second is to not take recombination into account, which works but cant differentiate between the male and female mutation rate (which isn't a problem if we are willing to assume that the male's mutation rate is higher).

One way to solve the bias problem is to first estimate with the recombination rate, and use the resulting estimate as a starting value for the algorithm without recombination.

### Example

```
library(RohMut)

data("fs_cousins", package = "RohMut")
data("prob_model", package = "RohMut")

fs_cousins <- data.table(fs_cousins)

probs <- predict(model, fs_cousins[, .N, by = subject], type = "response")
probs <- cbind(probs, 1-probs)

fs_cousins[, mutation := rpois(nrow(fs_cousins), length * (1e-8 + 2.5e-9 * (Nf - male_meiosis) + 7.5e-9 * (Nm + male_meiosis)))]
SQUAREM::squarem(c(1e-8, 2e-8, 1e-7), 
          function(theta) RohMut::estimate_EM(fs_cousins, theta, probs, maxit = 1, estimate_intercept = T), 
      control = list(tol=1e-15))
```

Estimates from 250 runs, with known genotyping error:

![](result.png)