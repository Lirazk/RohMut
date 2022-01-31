#' Estimate mutation parameters
#'
#' estimate_EM is used to estimate the male's and female's mutation rate based on ROH when the common ancestor of the parents is very close (the assumption here is offspring of either first or second cousin).
#' The model is based on a mixture of poissons.
#'
#' @param data a data.frame with (at least) the following cols: 
#' subject - an id of the subject from which the segment was taken from.
#' mutation - the number of mutation in the segment.
#' length - the physical length of the segment.
#' @param parms A vector of starting values for the parameters.
#' If male_map and female_map is given, the first parameter should be the male's mutation rate and the second the female's. Otherwise we can't promise which is which.
#' If estimate_intercept is true or there's a known genotyping error, a third parameter should be provided as well, and is the estimated (or known) genotyping error rate.
#' @param probs Matrix of probabilities, where the first column is the probability of being first cousin, and the second of second cousin.
#' @param maxit The maximum number of EM iterations
#' @param eps The algorithm stops when the parameters change is less than eps, |prev_mu - mu| < eps.
#' @param estimate_intercept Whether we should estimate the genotyping error.
#' @param male_map A vector with the male's recombination rate (in morgans) of each segment of specific length. Can be NULL.
#' @param female_map  A vector with the female's recombination rate (in morgans) of each segment of specific length.
#' @param male_start   A vector with the male's recombination rate (in morgans) around the starting point of each segment.
#' @param female_start  A vector with the female's recombination rate (in morgans) around the starting point of each segment.
#' @param male_end A vector with the male's recombination rate (in morgans) around the end point of each segment.
#' @param female_end A vector with the female's recombination rate (in morgans) around the end point of each segment.
#'
#' @return Returns the female(first parameter) and male(second parameter) mutation rate, as well as the (estimated or given) genotyping error.
#' The order of the first two is only if the female and male maps are provided, otherwise they depend on the starting values.
#' @examples 
#' library(data.table)
#' data("fs_cousins", package = "RohMut")
#' # Note: a somewhat larger file so it might take a few seconds to load
#' data("genetic_map", package = "RohMut")
#' fs_cousins <- data.table(fs_cousins)
#' genetic_map <- data.table(genetic_map)
#' # Calculate the male's and female's genetic distance for each segment by interpolation
#' male_dist <- vector(length = nrow(fs_cousins))
#' female_dist <- vector(length = nrow(fs_cousins))
#' for (i in 1:22) {
#'   index <- which(fs_cousins$chr == i)
#'   male_dist[index] <- approx(genetic_map[chr == i, pos],
#'                              genetic_map[chr == i, male_cM], xout = fs_cousins[index, end])$y -
#'     approx(genetic_map[chr == i, pos],
#'            genetic_map[chr == i, male_cM], xout = fs_cousins[index, start])$y
#'   
#'   female_dist[index] <- approx(genetic_map[chr == i, pos],
#'                                genetic_map[chr == i, female_cM], xout = fs_cousins[index, end])$y -
#'     approx(genetic_map[chr == i, pos],
#'            genetic_map[chr == i, female_cM], xout = fs_cousins[index, start])$y
#' }
#' # Same idea but for area around the start and end of each segment
#' male_start <- vector(length = nrow(fs_cousins))
#' male_end <- vector(length = nrow(fs_cousins))
#' female_start <- vector(length = nrow(fs_cousins))
#' female_end <- vector(length = nrow(fs_cousins))
#' range <- 25000
#' for (i in 1:22) {
#'   index <- which(fs_cousins$chr == i)
#'   male_end[index] <- approx(genetic_map[chr == i, pos],
#'                             genetic_map[chr == i, male_cM], xout = fs_cousins[index, end + range], rule = 2)$y -
#'     approx(genetic_map[chr == i, pos],
#'            genetic_map[chr == i, male_cM], xout = fs_cousins[index, end], rule = 2)$y
#'   
#'   male_start[index] <- approx(genetic_map[chr == i, pos],
#'                               genetic_map[chr == i, male_cM], xout = fs_cousins[index, start], rule = 2)$y -
#'     approx(genetic_map[chr == i, pos],
#'            genetic_map[chr == i, male_cM], xout = fs_cousins[index, start - range], rule = 2)$y
#'   
#'   female_end[index] <- approx(genetic_map[chr == i, pos],
#'                               genetic_map[chr == i, female_cM], xout = fs_cousins[index, end + range], rule = 2)$y -
#'     approx(genetic_map[chr == i, pos],
#'            genetic_map[chr == i, female_cM], xout = fs_cousins[index, end], rule = 2)$y
#'   female_start[index] <- approx(genetic_map[chr == i, pos],
#'                                 genetic_map[chr == i, female_cM], xout = fs_cousins[index, start], rule = 2)$y -
#'     approx(genetic_map[chr == i, pos],
#'            genetic_map[chr == i, female_cM], xout = fs_cousins[index, start - range], rule = 2)$y
#' }
#' 
#' # For the estimation of probability, we use a simple logistic regression which classify each 
#' # individual based on the number of segments he has, and trained on other simulated data. 
#' # A more sophisticated models can (and probably should) be used, since otherwise there's a small bias
#' # in the estimate.
#' data("prob_model", package = "RohMut")
#' probs <- predict(model, data.table::data.table(fs_cousins)[, .N, by = subject], type = "response")
#' RohMut::estimate_EM(fs_cousins, c(2e-8, 1e-8, 0.0), probs, maxit = 100, estimate_intercept = T, 
#' male_map = male_dist/100,  female_map = female_dist/100, male_start = male_start/100, 
#' female_start = female_start/100, male_end = male_end/100, female_end = female_end/100)
#' # It's also possible to use acceleration schemes, for example:
#' SQUAREM::squarem(c(2e-8, 1e-8), function(theta) RohMut::estimate_EM(fs_cousins, theta, probs, maxit = 1, 
#' male_map = male_dist/100, female_map = female_dist/100, male_start = male_start/100, 
#' female_start = female_start/100, male_end = male_end/100, female_end = female_end/100), 
#' control = list(tol=1e-15))
#' @export
estimate_EM <- function(data, parms, probs, maxit=100, eps=1e-15, estimate_intercept = F,
                        male_map=NULL, female_map=NULL, male_start=NULL, female_start=NULL,
                        male_end=NULL, female_end=NULL) {
  
  needed_cols <- c("length", "mutation", "subject")
  # Can probably do it in one line
  missing_cols <- c(!(needed_cols[1] %in% colnames(data)), 
                    !(needed_cols[2] %in% colnames(data)), 
                    !(needed_cols[3] %in% colnames(data)))
  if (any(missing_cols)) {
    stop(sprintf("The data frame 'data' is missing the columns: %s.", paste(needed_cols[missing_cols])))
  }
  if (length(parms) != 2 + estimate_intercept) {
    stop(sprintf("The number of starting parameter values given was %d, but %d are required when estimate_intercept is %s.", length(parms), 2 + estimate_intercept, ifelse(estimate_intercept, "true", "false")))
  }
  # 
  if (!estimate_intercept) {
    parms <- c(parms, 0.0)
  }
  # the cpp functions expect probability as a matrix, since it's easier to generalize it that way, 
  # so first make sure that it's a vector, and then turn it into matrix
  if (is.vector(probs)) {
    probs <- cbind(probs, 1-probs)
  }
  
  # Annoying type changes in R
  indexes <- split(as.integer(1:nrow(data)-1), data$subject)
  # if(is.null(male_map)) {
  #   cpp_em(parms[1], parms[2], 
  #              as.integer(data$mutation), as.integer(data$length), 
  #              indexes, as.integer(maxit), eps, probs,
  #              estimate_intercept = estimate_intercept, intercept = parms[3])
  # }
  # else {
  #   cpp_em_map(parms[1], parms[2], 
  #              as.integer(data$mutation), as.integer(data$length), 
  #              indexes, as.integer(maxit), eps, probs, male_map, female_map,
  #              male_start, female_start, male_end, female_end,
  #              estimate_intercept = estimate_intercept, intercept = parms[3])
  # }
  estimate <- cpp_em(parms[1], parms[2],
         as.integer(data$mutation), as.integer(data$length),
         indexes, as.integer(maxit), eps, probs,
         estimate_intercept = estimate_intercept, intercept = parms[3])
  if (!is.null(male_map) && is.null(female_map)) {
    if(!is.null(male_start && !is.null(male_end) && !is.null(female_start) && !is.null(female_end))) {
      ll1 <- cpp_likelihood_hier_boundary(estimate[1], estimate[2], 
                                 as.integer(data$mutation), as.integer(data$length),
                                 indexes, male_map, female_map, male_start, female_start, male_end, female_end,
                                 probs, estimate[3])
      ll2 <- cpp_likelihood_hier_boundary(estimate[2], estimate[1], 
                                 as.integer(data$mutation), as.integer(data$length),
                                 indexes, male_map, female_map, male_start, female_start, male_end, female_end,
                                 probs, estimate[3])
    }
    else {
      ll1 <- cpp_likelihood_hier(estimate[1], estimate[2], 
                                 as.integer(data$mutation), as.integer(data$length),
                                 indexes, male_map, female_map,
                                 probs, estimate[3])
      ll2 <- cpp_likelihood_hier(estimate[2], estimate[1], 
                                 as.integer(data$mutation), as.integer(data$length),
                                 indexes, male_map, female_map,
                                 probs, estimate[3])
    }
    
    if (ll2 < ll1) {
      if (estimate_intercept) {
        return(c(estimate[2], estimate[1], estimate[3]))
      }
      else {
        c(estimate[2], estimate[1])
      }
    }
  }
  if (estimate_intercept) {
    return(estimate)
  }
  else {
    return(estimate[1:2])
  }
}