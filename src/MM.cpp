#include <math.h>

#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// Too many versions of it.
double logsumexp(NumericMatrix &mat) {
  double maximum = mat(0, 0);
  double sum = 0.0;
  // for (size_t i = 0; i < mat.nrow(); i++) {
  //   for (size_t j = 0; j < mat.ncol(); j++) {
  for (size_t j = 0; j < mat.ncol(); j++) {
    for (size_t i = 0; i < mat.nrow(); i++) {
      if (mat(i, j) > maximum) {
        maximum = mat(i, j);
      }
    }
  }
  if (maximum < -1000000000) {
    return (maximum);
  }
  for (size_t j = 0; j < mat.ncol(); j++) {
    for (size_t i = 0; i < mat.nrow(); i++) {
      sum += exp(mat(i, j) - maximum);
    }
  }

  return (log(sum) + maximum);
}

double logsumexp(double a, double b) {
  double maximum = fmax(a, b);
  if (maximum < -1000000000) {
    return (maximum);
  }
  return (log(exp(a - maximum) + exp(b - maximum)) + maximum);
}

template <std::size_t N> double logsumexp(const double (&arr)[N]) {
  double maximum = arr[0];
  double sum = 0.0;
  for (size_t i = 0; i < N; i++) {
    if (arr[i] > maximum) {
      maximum = arr[i];
    }
  }
  // if (isinf(maximum)) {
  if (maximum < -1000000000) {
    return (maximum);
  }

  for (size_t i = 0; i < N; i++) {
    sum += exp(arr[i] - maximum);
  }
  return (log(sum) + maximum);
}

double logsumexp(const std::vector<double> &arr) {
  size_t N = arr.capacity();
  double maximum = arr[0];
  double sum = 0.0;
  for (size_t i = 0; i < N; i++) {
    if (arr[i] > maximum) {
      maximum = arr[i];
    }
  }
  // if (isinf(maximum)) {
  if (maximum < -1000000000) {
    return (maximum);
  }

  for (size_t i = 0; i < N; i++) {
    sum += exp(arr[i] - maximum);
  }
  return (log(sum) + maximum);
}

double logsumexp(const double a, const std::vector<double> &arr) {
  size_t N = arr.capacity();
  double maximum = a;
  double sum = 0.0;
  for (size_t i = 0; i < N; i++) {
    if (arr[i] > maximum) {
      maximum = arr[i];
    }
  }
  // if (isinf(maximum)) {
  if (maximum < -1000000000) {
    return (maximum);
  }

  for (size_t i = 0; i < N; i++) {
    sum += exp(arr[i] - maximum);
  }
  sum += exp(a - maximum);
  return (log(sum) + maximum);
}

// The likelihood of a specific segment, given that there are n1 female ancestors and n2 male ancestors
// d is the number of mutations, len is the lengths, mu_f and mu_m the female and male mutation rate correspondingly.
// The male/female lengths are in morgans.
NumericVector log_segment(int n1, int n2, const IntegerVector &d,
                           const IntegerVector &len, double mu_f, double mu_m, double intercept,
                           const NumericVector &male_map,
                           const NumericVector &female_map) {
  NumericVector result(d.size());
  // Optimization for small number of mutations... which should be fine for real mutation rate.
  static const double lgamma_lookup[] = {
      lgamma(1), lgamma(2), lgamma(3), lgamma(4),  lgamma(5),  lgamma(6),
      lgamma(7), lgamma(8), lgamma(9), lgamma(10), lgamma(11), lgamma(12)};
  for (size_t i = 0; i < d.size(); i++) {
    const double temp = lgamma_lookup[(int)d[i]];
    const double d1 = (n1 + 1) * female_map[i] + (n2 - 1) * male_map[i];
    const double d2 = (n1 - 1) * female_map[i] + (n2 + 1) * male_map[i];

    const double par1 = len[i] * (intercept + (n1 + 1) * mu_f + (n2 - 1) * mu_m);
    const double par2 = len[i] * (intercept + (n1 - 1) * mu_f + (n2 + 1) * mu_m);
    const double pois1 = log(0.5) + d[i] * log(par1) - par1 - temp - d1;
    const double pois2 = log(0.5) + d[i] * log(par2) - par2 - temp - d2;

    result[i] = logsumexp({pois1, pois2});
  }
  return (result);
}

// The likelihood of a specific segment, including the terms related to the start and end of each segment, given that there are n1 female ancestors and n2 male ancestors
// d is the number of mutations, len is the lengths, mu_f and mu_m the female and male mutation rate correspondingly.
// The male/female lengths are in morgans.
NumericVector log_segment(int n1, int n2, const IntegerVector &d, const IntegerVector &len,
                           double mu_f, double mu_m, double intercept, const NumericVector &male_map,
                           const NumericVector &female_map, const NumericVector &male_start,
                           const NumericVector &female_start, const NumericVector &male_end,
                           const NumericVector &female_end) {
  NumericVector result(d.size());
  // Optimization for small number of mutations... which should be fine for real mutation rate.
  static const double lgamma_lookup[] = {
      lgamma(1), lgamma(2), lgamma(3), lgamma(4),  lgamma(5),  lgamma(6),
      lgamma(7), lgamma(8), lgamma(9), lgamma(10), lgamma(11), lgamma(12)};
  for (size_t i = 0; i < d.size(); i++) {
    const double temp = lgamma_lookup[(int)d[i]];
    const double d1 = (n1 + 1) * female_map[i] + (n2 - 1) * male_map[i];
    const double d2 = (n1 - 1) * female_map[i] + (n2 + 1) * male_map[i];
    double w1 = log((n1 + 1) * female_start[i] + (n2 - 1) * male_start[i]) +
                log((n1 + 1) * female_end[i] + (n2 - 1) * male_end[i]);
    double w2 = log((n1 - 1) * female_start[i] + (n2 + 1) * male_start[i]) +
                log((n1 - 1) * female_end[i] + (n2 + 1) * male_end[i]);

    // Avoid taking log(0), start/end should be 0 when the points are around the boundary of the chromosome and so we shouldn't count them.
    if (male_end[i] == 0.0 && male_start[i] == 0.0) {
      w1 = 0.0;
      w2 = 0.0;
    } else if (male_end[i] == 0.0 || female_end[i] == 0.0) {
      w1 = log((n1 + 1) * female_start[i] + (n2 - 1) * male_start[i]);
      w2 = log((n1 - 1) * female_start[i] + (n2 + 1) * male_start[i]);
    } else if (male_start[i] == 0.0 || female_start[i] == 0.0) {
      w1 = log((n1 + 1) * female_end[i] + (n2 - 1) * male_end[i]);
      w2 = log((n1 - 1) * female_end[i] + (n2 + 1) * male_end[i]);
    }

    const double par1 = len[i] * (intercept + (n1 + 1) * mu_f + (n2 - 1) * mu_m);
    const double par2 = len[i] * (intercept + (n1 - 1) * mu_f + (n2 + 1) * mu_m);
    const double pois1 = log(0.5) + d[i] * log(par1) - par1 - temp - d1 + len[i]*w1;
    const double pois2 = log(0.5) + d[i] * log(par2) - par2 - temp - d2 + len[i]*w2;
    result[i] = logsumexp({pois1, pois2});
  }
  return (result);
}

// [[Rcpp::export]]
double cpp_likelihood_hier(double mu_f, double mu_m,
                           const IntegerVector &d, const IntegerVector &len,
                           const List &indexes,
                           const NumericVector &male_map,
                           const NumericVector &female_map,
                           const NumericMatrix &prob,
                           double intercept = 0.0)
{
  // mu_f = exp(mu_f);
  // mu_m = exp(mu_m);
  double ll = 0;
  
  // NumericVector N_ = {3, 5, 7};
  NumericVector N_ = {3, 5};

  // Not sure whether to preallocate everything.
  NumericMatrix segment (d.size(), sum(N_)); 
  
  size_t i = 0;
  // Calculate the likelihood for each segment.
  for (int N: N_) {
    for (size_t j = 0; j < N; j++) {
      segment(_, i) = log_segment(j+2, (N+3) - (j+2), d, len, mu_f, mu_m, intercept, male_map, female_map);
      i++;
    }
  }

  for (size_t i = 0; i < indexes.size(); ++i) {
    IntegerVector index = indexes[i];
    std::vector<double> sum_s (sum(N_), 0.0);
    for (size_t j = 0; j < index.size(); ++j) {
      for (size_t k = 0; k < sum(N_); k++) {
        sum_s[k] += segment(index[j], k);
      }
    }
    std::vector<double> segment_likelihood (sum(N_), 0.0);
    size_t m_all = 0;
    size_t n_counter = 0;
    for (int N: N_) {
      for (size_t j = 0; j < N; j++) {
        segment_likelihood[m_all] = prob(i, n_counter) + R::dbinom(j, N - 1, 0.5, 1) + sum_s[m_all];
        m_all++;
      }
      n_counter++;
    }
 
    ll += logsumexp(segment_likelihood);
  }
  return (-ll);
}

// [[Rcpp::export]]
double cpp_likelihood_hier_boundary(double mu_f, double mu_m,
                           const IntegerVector &d, const IntegerVector &len,
                           const List &indexes,
                           const NumericVector &male_map,
                           const NumericVector &female_map,
                           const NumericVector &male_start,
                           const NumericVector &female_start,
                           const NumericVector &male_end,
                           const NumericVector &female_end,
                           const NumericMatrix &prob,
                           double intercept = 0.0)
{
  // mu_f = exp(mu_f);
  // mu_m = exp(mu_m);
  double ll = 0;
  
  // NumericVector N_ = {3, 5, 7};
  NumericVector N_ = {3, 5};

  // Not sure whether to preallocate everything.
  NumericMatrix segment (d.size(), sum(N_)); 
  
  size_t i = 0;
  // Calculate the likelihood for each segment.
  for (int N: N_) {
    for (size_t j = 0; j < N; j++) {
      segment(_, i) = log_segment(j+2, (N+3) - (j+2), d, len, mu_f, mu_m, intercept, male_map, female_map, male_start, female_start, male_end, female_end);
      i++;
    }
  }

  for (size_t i = 0; i < indexes.size(); ++i) {
    IntegerVector index = indexes[i];
    std::vector<double> sum_s (sum(N_), 0.0);
    for (size_t j = 0; j < index.size(); ++j) {
      for (size_t k = 0; k < sum(N_); k++) {
        sum_s[k] += segment(index[j], k);
      }
    }
    std::vector<double> segment_likelihood (sum(N_), 0.0);
    size_t m_all = 0;
    size_t n_counter = 0;
    for (int N: N_) {
      for (size_t j = 0; j < N; j++) {
        segment_likelihood[m_all] = prob(i, n_counter) + R::dbinom(j, N - 1, 0.5, 1) + sum_s[m_all];
        m_all++;
      }
      n_counter++;
    }
 
    ll += logsumexp(segment_likelihood);
  }
  return (-ll);
}

// NumericVector cpp_em_map(
//     double mu1, double mu2, const NumericVector &mutation,
//     const NumericVector &len, List &indexes, int maxiter, double eps,
//     const NumericMatrix &prob, NumericVector male_map, NumericVector female_map,
//     NumericVector male_start = NA_REAL, NumericVector female_start = NA_REAL,
//     NumericVector male_end = NA_REAL, NumericVector female_end = NA_REAL,
//     bool stochastic = false, bool estimate_inner_pi = false,
//     bool estimate_outer_p = false,
//     bool estimate_intercept = false, double intercept = 0.0) 

//' Estimate the mutation rate using EM
//'
//' @param mu1 Starting value for the female's mutation parameter
//' @param mu2 Starting value for the male's mutation parameter
//' @param mutation Vector of the number of mutations in a segment
//' @param len Vector of the length of a segment
//' @param indexes List of vector, where each vector is the indexes of the segments for an individual
//' @param maxiter Maximum number of EM iterations
//' @param eps The algorithm stops when the parameters change is less than eps, |prev_mu - mu| < eps
//' @param prob Matrix of probabilities, where the first column is the probability of being first cousin, and the second of second cousin.
//' @param estimate_intercept Whether to estimate the intercept as well, intercept is used when there's genotyping error.
//' @param intercept Starting value for the intercept.
//' @param stochastic If true, an incremental version is used instead.
//' @return An estimate of the male mutation rate, female mutation rate and of the intercept.
// [[Rcpp::export]]
NumericVector cpp_em_map(
    double mu1, double mu2, const NumericVector &mutation,
    const NumericVector &len, List &indexes, int maxiter, double eps,
    const NumericMatrix &prob, IntegerVector male_map, IntegerVector female_map,
    NumericVector male_start = NA_REAL, NumericVector female_start = NA_REAL,
    NumericVector male_end = NA_REAL, NumericVector female_end = NA_REAL,
    bool stochastic = false,
    bool estimate_intercept = false, double intercept = 0.0) 
{
  // NumericVector N_ = {3, 5, 7};
NumericVector N_ = {3, 5};
  // Few tests
  if (prob.nrow() != indexes.size()) {
    REprintf("Number of rows in prob should be the same as the length of indexes.\n");
    // fprintf(
    //     stdout,
    //     "Number of rows in prob should be the same as the length of indexes.");
    return NA_REAL;
  }
  if (prob.ncol() != N_.size()) {
    REprintf("Number of cols in prob should be 2.\n");
      // fprintf(
      //   stdout,
      //   "Number of cols in prob should be 2.");
    return NA_REAL;
  }
  double prev_mu1 = mu1;
  double prev_mu2 = mu2;
  int n = mutation.size();

  int N_all = sum(N_);

  std::vector<double> pi1(N_all, 0.5);
  std::vector<double> pi2(N_all, 0.5);

  // std::vector<double> outer_p(N_all, 1/(double)N_all);

  // Precalcualte some of the possible things.
  std::vector<double> likelihood_precalc(n, 0.0);
  for (size_t k = 0; k < indexes.size(); k++) {
    IntegerVector index = indexes[k];
    for (size_t j = 0; j < index.size(); j++) {
      double temp = lgamma(mutation[index[j]] + 1);
      likelihood_precalc[index[j]] = mutation[index[j]] * log(len[index[j]]) - temp;
    }
  }

  // Another preallocation since it was a bit faster to do logsumexp in the end instead of in each iteration.
  // sm is the terms related to the mutation for each parameter
  // sl is the terms related to the lengths for each parameter
  // and si is the terms for the intercept.
  NumericMatrix sm1(n, N_all);
  NumericMatrix sm2(n, N_all);
  NumericMatrix sl1(n, N_all);
  NumericMatrix sl2(n, N_all);
  NumericMatrix si1(n, N_all);
  NumericMatrix si2(n, N_all);

  std::fill(sm1.begin(), sm1.end(), R_NegInf);
  std::fill(sm2.begin(), sm2.end(), R_NegInf);
  std::fill(sl1.begin(), sl1.end(), R_NegInf);
  std::fill(sl2.begin(), sl2.end(), R_NegInf);
  std::fill(si1.begin(), si1.end(), R_NegInf);
  std::fill(si2.begin(), si2.end(), R_NegInf);
  for (int i = 0; i < maxiter; i++) {
    // double sum_len[2] = {-INFINITY, -INFINITY};
    // double sum_mutation[2] = {-INFINITY, -INFINITY};

    // std::vector<double> sum_p1(N_all, -INFINITY);
    // std::vector<double> sum_p2(N_all, -INFINITY);

    // std::vector<double> sum_outer_p(N_all, -INFINITY);
    // double sum_intercept[2] = {-INFINITY, -INFINITY};

    // Calculate P(Z | Y_i)
    // for (size_t n_index = 0; n_index < N_.size(); n_index++) {
    // int N = N_[n_index];
    size_t p = indexes.size();

    for (size_t k = 0; k < p; k++) {
      IntegerVector index = indexes[k];
      if (stochastic) {
        index = sample(indexes.size(), 1, false, R_NilValue, false);
      }

      // Calculate wia
      std::vector<double> wi(N_all, 0.0);
      for (size_t j = 0; j < index.size(); j++) {
        double l = len[index[j]];
        size_t m_all = 0;
        for (int N : N_) {
          // Naming is hard.
          const int general_N = N + 3;
          for (size_t m = 0; m < N; m++) {
            const int n1 = m + 2;
            const int n2 = general_N - n1;

            const double par[2] = {
                (m + 1) * mu1 + (general_N - (m + 1)) * mu2 + intercept,
                (m + 3) * mu1 + (general_N - (m + 3)) * mu2 + intercept};
            const double d[2] = {
                (m + 1) * male_map[index[j]] +
                    (general_N - (m + 1)) * female_map[index[j]],
                (m + 3) * male_map[index[j]] +
                    (general_N - (m + 3)) * female_map[index[j]]};

            double w1 = log((n1 + 1) * female_start[i] + (n2 - 1) * male_start[i]) +
                log((n1 + 1) * female_end[i] + (n2 - 1) * male_end[i]);
            double w2 = log((n1 - 1) * female_start[i] + (n2 + 1) * male_start[i]) +
                        log((n1 - 1) * female_end[i] + (n2 + 1) * male_end[i]);

            if (male_end[i] == 0.0 && male_start[i] == 0.0) {
              w1 = 0.0;
              w2 = 0.0;
            } else if (male_end[i] == 0.0 || female_end[i] == 0.0) {
              w1 = log((n1 + 1) * female_start[i] + (n2 - 1) * male_start[i]);
              w2 = log((n1 - 1) * female_start[i] + (n2 + 1) * male_start[i]);
            } else if (male_start[i] == 0.0 || female_start[i] == 0.0) {
              w1 = log((n1 + 1) * female_end[i] + (n2 - 1) * male_end[i]);
              w2 = log((n1 - 1) * female_end[i] + (n2 + 1) * male_end[i]);
            }

            wi[m_all] +=
                logsumexp({mutation[index[j]] * log(par[0]) - l * par[0] +
                               likelihood_precalc[index[j]] - d[0] + w1,
                           mutation[index[j]] * log(par[1]) - l * par[1] +
                               likelihood_precalc[index[j]] - d[1] + w2});
            m_all++;
          }
        }
      }
      size_t m_all = 0;
      size_t n_counter = 0;
      for (int N : N_) {
        for (size_t m = 0; m < N; m++) {
          // Change 0.5 to prob for N=3, 5 etc
          // wi[m] += log(0.5) + R::dbinom(m, N-1, 0.5, 1);
        //   if (estimate_outer_p) {
        //       wi[m_all] += log(outer_p[m_all]);
        //   }
        //   else {
              wi[m_all] += log(prob(k, n_counter)) + R::dbinom(m, N - 1, 0.5, 1);
        //   }
          m_all++;
        }
        n_counter++;
      }
      double normalize = logsumexp(wi);
      for (size_t m = 0; m < N_all; m++) {
        wi[m] -= normalize;
      }

      for (size_t j = 0; j < index.size(); j++) {
        double l = len[index[j]];
        double log_mutation = log(mutation[index[j]]);
        double log_len = log(l);
        double a[N_all];
        double b[N_all];
        {
          size_t m_all = 0;
          for (int N : N_) {
            const int general_N = N + 3;
            for (size_t m = 0; m < N; m++) {
              const double par[2] = {
                  (m + 1) * mu1 + (general_N - (m + 1)) * mu2 + intercept,
                  (m + 3) * mu1 + (general_N - (m + 3)) * mu2 + intercept};
              const double log_par[2] = {log(par[0]), log(par[1])};
              const double name[2] = {-log_par[0] + log((m + 1) * mu1),
                                      -log_par[1] + log((m + 3) * mu1)};
              const double d[2] = {
                  (m + 1) * male_map[index[j]] +
                      (general_N - (m + 1)) * female_map[index[j]],
                  (m + 3) * male_map[index[j]] +
                      (general_N - (m + 3)) * female_map[index[j]]};

              const int n1 = m + 2;
              const int n2 = general_N - n1;

              double w1 = log((n1 + 1) * female_start[i] + (n2 - 1) * male_start[i]) +
                  log((n1 + 1) * female_end[i] + (n2 - 1) * male_end[i]);
              double w2 = log((n1 - 1) * female_start[i] + (n2 + 1) * male_start[i]) +
                          log((n1 - 1) * female_end[i] + (n2 + 1) * male_end[i]);

              if (male_end[i] == 0.0 && male_start[i] == 0.0) {
                w1 = 0.0;
                w2 = 0.0;
              } else if (male_end[i] == 0.0 || female_end[i] == 0.0) {
                w1 = log((n1 + 1) * female_start[i] + (n2 - 1) * male_start[i]);
                w2 = log((n1 - 1) * female_start[i] + (n2 + 1) * male_start[i]);
              } else if (male_start[i] == 0.0 || female_start[i] == 0.0) {
                w1 = log((n1 + 1) * female_end[i] + (n2 - 1) * male_end[i]);
                w2 = log((n1 - 1) * female_end[i] + (n2 + 1) * male_end[i]);
              }

              a[m_all] = mutation[index[j]] * log_par[0] - l * par[0] +
                         likelihood_precalc[index[j]] - d[0] + w1 + log(pi1[m_all]);
              b[m_all] = mutation[index[j]] * log_par[1] - l * par[1] +
                         likelihood_precalc[index[j]] - d[1] + w2 + log(pi2[m_all]);
              double normalize = logsumexp({a[m_all], b[m_all]});
              a[m_all] = a[m_all] - normalize;
              b[m_all] = b[m_all] - normalize;
              double temp_wi = logsumexp({a[m_all] + log_mutation + name[0],
                                          b[m_all] + log_mutation + name[1]}) +
                               wi[m_all];
              sm1(index[j], m_all) = temp_wi;

              const double name2[2] = {
                  -log_par[0] + log((general_N - (m + 1)) * mu2),
                  -log_par[1] + log((general_N - (m + 3)) * mu2)};

              temp_wi = logsumexp({a[m_all] + log_mutation + name2[0],
                                   b[m_all] + log_mutation + name2[1]}) +
                        wi[m_all];
              // sum_mutation[1] = logsumexp(sum_mutation[1], temp_wi);
              // sm2[index[j]] = temp_wi;
              sm2(index[j], m_all) = temp_wi;

              temp_wi = logsumexp({a[m_all] + log_len + log((m + 1)),
                                   b[m_all] + log_len + log((m + 3))}) +
                        wi[m_all];
              // sum_len[0] = logsumexp(sum_len[0], temp_wi);
              // sl1[index[j]] = temp_wi;
              sl1(index[j], m_all) = temp_wi;

              temp_wi =
                  logsumexp(
                      {a[m_all] + log_len + log(general_N - (m + 1)),
                       b[m_all] + log_len + log(general_N - (m + 3))}) +
                  wi[m_all];
              // sum_len[1] = logsumexp(sum_len[1], temp_wi);
              // sl2[index[j]] = temp_wi;
              sl2(index[j], m_all) = temp_wi;

              temp_wi =
                  logsumexp(
                      {a[m_all] + log_mutation - log_par[0] + log(intercept),
                       b[m_all] + log_mutation - log_par[1] + log(intercept)}) +
                  wi[m_all];
              // sum_intercept[0] = logsumexp(sum_intercept[0], temp_wi);
              // si1[index[j]] = temp_wi;
              si1(index[j], m_all) = temp_wi;

              temp_wi = logsumexp({a[m_all] + log_len, b[m_all] + log_len}) +
                        wi[m_all];
              // sum_intercept[1] = logsumexp(sum_intercept[1], temp_wi);
              // si2[index[j]] = temp_wi;
              si2(index[j], m_all) = temp_wi;

              m_all++;
            }
          }
        }

        // for (size_t m = 0; m < N_all; m++) {
        //   sum_p1[m] = logsumexp({sum_p1[m], a[m]});
        //   sum_p2[m] = logsumexp({sum_p2[m], b[m]});
        // }
        // for (size_t m = 0; m < N_all; m++) {
        //     sum_outer_p[m] = logsumexp({sum_outer_p[m], wi[m]});
        // }
      }
    }
    // mu1 = exp(sum_mutation[0] - sum_len[0]);
    // mu2 = exp(sum_mutation[1] - sum_len[1]);
    mu1 = exp(logsumexp(sm1) - logsumexp(sl1));
    mu2 = exp(logsumexp(sm2) - logsumexp(sl2));

    // if (estimate_inner_pi) {
    //   for (size_t m = 0; m < N_all; m++) {
    //     pi1[m] = exp(sum_p1[m]) / n;
    //     pi2[m] = exp(sum_p2[m]) / n;
    //   }
    // }
    if (estimate_intercept) {
      // intercept = exp(sum_intercept[0] - sum_intercept[1]);
      intercept = exp(logsumexp(si1) - logsumexp(si2));
    }
    // if (estimate_outer_p) {
    //     for (size_t m = 0; m < N_all; m++) {
    //         outer_p[m] = exp(sum_outer_p[m]) / n;
    //     }
    // }

    if (fabs(mu1 - prev_mu1) + fabs(mu2 - prev_mu2) < eps) {
      break;
    }
    prev_mu1 = mu1;
    prev_mu2 = mu2;
  }
  NumericVector result = NumericVector(3);
  result[0] = mu1;
  result[1] = mu2;
  result[2] = intercept;

//   if (estimate_inner_pi) {
//     fprintf(stderr, "Pi:\n");
//     for (size_t m = 0; m < N_all; m++) {
//       fprintf(stderr, "%d: (%f, %f) \n", m, pi1[m], pi2[m]);
//     }
//     fprintf(stderr, "\n");
//   }

//   if (estimate_outer_p) {
//     fprintf(stderr, "Pi:\n");
//     for (size_t m = 0; m < N_all; m++) {
//       fprintf(stderr, "%d: %f \n", m, outer_p[m]);
//     }
//     fprintf(stderr, "\n");
//   }
  // if (estimate_intercept) {
  //   fprintf(stderr, "Intercept: %f\n", intercept * 1e9);
  // }
  return (result);
}


//' Estimate the mutation rate using EM
//'
//' @param mu1 Starting value for the female's mutation parameter
//' @param mu2 Starting value for the male's mutation parameter
//' @param mutation Vector of the number of mutations in a segment
//' @param len Vector of the length of a segment
//' @param indexes List of vector, where each vector is the indexes of the segments for an individual
//' @param maxiter Maximum number of EM iterations
//' @param eps The algorithm stops when the parameters change is less than eps, |prev_mu - mu| < eps
//' @param prob Matrix of probabilities, where the first column is the probability of being first cousin, and the second of second cousin.
//' @param estimate_intercept Whether to estimate the intercept as well, intercept is used when there's genotyping error.
//' @param intercept Starting value for the intercept.
//' @param stochastic If true, an incremental version is used instead.
//' @return An estimate of the male mutation rate, female mutation rate and of the intercept.
// [[Rcpp::export]]
NumericVector cpp_em(
    double mu1, double mu2, const IntegerVector &mutation,
    const IntegerVector &len, List &indexes, int maxiter, double eps,
    const NumericMatrix &prob, bool stochastic = false,
    bool estimate_intercept = false, double intercept = 0.0) 
{
  // NumericVector N_ = {3, 5, 7};
  NumericVector N_ = {3, 5};
  // Few tests
  if (prob.nrow() != indexes.size()) {
    // fprintf(
    //     stdout,
    //     "Number of rows in prob should be the same as the length of indexes.");
    //     return NA_REAL;
    REprintf("Number of rows in prob should be the same as the length of indexes.\n");
    return NA_REAL;
  }
  if (prob.ncol() != N_.size()) {
    REprintf("Number of cols in prob should be 2.\n");
      // fprintf(
      //   stdout,
      //   "Number of cols in prob should be 3.");
    return NA_REAL;
  }
  double prev_mu1 = mu1;
  double prev_mu2 = mu2;
  int n = mutation.size();

  int N_all = sum(N_);

  std::vector<double> pi1(N_all, 0.5);
  std::vector<double> pi2(N_all, 0.5);

  // std::vector<double> outer_p(N_all, 1/(double)N_all);
  // double intercept = 0.0;
  if (estimate_intercept && intercept == 0.0) {
    intercept = 1e-8;
  }

  std::vector<double> likelihood_precalc(n, 0.0);
  for (size_t k = 0; k < indexes.size(); k++) {
    IntegerVector index = indexes[k];
    for (size_t j = 0; j < index.size(); j++) {
      double temp = lgamma(mutation[index[j]] + 1);
      likelihood_precalc[index[j]] = mutation[index[j]] * log(len[index[j]]) - temp;
    }
  }

  NumericMatrix sm1(n, N_all);
  NumericMatrix sm2(n, N_all);
  NumericMatrix sl1(n, N_all);
  NumericMatrix sl2(n, N_all);
  NumericMatrix si1(n, N_all);
  NumericMatrix si2(n, N_all);

  std::fill(sm1.begin(), sm1.end(), R_NegInf);
  std::fill(sm2.begin(), sm2.end(), R_NegInf);
  std::fill(sl1.begin(), sl1.end(), R_NegInf);
  std::fill(sl2.begin(), sl2.end(), R_NegInf);
  std::fill(si1.begin(), si1.end(), R_NegInf);
  std::fill(si2.begin(), si2.end(), R_NegInf);
  for (int i = 0; i < maxiter; i++) {
    // double sum_len[2] = {-INFINITY, -INFINITY};
    // double sum_mutation[2] = {-INFINITY, -INFINITY};

    // std::vector<double> sum_p1(N_all, -INFINITY);
    // std::vector<double> sum_p2(N_all, -INFINITY);

    // std::vector<double> sum_outer_p(N_all, -INFINITY);
    // double sum_intercept[2] = {-INFINITY, -INFINITY};

    // Calculate P(Z | Y_i)
    // for (size_t n_index = 0; n_index < N_.size(); n_index++) {
    // int N = N_[n_index];

    // When stochastic
    // double logsumexp_buffer[4] = {0.0};
    // double prev_buffer[4] = {0.0};
    size_t p = indexes.size();
    if (i == 1 && stochastic) {
      p = 1;
    }

    for (size_t k = 0; k < p; k++) {
      IntegerVector index = indexes[k];
      if (stochastic) {
        index = sample(indexes.size(), 1, false, R_NilValue, false);
        // prev_buffer[0] = 
      }

      // Calculate wia
      std::vector<double> wi(N_all, 0.0);
      for (size_t j = 0; j < index.size(); j++) {
        double l = len[index[j]];
        size_t m_all = 0;
        for (int N : N_) {
          for (size_t m = 0; m < N; m++) {
            const int general_N = N + 3;
            const double par[2] = {
            (m + 1) * mu1 + (general_N - (m + 1)) * mu2 + intercept,
            (m + 3) * mu1 + (general_N - (m + 3)) * mu2 + intercept};

            wi[m_all] +=
                logsumexp({mutation[index[j]] * log(par[0]) - l * par[0] +
                               likelihood_precalc[index[j]],
                           mutation[index[j]] * log(par[1]) - l * par[1] +
                               likelihood_precalc[index[j]]});
            m_all++;
          }
        }
      }
      size_t m_all = 0;
      size_t n_counter = 0;
      for (int N : N_) {
        for (size_t m = 0; m < N; m++) {
          // Change 0.5 to prob for N=3, 5 etc
          // wi[m] += log(0.5) + R::dbinom(m, N-1, 0.5, 1);
        //   if (estimate_outer_p) {
              // wi[m_all] += log(outer_p[m_all]);
        //   }
        //   else {
              wi[m_all] += log(prob(k, n_counter)) + R::dbinom(m, N - 1, 0.5, 1);
        //   }
          m_all++;
        }
        n_counter++;
      }
      double normalize = logsumexp(wi);
      for (size_t m = 0; m < N_all; m++) {
        wi[m] -= normalize;
      }

      for (size_t j = 0; j < index.size(); j++) {
        double l = len[index[j]];

        double log_mutation = log(mutation[index[j]]);
        double log_len = log(l);
        double a[N_all];
        double b[N_all];

        {
          size_t m_all = 0;
          for (int N : N_) {
            const int general_N = N + 3;
            for (size_t m = 0; m < N; m++) {
              const double par[2] = {
                  (m + 1) * mu1 + (general_N - (m + 1)) * mu2 + intercept,
                  (m + 3) * mu1 + (general_N - (m + 3)) * mu2 + intercept};
              const double log_par[2] = {log(par[0]), log(par[1])};
              const double name[2] = {-log_par[0] + log((m + 1) * mu1),
                                      -log_par[1] + log((m + 3) * mu1)};

              a[m_all] = mutation[index[j]] * log_par[0] - l * par[0] +
                         likelihood_precalc[index[j]] + log(pi1[m_all]);
              b[m_all] = mutation[index[j]] * log_par[1] - l * par[1] +
                         likelihood_precalc[index[j]] + log(pi2[m_all]);
              double normalize = logsumexp({a[m_all], b[m_all]});
              a[m_all] = a[m_all] - normalize;
              b[m_all] = b[m_all] - normalize;
              double temp_wi = logsumexp({a[m_all] + log_mutation + name[0],
                                          b[m_all] + log_mutation + name[1]}) +
                               wi[m_all];
              sm1(index[j], m_all) = temp_wi;

              const double name2[2] = {
                  -log_par[0] + log((general_N - (m + 1)) * mu2),
                  -log_par[1] + log((general_N - (m + 3)) * mu2)};

              temp_wi = logsumexp({a[m_all] + log_mutation + name2[0],
                                   b[m_all] + log_mutation + name2[1]}) +
                        wi[m_all];
              // sum_mutation[1] = logsumexp(sum_mutation[1], temp_wi);
              // sm2[index[j]] = temp_wi;
              sm2(index[j], m_all) = temp_wi;

              temp_wi = logsumexp({a[m_all] + log_len + log((m + 1)),
                                   b[m_all] + log_len + log((m + 3))}) +
                        wi[m_all];
              // sum_len[0] = logsumexp(sum_len[0], temp_wi);
              // sl1[index[j]] = temp_wi;
              sl1(index[j], m_all) = temp_wi;

              temp_wi = logsumexp({a[m_all] + log_len + log(general_N - (m + 1)),
                                   b[m_all] + log_len + log(general_N - (m + 3))}) +
                  wi[m_all];
              // sum_len[1] = logsumexp(sum_len[1], temp_wi);
              // sl2[index[j]] = temp_wi;
              sl2(index[j], m_all) = temp_wi;

              temp_wi = logsumexp({a[m_all] + log_mutation - log_par[0] + log(intercept),
                                   b[m_all] + log_mutation - log_par[1] + log(intercept)}) +
                  wi[m_all];
              // sum_intercept[0] = logsumexp(sum_intercept[0], temp_wi);
              // si1[index[j]] = temp_wi;
              si1(index[j], m_all) = temp_wi;

              temp_wi = logsumexp({a[m_all] + log_len, b[m_all] + log_len}) +
                        wi[m_all];
              // sum_intercept[1] = logsumexp(sum_intercept[1], temp_wi);
              // si2[index[j]] = temp_wi;
              si2(index[j], m_all) = temp_wi;

              m_all++;
            }
          }
        }

        // for (size_t m = 0; m < N_all; m++) {
        //   sum_p1[m] = logsumexp({sum_p1[m], a[m]});
        //   sum_p2[m] = logsumexp({sum_p2[m], b[m]});
        // }
        // for (size_t m = 0; m < N_all; m++) {
        //     sum_outer_p[m] = logsumexp({sum_outer_p[m], wi[m]});
        // }
      }
    }
    // mu1 = exp(sum_mutation[0] - sum_len[0]);
    // mu2 = exp(sum_mutation[1] - sum_len[1]);

    mu1 = exp(logsumexp(sm1) - logsumexp(sl1));
    mu2 = exp(logsumexp(sm2) - logsumexp(sl2));

    // if(!stochastic) {
    //   mu1 = exp(logsumexp(sm1) - logsumexp(sl1));
    //   mu2 = exp(logsumexp(sm2) - logsumexp(sl2));
    // }
    // else {
    //   if (i == 0) { // First iteration
    //     logsumexp_buffer[0] = logsumexp(sm1);
    //     logsumexp_buffer[1] = logsumexp(sl1);
    //     logsumexp_buffer[2] = logsumexp(sm2);
    //     logsumexp_buffer[3] = logsumexp(sl2);
    //   }
    // }
    

    // if (estimate_inner_pi) {
    //   for (size_t m = 0; m < N_all; m++) {
    //     pi1[m] = exp(sum_p1[m]) / n;
    //     pi2[m] = exp(sum_p2[m]) / n;
    //   }
    // }
    if (estimate_intercept) {
      // intercept = exp(sum_intercept[0] - sum_intercept[1]);
      intercept = exp(logsumexp(si1) - logsumexp(si2));
      // fprintf(stderr, "Intercept : %f\n", intercept * 1e9);
    }
    // if (estimate_outer_p) {
        // for (size_t m = 0; m < N_all; m++) {
        //     outer_p[m] = exp(sum_outer_p[m]) / n;
        // }
    // }

    if (fabs(mu1 - prev_mu1) + fabs(mu2 - prev_mu2) < eps) {
      break;
    }
    prev_mu1 = mu1;
    prev_mu2 = mu2;
  }
  // NumericVector result = NumericVector(3+outer_p.size());
  NumericVector result = NumericVector(3);
  result[0] = mu1;
  result[1] = mu2;
  result[2] = intercept;
  // for (size_t i = 0; i < outer_p.size(); i++) {
  //   result[i+3] = outer_p[i];
  // }

//   if (estimate_inner_pi) {
//     fprintf(stderr, "Pi:\n");
//     for (size_t m = 0; m < N_all; m++) {
//       fprintf(stderr, "%d: (%f, %f) \n", m, pi1[m], pi2[m]);
//     }
//     fprintf(stderr, "\n");
//   }

//   if (estimate_outer_p) {
//     fprintf(stderr, "Pi:\n");
//     for (size_t m = 0; m < N_all; m++) {
//       fprintf(stderr, "%d: %f \n", m, outer_p[m]);
//     }
//     fprintf(stderr, "\n");
//   }
  // if (estimate_intercept) {
  //   fprintf(stderr, "Intercept: %f\n", intercept * 1e9);
  // }
  return (result);
}