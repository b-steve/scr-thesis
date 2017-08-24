#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// ============================== //
//        eucdist (matrix)        //
// ============================== //
/*
 * Calculating the Euclidean distance between a point and each trap.
 * Returns a vector of distances.
 * - Not exported.
 */
// [[Rcpp::export]]
NumericMatrix eucdist_nll(NumericMatrix points,
                          NumericMatrix traplocations) {
  NumericMatrix dists(points.nrow(), traplocations.nrow());
  for(int i = 0;  i < points.nrow(); i++) {
    for(int j = 0; j < traplocations.nrow(); j++) {
      dists(i, j) = sqrt(pow(points(i, 0) - traplocations(j, 0), 2.0)
                         + pow(points(i, 1) - traplocations(j, 1), 2.0));
    }
  }
  return dists;
}

// =================================================================================== //
// =================================================================================== //
// ============================== //
//              distr             //
// ============================== //
/*
* Returns distribution function for a specified distribution
* - Not exported.
*/
//// [[Rcpp::export]]
/*
double distr(String distribution,
             double x,
             double lambda,
             double p,
             double n) {
  // Initialising the density value
  double dens;

  // Returning different densities based on `distr`
  if(Rcpp::tolower(distribution) == "pois") {
    dens = R::dpois(x, lambda, 1);
  } else if(::tolower(distribution) == "bin") {
    dens = R::dbinom(x, n, p, 1);
  }

  return dens;
}
*/
// ============================== //
//            scr_nll             //
// ============================== //
/*
* Calculating the Euclidean distance between a point and each trap.
* Returns a vector of distances.
*/
// [[Rcpp::export]]
double scr_nll(NumericVector pars,
               NumericMatrix caps,
               NumericMatrix traps,
               NumericMatrix mask,
               NumericMatrix maskDists) {
  // Storing/initialising (starting) parameter values.
  double D = exp(pars[0]);
  double g0 = 1 / (1 + exp(pars[1]));
  double sigma = exp(pars[2]);

  // Number of animals detected.
  int n = caps.nrow();
  // Number of traps. NB: NOT USED
  //int nTraps = traps.nrow();
  // Number of mask points.
  int nMask = mask.nrow();
  // Area of a single mask pixel.
  double area = mask.attr("area");

  /*
   * Constructing distance matrix.
   * - Element (i, j) gives dist. b/w ith mask pint and jth trap.
   */
  //NumericMatrix maskDists = eucdist_nll(mask, traps);

  /*
   * Constructing a detection probability matrix.
   * - Element (i, j) gives prob. of animal @ ith mask pt. being detected @ jth trap.
   */
  NumericMatrix maskProbs(maskDists.nrow(), maskDists.ncol());
  for(int i = 0; i < maskDists.nrow(); i++) {
    for(int j = 0; j < maskDists.ncol(); j++) {
      maskProbs(i, j) = g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0)));
    }
  }

  /*
   * Constructing a detection probability vector
   * - ith element = P(animal @ ith mask pt. is detected by >= 1 trap)
   */
  NumericVector pAvoid(maskProbs.nrow());
  for(int i = 0; i < maskProbs.nrow(); i++) {
    pAvoid[i] = 1 - maskProbs(i, 0);
    for(int j = 1; j < maskProbs.ncol(); j++) {
      pAvoid[i] *= 1 - maskProbs(i, j);
    }
  }

  // Vector of probability of detection at each mask point.
  NumericVector pDetected = 1 - pAvoid;

  /*
   * Calculating likelihood contribution
   * - Contribution from each detected animal's capt. hist.
   */
  NumericVector fCapt(n);
  NumericVector logfCapt_givenS(nMask);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < nMask; j++) {
      /*
       * Calculating log-probability of animal's capt. hist,
       * conditional on being at the jth mask point.
       * - Note that 'caps' and 'maskProbs' both have the same ncol();
       *    - i.e. # traps
       *  - Also note: R::dbinom(double x, double n, double p, int lg)
       *    - Where 'int lg' is 0 = F, 1 = T.
       */
      logfCapt_givenS[j] = R::dbinom(caps(i, 0), 1, maskProbs(j, 0), 1);
      for(int k = 1; k < caps.ncol(); k++) {
        logfCapt_givenS[j] += R::dbinom(caps(i, k), 1, maskProbs(j, k), 1);
      }
    }
    // Summing probabilities over all mask points.
    fCapt[i] = sum(exp(logfCapt_givenS));
  }

  /*
   * Log-likelihood contribution from all capture histories
   * - Calculated by log of sum of individual likelihood contributions.
   */
  double logfCapt = sum(log(fCapt + DBL_MIN));

  // Calculating effective survey area (unused in likelihood).
  double esa = area * sum(pDetected);

  // Log-likelihood contribution from number of animals detected.
  double logf_n = R::dpois(n, D * esa, 1);

  // Overall log-likelihood.
  double logLik = logf_n + logfCapt - n * log(sum(pDetected));


  return -logLik;
}

// =================================================================================== //
// =================================================================================== //

