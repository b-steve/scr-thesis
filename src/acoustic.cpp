#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// =================================================================================== //
// =================================================================================== //

// ============================== //
//          Acoustic NLL          //
// ============================== //
/*
* Calculating the Euclidean distance between a point and each trap.
* Returns a vector of distances.
*/
// [[Rcpp::export]]
double scr_nll_acoustic(NumericVector pars,
                        NumericMatrix caps,
			NumericMatrix toa_ssq,
                        NumericMatrix traps,
                        NumericMatrix mask,
                        NumericMatrix maskDists,
                        NumericVector nCalls) {
  /*
   *  Storing/initialising (starting) parameter values.
   *  - Note that parameters are back-transformed
   */
  double D = exp(pars[0]);
  double g0 = exp(pars[1]) / (1 + exp(pars[1]));
  double sigma = exp(pars[2]);
  double lambda_c = exp(pars[3]);
  double sigma_toa = exp(pars[4]);

  // Number of animals
  int nAnimals = nCalls.size();
  // Number of calls detected in total.
  //int n = caps.nrow();
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
  * - Line that fills in maskProbs(i, j) is the Hazard Half-Normal function (HHN)
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

  /*
  *  Probability of detecting on specific call emitted from s
  *
  */
  NumericVector pDetected = 1 - pAvoid;

  /*
  * Probability of detecting at least one call on at least one microphone
  *  from an individual located at s
  */
  NumericVector pAnimal = 1 - exp(-lambda_c * pDetected);

  // ========================================= //
  // ========================================= //
  NumericVector fCapt(nAnimals);
  // Probability for an animal being at each mask point
  NumericVector logfCapt_givenNS(nMask);
  NumericVector logfn_givenS(nMask);
  double subRow = 0;

  // Looping through all animals
  for (int i = 0; i < nAnimals; i++) {
    /*
     * Subsetting the capture matrix
     * - Sub-matrix: all cols; first row of sub-mat --- nCalls - 1
     */
    NumericMatrix subCaps = caps(Range(subRow, subRow + nCalls[i] - 1), _);
    NumericMatrix subtoas = toa_ssq(Range(subRow, subRow + nCalls[i] - 1), _);
    subRow += nCalls[i];

    // Looping through each mask point
    for (int j = 0; j < nMask; j++) {
      logfn_givenS[j] = log(R::dpois(nCalls[i], lambda_c * pDetected[j], 0) + DBL_MIN);

      // Looping through the calls (each sub-matrix)
      for (int k = 0; k < nCalls[i]; k++) {
        logfCapt_givenNS[j] = -log(pDetected[j] + DBL_MIN);
	if (use_toa){
	  logfCapt_givenNS[j] += (1 - n_dets(i))*log(sigma_toa) - (toa_ssq(i, j)/(2*pow(sigma_toa, 2)));
	}
        // Looping through each trap
        for (int m = 0; m < traps.nrow(); m++) {
          logfCapt_givenNS[j] += R::dbinom(subCaps(k, m), 1, maskProbs(j, m), 1);
        }
      }
    }
    fCapt[i] = sum(exp(logfn_givenS + logfCapt_givenNS));
  }

  // ========================================= //
  // ========================================= //

  /*
  * Log-likelihood contribution from all capture histories
  * - Calculated by log of sum of individual likelihood contributions.
  */
  double logfCapt = sum(log(fCapt + DBL_MIN));

  // Calculating effective survey area (unused in likelihood).
  double esa = area * sum(pAnimal);

  // Log-likelihood contribution from number of animals detected.
  double logf_n = R::dpois(nAnimals, D * esa, 1);

  // Overall log-likelihood.
  double logLik = logf_n + logfCapt - nAnimals * log(sum(pAnimal));

  // Returning log-likelihood
  return -logLik;
}
// =================================================================================== //
// =================================================================================== //


