#include <Rcpp.h>
using namespace Rcpp;

// BLICC model functions in C++

//' Calculate a double-sided normal selectivity curve for a length vector
//'
//' Slightly faster C++ version.
//' A double-sided normal takes two separate slope (sigma^2) parameters
//' around the mode parameter (mu) that change the slope independently
//' on either side of the mode. The selectivity varies from 0 to 1.0,
//' with 1.0 being at the mode. It is calculated for each position defined
//' in the length vector. Setting the second slope parameter (Ss2) to zero
//' creates a flat-top selectivity.
//'
//' @param  LMP  A vector of lengths (usually mid-points for length bins)
//' @param  Smx  Mode of the normal selectivity function
//' @param  Ss1  Left side slope, parameterized as 1/sigma^2
//' @param  Ss2  Right side slope, parameterized as 1/sigma^2. Zero implies
//' flat-topped selectivity
//' @return      A vector of selectivity values varying from 0.0 to 1.0
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::NumericVector Csel_dsnormal(Rcpp::NumericVector LMP,
                                  double Smx,
                                  double Ss1,
                                  double Ss2) {
  int n = LMP.size();
  Rcpp::NumericVector SL(n);
  for(int i = 0; i < n; ++i) {
    if (LMP[i] < Smx)
      SL[i] = exp(-Ss1 * pow((LMP[i]-Smx), 2.0));
    else
      SL[i] = exp(-Ss2 * pow((LMP[i]-Smx), 2.0));
  }
  return SL;
}


//' Calculate the survival of a fish cohort to sequential length boundaries
//'
//' This C++ version of the `RSurvival_Est()` function.
//' The model integrates over growth rate variability and sums mortality
//' piece-wise over intervening length intervals to derive the
//' survival to each length bin boundary taking into account variation in growth
//' and varying mortality-at-length. This is analogous to an age based catch
//' curve model, but for length. Survival can be used to estimate population
//' numbers-at-length and catch length frequency, among other things.
//' The model uses the Gauss-Laguerre quadrature rule for integration,
//' so nodes and weights for this must be provided.
//'
//' @param  node   Nodes for the Gauss-Laguerre quadrature rule
//' @param  wt     Weights for the Gauss-Laguerre quadrature rule
//' @param  Len    Vector of lower-bounds for length intervals
//' @param  Zki    Vector of total mortality (in time units of growth rate
//' parameter K) to be applied in each length interval.
//' This vector should be the same length as Len.
//' @param  Galpha Alpha parameter for the Gamma probability density function
//' that governs growth variability.
//' @param  Gbeta  Beta parameter for the Gamma probability density function
//' that governs growth variability.
//' @return A vector of proportion surviving to each length bin lower bound.
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::NumericVector CSurvival_Est(Rcpp::NumericVector node,
                                  Rcpp::NumericVector wt,
                                  Rcpp::NumericVector Len,
                                  Rcpp::NumericVector Zki,
                                  double Galpha,
                                  double Gbeta)  {
  int nv = node.size();
  int LN = Len.size();
  double lgamma_Galpha = lgamma(Galpha);
  double Galpha_1 = Galpha - 1.0;
  Rcpp::NumericVector surv(LN);
  Rcpp::NumericVector x_beta(nv);
  Rcpp::NumericVector log_x_beta(nv);
  Rcpp::NumericVector ss(nv);
  Rcpp::NumericVector Zii(LN);
  Zii[0] = -Zki[0];
  for(int i = 1; i < (LN-1); ++i) {
    Zii[i] = (Zki[i-1] - Zki[i]);
  }

  for (int i = 0; i < nv; ++i) {
    x_beta[i] = node[i] / Gbeta;
    log_x_beta[i] = log(x_beta[i]);
    surv[0] += exp(log(node[i] + Gbeta*Len[0]) * Galpha_1 -
      Gbeta*Len[0] -lgamma_Galpha) * wt[i];
    surv[1] += exp(log(x_beta[i] + Len[1] - Len[0]) * Zii[0] +
      log_x_beta[i] * Zki[0] +
      log(node[i] + Gbeta * Len[1]) * Galpha_1 -
      Gbeta * Len[1] - lgamma_Galpha) * wt[i];
  }

  for  (int Li = 2; Li < LN; ++Li) {
    double Ln = Len[Li];
    double v1 = Gbeta * Ln + lgamma_Galpha;
    Rcpp::NumericVector Lrange(LN);

    for (int Lii = 0; Lii < Li; ++Lii) {
      Lrange[Lii] = Ln - Len[Lii];
    }

    for (int i = 0; i < nv; ++i) {
      double lim = 0;
      for (int Lii = 0; Lii < Li; ++Lii) {
        // sum the previous log mortality for this length group
        lim += log(x_beta[i] + Lrange[Lii])*Zii[Lii];
      }
      surv[Li] += exp(lim + log_x_beta[i] * Zki[Li - 1] +
        log(node[i] + Gbeta * Ln) * Galpha_1 - v1) * wt[i];
    }

  }
  return surv;
}


//' Calculate the population relative numbers-at-length within length bins
//'
//' C++ Version of the `RNinInterval()` function.
//' The relative numbers in each length bin is based on the definite integral
//' over the bin growth transition time:
//'     N_i = (S_i - S_{i+1})/Z_{ki}
//'
//' @param  surv Survival to the lower boundary of each length bin
//' (see `RSurvival_est`)
//' @param  Zki  Total mortality within each bin (time in units of the
//' growth rate K)
//' @return A vector of relative numbers of fish within each length bin.
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::NumericVector CNinInterval(Rcpp::NumericVector surv,
                                 Rcpp::NumericVector Zki) {
  int n = surv.size();
  Rcpp::NumericVector NI(n);
  int n1 = n-1;

  for(int i = 0; i < n1; ++i) {
    NI[i] = (surv[i] - surv[i+1]) / Zki[i];
  }
  NI[n1] = surv[n1] / Zki[n1];
  return NI;
}



