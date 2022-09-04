
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fishblicc

<!-- badges: start -->
<!-- badges: end -->

The goal of fishblicc is to provide tools to fit a catch curve model to
individual length frequency samples, analogous to an age-based catch
curve. However, the model can account for length-based selectivity and
multiple gears. The model is fitted using MCMC in Stan (mc-stan.org) to
account for uncertainty. The approach could be useful for data-limited
stock assessments, such as assessing stock status at end of projects
where length data have been collected, or continuous monitoring of
low-data bycatch species as part of a harvest strategy. The model also
provides tools to assess how effective length sampling is in estimating
quantities of interest.

## Installation

You can install the development version of fishblicc like so:

``` r
if (require("remotes")) install.packages("remotes")
install_github("PaulAHMedley/fishblicc")
```

## Model Description

This is the Bayesian Length Interval Catch Curve project to estimate
mortality from length frequency samples that form a snapshot of a
fishery. Length samples are not treated as a time series. This is
explicitly a “data-limited” method designed to be used as a risk
assessment approach rather than provide definitive information on stock
status. Although the model does not require a data time series, it
assumes the data are a snapshot of the catch length composition from a
fishery in a stationary state. It is suitable for end-of-project
evaluation of length frequency data that may have been collected over a
short time period (e.g. a year) or to monitor status of many species
that may take up a small proportion of the catches.  
The assessment is carried out in R using Stan (mc-stan.org) to carry out
the MCMC.

In common with other such methods, there are significant assumptions
that will not be met in practice, therefore the focus of this method is
robustness and looking for ways to reduce the impact of model structural
errors or at least spot when such errors may invalidate the results.

The method can be used to evaluate stock status based on the spawning
potential ratio. In general, this uses the ratio of large (older,
mature) fish in a sample compared to small (immature young) fish. If the
ratio is too small compared to what might be expected in a sample of a
population that is unfished, the population will be at high risk of
overexploitation.

The model used here does not convert length to age, but instead accounts
for mortality as fish grow through length intervals. The estimation
method works by assuming that the von Bertalanffy growth model describes
the length at age and variation in individual fish growth is governed by
their maximum length which is drawn from a gamma probability density
function (constant CV). The model is flexible enough to include all
data, so, for example, observed lengths above $L_\infty$ are not
rejected.

The model fits the following 8 parameters: - Linf the asymptotic mean
length  
- Galpha the growth model inverse error parameter
($CV = Galpha^{-0.5}$) - Mk the natural mortality in time units of K,
the growth rate  
- Fk a scale parameter for the fishing mortality in time units of K, the
growth rate  
- Smx the maximum point for the double-sided normal selectivity model -
Ss1 the lower-side slope parameter for the double-sided normal
selectivity model (i.e. 1/variance) - Ss2 the upper-side slope parameter
for the double-sided normal selectivity model (equal to zero for a
flat-topped selectivity. - phi the over-dispersion parameter of the
counts in the length bins for the negative binomial.

The length catch curve is used to estimate the spawning potential ratio.
To do this, estimates of natural mortality, length-weight and maturity
are required. These cannot be estimated from a length frequency, so have
to be provided as priors from other sources. Being data-limited, only
the length-weight exponent parameter (b) and the logistic maturity at
length model parameters (length at 50% maturity and logistic steepness)
are required. The natural mortality is then inferred relying on life
history invariants as suggested by Prince et al. Hordyk et al. (2015?).
The additional parameters are: - $L_m$ and $L_s$ for the maturity at 50%
and maturity steepness for a logistic maturity curve.  
- $b$ parameter for the length weight relationship $W=aL^b$

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("fishblicc")
## Prepare some data in the required format
ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
## Fit the model to some data
slim <- blicc_fit(ld)
## Calculate reference points
res_df <- blicc_ref_pts(slim, ld)
## Get the model values by length
l_df <- blicc_expect_len(res_df, ld)
## Because the fit takes a little time, it might be best to save the results
save(ld, slim, l_df, file="fishblicc_analysis.rda") 
## Summarise the results
summary(res_df)
plot_expected_frequency(l_df, ld) #Plot the results to check the model fit
```
