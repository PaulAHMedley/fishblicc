#' Example small sample of length frequency data
#'
#' A data list created using simulated length frequency the function `blicc_dat`:
#' eg_ld <- blicc_dat(LLB = 10:35, a=1.0e-4, b=2.95, L50=24,
#'                    fq=c(0,1,2,5,26,70,72,66,36,34,25,24,20,10,12,5,3,5,6,4,2,0,2,0,1,0),
#'                    Linf=c(32, 3))
#'
#'@format A list of 26 data vectors used in the BLICC model:
#'\describe{
#'  \item{NK}{Number of nodes in the Gauss-Laguerre quadrature}
#'  \item{oBN}{Number of bins in the length frequency}
#'  \item{Len}{Vector of length bin lower boundaries}
#'  \item{LMP}{Length bin mid-points, used for selectivity and plotting}
#'  \item{fq}{Vector of length frequencies}
#'  \item{a}{Length-weight parameter in a L^b. Optional}
#'  \item{b}{Length-weight parameter in a L^b.}
#'  \item{Ls}{Logistic steepness parameter for the maturity ogive}
#'  \item{Lm}{Logistic parameter for 50% maturity}
#'  \item{Flat_top}{0 == Flat-topped selectivty; Otherwise selectivity is domed}
#'  \item{poLinfm}{Linf: Mean prior normal parameter}
#'  \item{poLinfs}{Linf: sd prior normal parameter}
#'  \item{polGam}{Galpha: Mean prior normal log parameter}
#'  \item{polGas}{Galpha: sd prior normal log parameter}
#'  \item{polMkm}{Mk: Mean prior normal log parameter}
#'  \item{polMks}{Mk: sd prior normal log parameter}
#'  \item{polFkm}{Fk: Mean prior normal log parameter}
#'  \item{polFks}{Fk: sd prior normal log parameter}
#'  \item{polSmxm}{Smx: Mean prior normal log parameter}
#'  \item{polSmxs}{Smx: sd prior normal log parameter}
#'  \item{polSs1m}{Ss1: Mean prior normal log parameter}
#'  \item{polSs1s}{Ss1: sd prior normal log parameter}
#'  \item{polSs2m}{Ss2: Mean prior normal log parameter}
#'  \item{polSs2s}{Ss2: sd prior normal log parameter}
#'  \item{polNB_phim}{NB_phi: Mean prior log normal parameter}
#'  \item{polNB_phis}{NB_phi: sd prior log normal parameter}
#'}
#'@source Data are simulated.
"eg_ld"

#'Example stanfit object from the BLICC model (Stan MCMC)
#'
#'A stanfit object from fitting the BLICC model using `BLICC_fit` to the "eg_ld"
#'data set. The draws can be extracted and reference points
#'calculated using function `blicc_ref_pts`. The MCMC simulations are also
#'accessible through a variety of `rstan` and `posterior` functions.
#'  eg_slim <- blicc_fit(eg_ld, thin=2)
#'
#'@format A `stanfit` object: variables can be obtained using `rstan::extract` and
#'other methods. see packages `rstan` and `posterior` for details.
#'
#'@source Data are simulated.
"eg_slim"

#'Example parameter draws and reference points from the BLICC model
#'
#'An example data table of parameter draws and reference points from the BLICC model
#'  eg_res_df <- blicc_ref_pts(eg_slim, eg_ld)
#'
#'@format A data frame of results of 21 variables from a BLICC model fit:
#'\describe{
#'  \item{Linf}{Mean maximum length from the von Bertalanffy growth model}
#'  \item{Galpha}{Gamma distribution parameter governing growth variability}
#'  \item{Mk}{Natural mortality (per unit K time)}
#'  \item{Fk}{Fishing mortality (per unit K time)}
#'  \item{Smx}{Selectivity mode (full selectivity)}
#'  \item{Ss1}{Left slope for the double normal selectivity (=1/sigma^2)}
#'  \item{Ss2}{Right slope for the double normal selectivity (=1/sigma^2)}
#'  \item{NB_phi}{Negative binomial parameter: Excess variance compared to the Poisson}
#'  \item{Gbeta}{Gamma distribution "rate" parameter: (=Galpha/Linf)}
#'  \item{SPR}{Spawning potential ratio}
#'  \item{lp__}{Log posterior probability for the draw}
#'  \item{.chain}{MCMC chain identifier}
#'  \item{.iteration}{MCMC iteration identifier}
#'  \item{.draw}{MCMC draw identifier}
#'  \item{F20}{Fishing mortality estimated to achieve 20% SPR}
#'  \item{F30}{Fishing mortality estimated to achieve 30% SPR}
#'  \item{F40}{Fishing mortality estimated to achieve 40% SPR}
#'  \item{F01}{Fishing mortality at 10% of the yield curve slope at the origin (F0.1)}
#'  \item{S20}{Selectivity mode estimated to achieve 20% SPR}
#'  \item{S40}{Selectivity mode estimated to achieve 40% SPR}
#'  \item{SMY}{Selectivity mode estimated to achieve maximum yield per recruit}
#'}
#'
#'@source Data are simulated.
"eg_rp"

