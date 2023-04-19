#' Example small sample of length frequency data
#'
#' A data list created using simulated length frequency data from the
#' function `blicc_dat()`:
#'
#' eg_ld <- blicc_dat(
#'   LLB = 10:35,
#'   fq=c(0,1,2,5,26,70,72,66,36,34,25,24,20,
#'        10,12,5,3,5,6,4,2,0,2,0,1,0),
#'   sel_fun = "dsnormal",
#'   gear_names="Trawl",
#'   Linf=c(32, 3),
#'   a=1.0e-4, b=2.95, L50=21)
#'
#'@format A list of 26 data vectors used in the BLICC model:
#'\describe{
#'  \item{mname}{Model and/or fit and/or species name}
#'  \item{NG}{Number of gear (selectivity) models}
#'  \item{NF}{Number of fishing mortality estimates}
#'  \item{NB}{Number of bins in the length frequency}
#'  \item{gname}{Fishing gear (selectivity model) names}
#'  \item{LLB}{Vector of length bin lower boundaries}
#'  \item{LMP}{Length bin mid-points, used for selectivity and plotting}
#'  \item{fq}{List of length frequency vectors one for each gear}
#'  \item{Catch}{Vector of proportional catch by gear of length NF}
#'  \item{Fkg}{Integer vector of index for gears link to fishing mortality}
#'  \item{fSel}{Index of selectivity model to use for each gear}
#'  \item{NP}{Number of selectivity parameters}
#'  \item{Pmx}{Maximum number of selectivity parameters for a single model}
#'  \item{spar}{Integer matrix of index linking gear models to selectivity parameters}
#'  \item{a}{Length-weight parameter in a L^b. Optional}
#'  \item{b}{Length-weight parameter in a L^b.}
#'  \item{L50}{Logistic parameter for 50% maturity}
#'  \item{Ls}{Logistic steepness parameter for the maturity ogive}
#'  \item{wt_L}{Vector of average weight in each length bin}
#'  \item{ma_L}{Vector of mature biomass in each length bin}
#'  \item{poLinfm}{Linf: Mean prior normal parameter}
#'  \item{poLinfs}{Linf: sd prior normal parameter}
#'  \item{polGam}{Galpha: Mean prior normal log parameter}
#'  \item{polGas}{Galpha: sd prior normal log parameter}
#'  \item{polMkm}{Mk: Mean prior normal log parameter}
#'  \item{polMks}{Mk: sd prior normal log parameter}
#'  \item{polFkm}{Fk: Mean prior normal log parameter}
#'  \item{polFks}{Fk: sd prior normal log parameter}
#'  \item{polSm}{Sm: Mean priors normal log selectivity parameter}
#'  \item{polSs}{Sm: sd priors normal log selectivity parameter}
#'  \item{polNB_phim}{NB_phi: Mean prior log normal parameter}
#'  \item{polNB_phis}{NB_phi: sd prior log normal parameter}
#'  \item{polCs}{Catch: sd prior log normal parameter on relative catch}
#'  \item{NK}{Number of nodes in the Gauss-Laguerre quadrature}
#'  \item{gl_nodes}{Gauss-Laguerre quadrature nodes}
#'  \item{gl_weights}{Gauss-Laguerre quadrature weights}
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

#'A list of results including reference points from the `blicc_ref_pts()` function
#'
#'An example data table of parameter draws and reference points from the BLICC model
#'  eg_rp <- blicc_ref_pts(eg_slim, eg_ld)
#'
#'@format A list of data frames of results from a BLICC model fit:
#'\describe{
#'  \item{vdir}{A vector describing the "direction" across gears for RP estimates}
#'  \item{rp_df}{A draws tibble of parameter estimates and reference points}
#'  \describe{
#'    \item{Linf}{Mean maximum length from the von Bertalanffy growth model}
#'    \item{Galpha}{Gamma distribution parameter governing growth variability}
#'    \item{Mk}{Natural mortality (per unit K time)}
#'    \item{Fk[.]}{A list vector of fishing mortality (per unit K time)}
#'    \item{Sm[.]}{A list vector of selectivity model parameters}
#'    \item{NB_phi}{Negative binomial parameter: Excess variance compared to the Poisson}
#'    \item{Gbeta}{Gamma distribution "rate" parameter: (=Galpha/Linf)}
#'    \item{SPR}{Spawning potential ratio}
#'    \item{lp__}{Log posterior probability for the draw}
#'    \item{.chain}{MCMC chain identifier}
#'    \item{.iteration}{MCMC iteration identifier}
#'    \item{.draw}{MCMC draw identifier}
#'    \item{F20}{Fishing mortality estimated to achieve 20% SPR}
#'    \item{F30}{Fishing mortality estimated to achieve 30% SPR}
#'    \item{F40}{Fishing mortality estimated to achieve 40% SPR}
#'    \item{F01}{Fishing mortality at 10% of the yield curve slope at the origin (F0.1)}
#'    \item{S20}{Selectivity mode estimated to achieve 20% SPR}
#'    \item{S40}{Selectivity mode estimated to achieve 40% SPR}
#'    \item{SMY}{Selectivity mode estimated to achieve maximum yield per recruit}
#'  }
#'  \item{lx_df}{A tibble of expected numbers of fish in each length bin}
#'  \item{ld}{The data object used to fit the model and calculate indicators and results }
#'}
#'
#'@source Data are simulated.
"eg_rp"

