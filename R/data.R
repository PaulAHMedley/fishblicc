#' Example small sample of length frequency data from gillnet fishery
#'
#' A data list created using simulated length frequency data from the
#' function [blicc_dat]:
#'
#'gillnet_ld <- fishblicc::blicc_dat(model_name = "Simulated Gillnet Data",
#'                                   LLB = seq(14, 62, by=2), 
#'                                   fq=list(Gillnet = c(0, 2, 7, 18, 32, 56, 115, 155, 227, 240, 246, 249, 
#'                                                       210, 158, 106, 73, 43, 28, 7, 6, 4, 3, 0, 1, 0)),
#'                                   sel_fun = "normal",
#'                                   gear_names = "Gillnet",
#'                                   Linf=c(64, 3),
#'                                   Mk=1.55,
#'                                   a=1.0e-4, b=2.95, L50=34, L95=43)
#'
#'@format A list of 26 data vectors used in the BLICC model:
#'\describe{
#'  \item{model_name}{Model and/or fit and/or species name}
#'  \item{NQ}{Number of length frequency data sets}
#'  \item{NG}{Number of gears}
#'  \item{NS}{Number of selectivity functions}
#'  \item{NT}{Number of time periods}
#'  \item{NF}{Number of fishing mortality estimates}
#'  \item{NB}{Number of bins in the length frequency}
#'  \item{NP}{Number of selectivity function parameters}
#'  \item{NM}{Number of mixture selectvity models}
#'  \item{fqname}{Length frequency data set names}
#'  \item{gname}{Fishing gear names}
#'  \item{tpname}{Time period names}
#'  \item{LLB}{Vector of length bin lower boundaries}
#'  \item{LMP}{Length bin mid-points, used for selectivity and plotting}
#'  \item{fq}{List of length frequency vectors, one for each gear}
#'  \item{Gi}{Integer vector linking length frequency data sets to gears}
#'  \item{Ti}{Integer vector linking length frequency data sets to time periods}
#'  \item{prop_catch}{Vector of proportional catch by gear of length NF}
#'  \item{Fkq}{Integer vector of length frequency data sets index linking to fishing mortality}
#'  \item{fSel}{Index of selectivity function to use for each selectvity model}
#'  \item{GSbase}{Index of based selectivity model for each gear}
#'  \item{GSmix1}{Index of start and end mixtures in GSmix2: 0 implies no mixture}
#'  \item{GSmix2}{Index of all selectivity mixture models referenced by gear in GSmix1}
#'  \item{poLinfm}{Linf: Mean prior normal parameter}
#'  \item{poLinfs}{Linf: sd prior normal parameter}
#'  \item{polGam}{Galpha: Mean prior normal log parameter}
#'  \item{polGas}{Galpha: sd prior normal log parameter}
#'  \item{a}{Length-weight parameter in a L^b. Optional}
#'  \item{b}{Length-weight parameter in a L^b.}
#'  \item{L50}{Logistic parameter for 50% maturity}
#'  \item{Ls}{Logistic steepness parameter for the maturity ogive}
#'  \item{wt_L}{Vector of average weight in each length bin}
#'  \item{ma_L}{Vector of mature biomass in each length bin}
#'  \item{polMkm}{Mk: Mean prior normal log parameter}
#'  \item{polMks}{Mk: sd prior normal log parameter}
#'  \item{polFkm}{Fk: Mean prior normal log parameter}
#'  \item{polFks}{Fk: sd prior normal log parameter}
#'  \item{sp_i}{Integer index of the parameter start in polSm for each selectivity function}
#'  \item{sp_e}{Integer index of the parameter end in polSm for each selectivity function}
#'  \item{polSm}{Sm: Mean priors normal log selectivity parameter including mixture weights at the end}
#'  \item{polSs}{Sm: sd priors normal log selectivity parameter including mixture weights at the end}
#'  \item{polNB_phim}{NB_phi: Mean prior log normal parameter}
#'  \item{polNB_phis}{NB_phi: sd prior log normal parameter}
#'  \item{polCs}{Catch: sd prior log normal parameter on relative catch}
#'  \item{NK}{Number of nodes in the Gauss-Laguerre quadrature}
#'  \item{gl_nodes}{Gauss-Laguerre quadrature nodes}
#'  \item{gl_weights}{Gauss-Laguerre quadrature weights}
#'}
#'@source Data are simulated.
"gillnet_ld"

#' Example small sample of length frequency data from trawl/gillnet fishery
#'
#' A data list created using simulated length frequency data from the
#' function [blicc_dat]:
#'
#'trgl_ld <- fishblicc::blicc_dat(model_name = "Simulated trawl/gillnet data",
#'                                LLB = seq(2, 76, by=2), 
#'                                fq=list(Trawl = c(0, 1, 0, 2, 11, 12, 39, 71, 111, 103, 116, 84, 91, 
#'                                                  69, 66, 69, 60, 56, 32, 33, 31, 27, 35, 28, 24, 19, 
#'                                                  21, 12, 10, 7, 5, 10, 3, 0, 0, 0, 1, 0),
#'                                        Gillnet  = c(0, 0, 0, 3, 5, 7, 11, 14, 22, 38, 49, 78, 106, 112, 
#'                                                     122, 98, 90, 78, 57, 49, 35, 22, 11, 4, 3, 1, 0, 2, 
#'                                                     0, 1, 0, 0, 0, 0, 0, 0, 0, 0)), 
#'                                sel_fun = c("logistic", "normal"),
#'                                gear_names = c("Trawl", "Gillnet"),
#'                                Catch = c(0.55, 0.45),
#'                                Linf=c(65, 3),
#'                                Mk=1.55,
#'                                a=1.0e-4, b=2.95, L50=50, L95=63)
#'
#'@format A list of 26 data vectors used in the BLICC model:
#'\describe{
#'  \item{model_name}{Model and/or fit and/or species name}
#'  \item{NQ}{Number of length frequency data sets}
#'  \item{NG}{Number of gears}
#'  \item{NS}{Number of selectivity functions}
#'  \item{NT}{Number of time periods}
#'  \item{NF}{Number of fishing mortality estimates}
#'  \item{NB}{Number of bins in the length frequency}
#'  \item{NP}{Number of selectivity function parameters}
#'  \item{NM}{Number of mixture selectvity models}
#'  \item{fqname}{Length frequency data set names}
#'  \item{gname}{Fishing gear names}
#'  \item{tpname}{Time period names}
#'  \item{LLB}{Vector of length bin lower boundaries}
#'  \item{LMP}{Length bin mid-points, used for selectivity and plotting}
#'  \item{fq}{List of length frequency vectors, one for each gear}
#'  \item{Gi}{Integer vector linking length frequency data sets to gears}
#'  \item{Ti}{Integer vector linking length frequency data sets to time periods}
#'  \item{prop_catch}{Vector of proportional catch by gear of length NF}
#'  \item{Fkq}{Integer vector of length frequency data sets index linking to fishing mortality}
#'  \item{fSel}{Index of selectivity function to use for each selectvity model}
#'  \item{GSbase}{Index of based selectivity model for each gear}
#'  \item{GSmix1}{Index of start and end mixtures in GSmix2: 0 implies no mixture}
#'  \item{GSmix2}{Index of all selectivity mixture models referenced by gear in GSmix1}
#'  \item{poLinfm}{Linf: Mean prior normal parameter}
#'  \item{poLinfs}{Linf: sd prior normal parameter}
#'  \item{polGam}{Galpha: Mean prior normal log parameter}
#'  \item{polGas}{Galpha: sd prior normal log parameter}
#'  \item{a}{Length-weight parameter in a L^b. Optional}
#'  \item{b}{Length-weight parameter in a L^b.}
#'  \item{L50}{Logistic parameter for 50% maturity}
#'  \item{Ls}{Logistic steepness parameter for the maturity ogive}
#'  \item{wt_L}{Vector of average weight in each length bin}
#'  \item{ma_L}{Vector of mature biomass in each length bin}
#'  \item{polMkm}{Mk: Mean prior normal log parameter}
#'  \item{polMks}{Mk: sd prior normal log parameter}
#'  \item{polFkm}{Fk: Mean prior normal log parameter}
#'  \item{polFks}{Fk: sd prior normal log parameter}
#'  \item{sp_i}{Integer index of the parameter start in polSm for each selectivity function}
#'  \item{sp_e}{Integer index of the parameter end in polSm for each selectivity function}
#'  \item{polSm}{Sm: Mean priors normal log selectivity parameter including mixture weights at the end}
#'  \item{polSs}{Sm: sd priors normal log selectivity parameter including mixture weights at the end}
#'  \item{polNB_phim}{NB_phi: Mean prior log normal parameter}
#'  \item{polNB_phis}{NB_phi: sd prior log normal parameter}
#'  \item{polCs}{Catch: sd prior log normal parameter on relative catch}
#'  \item{NK}{Number of nodes in the Gauss-Laguerre quadrature}
#'  \item{gl_nodes}{Gauss-Laguerre quadrature nodes}
#'  \item{gl_weights}{Gauss-Laguerre quadrature weights}
#'}
#'@source Data are simulated.
"trgl_ld"

#'Example stanfit object from the BLICC model (Stan MCMC)
#'
#'A stanfit object from fitting the BLICC model using [blicc_fit] to the "eg_ld"
#'data set. The draws can be extracted and reference points
#'calculated using function [blicc_ref_pts]. The MCMC simulations are also
#'accessible through a variety of `rstan` and `posterior` functions.
#'  trgl_slim <- blicc_fit(trgl_ld, nwarmup=500, thin=2)
#'
#'@format A `stanfit` object: variables can be obtained using `rstan::extract` and
#'other methods. see packages `rstan` and `posterior` for details.
#'
#'@source Data are simulated.
"trgl_slim"

#'A list of results including reference points from the [blicc_ref_pts] function
#'
#'An example data table of parameter draws and reference points from the BLICC model
#'  trgl_rp <- blicc_ref_pts(trgl_slim, trgl_ld)
#'
#'@format A list of data frames of results and input information for a BLICC model fit:
#'\describe{
#'  \item{dr_df}{A draws tibble of parameter estimates}
#'  \item{lx_df}{A tibble of expected numbers of fish in each length bin for each gear}
#'  \item{ld}{The data list from [blicc_dat] used to fit the model and calculate indicators and results}
#'  \item{scenario}{A list of time period, gears, a vector describing the "direction" across gears and data list used to calculate reference points}
#'  \item{rp_df}{A draws tibble of parameter estimates and reference points for the scenario}
#'}
#'@format rp_df is a tibble of parameter estimates from a BLICC model fit:
#'\describe{
#'  \item{Linf}{Mean maximum length from the von Bertalanffy growth model}
#'  \item{Galpha}{Gamma distribution parameter governing growth variability}
#'  \item{Mk}{Natural mortality (per unit K time)}
#'  \item{Fk}{A list vector of fishing mortality (per unit K time)}
#'  \item{Sm}{A list vector of selectivity model parameters}
#'  \item{NB_phi}{Negative binomial parameter: Excess variance compared to the Poisson}
#'  \item{Gbeta}{Gamma distribution "rate" parameter: (=Galpha/Linf)}
#'  \item{.chain}{MCMC chain identifier}
#'  \item{.iteration}{MCMC iteration identifier}
#'  \item{.draw}{MCMC draw identifier}
#'  \item{F20}{List vector of fishing mortality estimated to achieve 20% SPR}
#'  \item{F40}{List vector of fishing mortality estimated to achieve 40% SPR}
#'  \item{S40}{List vector of selectivity location estimated to achieve 40% SPR}
#'  \item{SMY}{List vector of selectivity location estimated to achieve maximum yield per recruit}
#'}
#'
#'@source Data are simulated.
"trgl_rp"


