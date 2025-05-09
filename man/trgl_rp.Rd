% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{trgl_rp}
\alias{trgl_rp}
\title{A list of results including reference points from the \link{blicc_ref_pts} function}
\format{
A list of data frames of results and input information for a BLICC model fit:
\describe{
\item{dr_df}{A draws tibble of parameter estimates}
\item{lx_df}{A tibble of expected numbers of fish in each length bin for each gear}
\item{ld}{The data list from \link{blicc_dat} used to fit the model and calculate indicators and results}
\item{scenario}{A list of time period, gears, a vector describing the "direction" across gears and data list used to calculate reference points}
\item{rp_df}{A draws tibble of parameter estimates and reference points for the scenario}
}

rp_df is a tibble of parameter estimates from a BLICC model fit:
\describe{
\item{Linf}{Mean maximum length from the von Bertalanffy growth model}
\item{Galpha}{Gamma distribution parameter governing growth variability}
\item{Mk}{Natural mortality (per unit K time)}
\item{Fk}{A list vector of fishing mortality (per unit K time)}
\item{Sm}{A list vector of selectivity model parameters}
\item{NB_phi}{Negative binomial parameter: Excess variance compared to the Poisson}
\item{Gbeta}{Gamma distribution "rate" parameter: (=Galpha/Linf)}
\item{.chain}{MCMC chain identifier}
\item{.iteration}{MCMC iteration identifier}
\item{.draw}{MCMC draw identifier}
\item{F20}{List vector of fishing mortality estimated to achieve 20\% SPR}
\item{F40}{List vector of fishing mortality estimated to achieve 40\% SPR}
\item{S40}{List vector of selectivity location estimated to achieve 40\% SPR}
\item{SMY}{List vector of selectivity location estimated to achieve maximum yield per recruit}
}
}
\source{
Data are simulated.
}
\usage{
trgl_rp
}
\description{
An example data table of parameter draws and reference points from the BLICC model
trgl_rp <- blicc_ref_pts(trgl_slim, trgl_ld)
}
\keyword{datasets}
