% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{eg_rp}
\alias{eg_rp}
\title{Example parameter draws and reference points from the BLICC model}
\format{
A data frame of results of 21 variables from a BLICC model fit:
\describe{
\item{Linf}{Mean maximum length from the von Bertalanffy growth model}
\item{Galpha}{Gamma distribution parameter governing growth variability}
\item{Mk}{Natural mortality (per unit K time)}
\item{Fk}{Fishing mortality (per unit K time)}
\item{Smx}{Selectivity mode (full selectivity)}
\item{Ss1}{Left slope for the double normal selectivity (=1/sigma^2)}
\item{Ss2}{Right slope for the double normal selectivity (=1/sigma^2)}
\item{NB_phi}{Negative binomial parameter: Excess variance comapred to the Poisson}
\item{Gbeta}{Gamma distribution "rate" parameter: (=Galpha/Linf)}
\item{SPR}{Spawning potential ratio}
\item{lp__}{Log posterior probability for the draw}
\item{.chain}{MCMC chain identifier}
\item{.iteration}{MCMC iteration identifier}
\item{.draw}{MCMC draw identifier}
\item{F20}{Fishing mortality estimated to achieve 20\% SPR}
\item{F30}{Fishing mortality estimated to achieve 30\% SPR}
\item{F40}{Fishing mortality estimated to achieve 40\% SPR}
\item{F01}{Fishing mortality at 10\% of the yield curve slope at the origin (F0.1)}
\item{S20}{Selectivity mode estimated to achieve 20\% SPR}
\item{S40}{Selectivity mode estimated to achieve 40\% SPR}
\item{SMY}{Selectivity mode estimated to achieve maximum yield per recruit}
}
}
\source{
Data are simulated.
}
\usage{
eg_rp
}
\description{
An example data table of parameter draws and reference points from the BLICC model
eg_res_df <- blicc_ref_pts(eg_slim, eg_ld)
}
\keyword{datasets}