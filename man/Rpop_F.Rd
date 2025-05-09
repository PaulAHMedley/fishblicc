% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_Model.R
\name{Rpop_F}
\alias{Rpop_F}
\title{Calculate the population size and fishing mortality applied within the length
bins for the BLICC model}
\usage{
Rpop_F(Galpha, Gbeta, Mk, Fk, FSel, blicc_ld)
}
\arguments{
\item{Galpha}{Alpha parameter for the Gamma probability density function that
governs growth variability.}

\item{Gbeta}{Rate parameter for the Gamma distribution growth variability
(=Galpha/Linf)}

\item{Mk}{Natural mortality divided by the growth rate K}

\item{Fk}{Fishing mortality divided by the growth rate K for each gear
making a contribution}

\item{FSel}{A list of all the selectivities for each length bin}

\item{blicc_ld}{A standard data list suitable for the model (see function
\link{blicc_dat})}
}
\value{
A list of the population size in each length bin for each time period
and a list of fishing mortalities at length for each gear / time period.
}
\description{
The model integrates over the growth rate variability and sums mortality
piece-wise over intervening length intervals to derive the relative
population number in each length bin taking into account variation in growth
and varying mortality-at-length. This is analogous to an age based catch
curve model, but for length. The model uses the Gauss-Laguerre quadrature
rule for integration, so nodes and weights for this must be provided. The
function returns the expected population size in each length bin, and the
fishing mortalities.
}
\examples{
Sel <- Rselectivities(exp(trgl_ld$polSm), trgl_ld)
S <- Rpop_F(100, 100/50, Mk=1.5, Fk=exp(trgl_ld$polFkm),
            FSel=Sel, blicc_ld=trgl_ld)
plot(y=S$N_L[[1]], x=trgl_ld$LMP, type="l")

}
