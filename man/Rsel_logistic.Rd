% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_Model.R
\name{Rsel_logistic}
\alias{Rsel_logistic}
\title{Calculate a logistic selectivity curve for a length vector}
\usage{
Rsel_logistic(Sp, LMP)
}
\arguments{
\item{Sp}{Vector (length 2) of the selectivity parameters (50\%, and slope)
for the logistic selectivity function.}

\item{LMP}{A vector of lengths (usually mid-points for length bins)}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{ A vector of selectivity values varying from 0.0 to 1.0
}\if{html}{\out{</div>}}
}
\description{
A logistic takes a 50\% selectivity and a slope parameter. The resulting
selectivity varies from 0 to 1.0, with 1.0 being the asymptote. It is
calculated for each position defined in the length vector.
}
\examples{
Sel <- Rsel_logistic(Sp = c(30, 0.5), LMP = seq(15.5, 55.5, by=1.0))
plot(y=Sel, x=seq(15.5, 55.5, by=1.0))

}
