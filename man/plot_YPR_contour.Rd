% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_Plots.R
\name{plot_YPR_contour}
\alias{plot_YPR_contour}
\title{Yield-per-recruit surface contour plot}
\usage{
plot_YPR_contour(blicc_rp, blicc_ld)
}
\arguments{
\item{blicc_rp}{A posterior draws and reference points tibble
from \code{blicc_ref_pts} function.}

\item{blicc_ld}{A data list from the \code{blicc_dat} function.}
}
\value{
ggplot geom object for plotting
}
\description{
The plot shows the YPR surface plotted against the fishing mortality (Fk) and
the full selectivity (Smx). The current state of the stock is marked as a
point and the F0.1 reference point is marked as a line. All values are
calculated as means and uncertainty is not represented to maintain
simplicity.
}
\examples{
plot_YPR_contour(eg_rp, eg_ld)

}
