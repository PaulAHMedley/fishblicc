# fishblicc Table Functions ------------------------------------------------


#' Returns a tibble containing a summary of the fishblicc priors being applied
#'
#' This function provides a summary of the priors contained in a
#' fishblicc data object in a tibble. This can be used, for example, to
#' create a table in a report.
#'
#' @export
#' @inheritParams blicc_mpd
#' @return A tibble describing the priors being applied.
#' @examples
#' table_prior(eg_ld)
#'
table_prior <- function(blicc_ld) {

  sel_fun <- Rsel_functions()

  gear_names <- c(NA, NA, NA, blicc_ld$gname[blicc_ld$Fkg])
  function_type <- c("Normal",
                     rep("Lognormal",
                         2+blicc_ld$NF))
  sel_par <- vector()
  for (i in 1:blicc_ld$NG) {
    npar <- sel_fun$npar[blicc_ld$fSel[i]]-1
    sel_par <- c(sel_par, sel_fun$par_names[[blicc_ld$fSel[i]]])
    gear_names <- c(gear_names,
                    blicc_ld$gname[i],
                    rep(NA, npar))
    function_type <- c(function_type,
                       "Lognormal",
                       rep(NA, npar))
    }

  gear_names <- c(gear_names, rep(NA, 4))
  function_type <- c(function_type, "Lognormal", rep(NA, 3))

  par_names <- c("Linf", "Galpha", "Mk", rep("Fk", blicc_ld$NF),
                 sel_par, "NB_phi", "b", "L50", "Ls")

  return(tibble::tibble(
    Gear = gear_names,
    Parameter = par_names,
    `Function Type` = function_type,
    Mean = with(blicc_ld, c(poLinfm, exp(c(polGam, polMkm, polFkm, polSm, polNB_phim)),
                            b, L50, Ls)),
    Mu = with(blicc_ld, c(poLinfm, polGam, polMkm, polFkm, polSm, polNB_phim,
                          b, L50, Ls)),
    SD = with(blicc_ld, c(poLinfs, polGas, polMks, rep(polFks, NF), polSs, polNB_phis,
                          NA, NA, NA)))
  )
}
