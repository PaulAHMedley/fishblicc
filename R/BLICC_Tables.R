# BLICC Table Functions -----------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

# To Do
# Function returning parameter names for tables
# Function table of reference point estimates


#' Return a tibble containing a summary of the fishblicc priors being applied
#'
#' This function provides a summary of the priors contained in a fishblicc data
#' object in a tibble. This can be used to inspect the model's assumptions or to
#' create a table in a report. (Weights have not been included yet).
#'
#' @export
#' @inheritParams blicc_mpd
#' @return A tibble describing the priors being applied.
#' @examples
#' blicc_prior(eg_ld)
#' 
blicc_prior <- function(blicc_ld) {

  sel_fun <- Rsel_functions()

  gear_names <- c(NA, NA, NA)
  function_type <- c("Normal", "Lognormal", "Lognormal")
  par_names <- c("Linf", "Galpha", "Mk")
  Mean <- with(blicc_ld, c(poLinfm, exp(c(polGam, polMkm))))
  Mu <- with(blicc_ld, c(poLinfm, polGam, polMkm))
  SD <- with(blicc_ld, c(poLinfs, polGas, polMks))

  if (blicc_ld$ref_length > 0) {
    gear_names <- c(gear_names, NA)
    function_type <- c(function_type, "Ref. length")
    par_names <- c(par_names, NA)
    Mean <- with(blicc_ld, c(Mean, ref_length))
    Mu <- c(Mu, NA)
    SD <- c(SD, NA)
  }
  
  gear_names <- c(gear_names, blicc_ld$gname[blicc_ld$Fkg])
  function_type <- c(function_type, rep("Lognormal", blicc_ld$NF))
  par_names <- c(par_names, rep("Fk", blicc_ld$NF))

  sel_par <- vector()
  for (i in 1:blicc_ld$NS) {
    npar <- sel_fun$npar[blicc_ld$fSel[i]]-1
    par_names <- c(par_names, sel_fun$par_names[[blicc_ld$fSel[i]]])
    
    gear_names <- c(gear_names,
                    #blicc_ld$gname[i],
                    rep(NA, npar+1))
    function_type <- c(function_type,
                       "Lognormal",
                       rep(NA, npar))
  }

  gear_names <- c(gear_names, rep(NA, 4))
  function_type <- c(function_type, "Lognormal", rep(NA, 3))
  par_names <- c(par_names, "NB_phi", "b", "L50", "Ls")
  Mean <- with(blicc_ld, c(Mean, exp(c(polFkm, polSm[1:blicc_ld$NP], polNB_phim)),
                          b, L50, Ls))
  Mu <- with(blicc_ld, c(Mu, polFkm, polSm[1:blicc_ld$NP], polNB_phim,
                        b, L50, Ls))
  SD <- with(blicc_ld, c(SD, rep(polFks, NF), polSs[1:blicc_ld$NP], polNB_phis,
                        NA, NA, NA))

  return(tibble::tibble(
    Gear = gear_names,
    Parameter = par_names,
    `Function Type` = function_type,
    Mean = Mean,
    Mu = Mu,
    SD = SD)
  )
}


#' Return a tibble containing a summary of the results from a fishblicc fit
#'
#' This function provides a summary of results from the fitted model in a form
#' suitable for displaying in a table. Note that other packages, such as rstan,
#' bayesplot and posterior can extract information from the stanfit object and
#' the [blicc_ref_pts] `rp_df` 'draws' object in more detailed form.
#'
#' @export
#' @param blicc_res Results from [blicc_fit] (stanfit object), [blicc_mpd] or
#'   [blicc_ref_pts]
#' @return A tibble summarising results of the model fit object.
#' @examples
#' blicc_results(eg_rp)
#' 
blicc_results <- function(blicc_res) {
  Fk=Linf=Parameter=Rhat=Sm=Value=lp__=median=mpd=n_eff=par=sd=se=slim=NULL
  `2.5%` = `97.5%` = NULL
  
  if (class(blicc_res)[1]=="stanfit") {
    NF <- blicc_res@par_dims$nFk
    NP <- blicc_res@par_dims$nSm
    params <- c(
      "Linf",
      "Galpha",
      "Mk",
      paste0("Fk[", as.character(1:NF), "]"),
      paste0("Sm[", as.character(1:NP), "]"),
      "NB_phi",
      "Gbeta",
      "SPR",
      "lp__"
    )
    par_value <- rstan::summary(blicc_res, pars = params, probs = 0.5)$summary
    res <- tibble::as_tibble(par_value) |>
      dplyr::mutate(Parameter = rownames(par_value)) |>
      dplyr::select(Parameter,
             Mean = mean,
             SD = sd,
             `2.5%`,
             `97.5%`,
             `N (eff)` = n_eff,
             Rhat)
  } else if (tibble::is_tibble(blicc_res) & paste(names(blicc_res), collapse=" ")=="par mpd se") {
    res <- blicc_res |>
      dplyr::rename(Parameter=par, `Max. Posterior` = mpd, `SE`=se)
  } else {
    if ("vdir rp_df lx_df ld" != paste(names(blicc_res), collapse=" ")) {
      stop(
        "The provided parameter must be a stanfit object, mpd fit or a results list from blicc_ref_pts() function."
      )
    }
    suppressWarnings(
      res <- blicc_res$rp_df |>
        dplyr::select(Linf:lp__) |>
        tidyr::unnest_wider(col=Fk, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("Fk")) |>
        tidyr::unnest_wider(col=Sm, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("Sm"))
    )
    p_order <- names(res)
    res <- res |>
      tidyr::pivot_longer(cols=tidyselect::everything(), names_to="Parameter",
                          values_to="Value") 
    if (nrow(blicc_res$rp_df) > 1) {
      res <- res |>
        dplyr::mutate(Parameter = factor(Parameter, levels=p_order)) |>
        dplyr::group_by(Parameter) |>
        dplyr::summarise(Mean=mean(Value), SD=sd(Value),
                         `10%`= stats::quantile(Value, probs = 0.1),
                         `50%`= median(Value),
                         `90%` = stats::quantile(Value, probs = 0.9)) |>
        dplyr::ungroup()
    }
    # dplyr::summarise(dplyr::across(Value, list(Mean=mean, SD=sd,
      #                                            `10%`= ~quantile(., probs = 0.1),
      #                                            `50%`=median, `Q90%` = ~quantile(., probs = 0.9)),
      #                                .names = "{.fn}"))
  }
  return(res)
}
