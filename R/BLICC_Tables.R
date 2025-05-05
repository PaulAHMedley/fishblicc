# BLICC Table Functions -----------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

# To Do
# parameter names for tables
# table of reference point estimates
# mixture weights in prior table


#' Return a tibble containing a summary of the fishblicc priors being applied
#'
#' This function provides a summary of the priors contained in a fishblicc data
#' list in a tibble. This can be used to inspect the model's assumptions or to
#' create a table in a report. (Mixture weights have not been included yet).
#'
#' @export
#' @inheritParams blicc_mpd
#' @return A tibble describing the priors being applied.
#' @examples
#' blicc_prior(trgl_ld)
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
  
  gear_names <- c(gear_names, blicc_ld$fqname[blicc_ld$Fkq>0])
  function_type <- c(function_type, rep("Lognormal", blicc_ld$NF))
  par_names <- c(par_names, rep("Fk", blicc_ld$NF))
  
  if (blicc_ld$NM > 0 | blicc_ld$NS != blicc_ld$NG)
    sel_comp <- paste("Sel", as.character(1:blicc_ld$NS))
  else
    sel_comp <- blicc_ld$gname
  
  
  sel_par <- vector()
  for (i in 1:blicc_ld$NS) {
    npar <- sel_fun$npar[blicc_ld$fSel[i]]-1
    par_names <- c(par_names, sel_fun$par_names[[blicc_ld$fSel[i]]])
    
    
    gear_names <- c(gear_names,
                    rep(sel_comp[i], npar+1))
    
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
#' blicc_results(trgl_slim)
#' 
blicc_results <- function(blicc_res) {
  Fk = Linf = Parameter = Rhat = Sm = Value = lp__ = median = mpd = NULL
  n_eff = par = sd = se = slim = NULL
  `2.5%` = `97.5%` = SPR = YPR = B_B0 = NULL
  
  if (class(blicc_res)[1]=="stanfit") {
    NF <- blicc_res@par_dims$nFk
    NP <- blicc_res@par_dims$nSm
    NT <- blicc_res@par_dims$SPR
    params <- c(
      "Linf",
      "Galpha",
      "Mk",
      paste0("Fk[", as.character(1:NF), "]"),
      paste0("Sm[", as.character(1:NP), "]"),
      "NB_phi",
      "Gbeta",
      paste0("SPR[", as.character(1:NT), "]"),
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
    if ("dr_df lx_df ld scenario rp_df" != paste(names(blicc_res), collapse=" ")) {
      stop(
        "The provided parameter must be a stanfit object, mpd fit or a results list from blicc_ref_pts() function."
      )
    }
    suppressWarnings(
      res <- blicc_res$dr_df |>
        dplyr::select(Linf:lp__, SPR, B_B0, YPR) |>
        tidyr::unnest_wider(col=Fk, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("Fk")) |>
        tidyr::unnest_wider(col=Sm, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("Sm")) |>
        tidyr::unnest_wider(col=SPR, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("SPR")) |>
        tidyr::unnest_wider(col=YPR, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("YPR")) |>
        tidyr::unnest_wider(col=B_B0, names_sep="[") |>
        dplyr::rename_with(~ paste0(.x, "]"), tidyselect::starts_with("B_B0"))
    )
    p_order <- names(res)
    res <- res |>
      tidyr::pivot_longer(cols=tidyselect::everything(), names_to="Parameter",
                          values_to="Value") 
    if (nrow(blicc_res$dr_df) > 1) {
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


#' Impact from the removal of each gear on the YPR and SPR in the specified
#' period.
#'
#' The change in YPR and SPR are calculated for each parameter set to produce a
#' tibble of the YPR and SPR change from the removal of each gear in sequence.
#' This is only useful if there is more than one gear.
#'
#' @export
#' @inheritParams blicc_expected_catches
#' @return A tibble with the changes in YPR and SPR from the removal of each
#'   gear for the scenario in blicc result list.
#' @examples
#' blicc_impact(blicc_ref_pts(blicc_mpd(trgl_ld), trgl_ld))
#'   
blicc_impact <- function(blicc_rp) {
  tp_ld <- blicc_rp$scenario$time_period_ld
  rp_df <- blicc_rp$rp_df
  
  suppressWarnings(
    dr_df <- blicc_rp$dr_df |>
      dplyr::select(.draw, YPR, SPR) |>
      dplyr::mutate(curValues = purrr::pmap(list(YPR, SPR), \(x, y) c(x, y))) |>
      dplyr::select(.draw, curValues)
  )
  suppressWarnings(
    rp_df <- rp_df |>
      dplyr::left_join(dr_df, by = ".draw") |>
      dplyr::mutate(
        Impact = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm, curValues),
          fPRImpact,
          blicc_ld = tp_ld,
          .progress = "Gear Impacts"
        )
      ) |>
      dplyr::select(.draw, Impact) |>
      tidyr::unnest(col = Impact)
  )
  return(rp_df)
}


