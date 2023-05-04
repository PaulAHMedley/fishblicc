# BLICC Plotting Functions ------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
#' Plot the observed and expected frequency data
#'
#' The graph shows the observed frequency and the expected frequency for the
#' specified gears. The plots include the 80% credible intervals for the
#' estimate (if an MCMC is provided) and 80% credible intervals for the
#' observations. The latter takes into account the negative binomial observation
#' error, so on average 8 out of ten observations would be expected within these
#' bounds.
#'
#' @export
#' @param blicc_rp  A list of posterior draws, reference points with associated
#'   direction, the data object and expected lengths from [blicc_ref_pts]
#'   function.
#' @param Gear Specifies the gear to plot
#' @return ggplot geom object plotting observed and expected frequency
#' @examples
#' plot_expected_frequency(eg_rp)
#'
plot_expected_frequency <-
  function(blicc_rp, Gear = NA) {
    .draw = Lgroup = LMP = fq = NB_phi = NULL  # Not necessary but stops CMD check notes
    efq = fq_lo = fq_hi = efq_m = efq_lo = efq_hi = dat_lo = dat_hi = NULL
    Sgroup = ofq = NULL

    rp_df <- blicc_rp$rp_df
    blicc_ld <- blicc_rp$ld
    blicc_lx <- blicc_rp$lx_df

    Gear <- parse_gear(Gear, blicc_ld)
    # if (any(is.na(Gear))) return()

    gear2plot <- blicc_ld$gname[Gear]
    blicc_lx <- blicc_lx |>
      dplyr::filter(Sgroup %in% gear2plot)

    dat_df <- with(blicc_ld,
                   tibble::tibble(
                     Sgroup = rep(gname, each = NB),
                     Lgroup = factor(rep(LLB, NG)),
                     LMP = rep(LMP, NG),
                     ofq = unlist(fq)
                   ))

    if (nrow(rp_df) == 1) {
      # MPD estimate
      MCMC_draws <- FALSE
      df1 <- rp_df |>
        dplyr::select(.draw, NB_phi) |>
        dplyr::right_join(blicc_lx, by = ".draw") |>
        dplyr::select(.draw, Sgroup, Lgroup, efq_m = efq, NB_phi) |>
        dplyr::mutate(
          dat_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq_m),
          dat_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq_m)
        ) |>
        dplyr::left_join(dat_df, by = c("Sgroup", "Lgroup"))
    } else {
      MCMC_draws <- TRUE
      suppressWarnings(
        df1 <- rp_df |>
          dplyr::select(.draw, NB_phi) |>
          dplyr::left_join(blicc_lx, by = ".draw") |>
          dplyr::select(.draw, Sgroup, Lgroup, efq, NB_phi) |>
          dplyr::mutate(
            fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
            fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
          ) |>
          dplyr::group_by(Sgroup, Lgroup) |>
          dplyr::summarise(
            efq_m = mean(efq),
            efq_lo = stats::quantile(efq, probs = c(0.1), names = F),
            efq_hi = stats::quantile(efq, probs = c(0.9), names = F),
            dat_lo = stats::quantile(fq_lo, probs = c(0.1), names = F),
            dat_hi = stats::quantile(fq_hi, probs = c(0.9), names = F),
          ) |>
          dplyr::ungroup() |>
          dplyr::left_join(dat_df, by = c("Sgroup", "Lgroup"))
      )
    }

    gp <- ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
      ggplot2::geom_col(ggplot2::aes(x = LMP, y = ofq),
                        fill = "lightblue",
                        alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(x = LMP, y = ofq),
                          shape = "-",
                          colour = "black") +
      ggplot2::geom_line(ggplot2::aes(y = efq_m))

    if (MCMC_draws) {
      gp <- gp + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = efq_lo, ymax = efq_hi),
        fill = "blue",
        alpha = 0.3
      )
    }

    gp <- gp +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
                           fill = "black",
                           alpha = 0.2) +
      ggplot2::labs(x = "Length", y = "Frequency") +
      ggplot2::geom_hline(yintercept = 0)

    if (length(Gear) > 1) {
      gp <- gp + ggplot2::facet_wrap(ggplot2::vars(Sgroup),
                            ncol = 2,
                            scales = "free_y")
    } else {
      gp <- gp + ggplot2::ggtitle(gear2plot)
    }
    return(gp)
  }

#' Plot standardised residuals against length
#'
#' The graph shows the standardised residuals plotted against length for all
#' gears combined, separated by colour. Due to the variance weighting, outliers
#' can be expected at the tails of the frequency.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object plotting standardised residuals
#' @examples
#' plot_residuals(eg_rp)
#'
plot_residuals <- function(blicc_rp) {
  .draw = NB_phi = Sgroup = Lgroup = efq = ofq = LMP = std_res = NULL

  blicc_ld <- blicc_rp$ld
  blicc_lx <- blicc_rp$lx_df
  dat_df <-
    with(blicc_ld,
         tibble::tibble(
           Sgroup = rep(gname, each = NB),
           Lgroup = factor(rep(LLB, NG)),
           LMP = rep(LMP, NG),
           ofq = unlist(fq)
         ))
  suppressWarnings(
    df1 <- blicc_rp$rp_df |>
      dplyr::select(.draw, NB_phi) |>
      dplyr::right_join(blicc_lx, by = ".draw") |>
      dplyr::select(.draw, Sgroup, Lgroup, efq, NB_phi) |>
      dplyr::left_join(dat_df, by = c("Sgroup", "Lgroup")) |>
      dplyr::mutate(std_res = (ofq - efq) / sqrt(efq + (efq ^ 2) / NB_phi))
  )
  gp <-
    ggplot2::ggplot(df1, ggplot2::aes(x = LMP, y = std_res, colour = Sgroup)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Length", y = "Standardised Residuals", colour = "Gear") +
    ggplot2::geom_hline(yintercept = 0)
  return(gp)
}


#' Plot selectivity and 80% credible interval, if available, by
#' length.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_selectivity(eg_rp)
#'
plot_selectivity <- function(blicc_rp) {
  Lgroup = sel = LMP = sel_01 = sel_90 = sel_25 = sel_75 = sel_m = NULL
  Sgroup = sel_10 = NULL

  blicc_ld <- blicc_rp$ld
  blicc_lx <- blicc_rp$lx_df
  gp <- blicc_lx |>
    dplyr::select(Sgroup, Lgroup, sel) |>
    dplyr::group_by(Sgroup, Lgroup) |>
    dplyr::summarise(
      sel_m = mean(sel),
      sel_10 = stats::quantile(sel, probs = c(0.10), names = F),
      sel_90 = stats::quantile(sel, probs = c(0.90), names = F),
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)]) |>
    ggplot2::ggplot(ggplot2::aes(x = LMP, fill = Sgroup, colour = Sgroup)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = sel_10, ymax = sel_90),
                         alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = sel_m)) +
    ggplot2::labs(
      x = "Length",
      y = "Selectivity",
      fill = "Gear",
      colour = "Gear"
    )
  return(gp)
}


#' Plot the estimated spawning potential ratio probability density, if available.
#'
#' This plot requires the MCMC to have been run.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_SPR_density(eg_rp)
#'
plot_SPR_density <- function(blicc_rp) {
  SPR = NULL
  rp_df <- blicc_rp$rp_df
  if (nrow(rp_df) <= 100) {
    stop("To obtain a density, you will need to obtain sufficient values (>100) from MCMC.")
  }
  ggplot2::ggplot(rp_df, ggplot2::aes(SPR)) +
    ggplot2::geom_density(fill = "lightblue") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0.2),
                        colour = "red",
                        linetype = 2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0.4),
                        colour = "green",
                        linetype = 2) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Spawning Potential Ratio", y = "Probability Density")
}



#' Plot fishing mortality relative to the SPR 40% reference point (Fk/F40)
#'
#' This plot requires the MCMC to have been run.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_FkF40_density(eg_rp)
#'
plot_FkF40_density <- function(blicc_rp) {
  Fk = F40 = `Fk/F40` = id = NULL # Not necessary but stops CMD check notes

  suppressWarnings(
    rp_df <- blicc_rp$rp_df |>
      dplyr::select(Fk, F40) |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidyr::unnest(cols = c(Fk, F40)) |>
      dplyr::mutate(`Fk/F40` = Fk / F40) |>
      dplyr::group_by(id) |>
      dplyr::summarise(`Fk/F40` = mean(`Fk/F40`)) |>
      dplyr::ungroup()
  )
  if (nrow(rp_df) <= 100) {
    stop("To obtain a density, you will need to obtain sufficient values (>100) from MCMC.")
  }

  StatAreaUnderDensity <- ggplot2::ggproto(
    "StatAreaUnderDensity",
    ggplot2::Stat,
    required_aes = "x",
    compute_group = function(data,
                             scales,
                             xlim = NULL,
                             n = 50) {
      fun <- stats::approxfun(stats::density(data$x))
      ggplot2::StatFunction$compute_group(data,
                                          scales,
                                          fun = fun,
                                          xlim = xlim,
                                          n = n)
    }
  )

  stat_aud <- function(mapping = NULL,
                       data = NULL,
                       geom = "area",
                       position = "identity",
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE,
                       n = 50,
                       xlim = NULL,
                       ...) {
    ggplot2::layer(
      stat = StatAreaUnderDensity,
      data = data,
      mapping = mapping,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(xlim = xlim, n = n, ...)
    )
  }

  fq_good <- sum(rp_df$`Fk/F40` <= 1.0)

  if (fq_good == nrow(rp_df)) {
    gp <- ggplot2::ggplot(rp_df, ggplot2::aes(`Fk/F40`)) +
      ggplot2::geom_density(fill = "#204030")
  } else if (fq_good == 0) {
    gp <- ggplot2::ggplot(rp_df, ggplot2::aes(`Fk/F40`)) +
      ggplot2::geom_density(fill = "#951826")
  } else {
    gp <- ggplot2::ggplot(rp_df, ggplot2::aes(`Fk/F40`)) +
      ggplot2::geom_density(fill = "#951826") +
      stat_aud(
        geom = "area",
        ggplot2::aes(fill = "good"),
        xlim = c(0, 1),
        alpha = 1.0,
        show.legend = FALSE
      )  +
      ggplot2::scale_fill_manual(values = c("good" = "#204030"))
  }

  gp <- gp +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 1.0),
                        colour = "green",
                        linetype = 1) +
    ggplot2::labs(y = "Probability Density")
  return(gp)
}


#' Plot expected length frequencies for a range of fishing mortality reference
#' points
#'
#' The facet plot covers the current estimated fishing mortality, together with
#' fishing mortalities required to obtain SPR 20%, SPR 30% and SPR 40%. The
#' graphs show the expected values and 80% credible intervals compared to
#' current observations, and can be used to assess whether length frequencies
#' should be able to detect changes in fishing mortality to these different
#' levels. Note that if reference points do not exist, they are not plotted, so
#' some graphs may be blank except for the data. The graphs show the number of
#' MCMC draws for which the reference points exist in each case.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting expected length frequencies by SPR
#' @examples
#' plot_efq_FRP(eg_rp)
#'
plot_efq_FRP <- function(blicc_rp, Gear = NA) {
  # Not necessary but stops CMD check notes
  Linf = Galpha = Mk = Fk = Sm = NB_phi = Sgroup = NULL
  .draw = Lgroup = Current = SPR20 = SPR40 = SMY = LMP = fq = Harvest_Level =
    NULL
  F20 = F30 = F40 = efq = fq_lo = fq_hi = N = efq_m = dat_lo = dat_hi =
    lab = x = y = NULL

  blicc_ld <- blicc_rp$ld
  blicc_lx <- blicc_rp$lx_df

  Gear <- parse_gear(Gear, blicc_ld)
  if (is.na(Gear))
    return()

  gear2plot <- blicc_ld$gname[Gear]
  blicc_lx <- blicc_lx |>
    dplyr::filter(Sgroup == gear2plot)

  df_dat <-
    tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq[[Gear]])

  suppressWarnings(
    df <- blicc_rp$rp_df |>
      dplyr::mutate(Lgroup = list(factor(blicc_ld$LLB))) |>
      dplyr::mutate(
        SPR20 = purrr::pmap(
          list(Linf, Galpha, Mk, F20, Sm),
          blicc_get_efq,
          Gear_i = Gear,
          blicc_ld = blicc_ld
        ),
        SPR30 = purrr::pmap(
          list(Linf, Galpha, Mk, F30, Sm),
          blicc_get_efq,
          Gear_i = Gear,
          blicc_ld = blicc_ld
        ),
        SPR40 = purrr::pmap(
          list(Linf, Galpha, Mk, F40, Sm),
          blicc_get_efq,
          Gear_i = Gear,
          blicc_ld = blicc_ld
        )
      ) |>
      dplyr::select(`.draw`, NB_phi, Lgroup, SPR20:SPR40) |>
      tidyr::unnest(Lgroup:SPR40) |>
      dplyr::left_join(
        dplyr::select(blicc_lx, .draw, Lgroup, Current = efq),
        by = c(".draw", "Lgroup")
      )
  )

  df1 <- df |>
    tidyr::pivot_longer(cols = SPR20:Current,
                        names_to = "Harvest_Level",
                        values_to = "efq") |>
    dplyr::mutate(
      fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
      fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
    ) |>
    dplyr::group_by(Harvest_Level, Lgroup) |>
    dplyr::summarise(
      efq_m = mean(efq),
      dat_lo = stats::quantile(
        fq_lo,
        probs = c(0.1),
        na.rm = T,
        names = F
      ),
      dat_hi = stats::quantile(
        fq_hi,
        probs = c(0.9),
        na.rm = T,
        names = F
      ),
      N = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])

  txt_df <- df1 |>
    dplyr::group_by(Harvest_Level) |>
    dplyr::summarise(N = dplyr::first(N)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      lab = paste0("N = ", format(N)),
      x = max(LMP = blicc_ld$LMP) * 0.95,
      y = max(df1$dat_hi)
    )

  ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
    ggplot2::geom_col(
      data = df_dat,
      ggplot2::aes(x = LMP, y = fq),
      fill = "lightblue",
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = df_dat,
      ggplot2::aes(x = LMP, y = fq),
      shape = 95,
      colour = "black"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = efq_m)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
                         fill = "black",
                         alpha = 0.2) +
    ggplot2::geom_text(data = txt_df,
                       ggplot2::aes(label = lab, x = x, y = y),
                       size = 2.5) +
    ggplot2::facet_wrap(ggplot2::vars(Harvest_Level), nrow = 2) +
    ggplot2::labs(x = "Length", y = "Frequency", title = gear2plot) +
    ggplot2::geom_hline(yintercept = 0)
}


#' Plot expected length frequencies for a range of selectivity reference points
#'
#' The facet plot covers the current estimated selectivity location parameters,
#' together with selectivity changes required to obtain SPR 20%, SPR 40% and the
#' maximum yield. The graphs show the expected values and 80% credible intervals
#' compared to current observations, and can be used to assess whether length
#' frequencies should be able to detect changes in selectivity to these
#' different levels. Note that if reference points do not exist, they are not
#' plotted, so some graphs may be blank except for data.  The graphs show the
#' number of MCMC draws for which the reference points exist in each case.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object plotting expected length frequencies by
#'   selectivity
#' @examples
#' plot_efq_SRP(eg_rp)
#'
plot_efq_SRP <- function(blicc_rp, Gear = NA) {
  Linf = Galpha = Mk = Fk = Sm = NB_phi = Sgroup = NULL
  .draw = Lgroup = Current = SPR20 = SPR40 = LMP = fq = Harvest_Level =
    NULL
  S20 = S40 = SMY = efq = fq_lo = fq_hi = N = efq_m = dat_lo = dat_hi =
    lab = x = y = NULL

  blicc_ld <- blicc_rp$ld
  blicc_lx <- blicc_rp$lx_df

  Max_Yield = Selectivity = NULL
  blicc_ld <- blicc_rp$ld
  Gear <- parse_gear(Gear, blicc_ld)
  if (is.na(Gear))
    return()

  gear2plot <- blicc_ld$gname[Gear]
  blicc_lx <- blicc_lx |>
    dplyr::filter(Sgroup == gear2plot)

  df_dat <-
    tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq[[Gear]])

  suppressWarnings(
    df <- blicc_rp$rp_df |>
      dplyr::mutate(Lgroup = list(factor(blicc_ld$LLB))) |>
      dplyr::mutate(
        SPR20 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, S20),
          blicc_get_efq,
          Gear_i = Gear,
          blicc_ld = blicc_ld
        ),
        SPR40 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, S40),
          blicc_get_efq,
          Gear_i = Gear,
          blicc_ld = blicc_ld
        ),
        Max_Yield = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, SMY),
          blicc_get_efq,
          Gear_i = Gear,
          blicc_ld = blicc_ld
        )
      ) |>
      dplyr::select(`.draw`, NB_phi, Lgroup, SPR20:Max_Yield) |>
      tidyr::unnest(Lgroup:Max_Yield) |>
      dplyr::left_join(
        dplyr::select(blicc_lx, .draw, Lgroup, Current = efq),
        by = c(".draw", "Lgroup")
      )
  )
  df1 <- df |>
    tidyr::pivot_longer(cols = SPR20:Current,
                        names_to = "Selectivity",
                        values_to = "efq") |>
    dplyr::filter(!is.na(efq)) |>
    dplyr::mutate(
      fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
      fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
    ) |>
    dplyr::group_by(Selectivity, Lgroup) |>
    dplyr::summarise(
      efq_m = mean(efq),
      dat_lo = stats::quantile(
        fq_lo,
        probs = c(0.1),
        na.rm = T,
        names = F
      ),
      dat_hi = stats::quantile(
        fq_hi,
        probs = c(0.9),
        na.rm = T,
        names = F
      ),
      N = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])

  txt_df <- df1 |>
    dplyr::group_by(Selectivity) |>
    dplyr::summarise(N = dplyr::first(N)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      lab = paste0("N = ", format(N)),
      x = max(LMP = blicc_ld$LMP) * 0.95,
      y = max(df1$dat_hi)
    )

  ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
    ggplot2::geom_col(
      data = df_dat,
      ggplot2::aes(x = LMP, y = fq),
      fill = "lightblue",
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = df_dat,
      ggplot2::aes(x = LMP, y = fq),
      shape = 95,
      colour = "black"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = efq_m)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
                         fill = "black",
                         alpha = 0.2) +
    ggplot2::geom_text(data = txt_df,
                       ggplot2::aes(label = lab, x = x, y = y),
                       size = 2.5) +
    ggplot2::facet_wrap(ggplot2::vars(Selectivity), nrow = 2) +
    ggplot2::labs(x = "Length", y = "Frequency", title = gear2plot) +
    ggplot2::geom_hline(yintercept = 0)
}


#' Spawning potential ratio surface contour plot
#'
#' The plot shows the SPR surface plotted against the fishing mortality (Fk) and
#' the 50% or maximum selectivity, dependent on the selectivity function. The
#' current state of the stock is marked as a point and the SPR 40% reference
#' point is marked as a line. All values are calculated as means and uncertainty
#' is not represented to maintain simplicity. The reference line is plotted if
#' it exists within a reasonable range of the current estimated status. For
#' multiple gears, the surface is calculated for a single reference gear.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_SPR_contour(eg_rp)
#'
plot_SPR_contour <- function(blicc_rp) {
  Smx_rp = spr = Lcur = Fcur = `..level..` = NULL  #R CMD Notes
  Smx = dF = level = NULL

  fSPR2 <- function(dF, dL) {
    vFk <- (1 + svdir * dF) * Fk
    vSm[indx] <- (1 + svdir * dL) * Sm[indx]
    Rsel <- Rselectivities(vSm, blicc_ld)
    res <- Rpop_F(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
    return(sum(res$N_L * blicc_ld$ma_L / SPR0))
  }

  blicc_ld <- blicc_rp$ld

  GridN <- 40
  Linf <- mean(blicc_rp$rp_df$Linf)
  Galpha <- mean(blicc_rp$rp_df$Galpha)
  Gbeta <- Galpha / Linf
  Mk <- mean(blicc_rp$rp_df$Mk)
  Fk <- as.vector(tapply(
    X = unlist(dplyr::pull(blicc_rp$rp_df, Fk)),
    INDEX = rep(1:blicc_ld$NF, nrow(blicc_rp$rp_df)),
    FUN = mean
  ))
  Sm <- as.vector(tapply(
    X = unlist(dplyr::pull(blicc_rp$rp_df, Sm)),
    INDEX = rep(1:blicc_ld$NP, nrow(blicc_rp$rp_df)),
    FUN = mean
  ))
  F20 <- tapply(
    X = unlist(dplyr::pull(blicc_rp$rp_df, F20)),
    INDEX = rep(1:blicc_ld$NF, nrow(blicc_rp$rp_df)),
    FUN = mean,
    na.rm = TRUE
  )
  F30 <- tapply(
    X = unlist(dplyr::pull(blicc_rp$rp_df, F30)),
    INDEX = rep(1:blicc_ld$NF, nrow(blicc_rp$rp_df)),
    FUN = mean,
    na.rm = TRUE
  )

  # Decide reference gear on highest fishing mortality to be controlled
  Gear_i <- which.max(Fk * unname(blicc_rp$vdir[blicc_ld$Fkg > 0]))
  gear2plot <- blicc_ld$gname[Gear_i]

  SPR0 <-
    RSPR_0(Galpha, Gbeta, Mk, blicc_ld) # Unexploited SPR

  MaxF <- 2 * max(Fk[Gear_i], F20[Gear_i], na.rm = T)
  Gear_Smx <- Sm[blicc_ld$spar[, 1]]

  # Grid Points
  mSmx <- seq(min(blicc_ld$LLB - 1), Linf,
              length.out = GridN)
  mFk <- seq(0, MaxF, length.out = GridN)
  dL <- (mSmx / Gear_Smx[Gear_i] - 1) / blicc_rp$vdir[Gear_i]

  vSm <- Sm
  indx <- with(blicc_ld, as.vector(t(spar[Fkg > 0, 1])))

  svdir <- blicc_rp$vdir[blicc_ld$Fkg > 0]

  F30 <- double(GridN)
  for (i in 1:GridN) {
    # Calculate the same proportional change to the selectivity parameters for all gears
    vSm[indx] <- (1 + svdir * dL[i]) * Sm[indx]
    F30_g <-
      FSPR_solve(Linf, Galpha, Mk, Fk, vSm, 0.3, blicc_rp$vdir, blicc_ld)
    F30[i] <- F30_g[blicc_ld$Fkg[Gear_i]]
  }
  SPRRP <- tibble::tibble(Smx_rp = mSmx, F30 = F30) |>
    dplyr::filter(F30 <= MaxF)

  SPRCurve <- tidyr::expand_grid(mFk = mFk, Smx = mSmx) |>
    dplyr::mutate(
      dF = (mFk / Fk[Gear_i] - 1) / blicc_rp$vdir[Gear_i],
      dL = (Smx / Gear_Smx[Gear_i] - 1) / blicc_rp$vdir[Gear_i],
      spr = purrr::map2_dbl(dF, dL, fSPR2)
    )

  All_Gear_Fk <- double(blicc_ld$NG)
  All_Gear_Fk[blicc_ld$Fkg] <- Fk  # With zero F's
  cp <-
    ggplot2::ggplot(SPRCurve, ggplot2::aes(x = Smx, y = mFk, z = spr)) +
    ggplot2::geom_point(
      tibble::tibble(
        Lcur = Gear_Smx[Gear_i],
        Fcur = All_Gear_Fk[Gear_i],
        Gear = blicc_ld$gname[Gear_i]
      ),
      mapping = ggplot2::aes(x = Lcur, y = Fcur),
      size = 3
    ) +
    ggplot2::labs(
      x = "Selectivity length parameter ",
      y = "Fishing Mortality (/k) ",
      colour =
        "SPR",
      title = paste0("SPR: ", blicc_ld$gname[Gear_i])
    ) +
    ggplot2::geom_contour(ggplot2::aes(colour = ggplot2::after_stat(level))) +
    ggplot2::geom_hline(yintercept = 0)

  if (nrow(SPRRP) > 0) {
    # Reference points exist / can be plotted
    cp <- cp +
      ggplot2::geom_line(
        data = SPRRP,
        ggplot2::aes(x = Smx_rp, y = F30),
        colour = "darkred",
        linewidth = 1.5,
        inherit.aes = F
      ) +
      ggplot2::annotate(
        "text",
        x = max(SPRRP$Smx_rp) + 0.5,
        y = max(SPRRP$F30) + 0.3,
        label = paste0("SPR30%"),
        colour = "darkred"
      )
  }
  return(cp)
}


#' Yield-per-recruit surface contour plot
#'
#' The plot shows the YPR surface plotted against the fishing mortality (Fk) and
#' the 50% or maximum selectivity, dependent on the selectivity function. The
#' current state of the stock is marked as a point and the F0.1 reference point
#' is marked as a line. All values are calculated as means and uncertainty is
#' not represented to maintain simplicity. The reference line is plotted if it
#' exists within a reasonable range of the current estimated status. For
#' multiple gears, the surface is calculated for a single reference gear.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_YPR_contour(eg_rp)
#'
plot_YPR_contour <- function(blicc_rp) {
  Smx_rp = ypr = F0.1 = Lcur = Fcur = `..level..` = NULL  # R CMD Notes
  Smx = dF = level = NULL

  fYPR2 <- function(dF, dL) {
    vFk <- (1 + svdir * dF) * Fk
    vSm[indx] <- (1 + svdir * dL) * Sm[indx]
    Rsel <- Rselectivities(vSm, blicc_ld)
    res <- Rpop_F(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
    Yield <- 0
    for (gi in 1:blicc_ld$NG) {
      if (blicc_ld$Fkg[gi] > 0)
        Yield <-
          Yield + sum(res$N_L * res$Fki[[gi]] * blicc_ld$wt_L)
    }
    return(Yield)
  }

  blicc_ld <- blicc_rp$ld
  GridN <- 40
  Linf <- mean(blicc_rp$rp_df$Linf)
  Galpha <- mean(blicc_rp$rp_df$Galpha)
  Gbeta <- Galpha / Linf
  Mk <- mean(blicc_rp$rp_df$Mk)
  Fk <-
    as.vector(tapply(
      X = unlist(dplyr::pull(blicc_rp$rp_df, Fk)),
      INDEX = rep(1:blicc_ld$NF, nrow(blicc_rp$rp_df)),
      FUN = mean
    ))
  Sm <-
    as.vector(tapply(
      X = unlist(dplyr::pull(blicc_rp$rp_df, Sm)),
      INDEX = rep(1:blicc_ld$NP, nrow(blicc_rp$rp_df)),
      FUN = mean
    ))
  F20 <-
    tapply(
      X = unlist(dplyr::pull(blicc_rp$rp_df, F20)),
      INDEX = rep(1:blicc_ld$NF, nrow(blicc_rp$rp_df)),
      FUN = mean,
      na.rm = TRUE
    )
  F30 <-
    tapply(
      X = unlist(dplyr::pull(blicc_rp$rp_df, F30)),
      INDEX = rep(1:blicc_ld$NF, nrow(blicc_rp$rp_df)),
      FUN = mean,
      na.rm = TRUE
    )

  # Decide reference gear on highest fishing mortality to be controlled
  Gear_i <- which.max(Fk * unname(blicc_rp$vdir[blicc_ld$Fkg > 0]))
  gear2plot <- blicc_ld$gname[Gear_i]

  MaxF <- 2 * max(Fk[Gear_i], F20[Gear_i], na.rm = T)
  Gear_Smx <- Sm[blicc_ld$spar[, 1]]

  # Grid Points
  mSmx <- seq(min(blicc_ld$LLB - 1), Linf,
              length.out = GridN)
  mFk <- seq(0, MaxF, length.out = GridN)
  dL <- (mSmx / Gear_Smx[Gear_i] - 1) / blicc_rp$vdir[Gear_i]

  vSm <- Sm
  indx <- with(blicc_ld, as.vector(t(spar[Fkg > 0, 1])))

  svdir <- blicc_rp$vdir[blicc_ld$Fkg > 0]

  F01 <- double(GridN)
  for (i in 1:GridN) {
    # Calculate the same proportional change to the selectivity parameters for all gears
    vSm[indx] <- (1 + svdir * dL[i]) * Sm[indx]
    F01_g <-
      F01_solve(Linf, Galpha, Mk, Fk, vSm, blicc_rp$vdir, blicc_ld)
    F01[i] <- F01_g[blicc_ld$Fkg[Gear_i]]
  }
  YPRRP <- tibble::tibble(Smx_rp = mSmx, F0.1 = F01) |>
    dplyr::filter(F0.1 <= MaxF)

  YPRCurve <- tidyr::expand_grid(mFk = mFk, Smx = mSmx) |>
    dplyr::mutate(
      dF = (mFk / Fk[Gear_i] - 1) / blicc_rp$vdir[Gear_i],
      dL = (Smx / Gear_Smx[Gear_i] - 1) / blicc_rp$vdir[Gear_i],
      ypr = purrr::map2_dbl(dF, dL, fYPR2)
    )

  All_Gear_Fk <- double(blicc_ld$NG)
  All_Gear_Fk[blicc_ld$Fkg] <- Fk  # With zero F's
  cp <-
    ggplot2::ggplot(YPRCurve, ggplot2::aes(x = Smx, y = mFk, z = ypr)) +
    ggplot2::geom_point(
      tibble::tibble(
        Lcur = Gear_Smx[Gear_i],
        Fcur = All_Gear_Fk[Gear_i],
        Gear = blicc_ld$gname[Gear_i]
      ),
      mapping = ggplot2::aes(x = Lcur, y = Fcur),
      size = 3
    ) +
    ggplot2::labs(
      x = "Selectivity location parameter ",
      y = "Fishing Mortality (/k) ",
      colour =
        "Relative Yield",
      title = paste0("Yield: ", blicc_ld$gname[Gear_i])
    ) +
    ggplot2::geom_contour(ggplot2::aes(colour = ggplot2::after_stat(level))) +
    ggplot2::geom_hline(yintercept = 0)

  if (nrow(YPRRP) > 0) {
    # Reference points exist / can be plotted
    cp <- cp +
      ggplot2::geom_line(
        data = YPRRP,
        ggplot2::aes(x = Smx_rp, y = F0.1),
        colour = "darkred",
        linewidth = 1.5,
        inherit.aes = F
      ) +
      ggplot2::annotate(
        "text",
        x = max(YPRRP$Smx_rp) + 0.5,
        y = max(YPRRP$F0.1) + 0.3,
        label = "F0.1",
        colour = "darkred"
      )
  }
  return(cp)
}
