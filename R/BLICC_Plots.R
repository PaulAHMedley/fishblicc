# BLICC Plotting Functions ------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


#' Plot the observed and expected frequency data based on the prior parameters
#'
#' The graph shows the observed frequency and the prior's predicted frequency
#' for all (or specified) gears. This allows a check that current prior
#' parameter settings are consistent with observations. The prior means are also
#' used as the start point for the [blicc_mpd] fits.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param gear Specifies the gear(s) to plot as an integer index or full name
#' @param time_period Specifies the period(s) to plot as an integer index or 
#'   full name
#' @return ggplot geom object plotting observed and the prior's expected
#'   frequency, separated by gear.
#' @examples
#' plot_prior(gillnet_ld)
#' 
plot_prior <- function(blicc_ld, gear = "All", time_period = "All") {
  Gear = LMP = fq = Name = NULL

  gear <- parse_gear(gear, blicc_ld)
  fq_used <- which(blicc_ld$Gi %in% gear)
  
  if (any(gear != "All" | time_period != "All")) { 
    gear <- parse_gear(gear, blicc_ld)
    time_period <- parse_period(time_period, blicc_ld)
    subset_fq <- with(blicc_ld, which(Gi %in% gear & Ti %in% time_period))
  } else {
    subset_fq <- seq(blicc_ld$NQ)
  }
  
  sel <- Rselectivities(exp(blicc_ld$polSm), blicc_ld)

  Zk <- with(blicc_ld, rep(list(exp(polMkm)*M_L), NT))
  for (qi in subset_fq) {
    if (blicc_ld$Fkq[qi] > 0)
      Zk[[blicc_ld$Ti[qi]]] <- with(blicc_ld, Zk[[Ti[qi]]] + 
                                      exp(polFkm[Fkq[qi]]) * sel[[Gi[qi]]])
  }
  pop <- list()
  ti_range <- unique(blicc_ld$Ti[subset_fq])
  
  for (ti in ti_range)
    pop[[ti]] <- with(blicc_ld, Rpop_len(gl_nodes, gl_weights, LLB, Zk[[ti]],
                                   exp(polGam), exp(polGam)/poLinfm))

  df <- tibble::tibble()
  for (qi in subset_fq) {
    efq <- with(blicc_ld, sel[[Gi[qi]]]*pop[[Ti[qi]]])
    efq <- with(blicc_ld, efq * sum(fq[[qi]]) / sum(efq))
    df <- with(blicc_ld, rbind(df, 
                               tibble::tibble(Name=fqname[qi], LMP = LMP,
                                fq=fq[[qi]], efq=efq)))
  }

  gp <- ggplot2::ggplot(data=df, ggplot2::aes(x = LMP)) +
    ggplot2::geom_line(ggplot2::aes(y=efq)) +
    ggplot2::geom_col(ggplot2::aes(y=fq), fill="darkgreen", alpha=0.4) +
    ggplot2::labs(x="Length", y="Frequency") +
    ggplot2::facet_wrap(ggplot2::vars(Name), scales="free_y")
  return(gp)
}

#' Plot the observed and expected length frequency data
#'
#' The graph shows the observed length frequency and the posterior predicted
#' frequency for the specified gears. The plots include the 80% credible
#' intervals for the estimate (if an MCMC is provided) and 80% credible
#' intervals for the observations. The latter takes into account the negative
#' binomial observation error, so on average 8 out of ten observations would be
#' expected within these bounds.
#'
#' @export
#' @inheritParams plot_prior
#' @inheritParams blicc_expected_catches
#' @return ggplot geom object plotting observed and expected frequency
#' @examples
#' plot_posterior(trgl_rp)
#' 
plot_posterior <-
  function(blicc_rp, gear = "All", time_period = "All") {
    .draw = Lgroup = LMP = fq = NB_phi = NULL  # Not necessary but stops CMD check notes
    efq = fq_lo = fq_hi = efq_m = efq_lo = efq_hi = dat_lo = dat_hi = NULL
    Qgroup = ofq = NULL

    dr_df <- blicc_rp$dr_df
    blicc_ld <- blicc_rp$ld
    blicc_lx <- blicc_rp$lx_df

    if (any(gear != "All" | time_period != "All")) { 
      gear <- parse_gear(gear, blicc_ld)
      time_period <- parse_period(time_period, blicc_ld)
      subset_fq <- with(blicc_ld, Gi %in% gear & Ti %in% time_period)
      fq2plot <- blicc_ld$fqname[subset_fq]
      ssfq <- blicc_ld$fq[subset_fq]
      
      blicc_lx <- blicc_lx |>
        dplyr::filter(Qgroup %in% fq2plot)
    } else {
      fq2plot <- blicc_ld$fqname
      ssfq <- blicc_ld$fq
    }
      
    Nplots <- length(fq2plot)

    dat_df <- with(blicc_ld,
                   tibble::tibble(
                     Qgroup = rep(fq2plot, each = NB),
                     Lgroup = factor(rep(LLB, Nplots)),
                     LMP = rep(LMP, Nplots),
                     ofq = unlist(ssfq)
                   ))

    if (nrow(dr_df) == 1) {
      # MPD estimate
      MCMC_draws <- FALSE
      df1 <- dr_df |>
        dplyr::select(.draw, NB_phi) |>
        dplyr::right_join(blicc_lx, by = ".draw") |>
        dplyr::select(.draw, Qgroup, Lgroup, efq_m = efq, NB_phi) |>
        dplyr::mutate(
          dat_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq_m),
          dat_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq_m)
        ) |>
        dplyr::left_join(dat_df, by = c("Qgroup", "Lgroup"))
    } else {
      MCMC_draws <- TRUE
      suppressWarnings(
        df1 <- dr_df |>
          dplyr::select(.draw, NB_phi) |>
          dplyr::left_join(blicc_lx, by = ".draw") |>
          dplyr::select(.draw, Qgroup, Lgroup, efq, NB_phi) |>
          dplyr::mutate(
            fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
            fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
          ) |>
          dplyr::group_by(Qgroup, Lgroup) |>
          dplyr::summarise(
            efq_m = mean(efq),
            efq_lo = stats::quantile(efq, probs = c(0.1), names = F),
            efq_hi = stats::quantile(efq, probs = c(0.9), names = F),
            dat_lo = stats::quantile(fq_lo, probs = c(0.1), names = F),
            dat_hi = stats::quantile(fq_hi, probs = c(0.9), names = F),
          ) |>
          dplyr::ungroup() |>
          dplyr::left_join(dat_df, by = c("Qgroup", "Lgroup"))
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

    if (length(fq2plot) > 1) {
      gp <- gp + ggplot2::facet_wrap(ggplot2::vars(Qgroup),
                            ncol = 2,
                            scales = "free_y")
    } else {
      gp <- gp + ggplot2::ggtitle(fq2plot)
    }
    return(gp)
  }

#' Plot standardised residuals against length
#'
#' The graph shows the standardised residuals plotted against length for all
#' gears combined, separated by colour. Due to the variance weighting, outliers
#' can most be expected at the tails of the frequency. For the MCMC output, the
#' residuals are plotted as 99% quartiles using error bars, together with the
#' outliers as points outside this range.
#'
#' If you want to alter the graph, the group names are `Sgroup` for the
#' selectivities and `Lgroup` for the lengths.
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object plotting standardised residuals
#' @examples
#' plot_residuals(trgl_rp)
#' 
plot_residuals <- function(blicc_rp, gear = "All", time_period = "All") {
  .draw = NB_phi = Qgroup = Lgroup = efq = ofq = LMP = std_res = NULL
  std_res_m = std_res_005 = std_res_995 = NULL

  dr_df <- blicc_rp$dr_df
  blicc_ld <- blicc_rp$ld
  blicc_lx <- blicc_rp$lx_df
  
  if (any(gear != "All" | time_period != "All")) { 
    gear <- parse_gear(gear, blicc_ld)
    time_period <- parse_period(time_period, blicc_ld)
    subset_fq <- with(blicc_ld, Gi %in% gear & Ti %in% time_period)
    fq2plot <- blicc_ld$fqname[subset_fq]
    ssfq <- blicc_ld$fq[subset_fq]
    
    blicc_lx <- blicc_lx |>
      dplyr::filter(Qgroup %in% fq2plot)
  } else {
    fq2plot <- blicc_ld$fqname
    ssfq <- blicc_ld$fq
  }
  
  Nplots <- length(fq2plot)
  
  dat_df <- with(blicc_ld,
                 tibble::tibble(
                   Qgroup = rep(fq2plot, each = NB),
                   Lgroup = factor(rep(LLB, Nplots)),
                   LMP = rep(LMP, Nplots),
                   ofq = unlist(ssfq)
                 ))
  suppressWarnings(
    df1 <- dr_df |>
      dplyr::select(.draw, NB_phi) |>
      dplyr::right_join(blicc_lx, by = ".draw") |>
      dplyr::select(.draw, Qgroup, Lgroup, efq, NB_phi) |>
      dplyr::left_join(dat_df, by = c("Qgroup", "Lgroup")) |>
      dplyr::mutate(std_res = (ofq - efq) / sqrt(efq + (efq ^ 2) / NB_phi))
  )

  if (nrow(dr_df) == 1) {
    gp <-
      ggplot2::ggplot(df1, ggplot2::aes(x = LMP, y = std_res, colour = Qgroup)) +
      ggplot2::geom_point() +
      ggplot2::geom_line(alpha=0.4)
  } else {
    df2 <- df1 |>
      dplyr::select(Qgroup, Lgroup, std_res) |>
      dplyr::group_by(Qgroup, Lgroup) |>
      dplyr::summarise(
        std_res_m = stats::median(std_res),
        std_res_005 = stats::quantile(std_res, probs = c(0.005), names = F),
        std_res_995 = stats::quantile(std_res, probs = c(0.995), names = F),
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])

    outlier_df <- df1 |>
      dplyr::select(Qgroup, Lgroup, std_res) |>
      dplyr::left_join(df2, by = c("Qgroup", "Lgroup")) |>
      dplyr::filter(std_res < std_res_005 | std_res > std_res_995)

    gp <-
      ggplot2::ggplot(data = df2) +
      ggplot2::geom_errorbar(ggplot2::aes(x = LMP,
                                          ymin = std_res_005,
                                          ymax = std_res_995,
                                          colour = Qgroup), width=0.2) +
      ggplot2::geom_point(ggplot2::aes(x=LMP, y=std_res_m,
                                      colour = Qgroup), shape="+") +
      ggplot2::geom_point(data=outlier_df, ggplot2::aes(x=LMP, y=std_res,
                                                   colour = Qgroup), alpha=0.5)
  }
  gp <- gp +
    ggplot2::labs(x = "Length", y = "Standardised Residuals", colour = "Sample") +
    ggplot2::geom_hline(yintercept = 0)

  return(gp)
}


#' Plot selectivity and 80% credible interval, if available, by length.
#'
#' The selectivity at length is plotted for the specified gears. All gears, or
#' specified gears, are plotted on a single graph. These curves have to be
#' multiplied by the fishing mortality (Fk) parameters to obtain the relative
#' fishing mortality at length. The 80%CI is only available for the MCMC fit.
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object for plotting
#' @examples
#' plot_selectivity(trgl_rp)
#' 
plot_selectivity <- function(blicc_rp, gear = "All") {
  Lgroup = sel = LMP = sel_10 = sel_90 = sel_m = Qgroup = NULL

  blicc_ld <- blicc_rp$ld
  blicc_lx <- blicc_rp$lx_df

  gear_i <- parse_gear(gear, blicc_ld)
  gear2plot <- blicc_ld$gname[gear_i]
  fq2plot <- character(0)
  for (i in gear_i) {
    fq2plot <- with(blicc_ld, c(fq2plot, fqname[min(which(i == Gi))]))
  }
  blicc_lx <- blicc_lx |>
    dplyr::filter(Qgroup %in% fq2plot)|>
    dplyr::mutate(Qgroup = gear2plot[match(Qgroup, fq2plot)]) # filter only one selectivity for each gear
  
  gp <- blicc_lx |>
    dplyr::select(Qgroup, Lgroup, sel) |>
    dplyr::group_by(Qgroup, Lgroup) |>
    dplyr::summarise(
      sel_m = stats::median(sel),
      sel_10 = stats::quantile(sel, probs = c(0.10), names = F),
      sel_90 = stats::quantile(sel, probs = c(0.90), names = F),
      .groups = "drop"
    ) |>
    #dplyr::ungroup() |>
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)]) |>
    ggplot2::ggplot(ggplot2::aes(x = LMP, fill = Qgroup, colour = Qgroup)) +
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


#' Plot the estimated spawning potential ratio probability density, if
#' available.
#'
#' The SPR estimate is converted to a probability density plot with 20% and 40%
#' reference point. This plot requires the MCMC to have been run. For multiple 
#' time periods, densities are overlaid.
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object for plotting
#' @examples
#' plot_SPR_density(trgl_rp)
#' 
plot_SPR_density <- function(blicc_rp) {
  SPR = .draw = Period = NULL
  if (nrow(blicc_rp$dr_df)/blicc_rp$ld$NT <= 500) {
    stop("To obtain a density, you will need to obtain sufficient values (>500) from MCMC.")
  }
  suppressWarnings({
    dr_df <- blicc_rp$dr_df |>
      dplyr::select(.draw, SPR) |>
      dplyr::mutate(Period = list(blicc_rp$ld$tpname)) |>
      tidyr::unnest(c(SPR, Period))
  })
  
  if (blicc_rp$ld$NT==1) {
    gp <- ggplot2::ggplot(dr_df, ggplot2::aes(SPR)) +
      ggplot2::geom_density(fill = "lightblue") 
  } else {
    gp <- ggplot2::ggplot(dr_df, ggplot2::aes(SPR, fill=Period)) +
      ggplot2::geom_density(alpha=0.5)
  }
  gp <- gp +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0.2),
                        colour = "red",
                        linetype = 2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0.4),
                        colour = "green",
                        linetype = 2) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Spawning Potential Ratio", y = "Probability Density")
  return(gp)
}



#' Plot fishing mortality relative to the SPR 40% reference point (`Fk/F40`)
#'
#' This plot requires the MCMC to have been run. The plot is for the time period 
#' for which the reference points were calculated (see [blicc_ref_pts]).
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object for plotting
#' @examples
#' plot_FkF40_density(trgl_rp)
#'
plot_FkF40_density <- function(blicc_rp) {
  Fk = F40 = `Fk/F40` = id = NULL 

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

  fq_good <- sum(rp_df$`Fk/F40` <= 1.0, na.rm=TRUE)

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


#' Plot expected length frequencies for a range of fishing mortality and
#' selectivity reference points
#'
#' The facet plot estimates the length frequency data expected for a range of
#' fishing mortality and selectivity to show how the length frequency data might
#' respond to changes in mortality and/or selectivity. Fishing mortalities
#' tested are SPR 20% and SPR 40% and selectivity to obtain SPR 40% and the
#' maximum yield.
#'
#' @details The graphs show the expected values and 80% credible intervals
#'   compared to the observations. They can be used to assess whether length
#'   frequencies should be able to detect changes in fishing mortality to these
#'   different levels. In many cases, differences between frequencies can be
#'   subtle and the precision of data collection may not be sufficient to detect
#'   such changes reliably.  Dependent on the selectivity, reference points may
#'   not exist. It may not be possible to achieve an SPR level without changing
#'   multiple parameters, whereas these scenarios are restricted to only
#'   changing either the fishing mortality or the selectivity function central
#'   locations. If reference points do not exist, they are not plotted, so some
#'   graphs may be blank except for the data. The graphs show the number of MCMC
#'   draws for which the reference points exist in each case.
#'
#'   If the apparent changes in length frequencies between the plots is small,
#'   relying on length frequency data alone may not be advisable.
#'
#'   If selectivity mixtures are used, all location parameters for the mixtures
#'   will be shifted equally by default. This may not make a lot of sense in
#'   practice and more realistic scenarios may need to be constructed to test
#'   possible management actions.
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object for plotting expected length frequencies by SPR
#'   level
#' @examples
#' plot_efq(trgl_rp)
#' 
plot_efq <- function(blicc_rp, gear = NULL) {
  Linf = Galpha = Mk = Fk = Sm = NB_phi = Sgroup = NULL
  .draw = Lgroup = Current = SPR20 = SPR40 = SMY = NULL
  F_SPR40 = S40 = Sel_MY = LMP = fq = Harvest_Level = NULL # Not necessary but stops CMD check notes
  F_SPR20 = F20 = F30 = F40 = efq = fq_lo = fq_hi = N = NULL
  efq_m = dat_lo = dat_hi = lab = x = y = NULL

  blicc_ld <- blicc_rp$scenario$time_period_ld
  svdir <-blicc_rp$scenario$vdir
  #  blicc_lx <- blicc_rp$lx_df

  if (is.null(gear)) {
    elig_gears <- seq(blicc_ld$NG)[svdir>0]
    gear_i <- elig_gears[which.max(blicc_ld$prop_catch[svdir>0])]
  } else {
    gear_i = parse_gear(gear, blicc_rp$ld)[1]
    if (! (gear_i %in% blicc_rp$scenario$gears))
      stop(paste("Error: The specified gear needs to be one of the following in the relevant time period:",
                 paste(blicc_rp$ld$gname[blicc_rp$scenario$gears], collapse=" ")))
    gear_i <- match(gear_i, blicc_rp$scenario$gears)  # translate to scenario
    if (svdir[gear_i] <= 0) 
      stop("Error: The specified gear management change (vdir) must be greater than zero.")
  }
  
  gear2plot <- blicc_ld$gname[gear_i]
  # blicc_lx <- blicc_lx |>
  #   dplyr::filter(Qgroup == gear2plot)

  df_dat <-
    tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq[[gear_i]])

  suppressWarnings(
    df <- blicc_rp$rp_df |>
      dplyr::mutate(Lgroup = list(factor(blicc_ld$LLB))) |>
      dplyr::mutate(
        F_SPR20 = purrr::pmap(
          list(Linf, Galpha, Mk, F20, Sm),
          blicc_get_efq,
          gear_i = gear_i,
          blicc_ld = blicc_ld
        ),
        # SPR30 = purrr::pmap(
        #   list(Linf, Galpha, Mk, F30, Sm),
        #   blicc_get_efq,
        #   gear_i = gear_i,
        #   blicc_ld = blicc_ld
        # ),
        F_SPR40 = purrr::pmap(
          list(Linf, Galpha, Mk, F40, Sm),
          blicc_get_efq,
          gear_i = gear_i,
          blicc_ld = blicc_ld
        ),
        Sel_SPR40 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, S40),
          blicc_get_efq,
          gear_i = gear_i,
          blicc_ld = blicc_ld
        ),
        Sel_MY = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, SMY),
          blicc_get_efq,
          gear_i = gear_i,
          blicc_ld = blicc_ld
        )
      ) |>
      dplyr::select(`.draw`, NB_phi, Lgroup, F_SPR20:Sel_MY) |>
      tidyr::unnest(Lgroup:Sel_MY) #|>
      # dplyr::left_join(
      #   dplyr::select(blicc_lx, .draw, Lgroup, Current = efq),
      #   by = c(".draw", "Lgroup")
      # )
  )

  df1 <- df |>
    tidyr::pivot_longer(cols = F_SPR20:Sel_MY,
                        names_to = "Harvest_Level",
                        values_to = "efq") |>
    dplyr::filter(!is.na(efq)) |>
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
#        na.rm = T,
        names = F
      ),
      dat_hi = stats::quantile(
        fq_hi,
        probs = c(0.9),
#        na.rm = T,
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


#' Spawning potential ratio surface contour plot
#'
#' The plot shows the SPR surface plotted against the fishing mortality (`Fk`)
#' and the 50% or maximum (location) selectivity, dependent on the selectivity
#' functions of the reference gear. The current state of the stock is marked as
#' a point and the SPR 40% reference point is marked as a line. All values are
#' calculated as means and uncertainty is not represented to maintain
#' simplicity. The reference line is plotted if it exists within a reasonable
#' range of the current estimated status. For multiple gears, the surface is
#' calculated for a single reference gear. For multiple time series, the plot is
#' for the time period used to estimate reference points (see [blicc_ref_pts]).
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object for plotting
#' @examples
#' plot_SPR_contour(trgl_rp)
#' 
plot_SPR_contour <- function(blicc_rp, gear = NULL) {
  Smx_rp = spr = Lcur = Fcur = `..level..` = NULL  #R CMD Notes
  Smx = dF = level = NULL

  fSPR2 <- function(dF, dL) {
    # Returns SPR calculation for delta- F and delta- selectivity location
    # parameter
    vFk <- (1 + svdir * dF) * Fk
    vSm[indx] <- (1 + svdir * dL) * Sm[indx]
    Rsel <- Rselectivities(vSm, blicc_ld)
    res <- Rpop_F(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
    return(sum(res$N_L[[1]] * blicc_ld$ma_L / SPR0))
  }

  blicc_ld <- blicc_rp$scenario$time_period_ld
  rp_df <- blicc_rp$rp_df
  svdir <-blicc_rp$scenario$vdir
  
  # If no reference gear defined, use gear with highest fishing
  if (is.null(gear)) {
    elig_gears <- seq(blicc_ld$NG)[svdir>0]
    gear_i <- elig_gears[which.max(blicc_ld$prop_catch[svdir>0])]
  } else {
    gear_i = parse_gear(gear, blicc_rp$ld)[1]
    if (! (gear_i %in% blicc_rp$scenario$gears))
      stop(paste("Error: The specified gear needs to be one of the following in the relevant time period:",
                 paste(blicc_rp$ld$gname[blicc_rp$scenario$gears], collapse=" ")))
    gear_i <- match(gear_i, blicc_rp$scenario$gears)  # translate to scenario
    if (svdir[gear_i] <= 0) 
      stop("Error: The specified gear management change (vdir) must be greater than zero.")
  }
  
  gear2plot <- blicc_ld$gname[gear_i]
  
  GridN <- 40
  Linf <- mean(rp_df$Linf)
  Galpha <- mean(rp_df$Galpha)
  Gbeta <- Galpha / Linf
  Mk <- mean(rp_df$Mk)
  Fk <- as.vector(tapply(
    X = unlist(dplyr::pull(rp_df, Fk)),
    INDEX = rep(1:blicc_ld$NF, nrow(rp_df)),
    FUN = mean
  ))
  Sm <- as.vector(tapply(
    X = unlist(dplyr::pull(rp_df, Sm)),
    INDEX = rep(1:(blicc_ld$NP+blicc_ld$NM), nrow(rp_df)),
    FUN = mean
  ))
  # Include checks for NA
  vF20 <- dplyr::pull(rp_df, F20)
  vF20 <- vF20[!is.na(vF20)]
  if (length(vF20) == 0)
    F20 <- NA
  else
    F20 <- tapply(
      X = unlist(vF20),
      INDEX = rep(1:blicc_ld$NF, length(vF20)),
      FUN = mean,
      na.rm = TRUE
    )
  vF40 <- dplyr::pull(rp_df, F40)
  vF40 <- vF40[!is.na(vF40)]
  if (length(vF40) == 0)
    F40 <- NA
  else
    F40 <- tapply(
      X = unlist(vF40),
      INDEX = rep(1:blicc_ld$NF, length(vF40)),
      FUN = mean,
      na.rm = TRUE
    )
  
  SPR0 <-
    RSPR_0(Galpha, Gbeta, Mk, blicc_ld) # Unexploited SPR

  MaxF <- 2 * max(Fk[gear_i], F20[gear_i], na.rm = TRUE)
  Gear_Smx <- Sm[blicc_ld$sp_i]

  # Grid Points
  mSmx <- seq(min(blicc_ld$LLB - 1), Linf,
              length.out = GridN)
  mFk <- seq(0, MaxF, length.out = GridN)
  dL <- (mSmx / Gear_Smx[gear_i] - 1) / svdir[gear_i]

  vSm <- Sm
  indx <- blicc_ld$sp_i

  F40 <- double(GridN)
  for (i in 1:GridN) {
    # Calculate the same proportional change to the selectivity parameters for all gears
    vSm[indx] <- (1 + svdir * dL[i]) * Sm[indx]
    F40_g <-
      FSPR_solve(Linf, Galpha, Mk, Fk, vSm, 0.3, svdir, blicc_ld)
    F40[i] <- F40_g[blicc_ld$Fkq[gear_i]]
  }
  SPRRP <- tibble::tibble(Smx_rp = mSmx, F40 = F40) |>
    dplyr::filter(F40 <= MaxF)

  SPRCurve <- tidyr::expand_grid(mFk = mFk, Smx = mSmx) |>
    dplyr::mutate(
      dF = (mFk / Fk[gear_i] - 1) / svdir[gear_i],
      dL = (Smx / Gear_Smx[gear_i] - 1) / svdir[gear_i],
      spr = purrr::map2_dbl(dF, dL, fSPR2)
    )

  All_Gear_Fk <- double(blicc_ld$NG)
  All_Gear_Fk[blicc_ld$Fkq] <- Fk  # With zero F's
  cp <-
    ggplot2::ggplot(SPRCurve, ggplot2::aes(x = Smx, y = mFk, z = spr)) +
    ggplot2::geom_point(
      tibble::tibble(
        Lcur = Gear_Smx[gear_i],
        Fcur = All_Gear_Fk[gear_i],
        Gear = gear2plot
      ),
      mapping = ggplot2::aes(x = Lcur, y = Fcur),
      size = 3
    ) +
    ggplot2::labs(
      x = "Selectivity length parameter ",
      y = "Fishing Mortality (/k) ",
      colour =
        "SPR",
      title = paste0("SPR: ", gear2plot)
    ) +
    ggplot2::geom_contour(ggplot2::aes(colour = ggplot2::after_stat(level))) +
    ggplot2::geom_hline(yintercept = 0)

  if (nrow(SPRRP) > 0) {
    # Reference points exist / can be plotted
    cp <- cp +
      ggplot2::geom_line(
        data = SPRRP,
        ggplot2::aes(x = Smx_rp, y = F40),
        colour = "darkred",
        linewidth = 1.5,
        inherit.aes = F
      ) +
      ggplot2::annotate(
        "text",
        x = max(SPRRP$Smx_rp) + 0.5,
        y = max(SPRRP$F40) + 0.3,
        label = paste0("SPR40%"),
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
#' multiple gears, the surface is calculated for a single reference gear. For
#' multiple time series, the plot is for the time period used to estimate
#' reference points (see [blicc_ref_pts])
#'
#' @export
#' @inheritParams plot_posterior
#' @return ggplot geom object for plotting
#' @examples
#' plot_YPR_contour(trgl_rp)
#' 
plot_YPR_contour <- function(blicc_rp, gear = NULL) {
  Smx_rp = ypr = F0.1 = Lcur = Fcur = `..level..` = NULL  # R CMD Notes
  Smx = dF = level = NULL

  fYPR2 <- function(dF, dL) {
    # Returns YPR calculation for delta F and delta selectivity location
    # parameter
    vFk <- (1 + svdir * dF) * Fk
    vSm[indx] <- (1 + svdir * dL) * Sm[indx]
    Rsel <- Rselectivities(vSm, blicc_ld)
    res <- Rpop_F(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
    Yield <- 0
    for (gi in 1:blicc_ld$NG) { # NG == NQ !!
      if (blicc_ld$Fkq[gi] > 0)
        Yield <-
          Yield + sum(res$N_L[[1]] * res$Fki[[gi]] * blicc_ld$wt_L)
    }
    return(Yield)
  }

  blicc_ld <- blicc_rp$scenario$time_period_ld
  rp_df <- blicc_rp$rp_df
  svdir <-blicc_rp$scenario$vdir
  
  # If no reference gear defined, use gear with highest fishing
  if (is.null(gear)) {
    elig_gears <- seq(blicc_ld$NG)[svdir>0]
    gear_i <- elig_gears[which.max(blicc_ld$prop_catch[svdir>0])]
  } else {
    gear_i = parse_gear(gear, blicc_rp$ld)[1]
    if (! (gear_i %in% blicc_rp$scenario$gears))
      stop(paste("Error: The specified gear needs to be one of the following in the relevant time period:",
           paste(blicc_rp$ld$gname[blicc_rp$scenario$gears], collapse=" ")))
    gear_i <- match(gear_i, blicc_rp$scenario$gears)  # translate to scenario
    if (svdir[gear_i] <= 0) 
      stop("Error: The specified gear management change (vdir) must be greater than zero.")
  }
  
  gear2plot <- blicc_ld$gname[gear_i]

  GridN <- 40
  Linf <- mean(rp_df$Linf)
  Galpha <- mean(rp_df$Galpha)
  Gbeta <- Galpha / Linf
  Mk <- mean(rp_df$Mk)
  Fk <-
    as.vector(tapply(
      X = unlist(dplyr::pull(rp_df, Fk)),
      INDEX = rep(1:blicc_ld$NF, nrow(rp_df)),
      FUN = mean
    ))
  Sm <-
    as.vector(tapply(
      X = unlist(dplyr::pull(rp_df, Sm)),
      INDEX = rep(1:(blicc_ld$NP+blicc_ld$NM), nrow(rp_df)),
      FUN = mean
    ))
  vF20 <- dplyr::pull(rp_df, F20)
  vF20 <- vF20[!is.na(vF20)]
  if (length(vF20) == 0)
    F20 <- NA
  else
    F20 <-
    tapply(
      X = unlist(vF20),
      INDEX = rep(1:blicc_ld$NF, length(vF20)),
      FUN = mean,
      na.rm = TRUE
    )
  vF40 <- dplyr::pull(rp_df, F40)
  vF40 <- vF40[!is.na(vF40)]
  if (length(vF40) == 0)
    F40 <- NA
  else
    F40 <-
    tapply(
      X = unlist(vF40),
      INDEX = rep(1:blicc_ld$NF, length(vF40)),
      FUN = mean,
      na.rm = TRUE
    )
  
  MaxF <- 2 * max(Fk[gear_i], F20[gear_i], na.rm = T)
  indx <- blicc_ld$sp_i
  Gear_Smx <- Sm[indx]

  # Grid Points
  mSmx <- seq(min(blicc_ld$LLB - 1), Linf,
              length.out = GridN)
  mFk <- seq(0, MaxF, length.out = GridN)
  svdir <- blicc_rp$scenario$vdir
  dL <- (mSmx / Gear_Smx[gear_i] - 1) / svdir[gear_i]

  vSm <- Sm

  F01 <- double(GridN)
  for (i in 1:GridN) {
    # Calculate the same proportional change to the selectivity parameters for all gears
    vSm[indx] <- (1 + svdir * dL[i]) * Sm[indx]
    F01_g <-
      F01_solve(Linf, Galpha, Mk, Fk, vSm, svdir, blicc_ld)
    F01[i] <- F01_g[blicc_ld$Fkq[gear_i]]
  }
  YPRRP <- tibble::tibble(Smx_rp = mSmx, `F0.1` = F01) |>
    dplyr::filter(F0.1 <= MaxF)

  YPRCurve <- tidyr::expand_grid(mFk = mFk, Smx = mSmx) |>
    dplyr::mutate(
      dF = (mFk / Fk[gear_i] - 1) / svdir[gear_i],
      dL = (Smx / Gear_Smx[gear_i] - 1) / svdir[gear_i],
      ypr = purrr::map2_dbl(dF, dL, fYPR2)
    )

  All_Gear_Fk <- double(blicc_ld$NG)
  All_Gear_Fk[blicc_ld$Fkq] <- Fk  # With zero F's
  cp <-
    ggplot2::ggplot(YPRCurve, ggplot2::aes(x = Smx, y = mFk, z = ypr)) +
    ggplot2::geom_point(
      tibble::tibble(
        Lcur = Gear_Smx[gear_i],
        Fcur = All_Gear_Fk[gear_i],
        Gear = blicc_ld$gname[gear_i]
      ),
      mapping = ggplot2::aes(x = Lcur, y = Fcur),
      size = 3
    ) +
    ggplot2::labs(
      x = "Selectivity location parameter ",
      y = "Fishing Mortality (/k) ",
      colour =
        "Relative Yield",
      title = paste0("Yield: ", blicc_ld$gname[gear_i])
    ) +
    ggplot2::geom_contour(ggplot2::aes(colour = ggplot2::after_stat(level))) +
    ggplot2::geom_hline(yintercept = 0)

  if (nrow(YPRRP) > 0) {
    # Reference points exist / can be plotted
    cp <- cp +
      ggplot2::geom_line(
        data = YPRRP,
        ggplot2::aes(x = Smx_rp, y = `F0.1`),
        colour = "darkred",
        linewidth = 1.5,
        inherit.aes = F
      ) +
      ggplot2::annotate(
        "text",
        x = max(YPRRP$Smx_rp) + 0.5,
        y = max(YPRRP$`F0.1`) + 0.3,
        label = "F0.1",
        colour = "darkred"
      )
  }
  return(cp)
}
