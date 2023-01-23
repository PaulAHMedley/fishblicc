# BLICC Plotting Functions ------------------------------------------------


#' Plot the observed and expected frequency data
#'
#' The graph shows the observed frequency and the expected frequency together
#' with the 80% credible intervals for the estimate as well as 80% credible
#' intervals for the observations. The latter takes into account the negative
#' binomial observation error, so on average 8 out of ten observations would be
#' expected within these bounds.
#'
#' @export
#' @param blicc_rp  A posterior draws and reference points tibble
#' from `blicc_ref_pts` function.
#' @param blicc_lx  A tibble from `blicc_expect_len` function.
#' @param blicc_ld  A data list from the `blicc_dat` function.
#' @return ggplot geom object plotting observed and expected frequency
#' @examples
#' eg_lx <- blicc_expect_len(eg_rp, eg_ld)
#' plot_expected_frequency(eg_rp, eg_lx, eg_ld)
#'
plot_expected_frequency <- function(blicc_rp, blicc_lx, blicc_ld) {
  .draw=Lgroup=LMP=fq=NB_phi=NULL  # Not necessary but stops CMD check notes
  efq=fq_lo=fq_hi=efq_m=efq_lo=efq_hi=dat_lo=dat_hi=NULL

  if (nrow(blicc_rp)==1) {
    df1 <- blicc_rp |>
      dplyr::select(.draw, NB_phi) |>
      dplyr::left_join(blicc_lx, by = ".draw") |>
      dplyr::select(.draw, Lgroup, efq, NB_phi) |>
      dplyr::mutate(
        dat_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
        dat_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
      ) |>
      dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])

    df2 <- tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq)

    ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
      ggplot2::geom_col(
        data = df2,
        ggplot2::aes(x = LMP, y = fq),
        fill = "lightblue",
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = df2,
        ggplot2::aes(x = LMP, y = fq),
        shape = 95,
        colour = "black"
      ) +
      ggplot2::geom_line(ggplot2::aes(y = efq)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
                           fill = "black",
                           alpha = 0.2) +
      ggplot2::labs(x = "Length", y = "Frequency") +
      ggplot2::geom_hline(yintercept = 0)

  } else {
    df1 <- blicc_rp |>
      dplyr::select(.draw, NB_phi) |>
      dplyr::left_join(blicc_lx, by = ".draw") |>
      dplyr::select(.draw, Lgroup, efq, NB_phi) |>
      dplyr::mutate(
        fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
        fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
      ) |>
      dplyr::group_by(Lgroup) |>
      dplyr::summarise(
        efq_m = mean(efq),
        efq_lo = stats::quantile(efq, probs = c(0.1), names = F),
        efq_hi = stats::quantile(efq, probs = c(0.9), names = F),
        dat_lo = stats::quantile(fq_lo, probs = c(0.1), names = F),
        dat_hi = stats::quantile(fq_hi, probs = c(0.9), names = F),
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])

    df2 <- tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq)

    ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
      ggplot2::geom_col(
        data = df2,
        ggplot2::aes(x = LMP, y = fq),
        fill = "lightblue",
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = df2,
        ggplot2::aes(x = LMP, y = fq),
        shape = 95,
        colour = "black"
      ) +
      ggplot2::geom_line(ggplot2::aes(y = efq_m)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = efq_lo, ymax = efq_hi),
                           fill = "blue",
                           alpha = 0.3) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
                           fill = "black",
                           alpha = 0.2) +
      ggplot2::labs(x = "Length", y = "Frequency") +
      ggplot2::geom_hline(yintercept = 0)
  }

  }


#' Plot selectivity and 50% and 80% credible intervals by length
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' eg_lx <- blicc_expect_len(eg_rp, eg_ld)
#' plot_selectivity(eg_lx, eg_ld)
#'
plot_selectivity <- function(blicc_lx, blicc_ld) {
  Lgroup=sel=LMP=sel_01=sel_90=sel_25=sel_75=sel_m=NULL
  blicc_lx |>
    dplyr::select(Lgroup, sel) |>
    dplyr::group_by(Lgroup) |>
    dplyr::summarise(
      sel_m = mean(sel),
      sel_01 = stats::quantile(sel, probs = c(0.10), names = F),
      sel_25 = stats::quantile(sel, probs = c(0.25), names = F),
      sel_75 = stats::quantile(sel, probs = c(0.75), names = F),
      sel_90 = stats::quantile(sel, probs = c(0.90), names = F),
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)]) |>
    ggplot2::ggplot(ggplot2::aes(x = LMP)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = sel_01, ymax = sel_90),
      fill = "#546E8A",
      alpha = 0.7
    ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = sel_25, ymax = sel_75),
      fill = "#B0BED5",
      alpha = 0.7
    ) +
    ggplot2::geom_line(ggplot2::aes(y = sel_m)) +
    ggplot2::labs(x = "Length", y = "Selectivity")
}


#' Plot of the estimated spawning potential ratio probability density
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_SPR_density(eg_rp)
#'
plot_SPR_density <- function(blicc_rp) {
  SPR=NULL
  if (nrow(blicc_rp) <= 100) {
    return("To obtain a density, you will need to obtain sufficient values (>100) from MCMC.")
  }
  ggplot2::ggplot(blicc_rp, ggplot2::aes(SPR)) +
    ggplot2::geom_density(fill = "lightblue") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0.2),
      colour = "red",
      linetype = 2
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0.4),
      colour = "green",
      linetype = 2
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Spawning Potential Ratio", y = "Probability Density")
}



#' Plot fishing mortality relative to the SPR 40% reference point (Fk/F40)
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_FkF40_density(eg_rp)
#'
plot_FkF40_density <- function(blicc_rp) {
  Fk=F40=NULL # Not necessary but stops CMD check notes
  if (nrow(blicc_rp) <= 100) {
    return("To obtain a density, you will need to obtain sufficient values (>100) from MCMC.")
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
        n = n
      )
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
  ggplot2::ggplot(blicc_rp, ggplot2::aes(Fk / F40)) +
    ggplot2::geom_density(fill = "#A11826") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 1.0),
      colour = "green",
      linetype = 1
    ) +
    stat_aud(
      geom = "area",
      ggplot2::aes(fill = "good"),
      xlim = c(0, 1),
      alpha = 1.0,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = c("good" = "#204030")) +
    ggplot2::labs(y = "Probability Density")
}


#' Plot expected length frequencies for a range of fishing mortality
#' reference points
#'
#' The facet plot covers the current estimated fishing mortality, together with
#' fishing mortalities required to obtain SPR 20%, SPR 30% and SPR 40%. The
#' graphs show the expected values and 80% credible intervals compared to
#' current observations, and can be used to assess whether length frequencies
#' should be able to detect changes in fishing mortality to these different
#' levels. Note that if reference points do not exist, they are not
#' plotted, so some graphs may be blank except for data.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting expected length frequencies by SPR
#' @examples
#' eg_lx <- blicc_expect_len(eg_rp, eg_ld)
#' plot_efq_FRP(eg_rp, eg_lx, eg_ld)
#'
plot_efq_FRP <- function(blicc_rp, blicc_lx, blicc_ld) {
  # Not necessary but stops CMD check notes
  Linf=Galpha=Mk=Fk=Smx=Ss1=Ss2=NB_phi=NULL
  .draw=Lgroup=Current=SPR20=SPR40=SMY=LMP=fq=Harvest_Level=NULL
  F20=F30=F40=efq=fq_lo=fq_hi=N=efq_m=dat_lo=dat_hi=lab=x=y=NULL

  glq <- statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  df <- blicc_rp |>
    dplyr::mutate(Lgroup = list(factor(blicc_ld$Len))) |>
    dplyr::mutate(
      SPR20 = purrr::pmap(
        list(Linf, Galpha, Mk, F20, Smx, Ss1, Ss2),
        blicc_get_efq,
        blicc_ld = blicc_ld,
        glq
      ),
      SPR30 = purrr::pmap(
        list(Linf, Galpha, Mk, F30, Smx, Ss1, Ss2),
        blicc_get_efq,
        blicc_ld = blicc_ld,
        glq
      ),
      SPR40 = purrr::pmap(
        list(Linf, Galpha, Mk, F40, Smx, Ss1, Ss2),
        blicc_get_efq,
        blicc_ld = blicc_ld,
        glq
      )
    ) |>
    dplyr::select(`.draw`, NB_phi, Lgroup, SPR20:SPR40) |>
    tidyr::unnest(Lgroup:SPR40) |>
    dplyr::left_join(dplyr::select(blicc_lx, .draw, Lgroup, Current = efq),
      by = c(".draw", "Lgroup")
    )

  df1 <- df |>
    tidyr::pivot_longer(
      cols = SPR20:Current,
      names_to = "Harvest_Level",
      values_to = "efq"
    ) |>
    dplyr::mutate(
      fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
      fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
    ) |>
    dplyr::group_by(Harvest_Level, Lgroup) |>
    dplyr::summarise(
      efq_m = mean(efq),
      dat_lo = stats::quantile(fq_lo, probs = c(0.1), na.rm = T, names = F),
      dat_hi = stats::quantile(fq_hi, probs = c(0.9), na.rm = T, names = F),
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

  df2 <- tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq)

  ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
    ggplot2::geom_col(
      data = df2,
      ggplot2::aes(x = LMP, y = fq),
      fill = "lightblue",
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = df2,
      ggplot2::aes(x = LMP, y = fq),
      shape = 95,
      colour = "black"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = efq_m)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
      fill = "black",
      alpha = 0.2
    ) +
    ggplot2::geom_text(
      data = txt_df,
      ggplot2::aes(label = lab, x = x, y = y),
      size = 2.5
    ) +
    ggplot2::facet_wrap(ggplot2::vars(Harvest_Level), nrow = 2) +
    ggplot2::labs(x = "Length", y = "Frequency") +
    ggplot2::geom_hline(yintercept = 0)
}


#' Plot expected length frequencies for a range of selectivity reference points
#'
#' The facet plot covers the current estimated selectivity mode (Smx), together
#' with selectivity changes required to obtain SPR 20%, SPR 40% and the maximum
#' yield. The graphs show the expected values and 80% credible intervals
#' compared to current observations, and can be used to assess whether length
#' frequencies should be able to detect changes in selectivity to these
#' different levels. Note that if reference points do not exist, they are not
#' plotted, so some graphs may be blank except for data.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object plotting expected length frequencies by
#' selectivity
#' @examples
#' eg_lx <- blicc_expect_len(eg_rp, eg_ld)
#' plot_efq_SRP(eg_rp, eg_lx, eg_ld)
#'
plot_efq_SRP <- function(blicc_rp, blicc_lx, blicc_ld) {
  Linf=Galpha=Mk=Fk=Smx=Ss1=Ss2=NB_phi=NULL
  .draw=Lgroup=Current=SPR20=SPR40=LMP=fq=Harvest_Level=NULL
  S20=S40=SMY=efq=fq_lo=fq_hi=N=efq_m=dat_lo=dat_hi=lab=x=y=NULL

  Max_Yield=Selectivity=NULL
  glq <- statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  df <- blicc_rp |>
    dplyr::mutate(Lgroup = list(factor(blicc_ld$Len))) |>
    dplyr::mutate(
      SPR20 = purrr::pmap(
        list(Linf, Galpha, Mk, Fk, S20, Ss1, Ss2),
        blicc_get_efq,
        blicc_ld = blicc_ld,
        glq
      ),
      SPR40 = purrr::pmap(
        list(Linf, Galpha, Mk, Fk, S40, Ss1, Ss2),
        blicc_get_efq,
        blicc_ld = blicc_ld,
        glq
      ),
      Max_Yield = purrr::pmap(
        list(Linf, Galpha, Mk, Fk, SMY, Ss1, Ss2),
        blicc_get_efq,
        blicc_ld = blicc_ld,
        glq
      )
    ) |>
    dplyr::select(`.draw`, NB_phi, Lgroup, SPR20:Max_Yield) |>
    tidyr::unnest(Lgroup:Max_Yield) |>
    dplyr::left_join(dplyr::select(blicc_lx, .draw, Lgroup, Current = efq),
      by = c(".draw", "Lgroup")
    )

  df1 <- df |>
    tidyr::pivot_longer(
      cols = SPR20:Current,
      names_to = "Selectivity",
      values_to = "efq"
    ) |>
    dplyr::filter(!is.na(efq)) |>
    dplyr::mutate(
      fq_lo = stats::qnbinom(0.1, size = NB_phi, mu = efq),
      fq_hi = stats::qnbinom(0.9, size = NB_phi, mu = efq)
    ) |>
    dplyr::group_by(Selectivity, Lgroup) |>
    dplyr::summarise(
      efq_m = mean(efq),
      dat_lo = stats::quantile(fq_lo, probs = c(0.1), na.rm = T, names = F),
      dat_hi = stats::quantile(fq_hi, probs = c(0.9), na.rm = T, names = F),
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

  df2 <- tibble::tibble(LMP = blicc_ld$LMP, fq = blicc_ld$fq)

  ggplot2::ggplot(df1, ggplot2::aes(x = LMP)) +
    ggplot2::geom_col(
      data = df2,
      ggplot2::aes(x = LMP, y = fq),
      fill = "lightblue",
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = df2,
      ggplot2::aes(x = LMP, y = fq),
      shape = 95,
      colour = "black"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = efq_m)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = dat_lo, ymax = dat_hi),
      fill = "black",
      alpha = 0.2
    ) +
    ggplot2::geom_text(
      data = txt_df,
      ggplot2::aes(label = lab, x = x, y = y),
      size = 2.5
    ) +
    ggplot2::facet_wrap(ggplot2::vars(Selectivity), nrow = 2) +
    ggplot2::labs(x = "Length", y = "Frequency") +
    ggplot2::geom_hline(yintercept = 0)
}


#' Spawning potential ratio surface contour plot
#'
#' The plot shows the SPR surface plotted against the fishing mortality (Fk) and
#' the full selectivity (Smx). The current state of the stock is marked as a
#' point and the SPR 40% reference point is marked as a line.
#' All values are calculated as means and uncertainty is not represented to
#' maintain simplicity. The reference line is plotted if it exists within a
#' reasonable range of the current estimated status.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_SPR_contour(eg_rp, eg_ld)
#'
plot_SPR_contour <- function(blicc_rp, blicc_ld) {
  Smx_rp=spr=Lcur=Fcur=`..level..`=NULL  #R CMD Notes

  GridN <- 40
  Linf <- mean(blicc_rp$Linf)
  Galpha <- mean(blicc_rp$Galpha)
  Gbeta <- Galpha/Linf
  Mk <- mean(blicc_rp$Mk)
  Fk <- mean(blicc_rp$Fk)
  Smx <- mean(blicc_rp$Smx)
  Ss1 <- mean(blicc_rp$Ss1)
  Ss2 <- mean(blicc_rp$Ss2)
  F20 <- mean(blicc_rp$F20, na.rm = TRUE)
  F30 <- mean(blicc_rp$F30, na.rm = TRUE)

  glq <-
    statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  Len <- blicc_ld$Len
  LMP <- blicc_ld$LMP
  mb <- mature_biomass_at_length(blicc_ld)
  SPR0 <-
    fSPR(0, Galpha, Gbeta, Mk, RSel = rep(0, length(Len)),
         Len, mb, glq) # Unexploited SPR

  MaxF <- 2 * max(Fk, F20, na.rm = T)
  mSmx <- seq(min(Len - 1), Linf,
    length.out = GridN
  )
  mF <- seq(0, MaxF, length.out = 40)

  SPRRP <- tibble::tibble(Smx_rp = mSmx) |>
    dplyr::mutate(F30 = purrr::map_dbl(
      Smx_rp,
      FSPR_solve,
      Linf,
      Galpha,
      Mk,
      Ss1,
      Ss2,
      0.3,
      blicc_ld,
      glq
    )) |>
    dplyr::filter(F30 <= MaxF)

  SPRCurve <- tibble::as_tibble(expand.grid(mF = mF, Smx = mSmx)) |>
    dplyr::mutate(spr = purrr::map2_dbl(Smx, mF, fSPR2, Galpha,
                                        Gbeta, Mk, Ss1, Ss2, Len, LMP, mb, glq)
                                  / SPR0)

  cp <- ggplot2::ggplot(SPRCurve, ggplot2::aes(x = Smx, y = mF, z = spr)) +
    ggplot2::geom_point(
      tibble::tibble(
        Lcur = Smx,
        Fcur = Fk,
        spr = 0
      ),
      mapping = ggplot2::aes(x = Lcur, y = Fcur),
      size = 3
    ) +
    ggplot2::labs(
      x = "Length at full selectivity", y = "Fishing Mortality (/k)", colour =
        "SPR"
    ) +
    ggplot2::geom_contour(ggplot2::aes(colour = ..level..))

  if (nrow(SPRRP) > 0) {   # Reference points exist / can be plotted
    cp <- cp +
      ggplot2::geom_line(
        data = SPRRP,
        ggplot2::aes(x = Smx_rp, y = F30),
        colour = "darkred",
        size = 2,
        inherit.aes = F
      ) +
      ggplot2::annotate(
        "text",
        x = max(SPRRP$Smx_rp) + 0.5,
        y = max(SPRRP$F30) + 0.3,
        label = "SPR30",
        colour = "darkred"
      )
  }
  return(cp)
}


#' Yield-per-recruit surface contour plot
#'
#' The plot shows the YPR surface plotted against the fishing mortality (Fk) and
#' the full selectivity (Smx). The current state of the stock is marked as a
#' point and the F0.1 reference point is marked as a line. All values are
#' calculated as means and uncertainty is not represented to maintain
#' simplicity. The reference line is plotted if it exists within a
#' reasonable range of the current estimated status.
#'
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' plot_YPR_contour(eg_rp, eg_ld)
#'
plot_YPR_contour <- function(blicc_rp, blicc_ld) {
  Smx_rp=ypr=F0.1=Lcur=Fcur=`..level..`=NULL  # R CMD Notes
  GridN <- 40
  Linf <- mean(blicc_rp$Linf)
  Galpha <- mean(blicc_rp$Galpha)
  Gbeta <- Galpha/Linf
  Mk <- mean(blicc_rp$Mk)
  Fk <- mean(blicc_rp$Fk)
  Smx <- mean(blicc_rp$Smx)
  Ss1 <- mean(blicc_rp$Ss1)
  Ss2 <- mean(blicc_rp$Ss2)
  F20 <- mean(blicc_rp$F20, na.rm = TRUE)
  F30 <- mean(blicc_rp$F30, na.rm = TRUE)

  glq <-
    statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  Len <- blicc_ld$Len
  LMP <- blicc_ld$LMP
  weight <- weight_at_length(blicc_ld)

  MaxF <- 2 * max(Fk, F20, na.rm = T)
  mSmx <- seq(min(Len - 1), Linf,
    length.out = GridN
  )
  mF <- seq(0, MaxF, length.out = GridN)

  YPRRP <- tibble::tibble(Smx_rp = mSmx) |>
    dplyr::mutate(F0.1 = purrr::map_dbl(Smx_rp, F01_solve, Linf,
                                        Galpha, Mk, Ss1, Ss2, blicc_ld, glq)) |>
    dplyr::filter(F0.1 <= MaxF)

  YPRCurve <- tibble::as_tibble(expand.grid(mF = mF, Smx = mSmx)) |>
    dplyr::mutate(ypr = purrr::map2_dbl(Smx, mF, fYPR2, Galpha, Gbeta, Mk,
                                        Ss1, Ss2, Len, LMP, weight, glq))

  cp <- ggplot2::ggplot(YPRCurve, ggplot2::aes(x = Smx, y = mF, z = ypr)) +
    ggplot2::geom_point(
      tibble::tibble(Lcur = Smx, Fcur = Fk),
      mapping = ggplot2::aes(x = Lcur, y = Fcur),
      size = 3,
      inherit.aes = F
    ) +
    ggplot2::labs(
      x = "Length at full selectivity", y = "Fishing Mortality (/k)", colour =
        "Relative Yield"
    ) +
    ggplot2::geom_contour(ggplot2::aes(colour = ..level..))

  if (nrow(YPRRP) > 0) {   # Reference points exist / can be plotted
    cp <- cp +
      ggplot2::geom_line(
        data = YPRRP,
        ggplot2::aes(x = Smx_rp, y = F0.1),
        colour = "darkred",
        size = 2,
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
