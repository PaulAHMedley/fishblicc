# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# PLOTTING FUNCTIONS
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

# BLICC Plotting Functions ------------------------------------------------

# Produces summary plots and tables based on yield-per-recruit
# SPR density plot
# SPR surface plot Smx vs Fk: current position +30% SPR reference point uncertainty as size for 80% credible interval
# YPR surface plot Smx vs Fk: current position +F0.1 reference point uncertainty as size for 80% credible interval


#' Plot the observed and expected frequency data
#' 
#' The graph shows the observed frequency and the expected frequency together with
#' the 80% credible intervals for the estimate as well as 80% credible intervals 
#' for the observations. The latter takes into account the negative binomial 
#' observation error, so on average 8 out of ten observations would be expected 
#' within these bounds. 
#'  
#' @export
#' @param blicc_res  A posterior draws and reference points tibble from `blicc_ref_pts` function.
#' @param blicc_lex  A tibble from `blicc_expect_len` function.
#' @param blicc_ld   A data list from the `blicc_dat` function.
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' l_df <- blicc_expect_len(res_df, ld)
#' plot_expected_frequency(l_df, ld)
#' }
#' 
plot_expected_frequency <- function(blicc_res, blicc_lex, blicc_ld) {
  df1 <- blicc_res %>%
    dplyr::select(.draw, NB_phi) %>%
    dplyr::left_join(blicc_lex, by=".draw") %>%
    dplyr::select(.draw, Lgroup, efq, NB_phi) %>%
    dplyr::mutate(fq_lo = qnbinom(0.1, size = NB_phi, mu = efq),
           fq_hi = qnbinom(0.9, size = NB_phi, mu = efq)) %>%
    dplyr::group_by(Lgroup) %>%
    dplyr::summarise(
      efq_m = mean(efq),
      efq_lo = quantile(efq, probs = c(0.1), names = F),
      efq_hi = quantile(efq, probs = c(0.9), names = F),
      dat_lo = quantile(fq_lo, probs = c(0.1), names = F),
      dat_hi = quantile(fq_hi, probs = c(0.9), names = F),
    )  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])
  
  df2 <- tibble::tibble(LMP=blicc_ld$LMP, fq=blicc_ld$fq)
  
  ggplot2::ggplot(df1, aes(x = LMP)) +
    geom_col(data=df2, aes(x=LMP, y=fq), fill="lightblue", alpha=0.5) +
    geom_point(data=df2, aes(x=LMP, y=fq), shape=95, size=5, colour="black") +
    geom_line(aes(y = efq_m)) +
    geom_ribbon(aes(ymin = efq_lo, ymax = efq_hi), fill="blue", alpha = 0.3) +
    geom_ribbon(aes(ymin = dat_lo, ymax = dat_hi), fill="black", alpha = 0.2) +
    labs(x = "Length", y = "Frequency") +
    geom_hline(yintercept=0)
} 


#' Plot selectivity and 50% and 80% credible intervals by length
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' l_df <- blicc_expect_len(res_df, ld)
#' plot_selectivity(l_df, ld)
#' }
#' 
plot_selectivity <- function(blicc_lex, blicc_ld) {
  blicc_lex %>%
    select(Lgroup, sel) %>%
    mutate(LMP=blicc_ld$LMP[as.integer(Lgroup)]) %>%
    group_by(Lgroup) %>%
  summarise(
    sel_m = mean(sel),
    sel_01 = quantile(sel, probs = c(0.10), names = F),
    sel_25 = quantile(sel, probs = c(0.25), names = F),
    sel_75 = quantile(sel, probs = c(0.75), names = F),
    sel_90 = quantile(sel, probs = c(0.90), names = F),
  )  %>%
  ungroup() %>%
  ggplot(df1, aes(x = LMP)) +
    geom_ribbon(aes(ymin = sel_01, ymax = sel_90), fill="#546E8A", alpha = 0.7) +
    geom_ribbon(aes(ymin = sel_25, ymax = sel_75), fill="#B0BED5", alpha = 0.7) +
    geom_line(aes(y = sel_m)) +
    labs(x = "Length (cm)", y = "Selectivity")
} 


#' Plot of the estimated spawning potential ratio probability density
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' l_df <- blicc_expect_len(res_df, ld)
#' plot_SPR_density(res_df)
#' }
#' 
plot_SPR_density <- function(blicc_res) {
  ggplot(blicc_res, aes(SPR)) +
    geom_density(fill = "lightblue") +
    geom_vline(aes(xintercept = 0.2),
               colour = "red",
               linetype = 2) +
    geom_vline(aes(xintercept = 0.4),
               colour = "green",
               linetype = 2) +
    coord_cartesian(xlim=c(0, 1)) +
    labs(x="Spawning Potential Ratio", y = "Probability Density") +
    theme_pubr()
}



#' Plot fishing mortality relative to the SPR 40% reference point (Fk/F40)
#' 
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' l_df <- blicc_expect_len(res_df, ld)
#' plot_FkF40_density(res_df)
#' }
#' 
plot_FkF40_density <- function(blicc_res) {
  StatAreaUnderDensity <- ggplot2::ggproto(
    "StatAreaUnderDensity", Stat,
    required_aes = "x",
    compute_group = function(data, scales, xlim = NULL, n = 50) {
      fun <- approxfun(density(data$x))
      StatFunction$compute_group(data, scales, fun = fun, xlim = xlim, n = n)
    }
  )
  
  stat_aud <- function(mapping = NULL, data = NULL, geom = "area",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, n = 50, xlim=NULL,  
                       ...) {
    layer(
      stat = StatAreaUnderDensity, data = data, mapping = mapping, geom = geom, 
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(xlim = xlim, n = n, ...))
  }
  
  ggplot(blicc_res, aes(Fk/F40)) +
    geom_density(fill = "#A11826") +
    geom_vline(aes(xintercept = 1.0),
               colour = "green",
               linetype = 1) +
    stat_aud(geom="area",
             aes(fill="good"),
             xlim = c(0, 1), 
             alpha = 1.0, show.legend=FALSE) +
    scale_fill_manual(values=c("good" = "#204030")) +
    labs(y = "Probability Density") +
    theme_pubr()
}



#' Plot expected length frequencies for a range of fishing mortality reference points
#' 
#' The facet plot covers the current estimated fishing mortality, together with
#' fishing mortalities required to obtain SPR 20%, SPR 30% and SPR 40%. The 
#' graphs show the expected values and 80% credible intervals compared to current
#' observations, and can be used to assess whether length frequencies should
#' be able to detect changes in fishing mortality to these different levels.
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' l_df <- blicc_expect_len(res_df, ld)
#' plot_efq_FRP(res_df, l_df, ld)
#' }
#' 
plot_efq_FRP <- function(blicc_res, blicc_lex, blicc_ld) {
  
  df <- blicc_res %>%
    dplyr::mutate(Lgroup=list(factor(blicc_ld$Len))) %>%
    dplyr::mutate(SPR20 = purrr::pmap(list(Linf, Galpha, Mk, F20, Smx, Ss1, Ss2), 
                                      blicc_get_efq, blicc_ld=blicc_ld),
                  SPR30 = purrr::pmap(list(Linf, Galpha, Mk, F30, Smx, Ss1, Ss2), 
                                      blicc_get_efq, blicc_ld=blicc_ld),
                  SPR40 = purrr::pmap(list(Linf, Galpha, Mk, F40, Smx, Ss1, Ss2), 
                                      blicc_get_efq, blicc_ld=blicc_ld)) %>%
    dplyr::select(`.draw`, NB_phi, Lgroup, SPR20:SPR40) %>%
    tidyr::unnest(Lgroup:SPR40) %>%
    left_join(select(blicc_lex, .draw, Lgroup, Current=efq), by=c(".draw", "Lgroup"))
  
  df1 <- df %>%
    pivot_longer(cols=SPR20:Current, names_to="Harvest_Level", values_to="efq") %>%
    mutate(fq_lo = qnbinom(0.1, size = NB_phi, mu = efq),
           fq_hi = qnbinom(0.9, size = NB_phi, mu = efq)) %>%
    group_by(Harvest_Level, Lgroup) %>%
    summarise(
      efq_m = mean(efq),
      dat_lo = quantile(fq_lo, probs = c(0.1), names = F),
      dat_hi = quantile(fq_hi, probs = c(0.9), names = F),
    )  %>%
    ungroup() %>%
    mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])
  
  df2 <- tibble(LMP=blicc_ld$LMP, fq=blicc_ld$fq)
  
  ggplot(df1, aes(x = LMP)) +
    geom_col(data=df2, aes(x=LMP, y=fq), fill="lightblue", alpha=0.5) +
    geom_point(data=df2, aes(x=LMP, y=fq), shape=95, size=5, colour="black") +
    geom_line(aes(y = efq_m)) +
    geom_ribbon(aes(ymin = dat_lo, ymax = dat_hi), fill="black", alpha = 0.2) +
    facet_wrap(vars(Harvest_Level), nrow=2) +
    labs(x = "Length", y = "Frequency") +
    geom_hline(yintercept=0)
} 



#' Plot expected length frequencies for a range of selectivity reference points
#' 
#' The facet plot covers the current estimated selectivity mode (Smx), together
#' with selectivity changes required to obtain SPR 20%, SPR 40% and the maximum 
#' yield. The graphs show the expected values and 80% credible intervals 
#' compared to current observations, and can be used to assess whether length 
#' frequencies should be able to detect changes in selectivity to these 
#' different levels.
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' l_df <- blicc_expect_len(res_df, ld)
#' plot_efq_SRP(res_df, l_df, ld)
#' }
#' 
plot_efq_SRP <- function(blicc_res, blicc_lex, blicc_ld) {
  
  df <- blicc_res %>%
    dplyr::mutate(Lgroup=list(factor(blicc_ld$Len))) %>%
    dplyr::mutate(SPR20 = purrr::pmap(list(Linf, Galpha, Mk, Fk, S20, Ss1, Ss2), 
                                      blicc_get_efq, blicc_ld=blicc_ld),
                  SPR40 = purrr::pmap(list(Linf, Galpha, Mk, Fk, S40, Ss1, Ss2), 
                                      blicc_get_efq, blicc_ld=blicc_ld),
                  Max_Yield = purrr::pmap(list(Linf, Galpha, Mk, Fk, SMY, Ss1, Ss2), 
                                          blicc_get_efq, blicc_ld=blicc_ld)) %>%
    dplyr::select(`.draw`, NB_phi, Lgroup, SPR20:Max_Yield) %>%
    tidyr::unnest(Lgroup:Max_Yield) %>%
    dplyr::left_join(select(blicc_lex, .draw, Lgroup, Current=efq), 
                     by=c(".draw", "Lgroup"))
  
  df1 <- df %>%
    dplyr::pivot_longer(cols=SPR20:Current, names_to="Selectivity", values_to="efq") %>%
    dplyr::filter(!is.na(efq)) %>%
    dplyr::mutate(fq_lo = qnbinom(0.1, size = NB_phi, mu = efq),
                  fq_hi = qnbinom(0.9, size = NB_phi, mu = efq)) %>%
    dplyr::group_by(Selectivity, Lgroup) %>%
    dplyr::summarise(
             efq_m = mean(efq),
             dat_lo = quantile(fq_lo, probs = c(0.1), names = F),
             dat_hi = quantile(fq_hi, probs = c(0.9), names = F),
             N = n()
             )  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(LMP = blicc_ld$LMP[as.integer(Lgroup)])
  
  txt_df <- df1 %>%
    dplyr::group_by(Selectivity) %>%
    dplyr::summarise (N = first(N), y=max(efq_m)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lab=paste0("N = ", format(N)), x=max(LMP=blicc_ld$LMP)*0.9)  
  
  df2 <- tibble::tibble(LMP=blicc_ld$LMP, fq=blicc_ld$fq)
  
  ggplot2::ggplot(df1, aes(x = LMP)) +
    geom_col(data=df2, aes(x=LMP, y=fq), fill="lightblue", alpha=0.5) +
    geom_point(data=df2, aes(x=LMP, y=fq), shape=95, size=5, colour="black") +
    geom_line(aes(y = efq_m)) +
    geom_ribbon(aes(ymin = dat_lo, ymax = dat_hi), fill="black", alpha = 0.2) +
    geom_text(data=txt_df, aes(label=lab, x=x, y=y)) +
    facet_wrap(vars(Selectivity), nrow=2) +
    labs(x = "Length", y = "Frequency") +
    geom_hline(yintercept=0)
} 


#' Spawning potential ratio surface contour plot
#' 
#' The plot shows the SPR surface plotted against the fishing mortality (Fk) and 
#' the full selectivity (Smx). The current state of the stock is marked as a point 
#' and the SPR 40% reference point is marked as a line. All values are calculated as 
#' means and uncertainty is not represented to maintain simplicity.  
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' blicc_SPR_contour_plot(res_df, ld)
#' }
#'
SPR_contour_plot <- function(blicc_res, blicc_ld)  {
  # SPRContour
  GridN = 40
  Linf <- mean(blicc_res$Linf)
  Galpha <- mean(blicc_res$Galpha)
  Mk <- mean(blicc_res$Mk)
  Fk <- mean(blicc_res$Fk)
  Smx <- mean(blicc_res$Smx)
  Ss1 <- mean(blicc_res$Ss1)
  Ss2 <- mean(blicc_res$Ss2)
  F20 <- mean(blicc_res$F20)
  F30 <- mean(blicc_res$F30)
  
  glq <- statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  Gbeta <- Galpha/Linf
  Len <- blicc_ld$Len
  LMP <- blicc_ld$LMP
  mb <- mature_biomass_at_length(blicc_ld)
  SPR0 <- fSPR(0, Galpha, Gbeta, Mk, RSel=rep(0, length(Len)), Len, mb, glq)   #Unexploited SPR
  
  MaxF <- 2 * max(Fk, F20)
  mSmx <- seq(min(Len-1), Linf,
              length.out = GridN)
  mF <- seq(0, MaxF, length.out = 40)
  
  SPRRP <- tibble(Smx_rp = mSmx) %>%
    mutate(F30 = map_dbl(Smx_rp, FSPR_solve, Linf, Galpha, Mk, Ss1, Ss2, 0.3, blicc_ld, glq)) %>%
    filter(F30 <= MaxF)
  
  # fSPR2 <- function(Smx, Fk, Galpha, Gbeta, Mk, Ss1, Ss2, Len, LMP, mb, glq)
  SPRCurve <- as_tibble(expand.grid(mF = mF, Smx = mSmx)) %>%
    mutate(spr =  map2_dbl(Smx, mF, fSPR2, Galpha, Gbeta, Mk, Ss1, Ss2, Len, LMP, mb, glq) / SPR0)
  
  ggplot(SPRCurve, aes(x = Smx, y = mF, z = spr)) +
    geom_point(
      tibble(
        Lcur = Smx,
        Fcur = Fk,
        spr = 0
      ),
      mapping = aes(x = Lcur, y = Fcur),
      size = 3
    ) +
    labs(x = "Length at 50% Capture", y = "Fishing Mortality (/k)", colour =
           "SPR") +
    geom_contour(aes(colour = ..level..)) +
    geom_line(
      data = SPRRP,
      aes(x = Smx_rp, y = F30),
      colour = "darkred",
      size = 2,
      inherit.aes = F
    ) +
    annotate(
      "text",
      x = max(SPRRP$Smx_rp) + 1,
      y = max(SPRRP$F30) + 0.1,
      label = "SPR30",
      colour = "darkred"
    )
}


#' Yield-per-recruit surface contour plot
#' 
#' The plot shows the YPR surface plotted against the fishing mortality (Fk) and 
#' the full selectivity (Smx). The current state of the stock is marked as a point 
#' and the F0.1 reference point is marked as a line. All values are calculated as 
#' means and uncertainty is not represented to maintain simplicity.  
#' 
#' @export
#' @inheritParams plot_expected_frequency
#' @return ggplot geom object for plotting
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#' slim <- blicc_fit(ld)
#' res_df <- blicc_ref_pts(slim, ld)
#' blicc_YPR_contour_plot(res_df, ld)
#' }
#'
YPR_contour_plot <- function(blicc_res, blicc_ld) {
  GridN = 40
  Linf <- mean(blicc_res$Linf)
  Galpha <- mean(blicc_res$Galpha)
  Mk <- mean(blicc_res$Mk)
  Fk <- mean(blicc_res$Fk)
  Smx <- mean(blicc_res$Smx)
  Ss1 <- mean(blicc_res$Ss1)
  Ss2 <- mean(blicc_res$Ss2)
  F20 <- mean(blicc_res$F20)
  F30 <- mean(blicc_res$F30)
  
  glq <- statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  Gbeta <- Galpha/Linf
  Len <- blicc_ld$Len
  LMP <- blicc_ld$LMP
  weight <- weight_at_length(blicc_ld)
  
  MaxF <- 2 * max(Fk, F20)
  mSmx <- seq(min(Len-1), Linf,
              length.out = GridN)
  mF <- seq(0, MaxF, length.out = GridN)
  
  YPRRP <- tibble(Smx_rp = mSmx) %>%
    mutate(F0.1 = map_dbl(Smx_rp, F01_solve, Linf, Galpha, Mk, Ss1, Ss2, blicc_ld, glq)) %>%
    filter(F0.1 <= MaxF)
  
  YPRCurve <- as_tibble(expand.grid(mF = mF, Smx = mSmx)) %>%
    mutate(ypr = map2_dbl(Smx, mF, fYPR2, Galpha, Gbeta, Mk, Ss1, Ss2, Len, LMP, weight, glq))
  
  ggplot(YPRCurve, aes(x = Smx, y = mF, z = ypr)) +
    geom_point(
      tibble(Lcur = Smx, Fcur = Fk),
      mapping = aes(x = Lcur, y = Fcur),
      size = 3,
      inherit.aes = F
    ) +
    labs(x = "Length at full selectivity", y = "Fishing Mortality (/k)", colour =
           "Relative Yield") +
    geom_contour(aes(colour = ..level..)) +
    geom_line(
      data = YPRRP,
      aes(x = Smx_rp, y = F0.1),
      colour = "darkred",
      size = 2,
      inherit.aes = F
    ) +
    annotate(
      "text",
      x = max(YPRRP$Smx_rp),
      y = max(YPRRP$F0.1) + 0.1,
      label = "F0.1",
      colour = "darkred"
    )
}


