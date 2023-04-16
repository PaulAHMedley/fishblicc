# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# Model functions
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

# BLICC Model Functions ---------------------------------------------------

# Model 1
#' Calculate a logistic selectivity curve for a length vector
#'
#' A logistic takes a 50% selectivity and a slope parameter. The resulting
#' selectivity varies from 0 to 1.0, with 1.0 being the asymptote. It is
#' calculated for each position defined in the length vector.
#'
#' @export
#' @param  Sp   Vector (length 2) of the selectivity parameters (50%,
#'              and slope) for the logistic selectivity function.
#' @param  LMP  A vector of lengths (usually mid-points for length bins)
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_logistic(Sp = c(30, 0.5), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_logistic <- function(Sp, LMP) {
  SL <- double(length(LMP))
  diff <- LMP - Sp[1L]
  SL <- 1.0/(1.0 + exp(-Sp[2L] * diff))
  return(SL)
}

# Model 2
#' Calculate a normal selectivity curve for a length vector
#'
#' A normal takes one slope (1/sigma^2) parameter for the
#' both sides of the mode parameter (mu). The resulting selectivity varies
#' from 0 to 1.0, with 1.0 being the mode. It is calculated for each position
#' defined in the length vector.
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 2) of the selectivity parameters (mode,
#'              and 1/sigma^2) for the normal selectivity function.
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_normal(Sp = c(35, 0.2), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_normal <- function(Sp, LMP) {
  SL <- double(length(LMP))
  diff <- (LMP - Sp[1L]) ^ 2
  SL <- exp(-Sp[2L] * diff)
  return(SL)
}

# Model 3
#' Calculate a single-sided normal selectivity curve for a length vector
#'
#' A single-sided normal takes one slope (1/sigma^2) parameter for the
#' left side of the mode parameter (mu). For the right side the slope is assumed
#' to be zero. The resulting selectivity varies from 0 to 1.0,
#' with 1.0 being the asymptote. It is calculated for each position defined
#' in the length vector.
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 2) of the selectivity parameters (mode,
#'              and left side 1/sigma^2) for the normal selectivity function.
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_ssnormal(Sp = c(35, 0.2), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_ssnormal <- function(Sp, LMP) {
  SL <- double(length(LMP))
  S1 <- LMP < Sp[1L]
  SL[S1] <- exp(-Sp[2L] * ((LMP[S1] - Sp[1L]) ^ 2))
  SL[!S1] <- 1.0
  return(SL)
}

# Model 4
#' Calculate a double-sided normal selectivity curve for a length vector
#'
#' A double-sided normal takes two separate slope (1/sigma^2) parameters
#' around the mode parameter (mu) that change the slope independently
#' on either side of the mode. The selectivity varies from 0 to 1.0,
#' with 1.0 being at the mode. It is calculated for each position defined
#' in the length vector. Setting the second slope parameter to zero
#' creates a flat-top selectivity (see `Rsel_ssnormal`).
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 3) of the selectivity parameters (mode,
#'              left  and right side 1/sigma^2) for the normal
#'              selectivity function
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_dsnormal(Sp = c(35, 0.1, 0.2), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_dsnormal <- function(Sp, LMP) {
  # Double sided normal
  SL <- double(length(LMP))
  S1 <- LMP < Sp[1L]
  SL[S1] <- exp(-Sp[2L] * (LMP[S1] - Sp[1L]) ^ 2)
  SL[!S1] <- exp(-Sp[3L] * (LMP[!S1] - Sp[1L]) ^ 2)
  return(SL)
}


#' Calculate the survival of a fish cohort to sequential length boundaries
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the
#' survival to each length bin boundary taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. Survival can be used to estimate population
#' numbers-at-length and catch length frequency, among other things.
#' The model uses the Gauss-Laguerre quadrature rule for integration,
#' so nodes and weights for this must be provided.
#'
#' @export
#' @param  node   Nodes for the Gauss-Laguerre quadrature rule
#' @param  wt     Weights for the Gauss-Laguerre quadrature rule
#' @param  Len    Vector of lower-bounds for length intervals
#' @param  Zki    Vector of total mortality (in time units of growth rate
#' parameter K) to be applied in each length interval.
#' This vector should be the same length as Len.
#' @param  Galpha Alpha parameter for the Gamma probability density function
#' that governs growth variability.
#' @param  Gbeta  Beta parameter for the Gamma probability density function
#' that governs growth variability.
#' @return A vector of proportion surviving to each length bin lower bound.
#' @examples
#' glq <- statmod::gauss.quad(90, kind = "laguerre", alpha = 0.0)
#' S <- RSurvival_Est(glq$nodes, glq$weights, Len=15:55,
#'                    Zki=c(rep(1.5, 10), rep(3, 31)), 100, 100/50)
#' plot(y=S, x=15:55, type="l")
RSurvival_Est <- function(node, wt, Len, Zki, Galpha, Gbeta)  {
  zsum <- function(x, Lr, Z) {
    # used to apply sequenced sum of mortality
    return(sum(log(x+Lr)*Z))
  }
  lgamma_Galpha <- lgamma(Galpha)
  nv <- length(node)
  LN <-  length(Len)
  surv <- double(LN)
  x_beta <- node / Gbeta
  log_x_beta <- log(x_beta)

  ss <- log(node + Gbeta * Len[1]) * (Galpha - 1.0) -
    Gbeta * Len[1] - lgamma_Galpha
  surv[1] <- sum(exp(ss) * wt)
  ss <- (log(x_beta + Len[2] - Len[1]) * (-Zki[1])) + log_x_beta * Zki[1] +
    log(node + Gbeta * Len[2]) * (Galpha - 1.0) -
    Gbeta * Len[2] - lgamma_Galpha
  surv[2] <- sum(exp(ss) * wt)

  for (Li in 3:LN) {
    Ln <- Len[Li]
    Lrange <- Ln - Len[1:(Li - 1)]

    Zii <- c(-Zki[1], Zki[1:(Li - 2)] - Zki[2:(Li - 1)])
    v2 <- log_x_beta * Zki[Li - 1]
    lim <- vapply(x_beta, zsum, FUN.VALUE=0, Lr=Lrange,
                  Z=Zii, USE.NAMES=FALSE)   # sapply not safe to use here
    ss <- lim + v2 + log(node + Gbeta * Ln) * (Galpha - 1.0) -
      Gbeta * Ln - lgamma_Galpha
    surv[Li] <- sum(exp(ss) * wt)
  }
  return(surv)
}


#' Calculate the population relative numbers-at-length within length bins
#'
#' The relative numbers in each length bin is based on the definite integral
#' over the bin growth transition time:
#'     N_i = (S_i - S_{i+1})/Z_{ki}
#'
#' @export
#' @param  surv Survival to the lower boundary of each length bin
#' (see `RSurvival_est()`)
#' @param  Zki  Total mortality within each bin (time in units of the
#' growth rate K)
#' @return A vector of relative numbers of fish within each length bin.
#' @examples
#' glq <- statmod::gauss.quad(90, kind = "laguerre", alpha = 0.0)
#' S <- RSurvival_Est(glq$nodes, glq$weights, Len=15:55, Zki=c(rep(1.5, 10),
#'                    rep(3, 31)), 100, 100/50)
#' P <- RNinInterval(S, Zki=c(rep(1.5, 10), rep(3, 31)))
#' plot(y=P, x=15:55, type="l")
#'
RNinInterval <- function (surv, Zki) {
  n <- length(surv)
  return(c(surv[-n] - surv[-1], surv[n]) / Zki)
}


#' Mature biomass at length
#'
#' Mature biomass is calculated based on the length-weight relationship and
#' the maturity ogive. In terms of calculating SPR and reference points,
#' only the exponent parameter b and maturity ogive have any effect.
#' The length-weight parameter a is only used for scaling purposes
#' and need not be provided.
#'
#' @export
#' @param  blicc_ld A standard data list created by `blicc_dat`
#' @examples
#' ld <- blicc_dat(LLB = 10:35, a=1.0e-4, b=2.95, L50=24,
#'                 fq=c(0,1,2,5,26,70,72,66,36,34,25,24,
#'                      20,10,12,5,3,5,6,4,2,0,2,0,1,0), Linf=c(32, 3))
#' mb <- mature_biomass_at_length(ld)
#' plot(y=mb, x=ld$LMP, type="l")
#'
mature_biomass_at_length <- function(blicc_ld) {
  return( (blicc_ld$a*blicc_ld$LMP^blicc_ld$b) /
            (1 + exp( -blicc_ld$Ls * (blicc_ld$LMP-blicc_ld$Lm) )) )
}


#' Weight at length
#'
#' The weight is calculated using the standard length-weight model W=a L^b.
#' In terms of reference points only the exponent parameter b is required,
#' and the parameter a is only used for scaling.
#'
#' @export
#' @inheritParams  mature_biomass_at_length
#' @examples
#' ld <- blicc_dat(LLB = 10:35, a=1.0e-4, b=2.95, L50=24,
#'                 fq=c(0,1,2,5,26,70,72,66,36,34,25,24,
#'                      20,10,12,5,3,5,6,4,2,0,2,0,1,0), Linf=c(32, 3))
#' wt <- weight_at_length(ld)
#' plot(y=wt, x=ld$LMP, type="l")
#'
weight_at_length <- function(blicc_ld) {
  return(blicc_ld$a*blicc_ld$LMP^blicc_ld$b)
}


#' Calculate the survival of a fish cohort to sequential length boundaries
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the
#' survival to each length bin boundary taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. Survival can be used to estimate population
#' numbers-at-length and catch length frequency, among other things.
#' The model uses the Gauss-Laguerre quadrature rule for integration,
#' so nodes and weights for this must be provided.
#'
#' @param Galpha Alpha parameter for the Gamma probability density function
#' that governs growth variability.
#' @param Gbeta  Rate parameter for the Gamma distribution growth variability
#' (=Galpha/Linf)
#' @param Mk     Natural mortality divided by the growth rate K
#' @param Fk     Fishing mortality divided by the growth rate K for each gear
#' making a contribution
#' @param FSel    A list of all the selectivities for each length bin
#' @param blicc_ld   A data list from the `blicc_dat` function.
#' @return A list of the population size in each length bin and a list of
#' fishing mortalities at length for each gear.
#' @examples
#' S <- RPop_F(100, 100, 0.2, exp(eg_ld$polFk), exp(eg_ld$polSm))
#' plot(y=S, x=15:55, type="l")
RPop_F <- function(Galpha, Gbeta, Mk, Fk, FSel, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L
  for (gi in 1:blicc_ld$NG) {
    if (blicc_ld$Fkg[gi] > 0) {
      FSel[[gi]] <- FSel[[gi]] * Fk[blicc_ld$Fkg[gi]]  # Fishing mortality
      Zki <- Zki + FSel[[gi]]                          # Total mortality
    }
  }

  N_L <- with(blicc_ld,
              CPop_Len(gl_nodes, gl_weights,
                       Len, Zki, Galpha, Gbeta) )
  return(list(N_L=N_L, Fki=FSel))
}


#' Calculate the survival of a fish cohort to sequential length boundaries
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the
#' survival to each length bin boundary taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. Survival can be used to estimate population
#' numbers-at-length and catch length frequency, among other things.
#' The model uses the Gauss-Laguerre quadrature rule for integration,
#' so nodes and weights for this must be provided.
#'
#' @param Sm  Vector of parameters for all the selectivity functions
#' flat-topped selectivity
#' @param blicc_ld  A data list from the `blicc_dat` function.
#' @return A list of fishing mortalities at length for each gear.
#' @examples
#' S <- RSelectivities(exp(eg_ld$polSm), eg_ld)
#' plot(y=S[[1]], x=eg_ld$LMP, type="l")
RSelectivities <- function(Sm, blicc_ld) {
  Ski <- list()
  for (gi in 1:blicc_ld$NG) {
    Indx <- blicc_ld$Spar[gi,]
    Indx <- Indx[Indx>0]
    #
    # Need to be changed to C functions or a C function
    #
    Ski[[gi]] <- switch(blicc_ld$fSel[gi],
                        Rsel_logistic(Sm[Indx], blicc_ld$LMP),
                        Rsel_normal(Sm[Indx], blicc_ld$LMP),
                        Rsel_ssnormal(Sm[Indx], blicc_ld$LMP),
                        Rsel_dsnormal(Sm[Indx], blicc_ld$LMP))
  }
  return(Ski)
}

#' Calculate the spawning potential with no fishing
#'
#' The model calculates the spawning potential by multiplying the proportion of
#' numbers in each length bin with no fishing by the mature biomass per recruit.
#'
#' @param Galpha Alpha parameter for the Gamma probability density function
#' that governs growth variability.
#' @param Gbeta  Rate parameter for the Gamma distribution growth variability
#' (=Galpha/Linf)
#' @param Mk     Natural mortality divided by the growth rate K
#' @param blicc_ld   A data list from the `blicc_dat` function.
#' @return A real of the spawning potential when F=0.
#' @examples
#' S <- RSPR_0(100, 100, 0.2, mature_biomass_at_length(eg_ld), eg_ld)
#' plot(y=S, x=15:55, type="l")
RSPR_0 <- function(Galpha, Gbeta, Mk, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L
  Rsurv <- with(blicc_ld,
                RSurvival_Est(gl_nodes, gl_weights,
                              Len, Zki, Galpha, Gbeta) )
  return(sum(RNinInterval(Rsurv, Zki) * blicc_ld$ma_L))
}
