# BLICC Model Functions ---------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


# Model 1

#' Calculate a logistic selectivity curve for a length vector
#'
#' A logistic takes a 50% selectivity and a slope parameter. The resulting
#' selectivity varies from 0 to 1.0, with 1.0 being the asymptote. It is
#' calculated for each position defined in the length vector.
#'
#' @export
#' @param  Sp   Vector (length 2) of the selectivity parameters (50%, and slope)
#'   for the logistic selectivity function.
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
#' A normal takes one slope (1/sigma^2) parameter for both sides of the location
#' parameter (mu). The resulting selectivity varies from 0 to 1.0, with 1.0
#' being the mode. It is calculated for each position defined in the length
#' vector.
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 2) of the selectivity parameters (mode, and
#'   1/sigma^2) for the normal selectivity function.
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
#' A single-sided normal takes one slope (1/sigma^2) parameter for the left side
#' of the location parameter (mu). For the right side the slope is assumed to be
#' zero. The resulting selectivity varies from 0 to 1.0, with 1.0 being the
#' asymptote. It is calculated for each position defined in the length vector.
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 2) of the selectivity parameters (mode, and left
#'   side 1/sigma^2) for the normal selectivity function.
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
#' A double-sided normal takes two separate slope (1/sigma^2) parameters around
#' the location parameter (mu) that change the slope independently on either
#' side of the mode. The selectivity varies from 0 to 1.0, with 1.0 being at the
#' mode. It is calculated for each position defined in the length vector.
#' Setting the second slope parameter to zero creates a flat-top selectivity
#' ([Rsel_ssnormal]).
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 3) of the selectivity parameters (mode, left  and
#'   right side 1/sigma^2) for the normal selectivity function
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


#' Calculate the population size and fishing mortality applied within the length
#' bins for the BLICC model
#'
#' The model integrates over the growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the relative
#' population number in each length bin taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. The model uses the Gauss-Laguerre quadrature
#' rule for integration, so nodes and weights for this must be provided. The
#' function returns the expected population size in each length bin, and the
#' fishing mortalities.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param Galpha Alpha parameter for the Gamma probability density function that
#'   governs growth variability.
#' @param Gbeta  Rate parameter for the Gamma distribution growth variability
#'   (=Galpha/Linf)
#' @param Mk     Natural mortality divided by the growth rate K
#' @param Fk     Fishing mortality divided by the growth rate K for each gear
#'   making a contribution
#' @param FSel    A list of all the selectivities for each length bin
#' @return A list of the population size in each length bin for each time period 
#'   and a list of fishing mortalities at length for each gear / time period.
#' @examples
#' Sel <- Rselectivities(exp(trgl_ld$polSm), trgl_ld)
#' S <- Rpop_F(100, 100/50, Mk=1.5, Fk=exp(trgl_ld$polFkm),
#'             FSel=Sel, blicc_ld=trgl_ld)
#' plot(y=S$N_L[[1]], x=trgl_ld$LMP, type="l")
#' 
Rpop_F <- function(Galpha, Gbeta, Mk, Fk, FSel, blicc_ld) {
  FFSel <- list()
  N_L <- list()
  Zki <- rep(list(Mk * blicc_ld$M_L), blicc_ld$NT)
  for (qi in seq(blicc_ld$NQ)) {
    if (blicc_ld$Fkq[qi] > 0) {
      ti <- blicc_ld$Ti[qi]
      gi <- blicc_ld$Gi[qi]
      FFSel[[qi]] <- FSel[[gi]] * Fk[blicc_ld$Fkq[qi]]  # Fishing mortality
      Zki[[ti]] <- Zki[[ti]] + FFSel[[qi]]
    }
  }
  
  for (ti in seq(blicc_ld$NT))
    N_L[[ti]] <- with(blicc_ld,
                 Cpop_len(gl_nodes, gl_weights,
                          LLB, Zki[[ti]], Galpha, Gbeta) )
  return(list(N_L=N_L, Fki=FFSel))
}


#' Calculate the population size within length bins
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the relative
#' population number in each length bin taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. The model uses the Gauss-Laguerre quadrature
#' rule for integration, so nodes and weights for this must be provided. This
#' function is the same as [Rpop_F], but does not return fishing mortality and
#' is written in R.
#'
#' @export
#' @inheritParams Rpop_F
#' @param  node   Nodes for the Gauss-Laguerre quadrature rule
#' @param  wt     Weights for the Gauss-Laguerre quadrature rule
#' @param  Len    Vector of lower-bounds for length intervals
#' @param  Zki    Vector of total mortality (in time units of growth rate
#'   parameter K) to be applied in each length interval. This vector should be
#'   the same length as Len
#' @return A vector of population size in each length bin
#' @examples
#' glq <- statmod::gauss.quad(90, kind = "laguerre", alpha = 0.0)
#' S <- Rpop_len(glq$nodes, glq$weights, Len=15:55,
#'               Zki=c(rep(1.5, 10), rep(3, 31)), 100, 100/50)
#' plot(y=S, x=15:55, type="l")
#' 
Rpop_len <- function(node, wt, Len, Zki, Galpha, Gbeta)  {
  lgamma_Galpha <- lgamma(Galpha)
  nv <- length(node)
  LN <-  length(Len)
  surv <- double(LN)
  x_beta <- node / Gbeta
  log_x_beta <- log(x_beta)
  
  ss <- log(node + Gbeta * Len[1]) * (Galpha - 1.0) -
    Gbeta * Len[1] - lgamma_Galpha
  surv[1] <- sum(exp(ss) * wt)
  ss <-
    (log(x_beta + Len[2] - Len[1]) * (-Zki[1])) +
    log_x_beta * Zki[1] +
    log(node + Gbeta * Len[2]) * (Galpha - 1.0) -
    Gbeta * Len[2] - lgamma_Galpha
  surv[2] <- sum(exp(ss) * wt)
  
  for (Li in 3:LN) {
    Ln <- Len[Li]
    Lrange <- Ln - Len[1:(Li - 1)]
    Zii <- c(-Zki[1], Zki[1:(Li - 2)] - Zki[2:(Li - 1)])
    lim <- log_x_beta * Zki[Li - 1] + 
      vapply(x_beta, function(x) sum(log(x + Lrange) * Zii), numeric(1))
    ss <- lim + log(node + Gbeta * Ln) * (Galpha - 1.0) - Gbeta * Ln - lgamma_Galpha
    surv[Li] <- sum(exp(ss) * wt)
  }
  pop <- c(surv[-LN] - surv[-1], surv[LN]) / Zki
  return(pop)
}


#' Calculate the selectivities for each length bin for each gear
#'
#' The function calculates the selectivities for each gear based on the relevant
#' selectivity functions and their parameters.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param  Sm  Vector of parameters for all the selectivity functions
#'   including mixture weights at the end of the vector
#' @return Selectivities as a list of vectors one for each gear's length bin.
#' @examples
#' S <- Rselectivities(exp(gillnet_ld$polSm), gillnet_ld)
#' plot(y=S[[1]], x=gillnet_ld$LMP, type="l")
Rselectivities <- function(Sm, blicc_ld) {
  Ski <- list()
  GSki <- list()
  for (si in 1:blicc_ld$NS) {
    Indx <- with(blicc_ld, sp_i[si] : sp_e[si])
    #
    # Need to be changed to C functions or a single C function if possible
    #
    Ski[[si]] <- with(blicc_ld,
                      switch(fSel[si],
                             Rsel_logistic(Sm[Indx], LMP),
                             Rsel_normal(Sm[Indx], LMP),
                             Rsel_ssnormal(Sm[Indx], LMP),
                             Rsel_dsnormal(Sm[Indx], LMP)))
  }
  
  for (gi in seq(blicc_ld$NG)) {
    GSki[[gi]] <- Ski[[blicc_ld$GSbase[gi]]]
    gii <- 1L + (gi-1L)*2L
    if (blicc_ld$GSmix1[gii] > 0) {
      for (si in blicc_ld$GSmix1[gii]:blicc_ld$GSmix1[gii+1L])
        GSki[[gi]] <- with(blicc_ld, GSki[[gi]] + Ski[[GSmix2[si]]] * Sm[NP+si])
    }
  }
  return(GSki)
}

#' Calculate the spawning potential with no fishing
#'
#' The model calculates the spawning potential by multiplying the proportion of
#' numbers in each length bin with no fishing by the mature biomass per recruit.
#'
#' @export
#' @inheritParams Rpop_F
#' @return The spawning potential (a double) when F=0.
#' @examples
#' RSPR_0(Galpha=100, Gbeta=100/50, Mk=1.5, blicc_ld=gillnet_ld)
#'
RSPR_0 <- function(Galpha, Gbeta, Mk, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L
  pop <- with(blicc_ld,
              Rpop_len(gl_nodes, gl_weights,
                       LLB, Zki, Galpha, Gbeta))
  return(sum(pop * blicc_ld$ma_L))
}

