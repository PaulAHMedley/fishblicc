# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# Per Recruit Reference point functions
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


# Although per recruit calculations could be done using age-based functions as is usually done, the per-recruit 
# calculations are carried out using the length-based model. This ensures consistency between the model used
# to estimate the parameters and that used to calculate the reference points.


# Solve For Reference Points ----------------------------------------------


#' Find the standard yield-per-recruit fishing mortality reference point F0.1
#' 
#' Estimates an F0.1 reference point consistent with parameters and YPR slope 
#' 10% of initial slope. The function solves for F0.1 using a simple 
#' linear approximation to slope at the origin and at each point on the curve.
#' 
#' @inheritParams blicc_get_expected
#' @param  glq  a list of `nodes` and `weights` for the Gauss-Laguerre rule (see `statmod::gauss.quad`) 
#' @return F0.1 per-recruit fishing mortality reference point
#' 
F01_solve <- function(Smx, Linf, Galpha, Mk, Ss1, Ss2, blicc_ld, glq) {
  maxval <- 30
  fYPRF01 <- function(mF) {
    #Crude approx: slope=(y2-y1)/(x2-x1). constant (x2-x1) so ignore
    slope <-
      (fYPR(mF + 0.01, Galpha, Gbeta, Mk, Sel, Len, weight, glq) - 
         fYPR(mF, Galpha, Gbeta, Mk, Sel, Len, weight, glq))
    return(slope - slope01)
  }
  
  Gbeta <- Galpha/Linf
  Len <- blicc_ld$Len
  weight <- weight_at_length(blicc_ld)
  Sel <- Rsel_dsnormal(blicc_ld$LMP, Smx, Ss1, Ss2)
  slope01 <- 0.1*( fYPR(0.01, Galpha, Gbeta, Mk, Sel, Len, weight, glq) - 
               fYPR(0, Galpha, Gbeta, Mk, Sel, Len, weight, glq) ) 
  if (fYPRF01(maxval) > 0)
    return(NA)  
  
  return(uniroot(f = fYPRF01,
                 interval = c(0, maxval))$root)
}


#' Solve for a fishing mortality which produces the target SPR 
#' 
#' The function finds the fishing mortality that gives the `tarSPR` spawning potential ratio (SPR).
#' 
#' calculates spawning biomass per recruit (Spawning Potential Ratio / SPR). 
#' This should generally work for sensible reference point target. However, note that fishing mortality estimate 
#' may not be finite, dependent on the selectivity.
#' 
#' @inheritParams F01_solve
#' @param  tarSPR target SPR reference point: usually 0.2, 0.3 or 0.4
#' @return Spawning potential ratio per-recruit fishing mortality consistent with parameters and reference point
#' 
FSPR_solve <- function(Smx, Linf, Galpha, Mk, Ss1, Ss2, tarSPR, blicc_ld, glq) {
  maxval <- 30
  SRP_eval <- function(mF) {
    return(fSPR(mF, Galpha, Gbeta, Mk, Sel, Len, mb, glq) / SPR0 - tarSPR) 
  }
    
  Gbeta <- Galpha/Linf
  Len <- blicc_ld$Len
  Sel <- Rsel_dsnormal(blicc_ld$LMP, Smx, Ss1, Ss2)
  mb <- mature_biomass_at_length(blicc_ld)
  SPR0 <- fSPR(0, Galpha, Gbeta, Mk, Sel, Len, mb, glq)   #Unexploited SPR
  max_SRP <- SRP_eval(maxval)
  
  if (is.na(max_SRP) | max_SRP > 0)
    return(NA)
  
  return(uniroot(
    f = SRP_eval,
    interval = c(0, maxval)
  )$root)
} 


#' Solve for a selectivity mode which produces the target SPR 
#' 
#' The function finds the selectivity mode (Smx) that gives the `tarSPR` spawning potential ratio (SPR).
#' 
#' Calculates spawning biomass per recruit. 
#' This should often work for sensible reference point target. However, note 
#' that selectivity may not achieve any particular reference point if the fishing mortality is
#' too low. In these cases, NA is returned.
#' 
#' @inheritParams blicc_get_expected
#' @param  tarSPR target SPR reference point: usually 0.2, 0.3 or 0.4
#' @param  glq  a list of `nodes` and `weights` for the Gauss-Laguerre rule (see `statmod::gauss.quad`) 
#' @return Selectivity mode (full selectivity) achieving the target SPR. NA indicates it cannot be achieved.
#'
SSPR_solve <- function(Fk, Linf, Galpha, Mk, Ss1, Ss2, tarSPR, blicc_ld, glq) {

  SRP_eval <- function(Smx) {
    return(fSPR2(Smx, Fk, Galpha, Gbeta, Mk, Ss1, Ss2, Len, LMP, mb, glq) / SPR0 - tarSPR) 
  }
  Gbeta <- Galpha/Linf
  Len <- blicc_ld$Len
  LMP <- blicc_ld$LMP
  mb <- mature_biomass_at_length(blicc_ld)
  SPR0 <- fSPR(0, Galpha, Gbeta, Mk, RSel=rep(0, length(Len)), Len, mb, glq)   #Unexploited SPR
    
  # First need to bracket tarSPR
  V2 <- SRP_eval(Linf)   
  if (V2 < 0) {
    S1 <- Linf
    V1 <- V2
    S2 <- Linf*(1 + 5/sqrt(Galpha))  # Should be a length that fish do not grow to
    V2 <- SRP_eval(S2)
    if (V2 < 0) return(NA)
  } else {   # V2 is positive
    S1 <- Linf - 0.5
    V1 <- SRP_eval(S1)
    while (V1 > 0 & S1 > 0) {
      V1 <- SRP_eval(S1)
      S1 <- S1 - 0.5
    }
    if (V1 > 0) {
      return(NA)
    }
    S1 <- S1 + 0.5
    S2 <- S1 + 0.5
  }
  return(uniroot(
    f = SRP_eval,
    interval = c(S1, S2)
  )$root)
} 


#' Solve for a selectivity mode which produces maximum yield 
#' 
#' The function finds the selectivity mode (Smx) that maximises the yield per
#' recruit (YPR), keeping all other parameters at their current value.
#' 
#' Calculates yield per recruit. 
#' This should generally work because there must be a maximum yield between
#' the extreme lengths (0 and Linf). 
#' 
#' @inheritParams SSPR_solve
#' @return maximum yield point for the selectivity mode
#' 
SMY_solve <- function(Fk, Linf, Galpha, Mk, Ss1, Ss2, blicc_ld, glq) {
  Gbeta <- Galpha/Linf
  Len <- blicc_ld$Len
  LMP <- blicc_ld$LMP
  weight <- weight_at_length(blicc_ld)
  res <- optimize(fYPR2, lower=0, upper=Linf, maximum=TRUE, tol=1.0e-5,
                  Fk=Fk, Galpha=Galpha, Gbeta=Gbeta, Mk=Mk, Ss1=Ss1, Ss2=Ss2, Len=Len, LMP=LMP, weight=weight, glq=glq)
  return(res$maximum)
} 


#' Yield per recruit with variable fishing mortality
#' 
#' Used to evaluate the YPR for different values of Fk. 
#' 
#' @inheritParams F01_solve
#' @param RSel   Vector of selectivity from function `Rsel_dsnormal` or another source
#' @param Len    Vector of lower length boundaries for the length bins
#' @param weight Vector of weights from `weight_at_length` function
#' @return yield-per-recruit 
#' 
fYPR <- function(Fk, Galpha, Gbeta, Mk, RSel, Len, weight, glq) {
  LN <- length(Len)
  Fki <- RSel * Fk
  Zki <- Fki + Mk
  Rsurv <- RSurvival_Est(glq$nodes, glq$weights, Len, Zki, Galpha, Gbeta)
  N_L <- c(Rsurv[1:(LN - 1)] - Rsurv[2:LN], Rsurv[LN]) / Zki
  efq <- N_L * Fki                     # Catch
  Yield <- sum(efq * weight)
  return(Yield)
}


#' Yield per recruit with variable selectivity
#' 
#' Used to evaluate the YPR for different values of Fk and Smx. 
#' 
#' @inheritParams fYPR
#' @param  Smx      Mode of the normal selectivity function (full selectivity)
#' @param  Ss1      Left side slope, parameterized as 1/\sigma^2
#' @param  Ss2      Right side slope, parameterized as 1/\sigma^2. Zero implies flat-topped selectivity
#' @param  LMP    Vector of length bin mid-points to calculate selectivity
#' @return yield-per-recruit 
#' 
fYPR2 <- function(Smx, Fk, Galpha, Gbeta, Mk, Ss1, Ss2, Len, LMP, weight, glq) {
  LN <- length(Len)
  RSel <- Rsel_dsnormal(LMP, Smx, Ss1, Ss2)
  Fki <- RSel * Fk
  Zki <- Fki + Mk
  Rsurv <- RSurvival_Est(glq$nodes, glq$weights, Len, Zki, Galpha, Gbeta)
  N_L <- c(Rsurv[1:(LN - 1)] - Rsurv[2:LN], Rsurv[LN]) / Zki
  efq <- N_L * Fki                     # Catch
  Yield <- sum(efq * weight)
  return(Yield)
}


#' Spawning biomass per recruit with variable fishing mortality
#' 
#' Used to evaluate the SPR for different values of Fk. 
#' 
#' @inheritParams fYPR
#' @param mb   Vector of mature biomass from `mature_biomass_at_length` function
#' @return spawning biomass per recruit
#' 
fSPR <- function(Fk, Galpha, Gbeta, Mk, RSel, Len, mb, glq) {
  # Spawning potential ratio
  LN <- length(Len)
  Zki <- RSel * Fk + Mk
  Rsurv <- RSurvival_Est(glq$nodes, glq$weights, Len, Zki, Galpha, Gbeta)
  N_L <- c(Rsurv[1:(LN - 1)] - Rsurv[2:LN], Rsurv[LN]) / Zki
  return(sum(N_L * mb))
}


#' Spawning biomass per recruit with variable selectivity
#' 
#' Used to evaluate the SPR for different values of Fk and Smx. 
#' 
#' @inheritParams fYPR2
#' @param mb   Vector of mature biomass from `mature_biomass_at_length` function
#' @return spawning biomass per recruit
#' 
fSPR2 <- function(Smx, Fk, Galpha, Gbeta, Mk, Ss1, Ss2, Len, LMP, mb, glq) {
  # Spawning potential ratio
  LN <- length(Len)
  RSel <- Rsel_dsnormal(LMP, Smx, Ss1, Ss2)
  Zki <- RSel * Fk + Mk
  Rsurv <- RSurvival_Est(glq$nodes, glq$weights, Len, Zki, Galpha, Gbeta)
  N_L <- c(Rsurv[1:(LN - 1)] - Rsurv[2:LN], Rsurv[LN]) / Zki
  return(sum(N_L * mb))
}


