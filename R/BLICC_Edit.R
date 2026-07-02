# BLICC Model Data Edit Functions -----------------------------------------------------


#' Returns data list from [blicc_dat] with a subset of model components
#'
#' Subsets the data list with only selectivity components, gears, frequencies
#' and parameters relevant indicated are retained. There are 
#' no checks on consistency of the resulting data list.
#'
#' @NoRd
#' @param blicc_ld Data list to be edited
#' @param NewName New name for model
#' @param Nindx Index for populations to include
#' @param Qindx Index for frequencies to include
#' @param Tindx Index for time periods to include
#' @param Findx Index for fishing mortalities to include
#' @param Gindx Index for gears to include
#' @param Sindx Index for selectivities to include
#' @param GTGindx Index for growth groups to include
#' @return The subset data object from blicc_ld  
#' 
blicc_subset <- function(blicc_ld,
                         NewName,
                         Nindx,
                         Qindx,
                         Tindx,
                         Gindx,
                         Sindx,
                         GTGindx) {
  
  ld <- blicc_ld
  Pindx <- integer(0)
  for (si in Sindx) {
    Pindx <- with(blicc_ld, c(Pindx, sp_i[si]:sp_e[si]))
  }
  Findx <- Qindx & blicc_ld$Fkq > 0

  ld$model_name <- NewName
  ld$NQ <- sum(Qindx)
  ld$NG <- length(Gindx)
  ld$NS <- length(Sindx)
  ld$fSel <- as.array(blicc_ld$fSel[Sindx])
  
  ld$NN <- length(Nindx)
  ld$NT <- length(Tindx)
  ld$NP <- length(Pindx)
  ld$NX <- length(GTGindx)
  ld$NF <- sum(Findx)
  
  ld$fqname <- as.array(blicc_ld$fqname[Qindx])
  ld$gname <- as.array(blicc_ld$gname[Gindx])
  ld$tpname <- as.array(blicc_ld$tpname[Tindx])
  
  ld$fq <- blicc_ld$fq[Qindx]
  ld$Ni <- as.array(match(blicc_ld$Ni[Qindx], Nindx))
  ld$Gi <- as.array(match(blicc_ld$Gi[Qindx], Gindx))
  ld$Ti <- as.array(match(blicc_ld$Ti[Qindx], Tindx))
  ld$Xi <- as.array(match(blicc_ld$Xi[Nindx], GTGindx))
  ld$poLinfm <- as.array(blicc_ld$poLinfm[GTGindx])
  ld$poLinfs <- as.array(blicc_ld$poLinfs[GTGindx])
  ld$polMkm <- as.array(blicc_ld$polMkm[GTGindx])
  ld$polMks <- as.array(blicc_ld$polMks[GTGindx])
  ld$a <- as.array(blicc_ld$a[GTGindx])
  ld$b <- as.array(blicc_ld$b[GTGindx])
  ld$L50 <- as.array(blicc_ld$L50[GTGindx])
  ld$Ls <- as.array(blicc_ld$Ls[GTGindx])
  ld$wt_L <- blicc_ld$wt_L[, GTGindx, drop = FALSE]
  ld$ma_L <- blicc_ld$ma_L[GTGindx, , drop = FALSE]
  
  pc <- double(blicc_ld$NQ)
  pc[blicc_ld$Fkq > 0] <- blicc_ld$prop_catch
  ld$prop_catch <- pc[Qindx]
  if (!any(ld$prop_catch > 0))
    stop("Error: at least one population must have catches greater than zero.")
  ld$prop_catch <- as.array(ld$prop_catch[ld$prop_catch > 0])
  
  ld$Fkq <- match(blicc_ld$Fkq[Qindx], blicc_ld$Fkq[Findx])
  ld$Fkq[is.na(ld$Fkq)] <- 0
  sums <- with(ld, tapply(prop_catch[prop_catch > 0], Ni[Fkq > 0], sum))
  names(sums) <- NULL
  ld$prop_catch <- with(ld, prop_catch / sums[Ni[Fkq > 0]]) # Normalise
  
  #seq_along(Gindx)
  ld$GSbase <- as.array(match(blicc_ld$GSbase[Gindx], Sindx)) # resequence
  mxn <- integer(0)
  mxpar <- mxpars <- double(0)
  
  ld$GSmix1 <- integer(2 * length(Gindx))
  mix_1 <- 1L
  mix_2 <- 1L
  for (gi in Gindx) {
    mgi <- (gi - 1L) * 2L + 1L
    if (blicc_ld$GSmix1[mgi] > 0) {
      mxindx <- with(blicc_ld, GSmix2[GSmix1[gi]:GSmix1[gi + 1L]])
      ld$GSmix1[mix_1] <- mix_2
      mix_1 <- mix_1 + 1L
      mix_2 <- mix_2 + length(mxindx) - 1L
      ld$GSmix1[mix_1] <- mix_2
      mix_1 <- mix_1 + 1L
      mix_2 <- mix_2 + 1L
      mxn <- c(mxn, match(mxindx, Sindx))
      mxpar <- with(blicc_ld, c(mxpar, polSm[NP + GSmix1[mgi]:GSmix1[mgi + 1L]]))
      mxpars <- with(blicc_ld, c(mxpars, polSs[NP + GSmix1[mgi]:GSmix1[mgi + 1L]]))
    } else {
      mix_1 <- mix_1 + 2L
    }
  }
  ld$GSmix2 <- mxn
  ld$NM <- length(mxn)
  ld$polFkm <- as.array(blicc_ld$polFkm[blicc_ld$Fkq[Findx]])
  
  npar <- Rsel_functions()$npar[ld$fSel]
  spare <- spari <- integer(ld$NS)
  np <- 1L
  for (i in seq_along(ld$fSel)) {
    spari[i] <- np
    np <- np + npar[i]
    spare[i] <- np - 1L
  }
  # Replace function list
  ld$sp_i <- as.array(spari)    #start
  ld$sp_e <- as.array(spare)    #end
  
  ld$polSm <- c(blicc_ld$polSm[Pindx], mxpar)
  ld$polSs <- c(blicc_ld$polSs[Pindx], mxpars)
  return(ld)
}


#' Returns data list from [blicc_dat] with only the selected populations
#' included
#'
#' Subsets the data list with only selectivity components, gears, frequencies
#' and parameters relevant to the defined populations.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param population  The names or an integer vector of the populations (time 
#'   period/growth groups) being selected.
#' @return The data object blicc_ld subset for the new populations 
#'   (growth groups / time period).
#' 
blicc_population_filter <- function(blicc_ld, population) {
  Nindx <- parse_population(population, blicc_ld)
  
  if (blicc_ld$NN==1L) return(blicc_ld)
  
  Qindx <- blicc_ld$Ni %in% Nindx
  Tindx <- unique(blicc_ld$Ti[Qindx])
  Gindx <- unique(blicc_ld$Gi[Qindx])
  Sindx <- get_selectivities(Gindx, blicc_ld)
  GTGindx <- unique(blicc_ld$Xi[Nindx])
  
  new_name <- paste(blicc_ld$model_name, "Subset Populations", 
                    paste(as.character(population), collapse=" "))
  
  ld <- blicc_subset(blicc_ld, new_name,
                     Nindx, Qindx, Tindx,
                     Gindx, Sindx, GTGindx)
  return(ld)    
}




#' Returns data list from [blicc_dat] with the populations/gears with zero 
#' fishing mortality removed
#'
#' Subsets the data list with only selectivity components, gears, frequencies
#' and parameters relevant to the frequencies associated with F > 0. 
#'
#' @export
#' @inheritParams blicc_mpd
#' @return The data object blicc_ld subset for F>0  
#' 
blicc_zeroF_filter <- function(blicc_ld) {

  Qindx <- blicc_ld$Fkq > 0
  if (blicc_ld$NQ==sum(Qindx)) return(blicc_ld) #nothing to remove
  
  Nindx <- unique(blicc_ld$Ni[Qindx])
  Tindx <- unique(blicc_ld$Ti[Qindx])
  Gindx <- unique(blicc_ld$Gi[Qindx])
  Sindx <- get_selectivities(Gindx, blicc_ld)
  GTGindx <- unique(blicc_ld$Xi[Nindx])
  new_name <- paste(blicc_ld$model_name, "Fishery Only", collapse=" ")

  ld <- blicc_subset(blicc_ld, new_name,
                     Nindx, Qindx, Tindx,
                     Gindx, Sindx, GTGindx)
  return(ld)    
}



#' Combines bins across a range into a single length frequency bin
#' 
#' The frequencies are summed over an interval range which is converted 
#' into a single bin. This may be useful for smoothing peak outliers across bins. 
#' The resulting frequency bins will have variable width.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param interval_range  The range of bins identified by their LLB that are to 
#'   be combined (lower bin and upper bin).
#' @return The same data object blicc_ld but with the bin adjustment.
#' 
blicc_combine_bins <- function(blicc_ld, interval_range) {
  collapse <- match(interval_range, blicc_ld$LLB)
  if (any(is.na(collapse))) stop("Error: Specified bin interval not found.")
  if (length(collapse) != 2 | collapse[1] >= collapse[2]) 
    stop("interval_range must be a vector matching lower length boundaries for 2 intervals.")
  
  cl_seq <- seq(collapse[1]+1L, collapse[2], by=1L)
  blicc_ld$LLB <- blicc_ld$LLB[-cl_seq]  
  blicc_ld$NB <- length(blicc_ld$LLB)
  blicc_ld$LMP <- with(blicc_ld, c((LLB[-NB] + LLB[-1]) * 0.5, 
                                   LLB[NB] + 0.5*(LLB[NB]-LLB[NB-1])))
  
  blicc_ld$fq <- lapply(blicc_ld$fq, 
                        FUN = \(x) { x[collapse[1]] <- x[collapse[1]]+sum(x[cl_seq])
                        return(x[-cl_seq])  })

  if (!any(is.na(blicc_ld$a) | is.na(blicc_ld$b) | is.na(blicc_ld$L50) | is.na(blicc_ld$Ls))) {
    L95 <- as.vector(blicc_ld$L50 - log(1 / 0.95 - 1) / blicc_ld$Ls)
    blicc_ld <- blip_LH(blicc_ld, a = as.vector(blicc_ld$a), b = as.vector(blicc_ld$b), 
                        L50 = as.vector(blicc_ld$L50), L95 = L95)    
  } else {
    warning("Weight/maturity at length bins removed.")
    blicc_ld$wt_L <- blicc_ld$wt_L[-cl_seq,, drop=FALSE]
    blicc_ld$ma_L <- blicc_ld$ma_L[-cl_seq,, drop=FALSE]
  }
  
  if (blicc_ld$ref_length > 0)
    blicc_ld$M_L <- with(blicc_ld, ref_length/LMP)
  else 
    blicc_ld$M_L <- rep(1, blicc_ld$NB)  # Fixed natural mortality
  
  return(blicc_ld)
}


