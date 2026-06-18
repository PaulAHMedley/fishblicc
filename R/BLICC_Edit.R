# BLICC Model Data Edit Functions -----------------------------------------------------

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
  population <- parse_population(population, blicc_ld)
  
  if (blicc_ld$NN==1L) return(blicc_ld)
  
  Xindx <- blicc_ld$Xi[population]
  Qindx <- blicc_ld$Ni %in% population
  Tindx <- unique(blicc_ld$Ti[Qindx])
  Findx <- Qindx & blicc_ld$Fkq > 0
  Gindx <- unique(blicc_ld$Gi[Qindx])
  Sindx <- get_selectivities(Gindx, blicc_ld)
  Pindx <- integer(0)
  for (si in Sindx) {
    Pindx <- c(Pindx, blicc_ld$sp_i[si]:blicc_ld$sp_e[si])
  }
  GTGindx <- sort(unique(Xindx))
  
  ld <- blicc_ld
  ld$model_name <- paste(blicc_ld$model_name, "Subset Populations", 
                         paste(as.character(population), collapse=" "))
  ld$NQ <- sum(Qindx)
  ld$NG <- length(Gindx)
  ld$NS <- length(Sindx)
  ld$fSel <- as.array(blicc_ld$fSel[Sindx])

  ld$NN <- length(population)
  ld$NT <- length(Tindx)
  ld$NF <- sum(Findx)
  ld$NP <- length(Pindx)
  ld$NX <- length(GTGindx)
  
  ld$fqname <- as.array(blicc_ld$fqname[Qindx])
  ld$gname <- as.array(blicc_ld$gname[Gindx])
  ld$tpname <- as.array(blicc_ld$tpname[Tindx])
  
  ld$fq <- blicc_ld$fq[Qindx]
  ld$Ni <- as.array(match(blicc_ld$Ni[Qindx], population))
  ld$Gi <- as.array(match(blicc_ld$Gi[Qindx], Gindx))
  ld$Ti <- as.array(match(blicc_ld$Ti[Qindx], Tindx))
  ld$Xi <- as.array(match(blicc_ld$Xi[Qindx], GTGindx))
  ld$poLinfm <- blicc_ld$poLinfm[GTGindx]
  ld$poLinfs <- blicc_ld$poLinfs[GTGindx]
  ld$polMkm <- blicc_ld$polMkm[GTGindx]
  ld$polMks <- blicc_ld$polMks[GTGindx]
  ld$a <- blicc_ld$a[GTGindx]
  ld$b <- blicc_ld$b[GTGindx]
  ld$L50 <- blicc_ld$L50[GTGindx]
  ld$Ls <- blicc_ld$Ls[GTGindx]
  ld$wt_L <- blicc_ld$wt_L[, GTGindx, drop=FALSE]
  ld$ma_L <- blicc_ld$ma_L[GTGindx, , drop=FALSE]
  
  ld$prop_catch <- as.array(blicc_ld$prop_catch[Qindx])
  ld$Fkq <- as.array(match(blicc_ld$Fkq[Qindx], which(Findx)))
  ld$Fkq[is.na(ld$Fkq)] <- 0
  sums <- with(ld, tapply(prop_catch, Ni[Fkq>0], sum))
  ld$prop_catch <- with(ld, prop_catch / sums[Ni[Fkq>0]]) # Normalise
  
  #seq_along(Gindx)
  ld$GSbase <- as.array(blicc_ld$GSbase[Gindx]) # resequence
  mxn <- integer(0)
  mxpar <- mxpars <- double(0)
  
  ld$GSmix1 <- integer(2*length(Gindx))
  mix_1 <- 1L
  mix_2 <- 1L
  for (gi in Gindx) {
    mgi <- (gi-1L)*2L + 1L
    if (blicc_ld$GSmix1[mgi] > 0) {
      mxindx <- blicc_ld$GSmix2[blicc_ld$GSmix1[gi]:blicc_ld$GSmix1[gi+1L]]
      ld$GSmix1[mix_1] <- mix_2
      mix_1 <- mix_1 + 1L
      mix_2 <- mix_2 + length(mxindx) - 1L
      ld$GSmix1[mix_1] <- mix_2
      mix_1 <- mix_1 + 1L
      mix_2 <- mix_2 + 1L
      mxn <- c(mxn, match(mxindx, Sindx))
      mxpar <- c(mxpar, blicc_ld$polSm[blicc_ld$NP+blicc_ld$GSmix1[mgi]:blicc_ld$GSmix1[mgi+1L]])
      mxpars <- c(mxpars, blicc_ld$polSs[blicc_ld$NP+blicc_ld$GSmix1[mgi]:blicc_ld$GSmix1[mgi+1L]])
    } else {
      mix_1 <- mix_1 + 2L
    }
  }
  ld$GSmix2 <- mxn
  ld$NM <- length(mxn)
  ld$polFkm <- as.array(blicc_ld$polFkm[Findx])

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


