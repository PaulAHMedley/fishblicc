 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 //
 /////////////    Bayes Length Interval Catch Curve Stock Assessment /////////////
 /////////////    Multiple Selectivity Fit                           /////////////
 /////////////    PAUL MEDLEY                                        /////////////
 /////////////    paulahmedley@gmail.com                             /////////////
 /////////////    July 2023                                        /////////////
 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>

 // Fits a length catch curve with with constant mortality within specified
 // length intervals (i.e. length bins).
 // Fits to multiple length frequency data
 // Selectivity models are uncoupled from gears to allow mixtures

 /////////////////////////////////////////////
 //////  FUNCTIONS   /////////////////////////
 /////////////////////////////////////////////


functions {

// Selectivity models

// MODEL 1
vector sel_logistic(vector LMP, vector par) {
  // Logistic selectivity model
  // par[1] = Sel50, par[2] = Ss1
  int nl = rows(LMP);
  vector[nl] Sel = inv_logit(par[2] * (LMP - par[1]));
  return Sel;
}  //sel_logistic

// MODEL 2
vector sel_normal(vector LMP, vector par) {
  // Normal selectivity model (i.e. same steepness on both sides of the mode)
  // par[1] = Smax, par[2] = Ss1
  int nl = rows(LMP);
  vector[nl] Sel = exp(-par[2] * (LMP - par[1])^2);
  return Sel;
} //sel_normal

// MODEL 3
vector sel_ssnormal(vector LMP, vector par) {
  // Single sided normal selectivity model (i.e. flat topped)
  // par[1] = Smax, par[2] = Ss1
  int nl = rows(LMP);
  int nv = 1;
  vector[nl] Sel;
  while (LMP[nv] < par[1])   //find bin above mode
    nv += 1;
  if (nv >= nl)
    Sel = exp(-par[2] * (LMP - par[1])^2);
  else {
    Sel = rep_vector(1.0, nl);
    if (nv > 1) {
      nv -= 1;
      Sel[1:nv] = exp(-par[2] * (LMP[1:nv] - par[1])^2);
    }
  }
  return Sel;
} //sel_ssnormal

// MODEL 4
vector sel_dsnormal(vector LMP, vector par) {
  // Double sided normal selectivity model
  // par[1] = Smax, par[2] = Ss1, par[3] = Ss2
  int nl = rows(LMP);
  vector[nl] Sel;
  vector[nl] pars;
  // copy correct scale parameter
  for (i in 1:nl)
    if (LMP[i] < par[1]) pars[i] = par[2]; else pars[i] = par[3];
  Sel = exp(-pars .* (LMP - par[1])^2);
  return Sel;
} //sel_dsnormal


// Population model: survival

vector Pop_L(vector gl_node, vector gl_wt, vector Len, vector Zki, real alpha, real beta) {
  // Expected numbers of fish within each length interval
  // survival function based on x to be integrated 0 - Inf using Gauss-Laguerre quadrature
  // (exp(-x) removed for Laguerre quadrature)
  int            nl = rows(Len);
  int            nv = rows(gl_node);
  real           lg_alpha = lgamma(alpha);
  vector[nv]     x_beta = gl_node / beta;
  vector[nv]     log_x_beta = log(x_beta);
  vector[nl-1]   Zin = append_row(-Zki[1], Zki[1:(nl-2)] - Zki[2:(nl-1)]);
  vector[nl]     surv;
  vector[nl]     pop;
  row_vector[nv] ss;

  ss = (log(gl_node + beta*Len[1]) * (alpha-1.0) - beta*Len[1] - lg_alpha)';
  surv[1] = exp(ss) * gl_wt;
  ss = (-log(x_beta + Len[2] - Len[1]) * Zki[1] + log_x_beta*Zki[1] +
       log(gl_node + beta*Len[2])*(alpha-1.0) - beta*Len[2] - lg_alpha)';
  surv[2] = exp(ss) * gl_wt;
  // Len >2
  for (n in 3:nl) {
    vector[n-1] Lrange = Len[1:(n-1)];
    vector[n-1] Zii = Zin[1:(n-1)];
    row_vector[nv] v2 = (log_x_beta*Zki[n-1] + log(gl_node + beta*Len[n])*(alpha-1.0))';
    real v3 = beta*Len[n] + lg_alpha;
    row_vector[nv] v1;
    for (i in 1:nv)
      v1[i] = sum(log(x_beta[i] + Len[n] - Lrange) .* Zii);
    ss =  v1 + v2 - v3;
    surv[n] = exp(ss) * gl_wt;
  }
  pop = append_row((surv[1:(nl-1)]-surv[2:nl]), surv[nl]) ./ Zki; // assumes no survival to the nl+1 interval
  return pop;
} //Pop_L


}  //FUNCTIONS



/////////////////////////////////////////////
//////  DATA   /////////////////////////
/////////////////////////////////////////////


data {
  //Likelihood are strictly applied to the bins supplied.
  //It is important to include all zero observations with non-negligible probability
  int<lower=1>            NG;           // Number of gears: separate length frequencies
  int<lower=1>            NS;           // Number of selectivity models
  int<lower=0, upper=NG>  NF;           // Number of F's: gears associated with non-zero catches
  int<lower=1>            NB;           // Number of length bins
  int<lower=2>            NP;           // Number of selectivity parameters
  int<lower=0>            NM;           // Number of mixture functions, & therefore mixture weights. Can be zero.

  vector[NB]              LLB;            // Bin lower bounds
  int                     fq[NG, NB];     // Frequency in each bin for each gear
  vector[NF]              prop_catch;     // Estimated total relative catch in numbers of fish, excluding zeroes for surveys etc.
  int<lower=0, upper=NG>  Fkg[NG];        // Index of the Fk associated with the gear gi. 0 implies catch negligible
  int<lower=1, upper=4>   fSel[NS];       // Selectivity function to use: 1 = logistic, 2 = normal, 3 = ss_normal, 4 = ds_normal
  int<lower=1, upper=NP>  sp_i[NS];       // Selectivity function parameter start index
  int<lower=1, upper=NS>  GSbase[NG];     // Integers linking gear to a selectivity function reference in fSel
  int<lower=0, upper=NM>  GSmix1[NG*2];   // Integer pairs linking gear to selectivity function reference in GSMix2
  int<lower=1, upper=NS>  GSmix2[NM];     // Integer references to fSel for mixtures. Could be zero length.
  row_vector[NB]          ma_L;           // Mature biomass at length
  //Constant Hyperparameters for priors - Linf has a normal prior; all other priors are log normal
  //see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  real            poLinfm;
  real            poLinfs;
  real            polGam;
  real            polGas;
  vector[NB]      M_L;          // Length-based adjustment for natural mortality
  real            polMkm;
  real            polMks;
  vector[NF]      polFkm;
  real            polFks;
  vector[NP+NM]   polSm;
  vector[NP+NM]   polSs;
  real            polNB_phim;
  real            polNB_phis;
  real            polCs;       // Catch sigma for lognormal
  int<lower=5>    NK;          // Number of Gauss-Laguerre nodes
  vector[NK]      gl_nodes;
  vector[NK]      gl_weights;
}


transformed data {
  real          eps;         // small number to avoid zeroes and assoicated numerical failure. It should have no influence on results.
  vector[NG]    NObs;
  vector[NF]    olC;         // Normalized observed log catch by gear
  vector[NB]    LMP;         // Length bin mid point
  array[NG] int GSmix0;  // Whether the gear selectivity includes mixtures > 0

  for (gi in 1:NG)
    NObs[gi] = sum(fq[gi]);

  {  //Check mixtures
    int Nmix = 0;
    if (NM == 0) {
      GSmix0 = rep_array(0, NG);
    } else {
      for (gi in 1:NG) {
        int si = 1+2*(gi-1);
        GSmix0[gi] = GSmix1[si];
        if (GSmix1[si] > 0)
          Nmix += GSmix1[si+1] - GSmix1[si] + 1;
      }
    }
    if (Nmix != NM)
      reject("Data object error: Number of mixtures / mixture references incorrect.");
  }

  {
    if (NF==1) {
      olC[1] = 0;
    } else {  // Identify gears contributing to catch
      for (gi in 1:NF) {
        olC[gi] = log(prop_catch[gi]);
        }
      olC -= log(sum(prop_catch));  //ensures catches proportional
    }
  }
  eps = 0.001/sum(NObs);

  // Length mid-point for bins: only used for selectivity models
  for (i in 1:(NB-1))
    LMP[i] = 0.5*(LLB[i] + LLB[i+1]);
  LMP[NB] = LLB[NB] + 0.5*(LLB[NB] - LLB[NB-1]);
}


/////////////////////////////////////////////
//////  PARAMETERS  /////////////////////////
/////////////////////////////////////////////


parameters {  //modelled param
  real<lower = -poLinfm/poLinfs>  nLinf;
  real                            nGalpha;
  real                            nMk;
  vector[NF]                      nFk;
  vector[NP+NM]                   nSm; //Includes selectivity mixture weights. NM may be zero.
  real                            nNB_phi;
}

transformed parameters {
  real Linf = poLinfm + nLinf*poLinfs;
  real Galpha = exp(polGam + nGalpha*polGas);
  real Mk = exp(polMkm + nMk*polMks);
  vector[NF] Fk = exp(polFkm + nFk * polFks);
  vector[NP+NM] Sm = exp(polSm + nSm .* polSs);
  real NB_phi = exp(polNB_phim + nNB_phi*polNB_phis);
  real Gbeta = Galpha / Linf;
}


/////////////////////////////////////////////
//////    MODEL     /////////////////////////
/////////////////////////////////////////////

model {
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  /////   PRIORS   ///  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>

  target += std_normal_lupdf(nLinf);
  target += std_normal_lupdf(nGalpha);
  target += std_normal_lupdf(nMk);
  target += std_normal_lupdf(nFk);
  target += std_normal_lupdf(nSm);
  target += std_normal_lupdf(nNB_phi);

  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  ///// Expecteds ///  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  {
    //calculate expected mortality for current parameter set
    vector[NB] efq;
    real       Total_Catch = 0;
    vector[NF] elC;
    real       eC_sum;
    vector[NB] eC;
    vector[NS] Seli[NS];                    // Selectivity models
    vector[NB] Fki[NG];                     // fishing mortality
    vector[NB] Zki = Mk * M_L;              // Total mortality
    vector[NB] Pop;                         // Population

    for (si in 1:NS) {
      if (fSel[si] == 1)
        Seli[si] = sel_logistic(LMP, segment(Sm, sp_i[si], 2));  // Logistic selectivity
      else if (fSel[si] == 2)
        Seli[si] = sel_normal(LMP, segment(Sm, sp_i[si], 2));    // Normal selectivity
      else if (fSel[si] == 3)
        Seli[si] = sel_ssnormal(LMP, segment(Sm, sp_i[si], 2));  // Normal flat selectivity
      else if (fSel[si] == 4)
        Seli[si] = sel_dsnormal(LMP, segment(Sm, sp_i[si], 3));  // Domed selectivity
    }

    for (gi in 1:NG) {
      Fki[gi] = Seli[GSbase[gi]];
      if (GSmix0[gi] != 0) {  // Add other selectivity if they exist
        int si = 1 + (gi-1)*2;
        for (i in GSmix1[si]:GSmix1[si+1])
          Fki[gi] += Sm[NP+i] * Seli[GSmix2[i]];
      }

      if (Fkg[gi] > 0) {
        Fki[gi] *= Fk[Fkg[gi]];                   // Fishing mortality
        Zki += Fki[gi];                           // Total mortality
      }
    }


    //calculate the expected survival at each length point integrating over age
    Pop = Pop_L(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta);

    //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
    ///// LIKELIHOOD ///  <><  ><>  <><  ><>  <><  ><>
    //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>

    for (gi in 1:NG) {
      eC = Fki[gi] .* Pop;
      eC_sum = sum(eC);
      efq = eC * NObs[gi]/eC_sum + eps;    // Normalise and raise to the expected numbers in the sample
      target += neg_binomial_2_lupmf(fq[gi] | efq, NB_phi);
      if (Fkg[gi] > 0) {
        Total_Catch += eC_sum;
        elC[Fkg[gi]] = log(eC_sum);
      }
    }
    elC -= log(Total_Catch);         // Normalise
    target += normal_lupdf(olC | elC, polCs);
  }


} //model


generated quantities {
  real SPR;
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  ///// Spawning Potential Ratio (SPR) ///  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>

  {
    real SPR0;
    real SPRF;
    vector[NS] Seli[NS];                    // Selectivity models
    vector[NB] Zki = Mk * M_L;             // Mortality with no fishing
    vector[NB] Sv;                         // Survival

    // spawning biomass for the unexploited stock
    SPR0 = ma_L * Pop_L(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta);

    // Add fishing mortality rate to each length bin
    for (gi in 1:NS) {
      if (fSel[gi] == 1)
        Seli[gi] = sel_logistic(LMP, segment(Sm, sp_i[gi], 2));  // Logistic selectivity
      else if (fSel[gi] == 2)
        Seli[gi] = sel_normal(LMP, segment(Sm, sp_i[gi], 2));    // Normal selectivity
      else if (fSel[gi] == 3)
        Seli[gi] = sel_ssnormal(LMP, segment(Sm, sp_i[gi], 2));  // Normal flat selectivity
      else if (fSel[gi] == 4)
        Seli[gi] = sel_dsnormal(LMP, segment(Sm, sp_i[gi], 3));  // Domed selectivity
    }

    for (gi in 1:NG) {
      if (Fkg[gi] > 0) {
        Zki += Seli[GSbase[gi]] * Fk[Fkg[gi]];
        if (GSmix0[gi] != 0) {  // Add other selectivity if they exist
          int si = 1 + (gi-1)*2;
          for (i in GSmix1[si]:GSmix1[si+1])
            Zki += Sm[NP+i] * Seli[GSmix2[i]] * Fk[Fkg[gi]];
        }
      }
    }
    //calculate the expected survival at each length point integrating over age
    SPRF = ma_L * Pop_L(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta);   // Spawning biomass calculation

    //SPR is the ratio of the spawning biomass with fishing to spawning biomass without fishing
    SPR = SPRF/SPR0;
  }
}

