 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 //
 /////////////    Bayes Length Interval Catch Curve Stock Assessment /////////////
 /////////////    Multiple Selectivity Fit                           /////////////
 /////////////    PAUL MEDLEY                                        /////////////
 /////////////    paulahmedley@gmail.com                             /////////////
 /////////////    May 2023                                        /////////////
 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>

 // Fits a length catch curve with with constant mortality within specified length interval (i.e. length bins).
 // Fits to multiple length frequency data

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
  // copy correct scale parameter
  while (LMP[nv] < par[1])
    nv += 1;
  if (nv >= nl)
    Sel = exp(-par[2] * (LMP - par[1])^2);
  else {
    Sel[1:nv] = exp(-par[2] * (LMP[1:nv] - par[1])^2);
    nv += 1;
    Sel[nv:nl] = rep_vector(1.0, nl-nv+1);
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

// MODEL 5

vector sel_dsnmix(vector LMP, vector par) {
  // 1 Location 1
  // 2 Location 2
  // 3 Height 1
  // 4 Slope 1
  // 5 Slope 2
  // 6 Slope 3
  // 7 Slope 4 (could be set to zero)
  int nl = rows(LMP);
  vector[nl] Sel = rep_vector(0, nl);
  vector[nl] Adiff2 = (LMP - par[1])^2;
  vector[nl] Bdiff2 = (LMP - par[2])^2;
  vector[nl] Apars;  //scale parameter
  vector[nl] Bpars;  //scale parameter
  real       lwt;

  if (par[1] < par[2])
    lwt = -log(1 + par[3]*exp(-((par[2]-par[1])^2)*par[6]));
  else
    lwt = -log(1 + par[3]*exp(-((par[2]-par[1])^2)*par[7]));

  for (i in 1:nl) {
    if (LMP[i] < par[1]) Apars[i] = par[4]; else Apars[i] = par[5];
    if (LMP[i] < par[2]) Bpars[i] = par[6]; else Bpars[i] = par[7];
  }
  Sel = exp(-Apars .* Adiff2 + lwt) + par[3]*exp(-Bpars .* Bdiff2 + lwt);
  return Sel;
}


// Population model: survival

vector Pop_L(vector gl_node, vector gl_wt, vector Len, vector Zki, real alpha, real beta) {
  // Expected numbers of fish within each length interval
  // survival function based on x to be integrated 0 - Inf using Gauss-Laguerre quadrature
  // (exp(-x) removed for Laguerre quadrature)
  // Due to numerical error, it is possible zero or, more likely, a very small negative number will be returned,
  // which will be rejected by the sampler.
  // Whether this is a problem in practice is not clear. Fixing it by setting minimum values might upset the sampler which uses gradients
  // and assumes continuous variables. Keep the length bins tight around observations (e.g. bracket the observations with a single zero).
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
  int<lower=1, upper=NG>  NF;           // Number of F's: gears associated with non-zero catches
  int<lower=1>            NB;           // Number of length bins

  vector[NB]              LLB;            // Bin lower bounds
  int                     fq[NG, NB];     // Frequency in each bin for each gear
  vector[NF]              prop_catch;     // Estimated total relative catch in numbers of fish, excluding zeroes for surveys etc.
  int<lower=0, upper=NG>  Fkg[NG];        // Index of the Fk associated with the gear gi. 0 implies catch negligible
  int<lower=1, upper=5>   fSel[NG];       // Selectivity function to use: 1 = logistic, 2 = normal, 3 = ss_normal, 4 = ds_normal
  int<lower=2>            NP;             // Number of selectivity parameters
  int<lower=1, upper=NP>  sp_i[NG];       // Selectivity function parameter start index
  row_vector[NB]          ma_L;           // Mature biomass at length
  //Constant Hyperparameters for priors - Linf has a normal prior; all other priors are log normal
  //see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  real         poLinfm;
  real         poLinfs;
  real         polGam;
  real         polGas;
  vector[NB]   M_L;          // Length-based adjustment for natural mortality
  real         polMkm;
  real         polMks;
  vector[NF]   polFkm;
  real         polFks;
  vector[NP]   polSm;
  vector[NP]   polSs;
  real         polNB_phim;
  real         polNB_phis;
  real         polCs;       // Catch sigma for lognormal
  int<lower=5> NK;          // Number of Gauss-Laguerre nodes
  vector[NK]   gl_nodes;
  vector[NK]   gl_weights;
}


transformed data {
  real          eps;     // small number to avoid zeroes and assoicated numerical failure. It should have no influence on results.
  vector[NG]    NObs;
  vector[NF]    olC;     // Normalized observed log catch by gear
  vector[NB]    LMP;     // Length bin mid point

  for (gi in 1:NG) {
    NObs[gi] = sum(fq[gi]);
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

  // LMP is only used for selectivity models
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
  vector[NP]                      nSm;   // 1-NG are Ss1, thereafter Ss2
  real                            nNB_phi;
}

transformed parameters {
  real Linf = poLinfm + nLinf*poLinfs;
  real Galpha = exp(polGam + nGalpha*polGas);
  real Mk = exp(polMkm + nMk*polMks);
  vector[NF] Fk = exp(polFkm + nFk*polFks);
  vector[NP] Sm = exp(polSm + nSm .* polSs);
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

  target += normal_lpdf(nLinf   | 0, 1);
  target += normal_lpdf(nGalpha | 0, 1);
  target += normal_lpdf(nMk     | 0, 1);
  target += normal_lpdf(nFk     | 0, 1);
  target += normal_lpdf(nSm     | 0, 1);
  target += normal_lpdf(nNB_phi | 0, 1);

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
    vector[NB] Fki[NG];                     // fishing mortality
    vector[NB] Zki = Mk * M_L;              // Total mortality
    vector[NB] Pop;                         // Population

    for (gi in 1:NG) {
      if (fSel[gi] == 1)
        Fki[gi] = sel_logistic(LMP, segment(Sm, sp_i[gi], 2));  // Logistic selectivity
      else if (fSel[gi] == 2)
        Fki[gi] = sel_normal(LMP, segment(Sm, sp_i[gi], 2));    // Normal selectivity
      else if (fSel[gi] == 3)
        Fki[gi] = sel_ssnormal(LMP, segment(Sm, sp_i[gi], 2));  // Normal flat selectivity
      else if (fSel[gi] == 4)
        Fki[gi] = sel_dsnormal(LMP, segment(Sm, sp_i[gi], 3));  // Domed selectivity
      else
        Fki[gi] = sel_dsnmix(LMP, segment(Sm, sp_i[gi], 7));    // Mixture selectivity

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
      target += neg_binomial_2_lpmf(fq[gi] | efq, NB_phi);
      if (Fkg[gi] > 0) {
        Total_Catch += eC_sum;
        elC[Fkg[gi]] = log(eC_sum);
      }
    }
    elC -= log(Total_Catch);         // Normalise
    target += normal_lpdf(olC | elC, polCs);
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
    vector[NB] Zki = Mk * M_L;             // Mortality with no fishing
    vector[NB] Sv;                         // Survival

    // spawning biomass for the unexploited stock
    SPR0 = ma_L * Pop_L(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta);

    // Add fishing mortality rate to each length bin
    for (gi in 1:NG) {
      if (Fkg[gi] > 0) {
        if (fSel[gi] == 1)
          Zki += sel_logistic(LMP, segment(Sm, sp_i[gi], 2))*Fk[Fkg[gi]];  // Flat top selectivity
        else if (fSel[gi] == 2)
          Zki += sel_normal(LMP, segment(Sm, sp_i[gi], 2))*Fk[Fkg[gi]];    // Normal selectivity
        else if (fSel[gi] == 3)
          Zki += sel_ssnormal(LMP, segment(Sm, sp_i[gi], 2))*Fk[Fkg[gi]];  // Flat top selectivity
        else if (fSel[gi] == 4)
          Zki += sel_dsnormal(LMP, segment(Sm, sp_i[gi], 3))*Fk[Fkg[gi]];  // Domed selectivity
        else
          Zki += sel_dsnmix(LMP, segment(Sm, sp_i[gi], 7))*Fk[Fkg[gi]];    // Mixed selectivity
      }
    }
    //calculate the expected survival at each length point integrating over age
    SPRF = ma_L * Pop_L(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta);   // Spawning biomass calculation

    //SPR is the ratio of the spawning biomass with fishing to spawning biomass without fishing
    SPR = SPRF/SPR0;
  }
}

