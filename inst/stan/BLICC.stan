 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 //  
 //###########    Bayes Length Interval Catch Curve Stock Assessment #############
 //###########    PAUL MEDLEY                                        #############
 //###########    paulahmedley@gmail.com                             #############
 //###########    August 2022                                        #############
 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>

// Fits a length catch curve with with constant mortality within specified length interval (i.e. length bins). 

//###########################################
//####  FUNCTIONS   #########################
//###########################################


functions {

  vector[] gauss_laguerre_quad(int nt, real alpha) {
    //  A Gauss-Laguerre rule for approximating
    //    Integral ( a <= x < +oo ) |x-a|^ALPHA exp(-B*x(x-a)) f(x) dx
    //  of order "nt".
    //  Order nt is the number of points in the rule:
    //  ALPHA is the exponent of |X|:
    //  Ported from C: Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    //                 C version by John Burkardt.https://people.math.sc.edu/Burkardt/c_src/laguerre_rule/laguerre_rule.c
    //                 distributed under the GNU LGPL License

    real r8_epsilon = 2.220446049250313E-016;
    int itn  = 30;
    vector[nt] aj;
    vector[nt] bj;
    vector[nt] wt = rep_vector(0, nt);
    real       zemu = tgamma(alpha+1.0);
    vector[nt] knot_wt[2];
    int  m;
    int  mml;
    int  k;
    real b;
    real c;
    real f;
    real g;
    real sign_g;
    real p;
    real r;
    real s;

    {
      vector[nt] i = cumulative_sum(rep_vector(1, nt));
      aj = 2.0*i - 1.0 + alpha;
      bj = sqrt(i .* (i+alpha));
    }
    wt[1] = sqrt(zemu);
    //  Diagonalize the Jacobi matrix.
    if (nt == 1) {
      knot_wt[1] = aj;
      knot_wt[2] = square(wt);
      return knot_wt;
    }
    
    bj[nt] = 0.0;
    
    for (ll in 1:nt) {
      int j = 0;
      int Do_Loop = 1;
      while (Do_Loop) {
        m = ll;
        while ((m<nt) && (fabs(bj[m]) > r8_epsilon * (fabs(aj[m]) + fabs(aj[m+1])))) { //set m
          m += 1;
        }
        if (m == ll) {  
          Do_Loop = 0;  //break; ends while(1) loop
        } else {
          p = aj[ll];
          if (itn <= j) {
            reject("gauss_laguerre_quad - Fatal error! Iteration limit exceeded");
          }
          j += 1;
          g = (aj[ll+1]-p) / (2.0*bj[ll]);
          if (g>=0) sign_g = 1.0; else sign_g = -1.0; 
          r =  sqrt(square(g) + 1.0);
          s = g + fabs(r)*sign_g;
          g = aj[m] - p + (bj[ll] / s);
          s = 1.0;
          c = 1.0;
          p = 0.0;
          mml = m - ll;
        
          for (ii in 1:mml) {
            int i = m - ii;
            f = s * bj[i];
            b = c * bj[i];
            if (fabs(g) <= fabs(f)) {
              c = g / f;
              r =  sqrt(square(c) + 1.0);
              bj[i+1] = f * r;
              s = 1.0 / r;
              c *= s;
            } else {
              s = f / g;
              r =  sqrt(square(s) + 1.0);
              bj[i+1] = g * r;
              c = 1.0 / r;
              s *= c;
            }
            g = aj[i+1] - p;
            r = (aj[i]-g)*s + 2.0 * c * b;
            p = s * r;
            aj[i+1] = g + p;
            g = c * r - b;
            f = wt[i+1];
            b = wt[i];
            wt[i+1] = s * b + c * f;
            wt[i] = c * b - s * f;
          }
          aj[ll] -= p;
          bj[ll] = g;
          bj[m] = 0.0;
        }
      }  //while
    }
    //  Sorting
    for (ii in 2:m)  {
      int i = ii - 1;
      k = i;
      p = aj[i];
      
      for (j in ii:nt) {
        if (aj[j] < p) {
          k = j;
          p = aj[j];
        }
      }
      
      if (k != i) {
        g = aj[i];
        aj[k] = g;
        aj[i] = p;
        p = wt[i];
        wt[i] = wt[k];
        wt[k] = p;
      }
    }
    knot_wt[1] = aj;
    knot_wt[2] = square(wt);
    return knot_wt;
  }  //gauss_laguerre_quad

// Selectivity models

vector sel_logistic(vector LMP, real Sel50, real Ss1) {
  // Logistic selectivity model
  int nl = rows(LMP);
  vector[nl] Sel = inv_logit(Ss1 * (LMP - Sel50));
  return Sel;
}  //sel_logistic


vector sel_dsnormal(vector LMP, real Smax, real Ss1, real Ss2) {
  // Double sided normal selectivity model
  int nl = rows(LMP);
  vector[nl] SL;
  for (i in 1:nl) {
    if (LMP[i] < Smax) {
      SL[i] = exp(-Ss1*(LMP[i]-Smax)^2);
    } else {
      SL[i] = exp(-Ss2*(LMP[i]-Smax)^2);
    }
  }
  return SL;
} //sel_dsnormal


vector sel_ssnormal(vector LMP, real Smax, real Ss1) {
  // Single sided normal selectivity model (i.e. flat topped)
  int nl = rows(LMP);
  vector[nl] SL;
  for (i in 1:nl) {
    if (LMP[i] < Smax) {
      SL[i] = exp(-Ss1*(LMP[i]-Smax)^2);
    } else {
      SL[i] = 1.0;
    }
  }
  return SL;
} //sel_ssnormal

// Population model: survival

vector Survival_Est(vector gl_node, vector gl_wt, vector Len, vector Zki, real alpha, real beta) {
  // survival function based on x to be integrated 0 - Inf using Gauss-Laguerre quadrature
  // (exp(-x) removed for Laguerre quadrature)
  int            nl = rows(Len);
  int            nv = rows(gl_node);
  real           lg_alpha = lgamma(alpha);
  vector[nv]     x_beta = gl_node / beta;
  vector[nv]     log_x_beta = log(x_beta);
  vector[nl-1]   Zin = append_row(-Zki[1], Zki[1:(nl-2)] - Zki[2:(nl-1)]); 
  vector[nl]     surv;
  row_vector[nv] ss; 

  ss = (log(gl_node + beta*Len[1]) * (alpha-1.0) - beta*Len[1] - lg_alpha)';
  surv[1] = exp(ss) * gl_wt;
  ss = (-log(x_beta + Len[2] - Len[1]) * Zki[1] + log_x_beta*Zki[1] + log(gl_node + beta*Len[2])*(alpha-1.0) - beta*Len[2] - lg_alpha)';
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
  return surv;
} //Survival_Est


vector NinInterval(vector Surv, vector Zki) {
  // Expected numbers of fish within each length interval
  int nl = rows(Surv);
  vector[nl] NinIntv = append_row((Surv[1:(nl-1)]-Surv[2:nl]), Surv[nl]) ./ Zki; // assumes no survival to the nl+1 interval
  return NinIntv;
}  // NinInterval



}  //FUNCTIONS



//###########################################
//####  DATA   #########################
//###########################################


data {
  //Likelihood are strictly applied to the bins supplied. 
  //It is important to include all zero observations with non-negligible probability
  int<lower=2>   NK;
  //Assume non-overlapping length bins
  int<lower=1>   oBN;          //Number of length bins

  vector[oBN]  Len;          //Bin lower bounds
  int          fq[oBN];      //Frequency in each bin
  real b;
  real Ls;
  real Lm;
  int  Flat_top;             // Whether selectivity is flattopped or not: 0 for flat topped selectivity, >0 implies double normal.
  //Constant Hyperparameters for priors - Linf has a normal prior; all other priors are log normal
  //see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  real poLinfm;
  real poLinfs;
  real polGam;
  real polGas;
  real polMkm;
  real polMks;
  real polFkm;
  real polFks;
  real polSmxm;
  real polSmxs;
  real polSs1m;
  real polSs1s;
  real polSs2m;
  real polSs2s;
  real polNB_phim; 
  real polNB_phis;
}

transformed data {
  real          eps = 1.0e-6;     // small number to avoid zeroes. Should be much less than 1.0. It should have no influence on results
  real          NObs = sum(fq);
  vector[NK]    kn_wt[2] = gauss_laguerre_quad(NK, 0);  //Prepare Gauss-Laguerre Quadrature
  vector[oBN]   LMP;
  int           N_Ss;
  if (Flat_top == 0)
    N_Ss = 1;
  else
    N_Ss = 2;
  
  // LMP is only used for selectivity model
  for (i in 1:(oBN-1)) 
    LMP[i] = 0.5*(Len[i]+Len[i+1]);
  LMP[oBN] = Len[oBN] + 0.5*(Len[oBN]-Len[oBN-1]);      // Final interval 
}


//###########################################
//####  PARAMETERS  #########################
//###########################################


parameters {  //modelled param
  real<lower = -poLinfm/poLinfs>  nLinf;  
  real                            nGalpha;
  real                            nMk;
  real                            nFk;
  real                            nSmx;
  vector[N_Ss]                    nSs;
  real                            nNB_phi;
}

transformed parameters {
  real Linf = poLinfm + nLinf*poLinfs;  
  real Galpha = exp(polGam + nGalpha*polGas);
  real Mk = exp(polMkm + nMk*polMks);
  real Fk = exp(polFkm + nFk*polFks);
  real Smx = exp(polSmxm + nSmx*polSmxs);
  real Ss1 = exp(polSs1m + nSs[1]*polSs1s);
  real Ss2;
  real NB_phi = exp(polNB_phim + nNB_phi*polNB_phis);
  real Gbeta = Galpha / Linf;
  
  if (Flat_top == 0)
    Ss2 = 0;
  else
    Ss2 = exp(polSs2m + nSs[2]*polSs2s);
}


//###########################################
//####    MODEL     #########################
//###########################################

model {
  vector[oBN] efq; 
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //###   PRIORS   ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  
  target += normal_lpdf(nLinf   | 0, 1);           
  target += normal_lpdf(nGalpha | 0, 1);
  target += normal_lpdf(nMk     | 0, 1);
  target += normal_lpdf(nFk     | 0, 1);
  target += normal_lpdf(nSmx    | 0, 1);
  target += normal_lpdf(nSs     | 0, 1);
  target += normal_lpdf(nNB_phi | 0, 1);

  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //### Expecteds ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  {
    //calculate expected mortality for current parameter set
    real        Sum_efq;
    vector[oBN] Fki;  // Fishing mortality with selectivity
    vector[oBN] Zki;                                 // Total mortality
    vector[oBN] Sv;
    if (Flat_top == 0)
      Fki = sel_ssnormal(LMP, Smx, Ss1) * Fk;          // Fishing mortality with selectivity
    else
      Fki = sel_dsnormal(LMP, Smx, Ss1, Ss2) * Fk;   // Fishing mortality with selectivity
    Zki = Fki+Mk;                                 // Total mortality
    //calculate the expected survival at each length point integrating over age
    Sv = Survival_Est(kn_wt[1], kn_wt[2], Len, Zki, Galpha, Gbeta);

    efq = Fki .* NinInterval(Sv, Zki);
    Sum_efq = sum(efq);
    efq *= NObs/Sum_efq;    // Normalise and raise to the expected numbers in the sample
    efq += eps;             // Add a small number to avoid an expected zero and increase numerical stability
    //efq = Fki .* append_row((Sv[1:(oBN-1)]-Sv[2:oBN]), Sv[oBN]) * NObs ./ Zki;
  }
  
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //### LIKELIHOOD ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  

  {
    target += neg_binomial_2_lpmf(fq | efq, NB_phi);
  }
} //model


generated quantities {
  real SPR;
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  //### Spawning Potential Ratio (SPR) ###  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>  
  
  {
    real SPR0;
    real SPRF;
    row_vector[oBN] MB = (exp(b*log(Len)) .* inv_logit(Ls*(Len - Lm)))';    // Mature biomass
    vector[oBN] Zki;           // Mortality with fishing
    vector[oBN] Sv;            // Survival
    // Calculte total mortality rate in each length bin
    if (Flat_top == 0)
      Zki = sel_ssnormal(LMP, Smx, Ss1) * Fk + Mk;          // Fishing mortality with selectivity
    else
      Zki = sel_dsnormal(LMP, Smx, Ss1, Ss2) * Fk + Mk;   // Fishing mortality with selectivity
    //calculate the expected survival at each length point integrating over age
    Sv = Survival_Est(kn_wt[1], kn_wt[2], Len, Zki, Galpha, Gbeta);
    SPRF = MB * NinInterval(Sv, Zki);   // Spawning biomass calculation
    
    // Same calculation of the spawning biomass but for the unexploited stock
    Zki = rep_vector(Mk, oBN);                                              // Only natural mortality for unexploited stock
    Sv = Survival_Est(kn_wt[1], kn_wt[2], Len, Zki, Galpha, Gbeta);
    SPR0 = MB * NinInterval(Sv, Zki);

    //SPR is the ratio of the spawning biomass with fishing to spawning biomass without fishing 
    SPR = SPRF/SPR0;
  }
}

