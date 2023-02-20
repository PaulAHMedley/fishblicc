 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 //
 //###########    Bayes Length Interval Catch Curve Stock Assessment #############
 //###########    Multiple Selectivity Fit                           #############
 //###########    PAUL MEDLEY                                        #############
 //###########    paulahmedley@gmail.com                             #############
 //###########    January 2023                                        #############
 //><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 // ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>

// Fits a length catch curve with with constant mortality within specified length interval (i.e. length bins).
// Fits to multiple length frequency data

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

// MODEL 1
vector sel_logistic(vector LMP, vector par) {
  // Logistic selectivity model
  // par[1] = Sel50, par[2] = Ss1
  int nl = rows(LMP);
  vector[nl] Sel = inv_logit(par[2] * (LMP - par[1]));
  return Sel;
}  //sel_logistic


// MODEL 2
vector sel_ssnormal(vector LMP, vector par) {
  // Single sided normal selectivity model (i.e. flat topped)
  // par[1] = Smax, par[2] = Ss1
  int nl = rows(LMP);
  vector[nl] Sel;
  for (i in 1:nl) {
    if (LMP[i] < par[1]) {
      Sel[i] = exp(-par[2] * (LMP[i] - par[1])^2);
    } else {
      Sel[i] = 1.0;
    }
  }
  return Sel;
} //sel_ssnormal


// MODEL 3
vector sel_dsnormal(vector LMP, vector par) {
  // Double sided normal selectivity model
  // par[1] = Smax, par[2] = Ss1, par[3] = Ss2
  int nl = rows(LMP);
  vector[nl] Sel;
  for (i in 1:nl) {
    if (LMP[i] < par[1]) {
      Sel[i] = exp(-par[2]*(LMP[i]-par[1])^2);
    } else {
      Sel[i] = exp(-par[3]*(LMP[i]-par[1])^2);
    }
  }
  return Sel;
} //sel_dsnormal



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
  return surv;
} //Survival_Est


vector NinInterval(vector Surv, vector Zki) {
  // Expected numbers of fish within each length interval
  // Due to numerical error, it is possible zero or, mre likely,  a very small negative number will be returned,
  // which will be rejected by the sampler.
  // Whether this is a problem in practice is not clear. Fixing it by setting minimum values might upset the sampler which uses gradients
  // and assumes continuous variables. Keeping the length bins tight around observations (e.g. bracket the observations with a single zero).
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
  int<lower=5>            NK;           // Number of Gauss-Laguerre nodes
  int<lower=1>            NG;           // Number of gears: separate length frequencies
  int<lower=1, upper=NG>  NF;           // Number of F's: gears associated with non-zero catches
  int<lower=1>            BN;           // Number of length bins
  int<lower=2>            NP;           // Number of selectivity parameters
  int<lower=2>            Pmx;          // Maximum parameters for available functions  (for matrix dimensions, currently = 3)

  int<lower=1, upper=3>         fSel[NG];       // Selectivity function to use: 1 = logistic, 2 = ss_normal, 3 = ds_normal
  int<lower=0, upper=3>         Spar[NG, Pmx];  // Selectivity function parameter index for each gear. Maximum number of parameters = 3
  vector[NF]                    Catch;          // Estimated total relative catch in numbers of fish, excluding zeroes for surveys etc.
  int<lower=1, upper=NG>        Fkg[NG];        // Index of the Fk associated with the gear gi. 0 implies catch negligible
  vector[BN]                    Len;            // Bin lower bounds
  int                           fq[NG, BN];     // Frequency in each bin for each gear
  real b;
  real Ls;
  real Lm;
  //Constant Hyperparameters for priors - Linf has a normal prior; all other priors are log normal
  //see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  real       poLinfm;
  real       poLinfs;
  real       polGam;
  real       polGas;
  real       polMkm;
  real       polMks;
  vector[NF] polFkm;
  real       polFks;
  vector[NP] polSm;
  vector[NP] polSs;
  real       polNB_phim;
  real       polNB_phis;
  real       polCs;          // Catch sigma for lognormal
}


transformed data {
  real          eps;     // small number to avoid zeroes and assoicated numerical failure. It should have no influence on results.
  vector[NK]    kn_wt[2] = gauss_laguerre_quad(NK, 0);  // Gauss-Laguerre quadrature nodes and weights
  vector[NG]    NObs;
  vector[NG]    olC;     // Normalized observed log catch by gear
  vector[BN]    LMP;     // Length bin mid point

  for (gi in 1:NG) {
    NObs[gi] = sum(fq[gi]);
  }

  {
    if (NF==1) {
      olC[1] = 0;
    } else {  // Identify gears contributing to catch
      for (gi in 1:NF) {
        olC[gi] = log(Catch[gi]);
        }
      olC -= log(sum(Catch));
    }
  }
  eps = 0.001/sum(NObs);

  // LMP is only used for selectivity models
  for (i in 1:(BN-1))
    LMP[i] = 0.5*(Len[i] + Len[i+1]);
  LMP[BN] = Len[BN] + 0.5*(Len[BN] - Len[BN-1]);
}


//###########################################
//####  PARAMETERS  #########################
//###########################################


parameters {  //modelled param
  real<lower = -poLinfm/poLinfs>  nLinf;
  real                            nGalpha;
  real                            nMk;
  vector[NF]                      nFk;
  vector[NP]                      nSel;   // 1-NG are Ss1, thereafter Ss2
  real                            nNB_phi;
}

transformed parameters {
  real Linf = poLinfm + nLinf*poLinfs;
  real Galpha = exp(polGam + nGalpha*polGas);
  real Mk = exp(polMkm + nMk*polMks);
  vector[NF] Fk = exp(polFkm + nFk*polFks);
  vector[NP] Sel = exp(polSm + nSel .* polSs);
  real NB_phi = exp(polNB_phim + nNB_phi*polNB_phis);
  real Gbeta = Galpha / Linf;
}


//###########################################
//####    MODEL     #########################
//###########################################

model {
  vector[BN] efq;
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  //###   PRIORS   ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>

  target += normal_lpdf(nLinf   | 0, 1);
  target += normal_lpdf(nGalpha | 0, 1);
  target += normal_lpdf(nMk     | 0, 1);
  target += normal_lpdf(nFk     | 0, 1);
  target += normal_lpdf(nSel    | 0, 1);
  target += normal_lpdf(nNB_phi | 0, 1);

  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  //### Expecteds ###  <><  ><>  <><  ><>  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
  {
    //calculate expected mortality for current parameter set
    real       Total_Catch = 0;
    vector[NF] elC;
    real       eC_sum;
    vector[BN] eC;
    vector[BN] Fki[NG];                     // fishing mortality
    vector[BN] Zki = rep_vector(Mk, BN);    // Total mortality
    vector[BN] Sv;                          // Survival
    vector[BN] Pop;                         // Population

    for (gi in 1:NG) {
      if (fSel[gi] == 1)
        Fki[gi] = sel_logistic(LMP, Sel[Spar[NG, 1:2]]);          // Flat top selectivity
      else if (fSel[gi] == 2)
        Fki[gi] = sel_ssnormal(LMP, Sel[Spar[NG, 1:2]]);          // Flat top selectivity
      else
        Fki[gi] = sel_dsnormal(LMP, Sel[Spar[NG, 1:3]]);          // Domed selectivity
      if (Fkg[gi] > 0) {
        Fki[gi] *= Fk[Fkg[gi]];                   // Fishing mortality
        Zki += Fki[gi];                           // Total mortality
      }
    }
    //calculate the expected survival at each length point integrating over age
    Sv = Survival_Est(kn_wt[1], kn_wt[2], Len, Zki, Galpha, Gbeta);
    Pop = NinInterval(Sv, Zki);

    //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>
    //### LIKELIHOOD ###  <><  ><>  <><  ><>  <><  ><>
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
  //### Spawning Potential Ratio (SPR) ###  <><  ><>
  //<><  ><>  <><  ><>  <><  ><>  <><  ><>  <><  ><>

  {
    real SPR0;
    real SPRF;
    row_vector[BN] MB = (exp(b*log(LMP)) .* inv_logit(Ls*(LMP - Lm)))';    // Mature biomass
    vector[BN] Zki = rep_vector(Mk, BN);           // Mortality with no fishing
    vector[BN] Sv;            // Survival

    // spawning biomass for the unexploited stock
    Sv = Survival_Est(kn_wt[1], kn_wt[2], Len, Zki, Galpha, Gbeta);
    SPR0 = MB * NinInterval(Sv, Zki);

    // Add fishing mortality rate to each length bin
    for (gi in 1:NG) {
      if (Fkg[gi] > 0) {
        if (fSel[gi] == 1)
          Zki += sel_logistic(LMP, Sel[Spar[NG, 1:2]])*Fk[Fkg[gi]];          // Flat top selectivity
        else if (fSel[gi] == 2)
          Zki += sel_ssnormal(LMP, Sel[Spar[NG, 1:2]])*Fk[Fkg[gi]];          // Flat top selectivity
        else
          Zki += sel_dsnormal(LMP, Sel[Spar[NG, 1:3]])*Fk[Fkg[gi]];          // Domed selectivity
      }
    }
    //calculate the expected survival at each length point integrating over age
    Sv = Survival_Est(kn_wt[1], kn_wt[2], Len, Zki, Galpha, Gbeta);
    SPRF = MB * NinInterval(Sv, Zki);   // Spawning biomass calculation

    //SPR is the ratio of the spawning biomass with fishing to spawning biomass without fishing
    SPR = SPRF/SPR0;
  }
}

