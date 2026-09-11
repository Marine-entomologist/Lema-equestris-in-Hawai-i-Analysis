//Assessing the impact of host plant on development - Model Script
//Contains two separate models (survival and length) as one call


data {
  int<lower=1> n_plant;                                  //Number of plants
  int<lower=1, upper=n_plant> control_plant;             //Which plant is the control
  int<lower=1> n_rep;                                    //Number of replicate dishes
  int<lower=1> n_outcome;                                //Number of survival steps (2: to 24 h, then to trial end)
  array[n_rep] int<lower=1, upper=n_plant> rep_plant;    //Grouping or "Rep | Plant"
  array[n_rep] int<lower=1> n_placed;                    //Larvae/rep
  array[n_outcome, n_rep] int<lower=0> successes;        //[1, ] survivors at 24 h; [2, ] survivors at trial end

  int<lower=1> n_larva;                                  //Number of larvae measured/plant
  array[n_larva] int<lower=1, upper=n_plant> larva_plant;//Grouping or "Length-measured larvae | Plant"
  vector<lower=0>[n_larva] length_mm;                    //Length of larvae

  real<lower=0> prior_bin_mu_sd;                         //Priors passed from R script
  real<lower=0> prior_bin_sigma_sd;
  real          prior_len_mu_loc;
  real<lower=0> prior_len_mu_sd;
  real<lower=0> prior_len_sigma_sd;
  real<lower=0> prior_len_resid_sd;

  int<lower=0, upper=1> overdispersed;                   //1 = beta-binomial dish-level overdispersion
  real<lower=0> prior_rho_a;                             //Beta prior on the within-dish correlation
  real<lower=0> prior_rho_b;

  int<lower=0, upper=1> prior_only;                      //1 = sample the prior only, used for check script
}

transformed data {
  array[n_outcome, n_rep] int n_at_risk;
  for (r in 1:n_rep) {
    n_at_risk[1, r] = n_placed[r];
    for (k in 2:n_outcome) n_at_risk[k, r] = successes[k - 1, r];
  }
}

parameters {
  vector[n_outcome] bin_mu;                              //Mean survival at each step (logit)
  vector<lower=0>[n_outcome] bin_sigma;                  //Survival spread across plants (logit)
  matrix[n_outcome, n_plant] bin_z;                      //Non-centred offset for survival (logit)
  array[overdispersed ? n_outcome : 0] real<lower=0, upper=1> rho; //Within-dish correlation
  real len_mu;                                           //Mean length (log)
  real<lower=0> len_sigma;                               //Length spread across plants (log)
  vector[n_plant] len_z;                                 //Non-centred offset for length (log)
  real<lower=0> len_resid;                               //Larval sd (log)
}

transformed parameters {
  matrix[n_outcome, n_plant] logit_prob;                 //Conditional survival (logit scale)
  vector[n_plant] log_len = len_mu + len_z * len_sigma;
  for (k in 1:n_outcome) logit_prob[k] = bin_mu[k] + bin_z[k] * bin_sigma[k];
}

model {
  bin_mu           ~ normal(0, prior_bin_mu_sd);
  bin_sigma        ~ normal(0, prior_bin_sigma_sd);            //Half-normal via <lower=0>
  to_vector(bin_z) ~ std_normal();
  len_mu           ~ normal(prior_len_mu_loc, prior_len_mu_sd);
  len_sigma        ~ normal(0, prior_len_sigma_sd);            //Half-normal
  len_z            ~ std_normal();
  len_resid        ~ normal(0, prior_len_resid_sd);            //Half-normal
  if (overdispersed) rho ~ beta(prior_rho_a, prior_rho_b);

  if (prior_only == 0) {
    for (k in 1:n_outcome) {                                   //Staged survival model
      if (overdispersed) {
        real kappa = (1 - rho[k]) / rho[k];                    //Beta-binomial precision
        for (r in 1:n_rep) {
          real q = inv_logit(logit_prob[k, rep_plant[r]]);
          if (n_at_risk[k, r] > 0)
            successes[k, r] ~ beta_binomial(n_at_risk[k, r], q * kappa, (1 - q) * kappa);
        }
      } else {
        successes[k] ~ binomial_logit(n_at_risk[k], to_vector(logit_prob[k])[rep_plant]);
      }
    }
    length_mm ~ lognormal(log_len[larva_plant], len_resid);     //Length model
  }
}

generated quantities {
  matrix[n_outcome, n_plant] prob_step;                  //Survival at each step, conditional on reaching it
  matrix[n_outcome, n_plant] prob;                       //Cumulative probability of surviving to each time point
  matrix[n_outcome, n_plant] prob_diff_vs_control;       //Cumulative survival - control cumulative survival
  vector[n_plant] length_mean = exp(log_len);            //Mean estimated length as mm (drops log scale)
  vector[n_plant] length_pct_of_control;                 //Length difference (estimated length - Control length) as a %
  array[n_outcome, n_rep] int successes_rep;             //Replicate data for the posterior predictive check
  vector[n_larva] length_rep;                            //Replicate data for the posterior predictive check
  matrix[n_outcome, n_rep] log_lik_surv;                 //Per-dish, per-step log likelihood

  for (k in 1:n_outcome)
    for (g in 1:n_plant)
      prob_step[k, g] = inv_logit(logit_prob[k, g]);

  for (g in 1:n_plant) {
    prob[1, g] = prob_step[1, g];
    for (k in 2:n_outcome) prob[k, g] = prob[k - 1, g] * prob_step[k, g];
  }
  for (k in 1:n_outcome) {
    real control_k = prob[k, control_plant];
    for (g in 1:n_plant) prob_diff_vs_control[k, g] = prob[k, g] - control_k;
  }

  {array[n_outcome, n_rep] int at_risk_rep;
    for (r in 1:n_rep) at_risk_rep[1, r] = n_placed[r];
    for (k in 1:n_outcome) {
      real kappa = 0;
      if (overdispersed) kappa = (1 - rho[k]) / rho[k];
      for (r in 1:n_rep) {
        real q = prob_step[k, rep_plant[r]];
        int  m = at_risk_rep[k, r];
        if (m == 0) {
          successes_rep[k, r] = 0;
        } else if (overdispersed) {
          successes_rep[k, r] = beta_binomial_rng(m, q * kappa, (1 - q) * kappa);
        } else {
          successes_rep[k, r] = binomial_rng(m, q);
        }
        if (n_at_risk[k, r] == 0) {
          log_lik_surv[k, r] = 0;
        } else if (overdispersed) {
          log_lik_surv[k, r] = beta_binomial_lpmf(successes[k, r] | n_at_risk[k, r],
                                                  q * kappa, (1 - q) * kappa);
        } else {
          log_lik_surv[k, r] = binomial_lpmf(successes[k, r] | n_at_risk[k, r], q);
        }
        if (k < n_outcome) at_risk_rep[k + 1, r] = successes_rep[k, r];
      }
    }
  }

  for (g in 1:n_plant)
    length_pct_of_control[g] = 100 * exp(log_len[g] - log_len[control_plant]);
  for (i in 1:n_larva)
    length_rep[i] = lognormal_rng(log_len[larva_plant[i]], len_resid);
}
