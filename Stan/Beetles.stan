//Assessing the impact of host plant on feeding - Model Script

data {
  int<lower=1> n;                                    //Trials
  int<lower=1> n_group;                              //Number of life stage and plant combinations
  int<lower=1> n_stage;                              //Number of life stages
  int<lower=1> n_plant;                              //Number of plant species
  array[n] int<lower=1, upper=n_group> group;        //Lifestage x Plant combos
  array[n_group] int<lower=1, upper=n_stage> stage;  //Life stage per group
  array[n_group] int<lower=1, upper=n_plant> group_plant; //Plant species per group
  array[n] int<lower=1> beetle_count;                //Beetles per trial
  vector<lower=0>[n] hours;                          //Trial length
  vector<lower=0>[n] area;                           //Leaf area consumed
  array[n] int<lower=0, upper=1> fed;                //Did they feed binomial, 1 = yes

  array[n_group] int<lower=1> ref_count;             //Reference # of beetles (to account for different lifestage/beetle # combos)
  vector<lower=0>[n_group] ref_hours;                //Reference # of hours (to account for different lifestage/hour combos)

  int<lower=1> n_mean_level;                         //# of Hyperparameters for mean (3)
  int<lower=1> n_scale_level;                        //# of Hyperparameters for scale (1)
  array[n_stage] int<lower=1, upper=n_mean_level>  mean_level;   //Hyperparameters for mean
  array[n_stage] int<lower=1, upper=n_scale_level> scale_level;  //Hyperparameters for scale

  int<lower=0, upper=1> plant_effect;                //1 = plant effect shared across life stages

  real<lower=0> prior_logit_p_mu_sd;                 //Priors passed from R script
  real<lower=0> prior_logit_p_sigma_sd;
  real          prior_log_mu_mu_loc;
  real<lower=0> prior_log_mu_mu_sd;
  real<lower=0> prior_log_mu_sigma_sd;
  real<lower=0> prior_phi_rate;
  real<lower=0> prior_plant_sigma_sd;                //Scale prior for the shared plant effect

  int<lower=0, upper=1> prior_only;                  //1 = sample the prior only, used for check script
}

transformed data {
  array[n_group] int g_mean  = mean_level[stage];
  array[n_group] int g_scale = scale_level[stage];
  int n_pl = plant_effect ? n_plant : 0;

  int n_fed = sum(fed);
  array[n_fed] int fed_idx = rep_array(0, n_fed);
  {
    int j = 0;
    for (i in 1:n) if (fed[i] == 1) { j += 1; fed_idx[j] = i; }
  }
  array[n_fed] int fed_group = group[fed_idx];
  array[n_fed] int fed_scale = g_scale[fed_group];
  vector[n_fed] fed_area     = area[fed_idx];
  vector[n_fed] fed_hours    = hours[fed_idx];
  vector[n_fed] fed_count    = to_vector(beetle_count[fed_idx]);
}

parameters {
  vector[n_mean_level]           logit_p_mu;         //Mean feeding probability (hurdle - logit)
  vector<lower=0>[n_scale_level] logit_p_sigma;      //Spread across combinations (hurdle - logit)
  vector[n_group]                logit_p_z;          //Non-centering offsets (hurdle - logit)
  vector[n_mean_level]           log_mu_mu;          //Mean consumption rate (gamma dist - log)
  vector<lower=0>[n_scale_level] log_mu_sigma;       //Spread across consumption rate (gamma dist - log)
  vector[n_group]                log_mu_z;           //Non-centering offset (gamma dist - log)
  vector<lower=0>[n_scale_level] phi;                //Shape per feeding beetle (gamma dist)

  vector[n_pl]                   logit_p_plant_z;    //Plant effect shared across stages (hurdle - logit)
  vector[n_pl]                   log_mu_plant_z;     //Plant effect shared across stages (gamma - log)
  vector<lower=0>[plant_effect ? 1 : 0] logit_p_plant_sigma;
  vector<lower=0>[plant_effect ? 1 : 0] log_mu_plant_sigma;
}

transformed parameters {
  vector[n_group] logit_p = logit_p_mu[g_mean] + logit_p_z .* logit_p_sigma[g_scale];
  vector[n_group] log_mu  = log_mu_mu[g_mean]  + log_mu_z  .* log_mu_sigma[g_scale];
  if (plant_effect) {
    logit_p += logit_p_plant_z[group_plant] * logit_p_plant_sigma[1];
    log_mu  += log_mu_plant_z[group_plant]  * log_mu_plant_sigma[1];
  }
}

model {
  logit_p_mu    ~ normal(0, prior_logit_p_mu_sd);
  logit_p_sigma ~ normal(0, prior_logit_p_sigma_sd);       //Half-normal via <lower=0>
  logit_p_z     ~ std_normal();
  log_mu_mu     ~ normal(prior_log_mu_mu_loc, prior_log_mu_mu_sd);
  log_mu_sigma  ~ normal(0, prior_log_mu_sigma_sd);        //Half-normal
  log_mu_z      ~ std_normal();
  phi           ~ exponential(prior_phi_rate);

  logit_p_plant_z     ~ std_normal();
  log_mu_plant_z      ~ std_normal();
  logit_p_plant_sigma ~ normal(0, prior_plant_sigma_sd);   //Half-normal
  log_mu_plant_sigma  ~ normal(0, prior_plant_sigma_sd);   //Half-normal

  if (prior_only == 0) {
    fed ~ bernoulli_logit(logit_p[group]);
    if (n_fed > 0)
      fed_area ~ gamma(phi[fed_scale] .* fed_count,
                       phi[fed_scale] ./ (exp(log_mu[fed_group]) .* fed_hours));
  }
}

generated quantities {
  vector[n_group] p           = inv_logit(logit_p);        //P(Feeding occurs)
  vector[n_group] mu          = exp(log_mu);               //Rate given feeding
  vector[n_group] rate_uncond = p .* mu;                   //Unconditional rate
  vector[n_group] y_ref_rate;                              //Predicted rate, new rep
  array[n] int fed_rep;                                    //Replicate data for the posterior predictive check
  vector[n] area_rep;                                      //Replicate data for the posterior predictive check
  vector[n] log_lik;                                       //For check, and LOO

  for (g in 1:n_group) {
    real phi_g = phi[g_scale[g]];
    y_ref_rate[g] = bernoulli_rng(p[g])
      ? gamma_rng(phi_g * ref_count[g], phi_g / (mu[g] * ref_hours[g]))
        / (ref_count[g] * ref_hours[g])
      : 0;
  }
  for (i in 1:n) {
    int g = group[i];
    real phi_i   = phi[g_scale[g]];
    real shape_i = phi_i * beetle_count[i];
    real rate_i  = phi_i / (mu[g] * hours[i]);

    fed_rep[i]  = bernoulli_logit_rng(logit_p[g]);
    area_rep[i] = fed_rep[i] ? gamma_rng(shape_i, rate_i) : 0;

    log_lik[i] = bernoulli_logit_lpmf(fed[i] | logit_p[g]);
    if (fed[i] == 1) log_lik[i] += gamma_lpdf(area[i] | shape_i, rate_i);
  }
}
