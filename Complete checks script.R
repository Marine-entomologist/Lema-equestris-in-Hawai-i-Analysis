#________________________####
#Assessing the impact of host plant on feeding - checks and balances ####
#________________________####

#1.0 Setup ####
#Note: This script REQUIRES the model to be run and loaded to work. It is intended as checks for the BAR guidelines

source(here::here("R", "Complete model script.R"))
library(bayesplot)
library(loo)
library(patchwork)
set.seed(SEED)

refit <- function(changes, ...) {
  model$sample(data = modifyList(stan_data, changes), chains = 4,
               parallel_chains = min(4, parallel::detectCores()),
               iter_warmup = 2000, iter_sampling = 2000,
               adapt_delta = 0.99, max_treedepth = 12, seed = SEED, refresh = 0, ...)
}

#1.1. Convergence (BARG 2B) ####
hyper <- c("logit_p_mu[1]", "logit_p_mu[2]", "logit_p_mu[3]", "logit_p_sigma[1]",
           "log_mu_mu[1]", "log_mu_mu[2]", "log_mu_mu[3]", "log_mu_sigma[1]", "phi[1]")
mcmc_rank_overlay(fit$draws(hyper))
ggsave(here("Plots", "feeding_check_rank.svg"), width = 24, height = 16, units = "cm")
p_pairs <- mcmc_pairs(fit$draws(hyper), np = nuts_params(fit)) #Too big to visualize in Rstudio - just look at the svg
ggsave(here("Plots", "feeding_check_pairs.svg"), p_pairs, width = 60, height = 30, units = "cm")

#1.2. Prior predictive check (BARG 1E) ####
#1.2.1 Checking priors ####
prior_fit <- refit(list(prior_only = 1L))
prior_fit %>%
  recover_types(select(beetles, group)) %>%
  spread_draws(y_ref_rate[group], p[group]) %>% ungroup() %>%
  pivot_longer(c(y_ref_rate, p)) %>%
  filter(name == "p" | value > 0) %>%
  mutate(value = if_else(name == "y_ref_rate", log10(value), value),
         name = if_else(name == "p", "P(feeding)",
                        "log10 predicted rate")) %>%
  ggplot(aes(value)) + geom_density() + facet_wrap(~ name, scales = "free") +
  theme_classic()
ggsave(here("Plots", "feeding_check_prior_predictive.svg"),
       width = 20, height = 9, units = "cm")

#1.2.2 Prior vs Posterior ####
all_pars <- c("logit_p_mu[1]", "logit_p_mu[2]", "logit_p_mu[3]",
              "log_mu_mu[1]",  "log_mu_mu[2]",  "log_mu_mu[3]",
              "logit_p_sigma[1]", "log_mu_sigma[1]", "phi[1]")
bind_rows(
  as_tibble(prior_fit$draws(all_pars, format = "df")) %>% mutate(source = "prior"),
  as_tibble(fit$draws(all_pars, format = "df")) %>% mutate(source = "posterior")) %>%
  pivot_longer(all_of(all_pars), names_to = "parameter") %>%
  ggplot(aes(value, fill = source)) +
  geom_density(alpha = 0.5, colour = NA) +
  facet_wrap(~ parameter, scales = "free") +
  scale_fill_manual(values = c(prior = "grey60", posterior = "firebrick")) +
  theme_classic() + theme(legend.title = element_blank())
ggsave(here("Plots", "feeding_check_prior_posterior.svg"),
       width = 22, height = 9, units = "cm")

#1.3. Posterior predictive check (BARG 3A) ####
#1.3.1 Hurdle Plot ####
rep_draws <- fit %>%
  spread_draws(fed_rep[trial], area_rep[trial], ndraws = 500) %>% ungroup() %>%
  left_join(beetles %>% mutate(trial = row_number()) %>%
              select(trial, group, Lifestage, Fed, AreaConsumed),
            by = "trial")

ppc_hurdle <- rep_draws %>%
  summarise(sim = mean(fed_rep), .by = c(group, .draw)) %>%
  summarise(lo = quantile(sim, 0.025), hi = quantile(sim, 0.975), .by = group) %>%
  left_join(summarise(beetles, obs = mean(Fed), .by = group), by = "group") %>%
  mutate(inside = obs >= lo & obs <= hi) %T>% print(n = Inf)
write_csv(ppc_hurdle, here("Output", "feeding_ppc_hurdle.csv"))

hurdle_plot <- ppc_hurdle %>%
  ggplot(aes(y = fct_rev(group))) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey60", linewidth = 3) +
  geom_point(aes(obs), colour = "firebrick", shape = 4, stroke = 1.2) +
  labs(x = "Proportion of trials with any feeding") +
  theme_classic() + theme(axis.title.y = element_blank())
hurdle_plot
ggsave(here("Plots", "feeding_check_ppc_hurdle.svg"),
       width = 16, height = 14, units = "cm")

#1.3.2 Gamma plot ####
gamma_plot <- rep_draws %>% filter(fed_rep == 1) %>%
  ggplot() +
  geom_density(aes(log10(area_rep), group = .draw), colour = "grey70", linewidth = 0.2) +
  geom_density(data = filter(beetles, Fed == 1), aes(log10(AreaConsumed)),
               colour = "firebrick", linewidth = 1) +
  facet_wrap(~ Lifestage, scales = "free_y", ncol = 1) +
  labs(x = "log10 leaf area consumed (mm^2)") + theme_classic()
gamma_plot
ggsave(here("Plots", "feeding_check_ppc_amount.svg"),
       width = 9, height = 20, units = "cm")

#1.3.3 Combo Plot ####
hurdle_plotc <- (hurdle_plot + labs(tag = "A") + theme(plot.tag = element_text(face = "bold")))
gamma_plotc <-  (gamma_plot + labs(tag = "B") + theme(plot.tag = element_text(face = "bold")))
Plot_checks_feeding <- wrap_elements(hurdle_plotc) | wrap_elements(gamma_plotc) +
  plot_layout(widths = c(2.5, 1))
Plot_checks_feeding
ggsave(here("Plots", "Plot_checks_feeding.svg"), Plot_checks_feeding, width = 30, height = 20, units = "cm")

#1.4. LOO ####
loo_fit <- fit$loo(cores = min(4, parallel::detectCores())) %T>% print()
tibble(row = seq_len(nrow(beetles)), group = beetles$group, fed = beetles$Fed,
       pareto_k = loo_fit$diagnostics$pareto_k) %>%
  filter(pareto_k > 0.7) %T>% print() %>%
  write_csv(here("Output", "feeding_loo.csv"))

#1.5. Prior sensitivity (BARG 5A, 5C) ####
watch <- c("Adult | S. americanum", "Adult | B. oleracea",
           "L1 | B. oleracea", "L1Naive | B. oleracea")
summarise_fit <- function(f, label) {
  f %>% recover_types(select(beetles, group)) %>%
    spread_draws(p[group], rate_uncond[group]) %>% ungroup() %>%
    filter(as.character(group) %in% watch) %>%
    summarise(across(c(p, rate_uncond), q_eti), .by = group) %>%
    mutate(setting = label,
           divergences = sum(f$diagnostic_summary()$num_divergent),
           rhat_max = max(f$summary()$rhat, na.rm = TRUE), .before = 1)
}

prior_grid <- list(
  "reported"                    = list(),
  "logit_p_mu sd 1"             = list(prior_logit_p_mu_sd = 1),
  "logit_p_mu sd 4"             = list(prior_logit_p_mu_sd = 4),
  "logit_p_sigma sd 0.5"        = list(prior_logit_p_sigma_sd = 0.5),
  "logit_p_sigma sd 2"          = list(prior_logit_p_sigma_sd = 2),
  "log_mu_mu centre 0.01"       = list(prior_log_mu_mu_loc = log(0.01)),
  "log_mu_mu centre 1"          = list(prior_log_mu_mu_loc = log(1)),
  "log_mu_mu sd 1"              = list(prior_log_mu_mu_sd = 1),
  "log_mu_mu sd 4"              = list(prior_log_mu_mu_sd = 4),
  "log_mu_sigma sd 0.5"         = list(prior_log_mu_sigma_sd = 0.5),
  "log_mu_sigma sd 2"           = list(prior_log_mu_sigma_sd = 2),
  "phi rate 0.5"                = list(prior_phi_rate = 0.5),
  "phi rate 2"                  = list(prior_phi_rate = 2)
)
sensitivity_priors <- imap(prior_grid, function(change, label)
  summarise_fit(refit(change), label)) %>% list_rbind() %T>% print(n = Inf)
write_csv(sensitivity_priors, here("Output", "feeding_sensitivity_priors.csv"))

#1.6. ROPE Sensitivity (BARG 4C, 5E) ####
sensitivity_rope <- map(c(0.5, 1, 2), function(k)
  draws %>%
    summarise(P_rope = mean(rate_uncond < k * rope), .by = c(Lifestage, Plant)) %>%
    mutate(rope_multiplier = k, .before = 1)) %>%
  list_rbind() %>%
  pivot_wider(names_from = rope_multiplier, values_from = P_rope,
              names_prefix = "P_rope_x") %T>% print(n = Inf)
write_csv(sensitivity_rope, here("Output", "feeding_sensitivity_rope.csv"))

#1.7. Structural sensitivity ####
structure_grid <- list(
  "reported: stage means, shared scale" = list(),
  "means pooled across stages"          = list(n_mean_level = 1L,
                                                mean_level = rep(1L, length(STAGES))),
  "scales and phi per stage"            = list(n_scale_level = length(STAGES),
                                                scale_level = seq_along(STAGES)),
  "plant effect shared across stages"   = list(plant_effect = 1L)
)
structure_fits <- imap(structure_grid, function(change, label) refit(change))
sensitivity_structure <- imap(structure_fits, summarise_fit) %>%
  list_rbind() %T>% print(n = Inf)
write_csv(sensitivity_structure, here("Output", "feeding_sensitivity_structure.csv"))

#1.7.1 What the shared plant effect does to combinations that never fed ####
crossed_vs_reported <- bind_rows(
  fit %>% recover_types(select(beetles, group)) %>%
    spread_draws(p[group], rate_uncond[group]) %>% ungroup() %>%
    mutate(structure = "reported"),
  structure_fits[["plant effect shared across stages"]] %>%
    recover_types(select(beetles, group)) %>%
    spread_draws(p[group], rate_uncond[group]) %>% ungroup() %>%
    mutate(structure = "plant effect")) %>%
  summarise(across(c(p, rate_uncond), q_eti), .by = c(group, structure)) %>%
  left_join(select(groups, group, Lifestage, Plant, n_trial, n_fed), by = "group") %>%
  arrange(Lifestage, Plant, structure) %T>% print(n = Inf)
write_csv(crossed_vs_reported, here("Output", "feeding_crossed_structure.csv"))

#________________________####
#Assessing the impact of host plant on development (larval growth and survival) - Checks ####
#________________________####
#2.0 Setup ####
#Note: This portion of the script utilizes some variables from the setup of 1.0 - run that section first
refit_dev <- function(changes, ...) {
  model_dev$sample(data = modifyList(stan_data_dev, changes), chains = 4,
               parallel_chains = min(4, parallel::detectCores()),
               iter_warmup = 2000, iter_sampling = 2000,
               adapt_delta = 0.99, max_treedepth = 12, seed = SEED, refresh = 0, ...)
}

#2.1. Convergence (BARG 2B) ####
hyper_dev <- c("bin_mu[1]", "bin_mu[2]", "bin_sigma[1]", "bin_sigma[2]",
           "len_mu", "len_sigma", "len_resid")
mcmc_rank_overlay(fit_dev$draws(hyper_dev))
ggsave(here("Plots", "development_check_rank.svg"), width = 24, height = 14, units = "cm")

p_pairs_dev <- mcmc_pairs(fit_dev$draws(hyper_dev), np = nuts_params(fit_dev))
p_pairs_dev
ggsave(here("Plots", "development_check_pairs.svg"), p_pairs_dev, width = 60, height = 30, units = "cm")

#2.2. Prior predictive check (BARG 1E) ####
#2.2.1 Checking priors ####
prior_fit_dev <- refit_dev(list(prior_only = 1L))
prior_fit_dev %>%
  spread_draws(prob[k, g], length_mean[g]) %>% ungroup() %>%
  filter(k == 2, g == stan_data_dev$control_plant) %>%
  pivot_longer(c(prob, length_mean)) %>%
  ggplot(aes(value)) + geom_density() +
  facet_wrap(~ name, scales = "free") + theme_classic()
ggsave(here("Plots", "development_check_prior_predictive.svg"),
       width = 20, height = 9, units = "cm")

#2.2.2 Prior vs Posterior ####
prior_posterior <- bind_rows(
  prior_fit_dev$draws(c("bin_sigma", "len_sigma", "len_resid"), format = "df") %>%
    as_tibble() %>% mutate(source = "prior"),
  fit_dev$draws(c("bin_sigma", "len_sigma", "len_resid"), format = "df") %>%
    as_tibble() %>% mutate(source = "posterior")) %>%
  pivot_longer(-c(source, starts_with(".")), names_to = "parameter") %>%
  summarise(across(value, q_eti), .by = c(source, parameter)) %T>% print(n = Inf)
write_csv(prior_posterior, here("Output", "development_prior_posterior.csv"))

#2.3. Posterior predictive check (BARG 3A) ####
#2.3.1 Survival Plot ####
ppc_survival <- fit_dev %>%
  spread_draws(successes_rep[k, r], ndraws = 500) %>% ungroup() %>%
  summarise(lo = quantile(successes_rep, 0.025),
            hi = quantile(successes_rep, 0.975),
            median = median(successes_rep), .by = c(k, r)) %>%
  mutate(Outcome = names(OUTCOMES)[k], Plant = survival$Plant[r],
         n_placed = survival$n_placed[r],
         observed = stan_data_dev$successes[cbind(k, r)]) %>%
  mutate(inside = observed >= lo & observed <= hi) %>%
  select(Outcome, Plant, r, n_placed, observed, lo, median, hi, inside) %T>%
  print(n = Inf)
write_csv(ppc_survival, here("Output", "development_ppc_survival.csv"))

ppcsc <- ppc_survival %>%
  mutate(dish = str_c(Plant, " #", r)) %>%
  ggplot(aes(y = fct_rev(dish))) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey60", linewidth = 3) +
  geom_point(aes(observed), colour = "firebrick", shape = 4, stroke = 1.2) +
  facet_wrap(~ Outcome, nrow = 1) +
  labs(x = "Larvae surviving per replicate dish") +
  theme_classic() + theme(axis.title.y = element_blank())
ppcsc
ggsave(here("Plots", "development_check_ppc_survival.svg"),
       width = 20, height = 14, units = "cm")
#2.3.2 Length Plot ####
lpc <- fit_dev %>%
  spread_draws(length_rep[i], ndraws = 200) %>% ungroup() %>%
  mutate(Plant = growth$Plant[i]) %>%
  ggplot() +
  geom_density(aes(length_rep, group = .draw), colour = "grey70", linewidth = 0.2) +
  geom_density(data = growth, aes(Length_mm), colour = "firebrick", linewidth = 1) +
  facet_wrap(~ Plant, scales = "free") +
  labs(x = "Final body length (mm)") + theme_classic()
lpc
ggsave(here("Plots", "development_check_ppc_length.svg"),
       width = 22, height = 14, units = "cm")

#2.3.3 Combo Plot ####
ppcscc <- (ppcsc + labs(tag = "A") + theme(plot.tag = element_text(face = "bold")))
lpcc <-  (lpc + labs(tag = "B") + theme(plot.tag = element_text(face = "bold")))
Plot_checks_dev <- wrap_elements(ppcscc) | wrap_elements(lpcc) +
  plot_layout(widths = c(2.5, 1))
Plot_checks_dev
ggsave(here("Plots", "Plot_checks_dev.svg"), Plot_checks_dev, width = 40, height = 20, units = "cm")

#2.4 Pooling check ####
pooling <- survival %>%
  summarise(placed = sum(n_placed), across(all_of(unname(OUTCOMES)), sum), .by = Plant) %>%
  pivot_longer(-c(Plant, placed), names_to = "column", values_to = "x") %>%
  mutate(Outcome = names(OUTCOMES)[match(column, OUTCOMES)],
         own_median = qbeta(0.500, 1 + x, 1 + placed - x),
         own_lower  = qbeta(0.025, 1 + x, 1 + placed - x),
         own_upper  = qbeta(0.975, 1 + x, 1 + placed - x)) %>%
  select(Outcome, Plant, placed, successes = x, starts_with("own_")) %>%
  left_join(select(survival_est, Outcome, Plant, starts_with("prob_")),
            by = c("Outcome", "Plant")) %T>% print(n = Inf)
write_csv(pooling, here("Output", "development_pooling.csv"))

#2.5. Prior sensitivity (BARG 5A, 5C) ####
summarise_fit_dev <- function(f, label) {
  f %>% spread_draws(prob[k, g], length_pct_of_control[g]) %>% ungroup() %>%
    filter(k == 2) %>%                                   # survival to trial end
    mutate(Plant = plants[g]) %>%
    summarise(across(c(prob, length_pct_of_control), q_eti), .by = Plant) %>%
    mutate(setting = label,
           divergences = sum(f$diagnostic_summary()$num_divergent),
           rhat_max = max(f$summary()$rhat, na.rm = TRUE), .before = 1)
}
prior_grid_dev <- list("reported"            = list(),
                   "prior_bin_mu sd 0.75"      = list(prior_bin_mu_sd = 0.75),
                   "prior_bin_mu sd 3"         = list(prior_bin_mu_sd = 3),
                   "prior_bin_sigma sd 0.75"   = list(prior_bin_sigma_sd = 0.75),
                   "prior_bin_sigma sd 3"      = list(prior_bin_sigma_sd = 3),
                   "prior_len_sigma sd 0.25"   = list(prior_len_sigma_sd = 0.25),
                   "prior_len_sigma sd 1"      = list(prior_len_sigma_sd = 1),
                   "prior_len_resid sd 0.125"  = list(prior_len_resid_sd = 0.125),
                   "prior_len_resid sd 0.5"    = list(prior_len_resid_sd = 0.5))
sensitivity_priors_dev <- imap(prior_grid_dev, function(change, label)
  summarise_fit_dev(refit_dev(change), label)) %>% list_rbind() %T>% print(n = Inf)
write_csv(sensitivity_priors_dev, here("Output", "development_sensitivity_priors.csv"))

#2.6. ROPE Sensitivity (BARG 4C, 5E) ####
sensitivity_rope_dev <- map(c(0.05, 0.10, 0.20), function(r)
  survival_draws %>%
    summarise(P_below_rope = mean(prob < r), .by = c(Outcome, Plant)) %>%
    mutate(rope = r, .before = 1)) %>%
  list_rbind() %>%
  pivot_wider(names_from = rope, values_from = P_below_rope,
              names_prefix = "P_below_") %T>% print(n = Inf)
write_csv(sensitivity_rope_dev, here("Output", "development_sensitivity_rope.csv"))

#2.7. Dish-level overdispersion ####
fit_dev_bb <- refit_dev(list(overdispersed = 1L))

dish_variation <- survival %>%
  mutate(across(all_of(unname(OUTCOMES)), ~ .x / n_placed, .names = "prop_{.col}")) %>%
  summarise(n_dish = n(), across(starts_with("prop_"), list(min = min, max = max)),
            .by = Plant) %T>% print(n = Inf)
write_csv(dish_variation, here("Output", "development_dish_variation.csv"))

overdispersion <- bind_rows(
  summarise_fit_dev(fit_dev, "binomial (reported)"),
  summarise_fit_dev(fit_dev_bb, "beta-binomial")) %>%
  arrange(Plant, setting) %T>% print(n = Inf)
write_csv(overdispersion, here("Output", "development_overdispersion.csv"))

rho_est <- fit_dev_bb %>% spread_draws(rho[k]) %>% ungroup() %>%
  mutate(Outcome = names(OUTCOMES)[k]) %>%
  summarise(across(rho, q_eti), .by = Outcome) %T>% print()
write_csv(rho_est, here("Output", "development_rho.csv"))

#2.8 LOO ####
dish_loglik <- function(f) {
  d   <- f$draws("log_lik_surv", format = "draws_matrix")
  rep_id <- as.integer(str_match(colnames(d), "log_lik_surv\\[\\d+,(\\d+)\\]")[, 2])
  sapply(split(seq_len(ncol(d)), rep_id), function(j) rowSums(d[, j, drop = FALSE]))
}
loo_binom <- loo(dish_loglik(fit_dev))
loo_betabin <- loo(dish_loglik(fit_dev_bb))
loo_survival <- loo_compare(list(binomial = loo_binom, `beta-binomial` = loo_betabin)) %T>% print()
write_csv(as_tibble(loo_survival, rownames = "model"),
          here("Output", "development_loo_survival.csv"))
