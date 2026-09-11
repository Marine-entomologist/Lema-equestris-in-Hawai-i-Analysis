#________________________####
#Assessing the impact of host plant on feeding ####
#________________________####

#1.1 Setting the stage (pun intended) ####
library(magrittr)
library(tidyverse)
library(here)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(ggdist)
library(patchwork)

SEED <- 20260902
set.seed(SEED)

STAGES <- c("Adult", "L1", "L1Naive")
STAGE_LABELS <- c(Adult = "Adult", L1 = "First Instar (L1)", L1Naive = "L1 - Naive")
PLANT_LATIN <- c(`American Nightshade`  = "S. americanum",
                 `Brazilian nightshade` = "S. seaforthianum",
                 Cabbage                = "B. oleracea",
                 Eggplant               = "S. melongena",
                 Gooseberry             = "P. peruviana",
                 Pepper                 = "C. annuum",
                 Potato                 = "S. tuberosum",
                 Tobacco                = "N. tabacum",
                 Tomatillo              = "P. philadelphica",
                 Tomato                 = "S. lycopersicum",
                 `Tree tobacco`         = "N. glauca")
ROPE_RATE <- c(Adult = 0.4, L1 = 0.02, L1Naive = 0.02)
ROPE_P    <- 0.10
q_eti <- list(median = median,
              lower  = function(x) unname(quantile(x, 0.025)),
              upper  = function(x) unname(quantile(x, 0.975)))
fmt <- function(m, l, u, dp = 3)
  ifelse(is.na(m), "", paste0(format(round(m, dp), nsmall = dp, trim = TRUE, scientific = FALSE), " [",
                              format(round(l, dp), nsmall = dp, trim = TRUE, scientific = FALSE), ", ",
                              format(round(u, dp), nsmall = dp, trim = TRUE, scientific = FALSE), "]"))

#1.2 Data import ####
beetles <- read_csv(here("Data", "SolaniDietFinal.csv"), show_col_types = FALSE) %>%
  mutate(stage = factor(Lifestage, STAGES),
         Plant_common = Plant,
         Plant = unname(PLANT_LATIN[Plant_common]),
         group = fct(str_c(Lifestage, Plant, sep = " | ")),
         rate_obs = AreaConsumed / (BeetleCount * Hours))

groups <- beetles %>%
  summarise(stage     = first(stage),
            ref_count = as.integer(round(median(BeetleCount))),
            ref_hours = median(Hours),
            n_trial   = n(),
            n_fed     = sum(Fed),
            fed_any   = any(Fed == 1),
            .by = c(group, Lifestage, Plant)) %>%
  arrange(as.integer(group))

plant_species <- sort(unique(beetles$Plant))

#1.3 Priors ####
priors <- list(prior_logit_p_mu_sd = 2, prior_logit_p_sigma_sd = 1,
               prior_log_mu_mu_loc = log(0.1), prior_log_mu_mu_sd = 2,
               prior_log_mu_sigma_sd = 1, prior_phi_rate = 1,
               prior_plant_sigma_sd = 1)

#1.4 Model ####
stan_data <- c(list(
  n = nrow(beetles), n_group = nrow(groups), n_stage = length(STAGES),
  n_plant = length(plant_species),
  group = as.integer(beetles$group), stage = as.integer(groups$stage),
  group_plant = match(groups$Plant, plant_species),
  beetle_count = beetles$BeetleCount, hours = beetles$Hours,
  area = beetles$AreaConsumed, fed = beetles$Fed,
  ref_count = groups$ref_count, ref_hours = groups$ref_hours,
  n_mean_level = length(STAGES), mean_level = seq_along(STAGES),
  n_scale_level = 1L, scale_level = rep(1L, length(STAGES)),
  plant_effect = 0L,   #1 adds a plant effect shared across life stages (see checks script)
  prior_only = 0L
), priors)

model <- cmdstan_model(here("Stan", "Beetles.stan"))
fit <- model$sample(data = stan_data, chains = 8,
                    parallel_chains = min(8, parallel::detectCores()),
                    iter_warmup = 2000, iter_sampling = 2000,
                    adapt_delta = 0.99, max_treedepth = 12,
                    seed = SEED, refresh = 0)
fit$save_object(here("Output", "feeding_fit.rds"))

#1.5 Convergence checks ####
fs <- fit$summary()
dg <- fit$diagnostic_summary()
write_csv(fs, here("Output", "feeding_parameter_summary.csv")) #Export for BARG (2B, 2C)
diagnostics <- tibble(rhat_max       = max(fs$rhat, na.rm = TRUE),
                      ess_bulk_min   = min(fs$ess_bulk, na.rm = TRUE),
                      ess_tail_min   = min(fs$ess_tail, na.rm = TRUE),
                      divergences    = sum(dg$num_divergent),
                      treedepth_hits = sum(dg$num_max_treedepth),
                      ebfmi_min      = min(dg$ebfmi)) %T>% print()
write_csv(diagnostics, here("Output", "feeding_diagnostics.csv"))

#1.6 Results ####
draws <- fit %>%
  recover_types(select(beetles, group)) %>%
  spread_draws(p[group], mu[group], rate_uncond[group], y_ref_rate[group]) %>%
  ungroup() %>%
  left_join(groups, by = "group") %>%
  mutate(rope = ROPE_RATE[as.character(Lifestage)])

plant_order <- draws %>%
  summarise(m = median(rate_uncond), .by = c(Lifestage, Plant)) %>%
  summarise(m = max(m), .by = Plant) %>%
  arrange(m) %>% pull(Plant)

draws  <- mutate(draws,  Plant = factor(Plant, plant_order),
                 Lifestage = factor(Lifestage, STAGES))
groups <- mutate(groups, Plant = factor(Plant, plant_order))
beetles <- mutate(beetles, Plant = factor(Plant, plant_order),
                  Lifestage = factor(Lifestage, STAGES))

#1.6.1 Estimates ####
estimates <- draws %>%
  summarise(across(c(p, mu, rate_uncond), q_eti),
            P_rope = mean(rate_uncond < rope),
            .by = c(Lifestage, Plant, n_trial, n_fed, fed_any)) %>%
  mutate(across(starts_with("mu_"), ~ if_else(fed_any, .x, NA_real_))) %>%
  arrange(Lifestage, desc(P_rope)) %T>% print(n = Inf)
write_csv(estimates, here("Output", "feeding_estimates.csv"))

estimates %>%
  transmute(Lifestage, Plant,
            `Trials (fed)` = str_c(n_trial, " (", n_fed, ")"),
            p     = fmt(p_median, p_lower, p_upper),
            mu    = fmt(mu_median, mu_lower, mu_upper),
            `p x mu` = fmt(rate_uncond_median, rate_uncond_lower, rate_uncond_upper),
            `P(ROPE)` = format(round(P_rope, 3), nsmall = 3, trim = TRUE, scientific = FALSE)) %>%
  print() %>%
  write_csv(here("Output", "feeding_table2.csv"))

#1.6.2 Contrasts ####
long <- draws %>%
  select(.draw, Lifestage, Plant, fed_any, p, mu, rate_uncond) %>%
  pivot_longer(c(p, mu, rate_uncond), names_to = "parameter") %>%
  filter(parameter != "mu" | fed_any) %>%
  select(-fed_any)
long <- bind_rows(long,
                  long %>% filter(Lifestage == "L1", Plant == "S. americanum") %>%
                    mutate(Lifestage = "L1Naive"))

plant_key <- bind_rows(groups,
                       groups %>% filter(Lifestage == "L1", Plant == "S. americanum") %>%
                         mutate(Lifestage = "L1Naive")) %>%
  select(Lifestage, Plant, fed_any)
pair_key <- plant_key %>%
  inner_join(plant_key, by = "Lifestage", suffix = c("_a", "_b"),
             relationship = "many-to-many") %>%
  filter(as.character(Plant_a) < as.character(Plant_b))

contrast_draws <- pair_key %>%
  inner_join(long, by = c("Lifestage", "Plant_a" = "Plant"),
             relationship = "many-to-many") %>%
  rename(value_a = value) %>%
  inner_join(long, by = c(".draw", "Lifestage", "parameter", "Plant_b" = "Plant")) %>%
  rename(value_b = value) %>%
  mutate(difference = value_a - value_b,
         ratio = case_when(
           parameter == "mu" & !(fed_any_a & fed_any_b) ~ NA_real_, .default = value_a / value_b),
         fold = case_when(
           parameter == "mu" & !(fed_any_a & fed_any_b) ~ NA_real_, .default = value_b / value_a),
         rope  = case_when(parameter == "p"           ~ ROPE_P,
                           parameter == "rate_uncond" ~ ROPE_RATE[as.character(Lifestage)],
                           .default = NA_real_),both_fed = fed_any_a & fed_any_b)

by_pair <- c("Lifestage", "parameter", "Plant_a", "Plant_b")
contrasts <- contrast_draws %>%
  summarise(across(difference, q_eti),
            P_direction = max(mean(difference > 0), mean(difference < 0)),
            P_rope      = mean(abs(difference) < rope),
            both_fed    = first(both_fed),
            .by = all_of(by_pair)) %>%
  left_join(contrast_draws %>% filter(!is.na(ratio)) %>%
              summarise(across(ratio, q_eti), .by = all_of(by_pair)), by = by_pair) %>%
  left_join(contrast_draws %>% filter(!is.na(fold)) %>%
              summarise(across(fold, q_eti), .by = all_of(by_pair)), by = by_pair) %>%
  arrange(Lifestage, parameter, desc(P_direction)) %T>% print(n = Inf)
write_csv(contrasts, here("Output", "feeding_contrasts.csv"))

#1.7 Plot it up ####
#1.7.1 Feeding estimates Combo ####
rope_band_rate <- tibble(Lifestage = factor(names(ROPE_RATE), STAGES),
                         rope = unname(ROPE_RATE))
obs_raw <- bind_rows(
  beetles %>% summarise(value = mean(Fed), .by = c(Lifestage, Plant)) %>%
    mutate(parameter = "p", zero = FALSE),
  beetles %>% filter(Fed == 1) %>%
    transmute(Lifestage, Plant, value = rate_obs, parameter = "mu", zero = FALSE),
  beetles %>%
    transmute(Lifestage, Plant, value = rate_obs, parameter = "rate_uncond",
              zero = AreaConsumed == 0))

param_draws <- draws %>%
  select(.draw, Lifestage, Plant, fed_any, p, mu, rate_uncond) %>%
  pivot_longer(c(p, mu, rate_uncond), names_to = "parameter") %>%
  filter(parameter != "mu" | fed_any)
PAR_LABEL <- list(
  p           = expression("Probability of feeding"),
  mu          = expression(atop("Rate given feeding",
                                "(mm"^2*" beetle"^-1*" h"^-1*")")),
  rate_uncond = expression(atop("Expected rate (p x mu)",
                                "(mm"^2*" beetle"^-1*" h"^-1*")")))

stage_plants <- groups %>%
  summarise(plants = list(levels(droplevels(Plant))), n = n(), .by = Lifestage) %>%
  mutate(Lifestage = factor(Lifestage, STAGES)) %>% arrange(Lifestage)

feeding_cell <- function(stg, par, xmax = NA) {
  lev <- stage_plants$plants[[match(stg, as.character(stage_plants$Lifestage))]]
  d   <- param_draws %>% filter(as.character(Lifestage) == stg, parameter == par) %>%
    mutate(Plant = factor(as.character(Plant), lev))
  raw <- obs_raw %>% filter(as.character(Lifestage) == stg, parameter == par) %>%
    mutate(Plant = factor(as.character(Plant), lev))
  g <- ggplot(d, aes(value, Plant))
  if (par == "p")
    g <- g + annotate("rect", xmin = 0, xmax = ROPE_P, ymin = -Inf, ymax = Inf,
                      fill = "firebrick", alpha = 0.12)
  if (par == "rate_uncond")
    g <- g + annotate("rect", xmin = 0, xmax = unname(ROPE_RATE[stg]),
                      ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.12)
  g +
    geom_jitter(data = raw, aes(shape = zero), height = 0.15, width = 0,
                size = 1.4, alpha = 0.7, colour = "grey25") +
    stat_pointinterval(.width = c(0.5, 0.95), point_size = 1.8) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.05))) +
    labs(x = PAR_LABEL[[par]], subtitle = STAGE_LABELS[[stg]]) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          plot.subtitle = element_text(face = "bold", size = 9))
}

mu_shared_xmax <- bind_rows(
  param_draws %>% filter(Lifestage %in% c("L1", "L1Naive"), parameter == "mu"),
  obs_raw     %>% filter(Lifestage %in% c("L1", "L1Naive"), parameter == "mu")
) %>%
  pull(value) %>% max(na.rm = TRUE)

cells <- expand_grid(Lifestage = STAGES, parameter = c("p", "mu", "rate_uncond")) %>%
  mutate(Lifestage = factor(Lifestage, rev(STAGES))) %>%
  arrange(Lifestage, factor(parameter, c("p", "mu", "rate_uncond"))) %>%
  mutate(xmax = if_else(parameter == "mu" & Lifestage %in% c("L1", "L1Naive"),
                        mu_shared_xmax-0.35, NA_real_), plot = pmap(list(as.character(Lifestage), parameter, xmax), feeding_cell))

cells$plot[[1]] <- cells$plot[[1]] + labs(tag = "A")
cells$plot[[2]] <- cells$plot[[2]] + labs(tag = "B")
cells$plot[[3]] <- cells$plot[[3]] + labs(tag = "C")

feeding_grid <- wrap_plots(cells$plot, ncol = 3, byrow = TRUE) +
  plot_layout(heights = stage_plants$n, axis_titles = "collect_x") &
  theme(plot.tag = element_text(face = "bold"))
feeding_grid
ggsave(here("Plots", "feeding_estimates.svg"), feeding_grid,
       width = 26, height = 21, units = "cm")

#1.7.2 Probability of feeding, all combinations on one axis ####
feeding_probability_plot <- draws %>%
  mutate(Lifestage = factor(Lifestage, levels = rev(STAGES))) %>%
  ggplot(aes(p, Plant)) +
  annotate("rect", xmin = 0, xmax = ROPE_P, ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.12) +
  geom_point(data = obs_raw %>% filter(parameter == "p") %>%
               mutate(Lifestage = factor(Lifestage, levels = rev(STAGES))),
             aes(value, Plant), colour = "grey30", shape = 4, stroke = 1) +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  facet_grid(Lifestage ~ ., scales = "free_y", space = "free_y",
             labeller = as_labeller(STAGE_LABELS)) +
  labs(x = "Probability of feeding") +
  theme_classic() + theme(axis.title.y = element_blank(), strip.text = element_text(face = "bold"))
feeding_probability_plot
ggsave(here("Plots", "feeding_probability.svg"), feeding_probability_plot,
       width = 16, height = 20, units = "cm")

#1.7.3 Posterior predictive feeding rate (supplemental) ####
VIS_CAP <- 0.99
draws_vis <- draws %>%
  mutate(cap = quantile(y_ref_rate, VIS_CAP), .by = Lifestage) %>%
  filter(y_ref_rate <= cap)

feeding_prediction_plot <- ggplot(draws_vis, aes(y_ref_rate, Plant)) +
  geom_rect(data = rope_band_rate, inherit.aes = FALSE,
            aes(xmin = 0, xmax = rope, ymin = -Inf, ymax = Inf),
            fill = "firebrick", alpha = 0.12) +
  geom_jitter(data = beetles, aes(rate_obs, Plant),
              height = 0.15, width = 0, size = 1.2, alpha = 0.5, colour = "grey30") +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  facet_wrap(~ Lifestage, ncol = 1, scales = "free",
             labeller = as_labeller(STAGE_LABELS)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression("Predicted feeding rate (mm"^2*" beetle"^-1*" h"^-1*")")) +
  theme_classic() + theme(axis.title.y = element_blank(), strip.text = element_text(face = "bold"))
feeding_prediction_plot
ggsave(here("Plots", "feeding_prediction.svg"), feeding_prediction_plot,
       width = 16, height = 22, units = "cm")

#1.7.4 Feeding Rate Contrasts ####
rope_band <- tibble(Lifestage = factor(names(ROPE_RATE), STAGES),
                    half_width = unname(ROPE_RATE))
contrast_draws %>%
  filter(parameter == "rate_uncond") %>%
  mutate(pair = fct_reorder(str_c(Plant_a, " - ", Plant_b), difference)) %>%
  ggplot(aes(difference, pair, shape = both_fed)) +
  geom_rect(data = rope_band, inherit.aes = FALSE,
            aes(xmin = -half_width, xmax = half_width, ymin = -Inf, ymax = Inf),
            fill = "firebrick", alpha = 0.12) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = "Both plants fed on") +
  facet_wrap(~ Lifestage, ncol = 1, scales = "free",
             labeller = as_labeller(STAGE_LABELS)) +
  labs(x = expression("Difference in expected rate (mm"^2*" beetle"^-1*" h"^-1*")")) +
  theme_classic() + theme(axis.title.y = element_blank(), legend.position = "bottom")
ggsave(here("Plots", "feeding_contrasts.svg"), width = 18, height = 34, units = "cm")

#1.7.5 Feeding Probability Contrasts ####
contrast_draws %>%
  filter(parameter == "p") %>%
  mutate(pair = fct_reorder(str_c(Plant_a, " - ", Plant_b), difference)) %>%
  ggplot(aes(difference, pair)) +
  annotate("rect", xmin = -ROPE_P, xmax = ROPE_P, ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.12) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  facet_wrap(~ Lifestage, ncol = 1, scales = "free_y",
             labeller = as_labeller(STAGE_LABELS)) +
  labs(x = "Difference in the feeding probability") +
  theme_classic() + theme(axis.title.y = element_blank())
ggsave(here("Plots", "feeding_contrasts_probability.svg"),
       width = 18, height = 40, units = "cm")

#________________________####
#Assessing the impact of host plant on development (larval growth and survival) ####
#________________________####
#2. Development and survival models ####
#2.1 Setup ####
CONTROL <- "S. americanum"
ROPE_SURVIVAL <- 0.10 #Not really relevant to the data (Turned out to be a strict 0 or 100), but keeping consistant with the feeding model
OUTCOMES <- c(survival_24h = "n_survived_24h", survival_end = "n_survived_end")
OUTCOME_LABELS <- c(survival_24h = "Survival to 24 hrs",
                    survival_end = "Survival to trial end (72 hrs)")
STEP_LABELS    <- c(survival_24h = "Survival to 24 hrs (of larvae placed)",
                    survival_end = "Survival to trial end (of larvae alive at 24 hrs)")

#2.2 Data import ####
survival <- read_csv(here("Data", "SolaniDevelopment_survival.csv"), show_col_types = FALSE)
growth   <- read_csv(here("Data", "SolaniDevelopment_growth.csv"),   show_col_types = FALSE)
plants <- c(CONTROL, sort(setdiff(union(survival$Plant, growth$Plant), CONTROL)))
pid <- function(x) match(x, plants)

#2.3 Priors ####
priors_dev <- list(prior_bin_mu_sd = 1.5, prior_bin_sigma_sd = 1.5,
               prior_len_mu_loc = log(1.5), prior_len_mu_sd = 0.5,
               prior_len_sigma_sd = 0.5, prior_len_resid_sd = 0.25,
               prior_rho_a = 1, prior_rho_b = 1)
#2.4 Model ####
stan_data_dev <- c(list(
  n_plant = length(plants), control_plant = pid(CONTROL),
  n_rep = nrow(survival), n_outcome = length(OUTCOMES),
  rep_plant = pid(survival$Plant), n_placed = survival$n_placed,
  successes = t(as.matrix(survival[unname(OUTCOMES)])),
  n_larva = nrow(growth), larva_plant = pid(growth$Plant),
  length_mm = growth$Length_mm,
  overdispersed = 0L,  #1 fits a beta-binomial effect (see checks script)
  prior_only = 0L
), priors_dev)

model_dev <- cmdstan_model(here("Stan", "Development.stan"))
fit_dev <- model_dev$sample(data = stan_data_dev, chains = 8,
                    parallel_chains = min(8, parallel::detectCores()),
                    iter_warmup = 2000, iter_sampling = 2000,
                    adapt_delta = 0.99, max_treedepth = 12,
                    seed = SEED, refresh = 0)
fit_dev$save_object(here("Output", "development_fit.rds"))

#2.5 Convergence checks ####
fs_dev <- fit_dev$summary()
dg_dev <- fit_dev$diagnostic_summary()
write_csv(fs_dev, here("Output", "development_parameter_summary.csv"))  #Export for BARG (2B, 2C)
diagnostics_dev <- tibble(rhat_max       = max(fs_dev$rhat, na.rm = TRUE),
                      ess_bulk_min   = min(fs_dev$ess_bulk, na.rm = TRUE),
                      ess_tail_min   = min(fs_dev$ess_tail, na.rm = TRUE),
                      divergences    = sum(dg_dev$num_divergent),
                      treedepth_hits = sum(dg_dev$num_max_treedepth),
                      ebfmi_min      = min(dg_dev$ebfmi)) %T>% print()
write_csv(diagnostics_dev, here("Output", "development_diagnostics.csv"))

#2.6 Results ####
#2.6.1 Survival Estimates ####
survival_draws <- fit_dev %>%
  spread_draws(prob[k, g], prob_step[k, g], prob_diff_vs_control[k, g]) %>% ungroup() %>%
  mutate(Outcome = factor(names(OUTCOMES)[k], names(OUTCOMES)), Plant = plants[g])

observed_survival <- survival %>%
  summarise(n_rep = n(), placed = sum(n_placed),
            across(all_of(unname(OUTCOMES)), sum), .by = Plant) %>%
  pivot_longer(all_of(unname(OUTCOMES)), names_to = "column", values_to = "survived") %>%
  mutate(Outcome = factor(names(OUTCOMES)[match(column, OUTCOMES)], names(OUTCOMES)),
         at_risk = if_else(Outcome == "survival_24h", placed, NA_integer_)) %>%
  select(-column) %>%
  arrange(Plant, Outcome) %>%
  mutate(at_risk = if_else(is.na(at_risk), lag(survived), at_risk), .by = Plant) %>%
  mutate(observed = survived / placed,
         label_cum  = str_c(survived, "/", placed),
         label_step = str_c(survived, "/", at_risk))

survival_est <- survival_draws %>%
  summarise(across(c(prob, prob_step, prob_diff_vs_control), q_eti),
            P_below_rope = mean(prob < ROPE_SURVIVAL),
            .by = c(Outcome, Plant)) %>%
  left_join(select(observed_survival, Outcome, Plant, n_rep, placed, at_risk, survived),
            by = c("Outcome", "Plant")) %>%
  arrange(Outcome, desc(prob_median)) %T>% print(n = Inf)
write_csv(survival_est, here("Output", "development_survival.csv"))

survival_est %>%
  transmute(Outcome, Plant,
            Observed = str_c(survived, "/", placed, " (", n_rep, " reps)"),
            `Survival (cumulative)` = fmt(prob_median, prob_lower, prob_upper),
            `Survival (this step)`  = fmt(prob_step_median, prob_step_lower, prob_step_upper),
            `Difference vs control` = fmt(prob_diff_vs_control_median,
                                          prob_diff_vs_control_lower,
                                          prob_diff_vs_control_upper),
            `P(< ROPE)` = format(round(P_below_rope, 3), nsmall = 3, trim = TRUE, scientific = FALSE)) %>%
  print() %>%
  write_csv(here("Output", "development_survival_formatted.csv"))

#2.6.2 Length Estimates ####
length_est <- fit_dev %>%
  spread_draws(length_mean[g], length_pct_of_control[g]) %>% ungroup() %>%
  mutate(Plant = plants[g]) %>%
  summarise(across(c(length_mean, length_pct_of_control), q_eti), .by = Plant) %>%
  left_join(summarise(growth, n_larva = n(), raw_mean = mean(Length_mm),
                      raw_sd = sd(Length_mm), .by = Plant), by = "Plant") %>%
  arrange(desc(length_pct_of_control_median)) %T>% print(n = Inf)
write_csv(length_est, here("Output", "development_length.csv"))

length_est %>%
  transmute(Plant, n_larva,
            `Raw mean (SD)` = paste0(format(round(raw_mean, 2), nsmall = 2, trim = TRUE, scientific = FALSE), " (",
                                     format(round(raw_sd, 2), nsmall = 2, trim = TRUE, scientific = FALSE), ")"),
            `Modelled length, mm` = fmt(length_mean_median, length_mean_lower,
                                        length_mean_upper, 2),
            `% of control` = fmt(length_pct_of_control_median,
                                 length_pct_of_control_lower,
                                 length_pct_of_control_upper, 1)) %>%
  print() %>%
  write_csv(here("Output", "development_length_formatted.csv"))

#2.7 Plot it up ####
#2.7.1 Survival plot ####
surv_order <- survival_draws %>%
  filter(Outcome == "survival_end") %>%
  summarise(m = median(prob), .by = Plant) %>% arrange(m) %>% pull(Plant)

survival_dev_plot <- survival_draws %>%
  mutate(Plant = factor(Plant, surv_order)) %>%
  ggplot(aes(prob, Plant)) +
  annotate("rect", xmin = 0, xmax = ROPE_SURVIVAL, ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.12) +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  geom_point(data = mutate(observed_survival, Plant = factor(Plant, surv_order)),
             aes(observed, Plant), colour = "firebrick", shape = 4, stroke = 1.1) +
  geom_text(data = mutate(observed_survival, Plant = factor(Plant, surv_order)),
            aes(x = 1.06, y = Plant, label = label_cum), hjust = 0, size = 2.7) +
  facet_wrap(~ Outcome, ncol = 1, labeller = as_labeller(OUTCOME_LABELS)) +
  coord_cartesian(xlim = c(0, 1.02), clip = "off") +
  scale_x_continuous(expand = expansion(0), breaks = seq(0, 1, 0.25)) +
  labs(x = "Cumulative survival probability per larva") +
  theme_classic() +
  theme(axis.title.y = element_blank(), strip.text = element_text(face = "bold"),
        plot.margin = margin(5.5, 34, 5.5, 5.5))
survival_dev_plot
ggsave(here("Plots", "development_survival.svg"), survival_dev_plot,
       width = 12, height = 15, units = "cm")

#2.7.2 Length Plot ####
len_order <- length_est %>% arrange(length_pct_of_control_median) %>% pull(Plant)
length_draws <- fit_dev %>%
  spread_draws(length_mean[g], length_pct_of_control[g]) %>% ungroup() %>%
  mutate(Plant = factor(plants[g], len_order))

length_abs_plot <- ggplot(length_draws, aes(length_mean, Plant)) +
  geom_jitter(data = mutate(growth, Plant = factor(Plant, len_order)),
              aes(Length_mm, Plant), height = 0.15, width = 0,
              size = 1.2, alpha = 0.5, colour = "grey30") +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  labs(x = "Body length (mm)") +
  theme_classic() + theme(axis.title.y = element_blank())

length_pct_plot <- ggplot(length_draws, aes(length_pct_of_control, Plant)) +
  geom_vline(xintercept = 100, linetype = 2, colour = "grey40") +
  stat_pointinterval(.width = c(0.5, 0.95)) +
  labs(x = "Length as % of control") +
  theme_classic() + theme(axis.title.y = element_blank())

development_length <- length_abs_plot + length_pct_plot +
  plot_layout(axes = "collect")
development_length
ggsave(here("Plots", "development_length.svg"), development_length,
       width = 20, height = 10, units = "cm")

#2.7.3 Combo plot ####
development_row <- (survival_dev_plot + labs(tag = "D")) +
  (length_abs_plot + labs(tag = "E")) +
  (length_pct_plot + labs(tag = "F")) +
  plot_layout(widths = c(1.15, 1, 1)) &
  theme(plot.tag = element_text(face = "bold"))
development_row 
ggsave(here("Plots", "development_plot.svg"), development_row,
       width = 26, height = 10, units = "cm")

Comboplot <- wrap_elements(feeding_grid) / wrap_elements(development_row) +
  plot_layout(heights = c(2.2, 1))
Comboplot
ggsave(here("Plots", "combo_plot.svg"), Comboplot, width = 26, height = 30, units = "cm")

#________________________####
#3. Reproducibility ####
#Per BARG (2A, 6B, 6H)
writeLines(c(capture.output(sessionInfo()),
             str_c("cmdstan version: ", cmdstan_version()),
             str_c("seed: ", SEED)),
           here("Output", "model_session_info.txt"))
