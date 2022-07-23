library(tidyverse)
library(rstatix)

library(sjPlot)
library(visreg)
library(MASS)
library(ggplot2)
library(car)
library(jtools)


# get this from thes.R, shitty solution i know
view(analysis_host_allele_counts_meta_last_100_by_config_id_summary)

analysis_host_allele_counts_meta_last_100_summary_grid_stat_shapiro_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  # scenario 1: group_by(pathogens.introgression_individuals_per_generation,hosts.species_n,pathogens.species_n,infection.merit_threshold) %>%
  filter(infection.merit_threshold == 4) %>%
  #group_by(derived_sim_mode, pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  group_by(id_same_config_different_mode, derived_sim_mode) %>%
  shapiro_test(allele_count_mean) %>%
  add_significance()

analysis_host_allele_counts_meta_last_100_summary_grid_stat_pw_t_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  # scenario 1: group_by(pathogens.introgression_individuals_per_generation,hosts.species_n,pathogens.species_n,infection.merit_threshold) %>%
  filter(infection.merit_threshold == 4) %>%
  group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  pairwise_t_test(allele_count_mean ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")


analysis_allele_counts_no_merit_change = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  filter(infection.merit_threshold == 4)


analysis_allele_counts_no_merit_change_neutrality = analysis_allele_counts_no_merit_change %>% filter(derived_sim_mode == "Neutrality")
analysis_allele_counts_no_merit_change_coevolution = analysis_allele_counts_no_merit_change %>% filter(derived_sim_mode == "Coevolution")
analysis_allele_counts_no_merit_change_no_coevolution = analysis_allele_counts_no_merit_change %>% filter(derived_sim_mode == "No-coevolution")


summary(analysis_allele_counts_no_merit_change_neutrality)

neutrality_lm <- lm(allele_count_mean ~ (pathogens.introgression_individuals_per_generation + pathogens.species_n + hosts.species_n)^2,
                    data = analysis_allele_counts_no_merit_change_neutrality)

anova(neutrality_lm)

par(mfrow=c(2,2))
plot(neutrality_lm)
par(mfrow=c(1,1))


coevolution_lm <- lm(allele_count_mean ~ (pathogens.introgression_individuals_per_generation + pathogens.species_n + hosts.species_n)^2,
                    data = analysis_allele_counts_no_merit_change_coevolution)

anova(coevolution_lm)

par(mfrow=c(2,2))
plot(coevolution_lm)
par(mfrow=c(1,1))



no_coevolution_lm <- lm(allele_count_mean ~ (pathogens.introgression_individuals_per_generation + pathogens.species_n + hosts.species_n)^2,
                    data = analysis_allele_counts_no_merit_change_no_coevolution)

anova(no_coevolution_lm)

par(mfrow=c(2,2))
plot(no_coevolution_lm)
par(mfrow=c(1,1))



neutrality_aov = aov(allele_count_mean ~ pathogens.introgression_individuals_per_generation + pathogens.species_n + hosts.mutation_rate_per_peptide + infection.merit_threshold,
                     data = neutrality_only)

anova(neutrality_aov)

no_thresh = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  filter(derived_sim_mode == "Coevolution" & infection.merit_threshold == 4)

test.lm <- lm(allele_count_mean ~ (pathogens.introgression_individuals_per_generation + pathogens.mutation_rate_per_peptide + pathogens.species_n + hosts.mutation_rate_per_peptide)^2,
              data = no_thresh)


anova(test.lm)


test_aov <- aov(allele_count_mean ~ derived_sim_mode + pathogens.introgression_individuals_per_generation + pathogens.species_n + hosts.mutation_rate_per_peptide, 
                data = analysis_host_allele_counts_meta_last_100_by_config_id_summary)


par(mfrow=c(2,2))
plot(test_aov)
par(mfrow=c(1,1))


anova(test_aov)
