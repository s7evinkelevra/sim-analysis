rm(list = ls())

library(tools)
library(rjson)
library(jsonlite)

library(goeveg)
library(digest)

library(tidyverse)
library(plotly)
library(Cairo)
library(ggpubr)
library(rstatix)

library(esquisse)

source("./help.R")



data_dir <- "./data"
output_dir <- "./output"

analysis_id <- "scenario_4"
output_path <- file.path(output_dir, analysis_id)
dir.create(output_path, showWarnings = FALSE)



run_ids = build_run_ids_same_date("2022-07-15",1,32)

analysis_configs_with_config_ids = read_analysis_configs(run_ids)

config_ids_sample = analysis_configs_with_config_ids %>% pull(config_id)
  
analysis_configs = analysis_configs_with_config_ids %>% select(!config_id)

# needs to be executed in one go, else hash might change???
analysis_configs$hash_no_sim_mode = analysis_configs %>% select((!run_id & !simulation_mode)) %>% apply(1, digest)

analysis_configs_merged_by_no_sim_hash = analysis_configs %>%
  group_by(hash_no_sim_mode) %>%
  summarise_all(merge_by_paste_if_unequal) 

analysis_configs_hash_run_ids_map = analysis_configs_merged_by_no_sim_hash %>%
  mutate(run_id_same_config_different_mode = run_id) %>%
  select(hash_no_sim_mode, run_id_same_config_different_mode)


analysis_configs_merged = analysis_configs %>%
  left_join(analysis_configs_hash_run_ids_map, by = "hash_no_sim_mode")

analysis_configs_merged_unique = analysis_configs_merged %>%
  select((!run_id & !simulation_mode)) %>%
  distinct(hash_no_sim_mode, .keep_all = TRUE) %>%
  mutate(id_same_config_different_mode = 1:n())



analysis_configs_merged_changed_mask = analysis_configs_merged %>% sapply(function(x) !length(unique(x)) == 1)
analysis_configs_merged_changed = analysis_configs_merged[, analysis_configs_merged_changed_mask]
analysis_configs_merged_common = analysis_configs_merged[1, !analysis_configs_merged_changed_mask]

analysis_configs_merged_unique_changed_mask = analysis_configs_merged_unique %>% sapply(function(x) !length(unique(x)) == 1)
analysis_configs_merged_unique_changed = analysis_configs_merged_unique[, analysis_configs_merged_unique_changed_mask]
analysis_configs_merged_unique_common = analysis_configs_merged_unique[1, !analysis_configs_merged_unique_changed_mask]


write.csv(analysis_configs_merged, file.path(output_path, "config_summary.csv"))
write.csv(analysis_configs_merged_changed, file.path(output_path, "config_changed.csv"))
write.csv(analysis_configs_merged_common, file.path(output_path, "config_common.csv"))

gc()
analysis_meta_data = read_analysis_data(run_ids, "meta_data.csv") %>% add_sim_mode()
gc()
analysis_host_allele_data = read_analysis_data(run_ids, "host_allele_data.csv")

analysis_host_allele_data_downsampled = analysis_host_allele_data %>%
  filter(config_id %in% config_ids_sample) %>% # only get one config from each run
  filter(species == 0) %>% # only use first host species
  add_meta_data_analysis(analysis_meta_data) %>%
  add_run_config_analysis(analysis_configs_merged) %>%
  filter(infection.merit_threshold == 4) %>% # filter out merit, do later separately
  group_by(run_id_same_config_different_mode, config_id, locus_id, derived_sim_mode) %>% 
  mutate(age = generation - created_at) # %>%
  # filter for scenario 3, where not all host mut rate and patho mut rate combinations are present
  # filter(hosts.mutation_rate_per_peptide == 5e-5)



analysis_meta_data_summary = analysis_meta_data %>%
  group_by(run_id, config_id, derived_sim_mode) %>%
  summarise(generation_max = max(generation), step_duration_min = min(step_duration), step_duration_max = max(step_duration), step_duration_mean = mean(step_duration), .groups = "keep")

analysis_host_allele_counts = analysis_host_allele_data %>%
  count(run_id, config_id, generation, species, locus_id)

analysis_host_allele_counts_meta = analysis_host_allele_counts %>%
  add_meta_data_analysis(analysis_meta_data) %>%
  add_run_config_analysis(analysis_configs_merged)



analysis_host_allele_counts_meta_last_100_by_sim_mode = analysis_host_allele_counts_meta %>%
  group_by(run_id_same_config_different_mode, config_id, locus_id, derived_sim_mode) %>%
  filter(generation > max(generation) - 100)



analysis_host_allele_counts_meta_last_100_by_config_id_summary = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
  summarise(allele_count_mean = mean(n), allele_count_sd = sd(n), allele_count_median = median(n), .groups = "keep") %>%
  add_merged_run_config_analysis(analysis_configs_merged_unique)


analysis_host_allele_counts_meta_last_100_by_config_id_summary_summary = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  group_by(run_id_same_config_different_mode, derived_sim_mode) %>%
  summarise(allele_count_mean_mean = mean(allele_count_mean), allele_count_mean_min = min(allele_count_mean), allele_count_mean_max = max(allele_count_mean), .groups = "keep")

#view(analysis_host_allele_counts_meta_last_100_by_config_id_summary_summary)

analysis_host_allele_counts_meta_last_100_summary = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
  group_by(run_id_same_config_different_mode, locus_id, derived_sim_mode) %>%
  summarise(allele_count_mean = mean(n), allele_count_sd = sd(n), allele_count_median = median(n), allele_count_min = min(n), allele_count_max = max(n), .groups = "keep") %>%
  add_merged_run_config_analysis(analysis_configs_merged_unique)

###### Testing ######
# Test for normality
# test = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
#   group_by(run_id_same_config_different_mode, locus_id, derived_sim_mode) %>%
#   group_split(.keep = TRUE)
# 
# test2 = as.numeric(unlist(test[[5]]["n"]))
# 
# ggqqplot(test2)
# 
# qqnorm(test2)
# qqline(test2)
# 
# ks.test(test2, 'pnorm')
# 
# shapiro.test(rnorm(100, mean = 5, sd = 3))
# shapiro.test(test2)
# 
# test4 = compare_means(n ~ derived_sim_mode, group.by = c("run_id_same_config_different_mode"), method = "wilcox.test", data = analysis_host_allele_counts_meta_last_100_by_sim_mode)
# 
# test5 = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
#   group_by(run_id_same_config_different_mode) %>%
#   wilcox_test(n ~ derived_sim_mode) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance() %>%
#   add_xy_position(x = "derived_sim_mode")
# 
# test6 = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
#   filter(pathogens.introgression_individuals_per_generation == 0) %>%
#   group_by(run_id_same_config_different_mode) %>%
#   wilcox_test(allele_count_mean ~ derived_sim_mode) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance() %>%
#   add_xy_position(x = "derived_sim_mode")
# 
# 
# test7 = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
#   filter(pathogens.introgression_individuals_per_generation == 0) %>%
#   group_by(run_id_same_config_different_mode) %>%
#   shapiro_test(allele_count_mean) %>%
#   add_significance()
#  

write.csv(analysis_host_allele_counts_meta_last_100_summary, file.path(output_path, "allele_counts_summary.csv"))
write.table(analysis_host_allele_counts_meta_last_100_summary, file.path(output_path, "allele_counts_summary_excel.csv"), sep = ";", dec = ",", row.names = FALSE)


analysis_host_allele_counts_meta_last_100_summary_no_introgression_stat_pw_t_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  filter(pathogens.introgression_individuals_per_generation == 0) %>%
  group_by(id_same_config_different_mode) %>%
  pairwise_t_test(allele_count_mean ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")

analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
                                                                     filter(pathogens.introgression_individuals_per_generation == 0)) +
  aes(
    x = derived_sim_mode,
    y = allele_count_mean,
    color = derived_sim_mode,
    facet = id_same_config_different_mode
  ) + 
  geom_boxplot() +
  theme_minimal() + 
  stat_pvalue_manual(analysis_host_allele_counts_meta_last_100_summary_no_introgression_stat_pw_t_test) +
  facet_wrap(
    ~ id_same_config_different_mode
  )

analysis_host_allele_counts_meta_last_100_summary_introgression_stat_pw_t_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  filter(pathogens.introgression_individuals_per_generation > 0) %>%
  group_by(id_same_config_different_mode) %>%
  pairwise_t_test(allele_count_mean ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")

analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
                                                                                      filter(pathogens.introgression_individuals_per_generation > 0)) +
  aes(
    x = derived_sim_mode,
    y = allele_count_mean,
    color = derived_sim_mode,
    facet = id_same_config_different_mode
  ) + 
  geom_boxplot() +
  theme_minimal() + 
  stat_pvalue_manual(analysis_host_allele_counts_meta_last_100_summary_introgression_stat_pw_t_test) +
  facet_wrap(
    vars(id_same_config_different_mode)
  )


# todo(ajn): add proper test for normality. Previously, not all groups may have been properly set.
analysis_host_allele_counts_meta_last_100_summary_stat_pw_t_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  group_by(id_same_config_different_mode) %>%
  pairwise_t_test(allele_count_mean ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")


analysis_host_allele_counts_meta_last_100_summary_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary) +
  aes(
    x = derived_sim_mode,
    y = allele_count_mean,
    color = derived_sim_mode,
    facet = id_same_config_different_mode
  ) + 
  geom_boxplot() +
  theme_minimal() + 
  stat_pvalue_manual(analysis_host_allele_counts_meta_last_100_summary_stat_pw_t_test) +
  facet_wrap(
    . ~ id_same_config_different_mode,
    labeller = labeller(id_same_config_different_mode = label_from_id )
  )


analysis_host_allele_counts_meta_last_100_summary_grid_stat_shapiro_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  # scenario 2 filter(infection.merit_threshold == 4) %>%
  # scenario 2 group_by(pathogens.introgression_individuals_per_generation,hosts.species_n,pathogens.species_n) %>%
  # scenario 2 with threshold group_by(pathogens.introgression_individuals_per_generation,hosts.species_n,pathogens.species_n,infection.merit_threshold) %>%
  group_by(derived_sim_mode, pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  #group_by(id_same_config_different_mode, derived_sim_mode) %>%
  shapiro_test(allele_count_mean) %>%
  add_significance()

analysis_host_allele_counts_meta_last_100_summary_grid_stat_pw_t_test = analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
  # scenario 2 
  # filter(infection.merit_threshold == 4) %>%
  # scenario 2
  # group_by(pathogens.introgression_individuals_per_generation,hosts.species_n,pathogens.species_n) %>%
  # scenario 2 with threshold group_by(pathogens.introgression_individuals_per_generation,hosts.species_n,pathogens.species_n,infection.merit_threshold) %>%
  # other scenarios 
  group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  pairwise_t_test(allele_count_mean ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")

analysis_host_allele_counts_meta_last_100_summary_grid_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary %>% filter(infection.merit_threshold == 4)) +
  aes(
    x = derived_sim_mode,
    y = allele_count_mean,
    color = derived_sim_mode,
    #facet = id_same_config_different_mode
  ) + 
  geom_boxplot() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  #theme_minimal() + 
  stat_pvalue_manual(analysis_host_allele_counts_meta_last_100_summary_grid_stat_pw_t_test) +
  facet_grid(
    # scenario 2 
    # pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
    # scenario 2 with threshold 
    # pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n + infection.merit_threshold,
    #other scenarios: 
    pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide  ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free_x",
    labeller = label_both
  )



analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box 
analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box
analysis_host_allele_counts_meta_last_100_summary_plt_box
analysis_host_allele_counts_meta_last_100_summary_grid_plt_box

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_grid_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_grid_plt_box, 10000, 5000 )


analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
                                                                                      filter(pathogens.introgression_individuals_per_generation == 0)) +
  aes(
    x = id_same_config_different_mode,
    y = allele_count_mean,
    group = id_same_config_different_mode,
    color = derived_sim_mode,
    facet = derived_sim_mode
  ) + 
  geom_boxplot() +
  theme_minimal() + 
  scale_x_continuous(breaks = seq(1,80,1)) +
  facet_wrap(
    ~ derived_sim_mode, scales = "free_x"
  )


analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary %>%
                                                                                   filter(pathogens.introgression_individuals_per_generation > 0)) +
  aes(
    x = id_same_config_different_mode,
    y = allele_count_mean,
    group = id_same_config_different_mode,
    color = derived_sim_mode,
    facet = derived_sim_mode
  ) + 
  geom_boxplot() +
  theme_minimal() + 
  scale_x_continuous(breaks = seq(1,80,1)) +
  facet_wrap(
    ~ derived_sim_mode, scales = "free_x"
  )

analysis_host_allele_counts_meta_last_100_summary_id_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary %>% filter(infection.merit_threshold == 4)) +
  aes(
    # scenario 2
    x = hosts.species_n,
    # other scenarios 
    # x = pathogens.mutation_rate_per_peptide,
    y = allele_count_mean,
    group = hosts.species_n,
    # other scenarios group = pathogens.mutation_rate_per_peptide,
    color = derived_sim_mode,
    facet = derived_sim_mode
  ) + 
  geom_boxplot() +
  # other scenarios theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
  facet_grid(
    # scenario 2
    . ~ derived_sim_mode + pathogens.introgression_individuals_per_generation + pathogens.species_n,
    # scenario 2 with threshold . ~ derived_sim_mode + pathogens.introgression_individuals_per_generation + infection.merit_threshold + pathogens.species_n,
    # other scenarios . ~ derived_sim_mode + hosts.mutation_rate_per_peptide +  pathogens.introgression_individuals_per_generation + pathogens.species_n, 
    # labeller = label_both,
    # scenario 1: scales = "free_x"
  )



analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box
analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box
analysis_host_allele_counts_meta_last_100_summary_id_plt_box

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_id_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_id_plt_box, 7600, 3000 )




summary(analysis_host_allele_counts_meta_last_100_by_sim_mode)


analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt = ggplot(analysis_host_allele_counts_meta_last_100_by_sim_mode %>% 
                                                              filter(pathogens.introgression_individuals_per_generation == 0)) +
  aes(
    n,
    group = derived_sim_mode,
    color = derived_sim_mode,
    fill = derived_sim_mode,
    facet = run_id_same_config_different_mode
  ) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.2) +
  theme_minimal() +
  facet_wrap(vars(
    run_id_same_config_different_mode
  ))

analysis_host_allele_counts_meta_last_100_dist_introgression_plt = ggplot(analysis_host_allele_counts_meta_last_100_by_sim_mode %>% 
                                                                               filter(pathogens.introgression_individuals_per_generation > 0)) +
  aes(
    n,
    group = derived_sim_mode,
    color = derived_sim_mode,
    fill = derived_sim_mode,
    facet = run_id_same_config_different_mode
  ) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.2) +
  theme_minimal() +
  facet_wrap(vars(
    run_id_same_config_different_mode
  ))

analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt
analysis_host_allele_counts_meta_last_100_dist_introgression_plt

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt.png"),analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_dist_introgression_plt.png"),analysis_host_allele_counts_meta_last_100_dist_introgression_plt, 3000, 3000 )


##### allele counts over time

analysis_host_allele_counts_meta_grouped = analysis_host_allele_counts_meta %>%
  filter(infection.merit_threshold == 4) %>%
  group_by(run_id_same_config_different_mode, config_id, locus_id, derived_sim_mode) %>%
  filter(generation >= max(generation) - 1500) %>%
  mutate(neg_generation = generation - max(generation))

analysis_host_allele_counts_meta_grouped_mean = analysis_host_allele_counts_meta_grouped %>%
  group_by(run_id_same_config_different_mode, generation, locus_id, derived_sim_mode) %>%
  mutate(allele_count_mean = mean(n)) %>%
  ungroup()

analysis_host_allele_counts_meta_grouped_plt_grid = ggplot(analysis_host_allele_counts_meta_grouped_mean) +
  aes(
    x = neg_generation,
    y = n,
    group = interaction(species, config_id, derived_sim_mode),
    color = derived_sim_mode
  ) +
  geom_line(size = 0.2) +
  # geom_line(aes(
  #   y = allele_count_mean,
  #   color = "Mean",
  # ),
  # size = 0.4) +
  facet_grid(
    # scenario 2 
    # pathogens.introgression_individuals_per_generation + hosts.species_n ~ pathogens.species_n,
    # scenario 2 with threshold 
    # pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n + infection.merit_threshold,
    # other scenarios: 
    pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free_y",
    labeller = label_both
  )

analysis_host_allele_counts_meta_grouped_plt_grid

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_grouped_plt_grid_test.png"),analysis_host_allele_counts_meta_grouped_plt_grid, 5000, 8000 )


## allele frequency
# layout as grid, same as allele counts by sim mode
# except sim mode is additional split in y axis


analysis_host_allele_data_first_n_last_n = analysis_host_allele_data_downsampled %>% 
  filter(generation <= min(generation) + 100 | generation >= max(generation) - 100) %>%
  mutate(neg_generation = generation - max(generation)) %>% 
  mutate(first_last_generation = case_when(
    generation <= min(generation) + 100 ~ generation - min(generation),
    generation >= max(generation) - 100 ~ neg_generation,
  )) %>%
  mutate(first_or_last = case_when(
    generation <= min(generation) + 100 ~ "First",
    generation >= max(generation) - 100 ~ "Last"
  ))

analysis_host_allele_data_last_n = analysis_host_allele_data_first_n_last_n %>% filter(generation >= max(generation) - 100)

# only show last n generations
analysis_host_allele_data_last_n_freq_plt_grid = ggplot(analysis_host_allele_data_last_n) +
  aes(
    # last n generations
    x = neg_generation,
    y = frequency,
    color = created_at,
    group = interaction(allele_id, config_id)
  ) +
  geom_line(size = 0.2) +
  scale_color_viridis_c(option = "turbo") +
  facet_grid(
    # scenario 2 last n
    pathogens.introgression_individuals_per_generation + hosts.species_n + derived_sim_mode ~ pathogens.species_n,
    # other scenarios last n:
    # pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide + derived_sim_mode ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free",
    labeller = label_both
  )

analysis_host_allele_data_first_n_last_n_freq_plt_grid = ggplot(analysis_host_allele_data_first_n_last_n) +
  aes(
    # first and last gen
    x = first_last_generation,
    y = frequency,
    color = created_at,
    group = interaction(allele_id, config_id)
  ) +
  geom_line(size = 0.2) +
  #theme_minimal() +
  scale_color_viridis_c(option = "turbo") +
  facet_grid(
    # scenario 2 first n last n
    pathogens.introgression_individuals_per_generation + hosts.species_n + derived_sim_mode ~ pathogens.species_n + first_or_last,
    # other scenarios first n last n: 
    # pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide + derived_sim_mode ~ pathogens.species_n + pathogens.mutation_rate_per_peptide + first_or_last,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free",
    labeller = label_both
  )

analysis_host_allele_data_last_n_freq_plt_grid
analysis_host_allele_data_first_n_last_n_freq_plt_grid

save_plot_defaults(file.path(output_path, "analysis_host_allele_data_last_n_freq_plt_grid.png"),analysis_host_allele_data_last_n_freq_plt_grid, 8000, 10000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_data_first_n_last_n_freq_plt_grid.png"),analysis_host_allele_data_first_n_last_n_freq_plt_grid, 10000, 10000 )


#### allele age and recruitment rate

analysis_host_allele_data_downsampled_alleles_ids = analysis_host_allele_data_downsampled %>%
  group_by(config_id, derived_sim_mode) %>%
  distinct(allele_id) %>%
  slice(sample(n(), min(20, n()))) %>%
  pull(allele_id)

analysis_host_allele_data_downsampled_alleles = analysis_host_allele_data_downsampled %>% 
  filter(allele_id %in% analysis_host_allele_data_downsampled_alleles_ids) %>%
  filter(age <= 100)


# frequency over allele age plot
analysis_host_allele_data_last_n_age_plt_grid = ggplot(analysis_host_allele_data_downsampled_alleles) +
  aes(
    x = age,
    y = frequency,
    group = interaction(allele_id,config_id),
    color = derived_sim_mode
  ) +
  geom_line() +
  # geom_line(data = analysis_host_allele_data_downsampled_alleles_counts, aes(
  #   y = alleles_retained_percent,
  #   group = c(),
  #   color = c("black")
  # )) + 
  facet_grid(
    # scenario 2 first n last n
    # pathogens.introgression_individuals_per_generation + hosts.species_n + derived_sim_mode ~ pathogens.species_n,
    # other scenarios first n last n: 
    pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide + derived_sim_mode ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free",
    labeller = label_both
  )

analysis_host_allele_data_last_n_age_plt_grid

save_plot_defaults(file.path(output_path, "analysis_host_allele_data_last_n_age_plt_grid_b.png"),analysis_host_allele_data_last_n_age_plt_grid, 8000, 10000 )


# calculate allele max age 
analysis_host_allele_data_downsampled_max_age = analysis_host_allele_data_downsampled %>%
  group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species, allele_id) %>%
  filter(age == max(age)) %>%
  group_by(run_id_same_config_different_mode, derived_sim_mode, species) %>%
  filter(created_at >= max(generation) - 100) %>%
  summarise(max_age_allele_n = n(), max_age_mean = mean(age), max_age_median = median(age), max_age_max = max(age), max_age_min = min(age), .groups = "keep")


# allele recruitment rate
analysis_host_allele_data_downsampled_alleles_counts = analysis_host_allele_data_downsampled %>%
  # group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species, allele_id) %>%
  # filter(age == max(age)) %>%
  group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species) %>%
  # filters alleles created in the last 100 generations 
  # problem: in a window only 100 generations wide, there is only one possible generation in which an allele could have emerged, that reaces 100 generations of age
  # so: very old alleles (those that are most likely to not be in the 100 generation window) are not captured
  #filter(created_at >= max(generation) - 100) %>%
  # Filter all alleles that are currently at or below 100 generations old.
  # Foregoes issue described above, but also includes alleles from non-equilibrium generations
  filter(age <= 100) %>%
  group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species, age) %>%
  summarise(frequency_mean = mean(frequency), alleles_retained = n(), .groups = "keep") %>%
  group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species) %>%
  mutate(alleles_retained_percent = alleles_retained/alleles_retained[1L]) %>%
  add_merged_run_config_analysis(analysis_configs_merged_unique)

analysis_host_allele_data_downsampled_alleles_counts_plt_grid = ggplot(analysis_host_allele_data_downsampled_alleles_counts) +
  aes(
    x = age,
    y = alleles_retained_percent,
    group = interaction(config_id, derived_sim_mode),
    color = derived_sim_mode
  ) + 
  geom_line(size = 0.4) +
  facet_grid(
    # scenario 2
    pathogens.introgression_individuals_per_generation + hosts.species_n ~ pathogens.species_n,
    # scenario 2 split by sim mode
    # pathogens.introgression_individuals_per_generation + hosts.species_n + derived_sim_mode ~ pathogens.species_n,
    # other scenarios: 
    # pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    #scales = "free",
    labeller = label_both
  )


analysis_host_allele_data_downsampled_alleles_counts_plt_grid
save_plot_defaults(file.path(output_path, "analysis_host_allele_data_downsampled_alleles_counts_plt_grid.png"),analysis_host_allele_data_downsampled_alleles_counts_plt_grid, 3000, 5000 )


## non downsampled recruitment analysis
gc()
analysis_host_allele_data_created_last_n = analysis_host_allele_data %>%
  # group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species, allele_id) %>%
  # filter(age == max(age)) %>%
  add_meta_sim_mode_analysis(analysis_meta_data) %>%
  group_by(config_id, derived_sim_mode, species) %>%
  mutate(age = generation - created_at) %>%
  # filters alleles created in the last 100 generations 
  # problem: in a window only 100 generations wide, there is only one possible generation in which an allele could have emerged, that reaces 100 generations of age
  # so: very old alleles (those that are most likely to not be in the 100 generation window) are not captured
  filter(created_at >= max(generation) - 100) %>%
  # Filter all alleles that are currently at or below 100 generations old.
  # Foregoes issue described above, but also includes alleles from non-equilibrium generations
  # filter(age <= 1000) %>%
  add_run_config_run_id_same_config_different_mode_analysis(analysis_configs_merged)
  
  
gc()
analysis_host_allele_data_alleles_counts_by_age = analysis_host_allele_data_created_last_n %>%
  group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species, age) %>%
  filter(species == 0) %>%
  summarise(frequency_mean = mean(frequency), alleles_retained = n(), .groups = "keep") %>%
  group_by(run_id_same_config_different_mode, config_id, derived_sim_mode, species) %>%
  mutate(alleles_retained_percent = alleles_retained/alleles_retained[1L])  %>% 
  add_merged_run_config_analysis(analysis_configs_merged_unique) %>%
  filter(infection.merit_threshold == 4)
  


# allele recruitment proportion
analysis_host_allele_data_alleles_retained_between_0_n = analysis_host_allele_data_alleles_counts_by_age %>%
  filter(age == 20)

# allele overall retention proportion
analysis_host_allele_data_alleles_retained_between_0_end = analysis_host_allele_data_alleles_counts_by_age %>%
  filter(age == 100)
  
# allele retention proportion
analysis_host_allele_data_alleles_retained_between_n_m = analysis_host_allele_data_alleles_counts_by_age %>%
  filter(age == 20 | age == 100) %>%
  mutate(allele_retained_n_m = alleles_retained/alleles_retained[1L]) %>%
  filter(age == 100)



# test recruitment proportion for differences
analysis_host_allele_data_alleles_retained_between_0_n_shapiro_test = analysis_host_allele_data_alleles_retained_between_0_n %>%
  # scenario 2 pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
  group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n) %>%
  # other scenarios
  # group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  shapiro_test(alleles_retained_percent) %>%
  add_significance()

analysis_host_allele_data_alleles_retained_between_0_n_pw_test = analysis_host_allele_data_alleles_retained_between_0_n %>%
  # scenario 2 pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
  group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n) %>%
  # other scenarios
  # group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  wilcox_test(alleles_retained_percent ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")


# test overall retention proportion
analysis_host_allele_data_alleles_retained_between_0_end_shapiro_test = analysis_host_allele_data_alleles_retained_between_0_end %>%
  # scenario 2 pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
  group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n) %>%
  # other scenarios
  # group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  shapiro_test(alleles_retained_percent) %>%
  add_significance()

analysis_host_allele_data_alleles_retained_between_0_end_pw_test = analysis_host_allele_data_alleles_retained_between_0_end %>%
  # scenario 2 pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
  group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n) %>%
  # other scenarios
  # group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  wilcox_test(alleles_retained_percent ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")

# test retainment proportion for differences
analysis_host_allele_data_alleles_retained_between_n_m_shapiro_test = analysis_host_allele_data_alleles_retained_between_n_m %>%
  group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  shapiro_test(alleles_retained_percent) %>%
  add_significance()

analysis_host_allele_data_alleles_retained_between_n_m_pw_test = analysis_host_allele_data_alleles_retained_between_n_m %>%
  # scenario 2 pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
  group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n) %>%
  # other scenarios
  # group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide) %>%
  wilcox_test(alleles_retained_percent ~ derived_sim_mode) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")


# alleles retained after n generation box plt with significances
analysis_host_allele_data_alleles_retained_between_0_n_plt_box = ggplot(analysis_host_allele_data_alleles_retained_between_0_n) +
  aes(
    x = derived_sim_mode,
    y = alleles_retained_percent,
    group = derived_sim_mode,
    color = derived_sim_mode,
  ) +
  geom_boxplot() +
  stat_pvalue_manual(analysis_host_allele_data_alleles_retained_between_0_n_pw_test) +
  facet_grid(
    # scenario 2 
    pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
    # scenario 2 with threshold 
    # pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n + infection.merit_threshold,
    #other scenarios: 
    # pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide  ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    #scales = "free_x",
    labeller = label_both
  )

# alleles retained after investigated total interval generations box plt with significances
analysis_host_allele_data_alleles_retained_between_0_end_plt_box = ggplot(analysis_host_allele_data_alleles_retained_between_0_end) +
  aes(
    x = derived_sim_mode,
    y = alleles_retained_percent,
    group = derived_sim_mode,
    color = derived_sim_mode,
  ) +
  geom_boxplot() +
  stat_pvalue_manual(analysis_host_allele_data_alleles_retained_between_0_end_pw_test) +
  facet_grid(
    # scenario 2 
    pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
    # scenario 2 with threshold 
    # pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n + infection.merit_threshold,
    #other scenarios: 
    #pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide  ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    #scales = "free_x",
    labeller = label_both
  )


analysis_host_allele_data_alleles_retained_between_0_n_plt_box
analysis_host_allele_data_alleles_retained_between_0_end_plt_box

save_plot_defaults(file.path(output_path, "analysis_host_allele_data_alleles_retained_between_0_n_plt_box.png"),analysis_host_allele_data_alleles_retained_between_0_n_plt_box, 5000, 5000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_data_alleles_retained_between_0_end_plt_box.png"),analysis_host_allele_data_alleles_retained_between_0_end_plt_box, 5000, 5000 )


analysis_host_allele_data_alleles_retained_per_config_all_generations_plt_grid = ggplot(analysis_host_allele_data_alleles_counts_by_age) +
  aes(
    x = age,
    y = alleles_retained_percent,
    group = interaction(config_id, derived_sim_mode, species),
    color = derived_sim_mode
  ) + 
  geom_line(size = 0.4) +
  facet_grid(
    # scenario 2 first n last n
    pathogens.introgression_individuals_per_generation + hosts.species_n + derived_sim_mode ~ pathogens.species_n,
    # other scenarios first n last n: 
    # pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    #scales = "free",
    labeller = label_both
  )

analysis_host_allele_data_alleles_retained_per_config_all_generations_plt_grid
save_plot_defaults(file.path(output_path, "analysis_host_allele_data_alleles_retained_per_config_all_generations_plt_grid.png"),analysis_host_allele_data_alleles_retained_per_config_all_generations_plt_grid, 5000, 5000 )



# mean over all configs 
analysis_host_allele_data_alleles_counts_by_age_summary = analysis_host_allele_data_alleles_counts_by_age %>%
  group_by(run_id_same_config_different_mode, derived_sim_mode, species, age) %>%
  summarise(frequency_mean_mean = mean(frequency_mean), alleles_retained_mean = mean(alleles_retained), alleles_retained_percent_mean = mean(alleles_retained_percent), .groups = "keep") %>%
  add_merged_run_config_analysis(analysis_configs_merged_unique)


analysis_host_allele_data_alleles_retained_summary_all_generations_plt_grid = ggplot(analysis_host_allele_data_alleles_counts_by_age_summary) +
  aes(
    x = age,
    y = alleles_retained_percent_mean,
    group = derived_sim_mode,
    color = derived_sim_mode
  ) + 
  geom_line(size = 0.4) +
  facet_grid(
    # scenario 2 first n last n
    pathogens.introgression_individuals_per_generation + hosts.species_n ~ pathogens.species_n,
    # other scenarios first n last n: 
    # pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide ~ pathogens.species_n + pathogens.mutation_rate_per_peptide,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    #scales = "free",
    labeller = label_both
  )


analysis_host_allele_data_alleles_retained_summary_all_generations_plt_grid
save_plot_defaults(file.path(output_path, "analysis_host_allele_data_alleles_retained_summary_all_generations_plt_grid.png"),analysis_host_allele_data_alleles_retained_summary_all_generations_plt_grid, 3000, 5000 )





gc()
analysis_host_data_neutrality = read_analysis_data(run_ids, "host_data.csv", min_generation = 1400, max_generation = 1500)
gc()
analysis_host_genome_data_neutrality = read_analysis_data(run_ids, "host_genome_data.csv", min_generation = 1400, max_generation = 1500)
gc()

gc()
analysis_host_data_post_neutrality = read_analysis_data(run_ids, "host_data.csv", min_generation = 2000)
gc()
analysis_host_genome_data_post_neutrality = read_analysis_data(run_ids, "host_genome_data.csv", min_generation = 2000)
gc()


analysis_host_data_meta_last_n = analysis_host_data_neutrality %>% 
  bind_rows(analysis_host_data_post_neutrality) %>%
  add_meta_sim_mode_analysis(analysis_meta_data) %>%
  group_by(config_id, derived_sim_mode) %>%
  filter(generation >= max(generation) - 100)


analysis_host_genome_data_meta_last_n = analysis_host_genome_data_neutrality %>% 
  bind_rows(analysis_host_genome_data_post_neutrality) %>%
  add_meta_sim_mode_analysis(analysis_meta_data) %>%
  group_by(config_id, derived_sim_mode) %>%
  filter(generation >= max(generation) - 100) %>%
  mutate(zygosity = case_when(
    allele_1_id == allele_2_id ~ "Homozygous",
    allele_1_id != allele_2_id ~ "Heterozygous"
  ))



gc()


analysis_host_genome_meta_combined = analysis_host_data_meta_last_n %>%
  left_join(analysis_host_genome_data_meta_last_n, by = c("config_id", "generation", "species", "id"), suffix = c("", ".y")) %>%
  select(-ends_with(".y"))
  


analysis_host_genome_meta_combined_zygosity_summary = analysis_host_genome_meta_combined %>%
  group_by(run_id, config_id, derived_sim_mode, species, locus_id, zygosity) %>%
  summarise(fitness_mean = mean(fitness), fitness_sd = sd(fitness), fitness_cv = cv(fitness), count = n(), .groups = "keep") %>%
  filter(derived_sim_mode != "Neutrality") %>%
  add_run_config_analysis(analysis_configs_merged) %>%
  filter(infection.merit_threshold == 4) %>%
  filter(species == 0)




## statistical tests 
# shapiro test for normality
analysis_host_genome_meta_combined_zygosity_summary_shapiro_test = analysis_host_genome_meta_combined_zygosity_summary %>%
  # scenario 2 pathogens.introgression_individuals_per_generation + hosts.species_n  ~ pathogens.species_n,
  # group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n, derived_sim_mode, zygosity) %>%
  # other scenarios
  group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide, derived_sim_mode, zygosity, infection.merit_threshold) %>%
  shapiro_test(fitness_mean) %>%
  add_significance()

# paired t test probs would also work, almost all groups normally distributed
analysis_host_genome_meta_combined_zygosity_summary_sim_mode_pw_test = analysis_host_genome_meta_combined_zygosity_summary %>%
  # scenario 2 
  # group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n, zygosity) %>%
  # other scenarios
  group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide, zygosity) %>%
  pairwise_wilcox_test(fitness_mean ~ derived_sim_mode, exact = TRUE, p.adjust.method = "holm") %>%
  #adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "derived_sim_mode")


analysis_host_genome_meta_combined_zygosity_summary_zygosity_pw_test = analysis_host_genome_meta_combined_zygosity_summary %>%
  # scenario 2 
  # group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n, derived_sim_mode) %>%
  # other scenarios
  group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide, derived_sim_mode) %>%
  pairwise_wilcox_test(fitness_mean ~ zygosity, exact = TRUE, p.adjust.method = "holm") %>%
  #adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "zygosity")

analysis_host_genome_meta_combined_zygosity_summary_anova = analysis_host_genome_meta_combined_zygosity_summary %>%
  # scenario 2 
  # group_by(pathogens.introgression_individuals_per_generation, hosts.species_n, pathogens.species_n, derived_sim_mode) %>%
  # other scenarios
  group_by(derived_sim_mode) %>%
  anova_test(fitness_mean ~ zygosity)

#### wilcox experiments
# 
# analysis_host_genome_meta_combined_zygosity_summary_stats_test = analysis_host_genome_meta_combined_zygosity_summary %>%
#   group_by(pathogens.introgression_individuals_per_generation, hosts.mutation_rate_per_peptide, pathogens.species_n, pathogens.mutation_rate_per_peptide, zygosity) %>%
#   group_split()
# 
# 
# first_group = analysis_host_genome_meta_combined_zygosity_summary_stats_test[[1]]
# hist(first_group$fitness_mean)
# wilcox.test(fitness_mean ~ derived_sim_mode, data = analysis_host_genome_meta_combined_zygosity_summary_stats_test[[5]], exact = TRUE, correct = FALSE)

# todo: fix colors to match other colors for the sim modes
analysis_host_genome_meta_combined_zygosity_fitness_sim_mode_plt_box = ggplot(analysis_host_genome_meta_combined_zygosity_summary) +
  aes(
    x = derived_sim_mode,
    y = fitness_mean,
    group = derived_sim_mode,
    color = derived_sim_mode
  ) + 
  geom_boxplot() +
  stat_pvalue_manual(analysis_host_genome_meta_combined_zygosity_summary_sim_mode_pw_test) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
  facet_grid(
    # scenario 2 first n last n
    # pathogens.introgression_individuals_per_generation + hosts.species_n ~ pathogens.species_n + zygosity,
    # other scenarios first n last n: 
    pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide ~ pathogens.species_n + pathogens.mutation_rate_per_peptide + zygosity,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free_x",
    labeller = label_both
  )

# todo: fix colors, figure out colors for zygosity
analysis_host_genome_meta_combined_zygosity_fitness_zygosity_plt_box = ggplot(analysis_host_genome_meta_combined_zygosity_summary) +
  aes(
    x = zygosity,
    y = fitness_mean,
    group = zygosity,
    color = zygosity
  ) + 
  geom_boxplot() +
  stat_pvalue_manual(analysis_host_genome_meta_combined_zygosity_summary_zygosity_pw_test) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
  facet_grid(
    # scenario 2 first n last n
    # pathogens.introgression_individuals_per_generation + hosts.species_n ~ pathogens.species_n + derived_sim_mode,
    # other scenarios first n last n: 
    pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide ~ pathogens.species_n + pathogens.mutation_rate_per_peptide + derived_sim_mode,
    #labeller = labeller(id_same_config_different_mode = label_from_id ),
    scales = "free_x",
    labeller = label_both
  )



analysis_host_genome_meta_combined_zygosity_fitness_sim_mode_plt_box
analysis_host_genome_meta_combined_zygosity_fitness_zygosity_plt_box
save_plot_defaults(file.path(output_path, "analysis_host_genome_meta_combined_zygosity_fitness_sim_mode_plt_box.png"),analysis_host_genome_meta_combined_zygosity_fitness_sim_mode_plt_box, 5000, 5000 )
save_plot_defaults(file.path(output_path, "analysis_host_genome_meta_combined_zygosity_fitness_zygosity_plt_box.png"),analysis_host_genome_meta_combined_zygosity_fitness_zygosity_plt_box, 5000, 5000 )


gc()


gc()
analysis_host_data = read_analysis_data(run_ids, "host_data.csv")
gc()
analysis_host_genome_data = read_analysis_data(run_ids, "host_genome_data.csv")
gc()
