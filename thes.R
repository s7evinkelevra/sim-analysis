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

analysis_id <- "scenario_3"
output_path <- file.path(output_dir, analysis_id)
dir.create(output_path, showWarnings = FALSE)


build_run_ids_same_date <- function(date_base, from, to){
  run_ids = c()
  for(id in from:to){
    run_ids = append(run_ids, paste0(date_base,"-",id))
  }
  
  return(run_ids)
}

read_analysis_configs <- function(run_ids) {
  analysis_configs = list()
  for(id in run_ids) {
    run_path = file.path(data_dir, id)
    configs = read_run_configs(run_path)
    run_configs_tibble = generate_run_config_summary(id, configs)
    analysis_configs[[id]] = run_configs_tibble[1,]
  }
  
  analysis_configs_tbl = bind_rows(analysis_configs, .id = "run_id") %>% replace_na(list("hosts.introgression_individuals_per_generation" = "0", "pathogens.introgression_individuals_per_generation" = "0" ))
  
  return(analysis_configs_tbl)
}

read_analysis_data <- function(run_ids, csv_name, min_generation = 0, max_generation = 1000000) {
  analysis_data = list()
  for(id in run_ids){
    gc()
    print(id)
    run_path = file.path(data_dir, id)
    data = read_run_data(run_path, csv_name)
    data_combined = bind_rows(data, .id = "config_id") %>% filter(generation > min_generation & generation < max_generation)
    analysis_data[[id]] = data_combined
  }
  analysis_data_combined = bind_rows(analysis_data, .id = "run_id")
  return(analysis_data_combined)
}

add_meta_data_analysis <- function(left, meta_data) {
  return(
    left_join(left, meta_data, by = c("run_id", "config_id", "generation"))
  )
}

add_run_config_analysis <- function(left, analysis_configs) {
  return(
    left_join(left, analysis_configs, by = c("run_id"))
  )
}

add_merged_run_config_analysis <- function(left, analysis_configs) {
  return(
    left_join(left, analysis_configs, by = c("run_id_same_config_different_mode"))
  )
}

merge_by_paste_if_unequal <- function(x) {
  if(length(unique(x)) > 1){
    paste(x, collapse = ".")
  }else{
    x[1]  
  }
}

run_ids = build_run_ids_same_date("2022-07-14",1,48)


analysis_configs = read_analysis_configs(run_ids)
write.csv(analysis_configs, file.path(output_path, "config_summary.csv"))

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





gc()
analysis_meta_data = read_analysis_data(run_ids, "meta_data.csv") %>% add_sim_mode()
gc()
analysis_host_allele_data = read_analysis_data(run_ids, "host_allele_data.csv")


analysis_host_allele_counts = analysis_host_allele_data %>%
  count(run_id, config_id, generation, species, locus_id)

analysis_host_allele_counts_meta = analysis_host_allele_counts %>%
  add_meta_data_analysis(analysis_meta_data) %>%
  add_run_config_analysis(analysis_configs_merged_changed)



analysis_host_allele_counts_meta_last_100_by_sim_mode = analysis_host_allele_counts_meta %>%
  group_by(run_id_same_config_different_mode, config_id, locus_id, derived_sim_mode) %>%
  filter(generation > max(generation) - 100)



analysis_host_allele_counts_meta_last_100_by_config_id_summary = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
  summarise(allele_count_mean = mean(n), allele_count_sd = sd(n), allele_count_median = median(n), .groups = "keep") %>%
  add_merged_run_config_analysis(analysis_configs_merged_unique_changed)


analysis_host_allele_counts_meta_last_100_summary = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
  group_by(run_id_same_config_different_mode, locus_id, derived_sim_mode) %>%
  summarise(allele_count_mean = mean(n), allele_count_sd = sd(n), allele_count_median = median(n), .groups = "keep") %>%
  add_merged_run_config_analysis(analysis_configs_merged_unique_changed)


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

analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box 
analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_no_introgression_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_introgression_plt_box, 3000, 3000 )


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

analysis_host_allele_counts_meta_last_100_summary_id_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_config_id_summary) +
  aes(
    x = pathogens.mutation_rate_per_peptide,
    y = allele_count_mean,
    group = pathogens.mutation_rate_per_peptide,
    color = derived_sim_mode,
    facet = derived_sim_mode
  ) + 
  geom_boxplot() +
  facet_grid(
    . ~ derived_sim_mode + pathogens.introgression_individuals_per_generation + hosts.mutation_rate_per_peptide + pathogens.species_n, scales = "free_x"
  )



analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box
analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box
analysis_host_allele_counts_meta_last_100_summary_id_plt_box

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_no_introgression_id_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_introgression_id_plt_box, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_summary_id_plt_box.png"),analysis_host_allele_counts_meta_last_100_summary_id_plt_box, 5000, 3000 )




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
  geom_histogram(binwidth = 10, position = "identity", alpha = 0.2) +
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
  geom_histogram(binwidth = 10, position = "identity", alpha = 0.2) +
  theme_minimal() +
  facet_wrap(vars(
    run_id_same_config_different_mode
  ))

analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt
analysis_host_allele_counts_meta_last_100_dist_introgression_plt

save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt.png"),analysis_host_allele_counts_meta_last_100_dist_no_introgression_plt, 3000, 3000 )
save_plot_defaults(file.path(output_path, "analysis_host_allele_counts_meta_last_100_dist_introgression_plt.png"),analysis_host_allele_counts_meta_last_100_dist_introgression_plt, 3000, 3000 )



gc()
analysis_host_data = read_analysis_data(run_ids, "host_data.csv", min_generation = 1500)
gc()
analysis_host_genome_data = read_analysis_data(run_ids, "host_genome_data.csv", min_generation = 1500)
gc()
