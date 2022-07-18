rm(list = ls())

library(tools)
library(rjson)
library(jsonlite)

library(goeveg)

library(tidyverse)
library(plotly)
library(Cairo)

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

read_analysis_data <- function(run_ids, csv_name) {
  analysis_data = list()
  for(id in run_ids){
    print(id)
    run_path = file.path(data_dir, id)
    data = read_run_data(run_path, csv_name)
    data_combined = bind_rows(data, .id = "config_id")
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


run_ids = build_run_ids_same_date("2022-07-14",1,56)


analysis_configs = read_analysis_configs(run_ids)
write.csv(analysis_configs, file.path(output_path, "config_summary.csv"))


analysis_meta_data = read_analysis_data(run_ids, "meta_data.csv")
analysis_meta_data = add_sim_mode(analysis_meta_data)

analysis_host_allele_data = read_analysis_data(run_ids, "host_allele_data.csv")


analysis_host_allele_counts = analysis_host_allele_data %>%
  count(run_id, config_id, generation, species, locus_id)

analysis_host_allele_counts_meta = add_meta_data_analysis(analysis_host_allele_counts, analysis_meta_data)


analysis_host_allele_counts_meta_last_100_by_sim_mode = analysis_host_allele_counts_meta %>%
  group_by(run_id, config_id, locus_id, derived_sim_mode) %>%
  filter(generation > max(generation) - 100) %>%
  add_run_config_analysis(analysis_configs)



analysis_host_allele_counts_meta_last_100_by_config_id_summary = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
  summarise(allele_count_mean = mean(n), allele_count_sd = sd(n), allele_count_median = median(n), .groups = "keep") %>%
  add_run_config_analysis(analysis_configs)


analysis_host_allele_counts_meta_last_100_summary = analysis_host_allele_counts_meta_last_100_by_sim_mode %>%
  group_by(run_id, locus_id, derived_sim_mode) %>%
  summarise(allele_count_mean = mean(n), allele_count_sd = sd(n), allele_count_median = median(n), .groups = "keep") %>%
  add_run_config_analysis(analysis_configs)

write.csv(analysis_host_allele_counts_meta_last_100_summary, file.path(output_path, "allele_counts_summary.csv"))

analysis_host_allele_counts_meta_last_100_summary_plt_box = ggplot(analysis_host_allele_counts_meta_last_100_by_sim_mode) +
  aes(
    x = derived_sim_mode,
    y = n,
    facet = run_id
  ) + 
  geom_boxplot() +
  theme_minimal() + 
  facet_wrap(
    vars(run_id)
  )

analysis_host_allele_counts_meta_last_100_summary_plt_box
