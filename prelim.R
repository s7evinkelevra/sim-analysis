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
current_run_id <- "2022-07-14-25"
current_run_path <- file.path(data_dir, current_run_id)

output_dir <- "./output"



process_run_data <- function(data_dir, run_id) {
  run_path = file.path(data_dir, run_id)

  

}


##### Read & combine run data ####
# read run data and combine into one tibble
# add meta data to most

# read confings
run_configs = read_run_configs(current_run_path)

# read run meta data
run_meta_data = read_run_data(current_run_path, "meta_data.csv")
run_meta_data_combined = bind_rows(run_meta_data, .id="config_id")

# build mode info from flags
run_meta_data_combined = add_sim_mode(run_meta_data_combined)

# read allele frequency data
run_host_allele_data = read_run_data(current_run_path, "host_allele_data.csv")
run_host_allele_data_combined = bind_rows(run_host_allele_data, .id = "config_id") %>% mutate(age = generation - created_at)
run_host_allele_data_combined_meta <- add_data(run_host_allele_data_combined, run_meta_data_combined)

# read host data (infection/fitness data + ancestry)
run_host_data = read_run_data(current_run_path, "host_data.csv")
run_host_data_combined = bind_rows(run_host_data, .id = "config_id")
run_host_data_combined_meta = add_data(run_host_data_combined, run_meta_data_combined) 

# read host genome data (allele ids of all host individuals)
run_host_genome_data = read_run_data(current_run_path, "host_genome_data.csv")
run_host_genome_data_combined = bind_rows(run_host_genome_data, .id = "config_id")
run_host_genome_data_combined_meta = add_data(run_host_genome_data_combined, run_meta_data_combined)

# read pathogen data
run_pathogen_data = read_run_data(current_run_path, "pathogen_data.csv")
run_pathogen_data_combined = bind_rows(run_pathogen_data, .id = "config_id")
run_pathogen_data_combined_meta = left_join(run_pathogen_data_combined, run_meta_data_combined, by = c("config_id", "generation" = "pathogen_generation"))
  


#### Housekeeping ####
# do housekeeping stuff like creating dirs and cleaning the workspace

# Create output folder for this run
dir.create(file.path(output_dir, current_run_id), showWarnings = FALSE)



#### Configs ####
# process the configs to generate config summaries etc.

# Write summary of the configs to csv
run_config_summary = generate_run_config_summary(current_run_id, run_configs)
write.csv(run_config_summary, file.path(output_dir, current_run_id, "config_summary.csv"))

# Get all config vars that have been changed throughout the run
# not really necessary anymore as a run only contains replicates
run_config_changed_mask = run_config_summary %>% sapply(function(x) !length(unique(x)) == 1)
run_config_changed = run_config_summary[, run_config_changed_mask]
run_config_common = run_config_summary[1, !run_config_changed_mask]

# write common config values in long format to file for better use later
run_config_key_value = run_config_common %>% pivot_longer(cols = everything(), names_to = "property", values_to = "values")
write.csv(run_config_key_value, file.path(output_dir, current_run_id, "config_common.csv"), row.names = FALSE)





#### Process host fitness data

## successful antigen presentations in the last generation
run_host_data_last_gen_successful_presentation <- run_host_data_combined_meta %>%
  group_by(config_id, species) %>%
  filter(generation == max(generation)) %>%
  mutate(successful_presentations_fac = as.factor(successful_presentations))
  

run_host_data_last_gen_successful_presentation_dist_plt <- ggplot(run_host_data_last_gen_successful_presentation) +
  aes(
    x = config_id,
    group = successful_presentations,
    fill = successful_presentations_fac
  ) +
  geom_bar(position = "dodge") +
  theme_minimal()

run_host_data_last_gen_successful_presentation_dist_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_presentation_counts_last_gen.png"), run_host_data_last_gen_successful_presentation_dist_plt, 3000, 1500)

## successful antigen presentations sum over all generations
run_host_data_successful_presentation <- run_host_data_combined_meta %>%
  group_by(config_id, species, successful_presentations) %>%
  filter(bBurninMode == FALSE)


run_host_data_successful_presentation_dist_plt <- ggplot(run_host_data_successful_presentation) +
  aes(
    x = config_id,
    group = successful_presentations,
    fill = as.factor(successful_presentations)
  ) +
  geom_bar(position = "dodge") +
  theme_minimal()

run_host_data_successful_presentation_dist_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_presentation_counts_all.png"), run_host_data_successful_presentation_dist_plt, 3000, 1500)

## host fitness distribution
run_host_data_combined_meta_no_burnin = run_host_data_combined_meta %>%
  group_by(config_id, species) %>%
  filter(bBurninMode == FALSE)

run_host_data_combined_meta_last_generation = run_host_data_combined_meta %>%
  group_by(config_id, species) %>%
  filter(generation == max(generation))

run_host_fitness_dist_last_gen_plt = ggplot(run_host_data_combined_meta_last_generation) +
  aes(
    x = fitness,
    group = config_id,
    color = fitness,
    fill = fitness,
    facet = config_id
  ) +
  geom_bar(aes(fill = fitness)) +
  theme_minimal() +
  facet_wrap(
    vars(config_id)
  )


run_host_fitness_dist_last_gen_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_dist_last_gen.png"), run_host_fitness_dist_last_gen_plt, 3000, 1500)


run_host_fitness_dist_no_burnin_plt = ggplot(run_host_data_combined_meta_no_burnin) +
  aes(
    x = fitness,
    group = config_id,
    color = fitness,
    fill = fitness,
    facet = config_id
  ) +
  geom_bar(aes(fill = fitness)) +
  theme_minimal() +
  facet_wrap(
    vars(config_id)
  )

run_host_fitness_dist_no_burnin_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_dist_no_burnin.png"), run_host_fitness_dist_no_burnin_plt, 3000, 1500)


## host fitness over time summary
run_host_data_fitness_summary <- run_host_data_combined_meta %>%
  group_by(config_id, generation, species, bBurninMode) %>%
  summarize(fitness_mean = mean(fitness), fitness_sd = sd(fitness), fitness_cv = cv(fitness), fitness_median = median(fitness), .groups = "keep")

run_host_data_fitness_summary_no_burnin = run_host_data_fitness_summary %>%
  filter(bBurninMode == FALSE)


# plt: fitness mean over time
#TODO(JAN): add CI whiskers
run_host_data_fitness_mean_over_time_plt_autoscale <- ggplot(run_host_data_fitness_summary_no_burnin) +
  aes(
    x = generation,
    y = fitness_mean,
    group = interaction(config_id, species),
    color = config_id
  ) +
  geom_hline(yintercept = as.numeric(run_config_common[1,"hosts.fitness_minimum"])) +
  geom_point() +
  geom_line() +
  theme_minimal()

run_host_data_fitness_mean_over_time_plt_autoscale
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_mean_autoscale.png"), run_host_data_fitness_mean_over_time_plt_autoscale, 3000, 1500)

run_host_data_fitness_mean_over_time_plt = run_host_data_fitness_mean_over_time_plt_autoscale + ylim(0,1)
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_mean.png"), run_host_data_fitness_mean_over_time_plt, 3000, 1500)


run_host_data_fitness_mean_over_time_plt_autoscale_facet <- ggplot(run_host_data_fitness_summary_no_burnin) +
  aes(
    x = generation,
    y = fitness_mean,
    group = interaction(config_id, species),
    color = config_id,
    facet = config_id
  ) +
  geom_hline(yintercept = as.numeric(run_config_common[1,"hosts.fitness_minimum"])) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  facet_wrap(
    vars(config_id)
  )

run_host_data_fitness_mean_over_time_plt_autoscale_facet
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_mean_autoscale_facet.png"), run_host_data_fitness_mean_over_time_plt_autoscale_facet, 3000, 1500)


# plt host fitness coefficient of variation (measure of dispersion) over time
# here: standin for selection
run_host_data_fitness_cv_no_burnin_plt <- ggplot(run_host_data_fitness_summary_no_burnin) +
  aes(
    x = generation,
    y = fitness_cv,
    group = interaction(config_id, species),
    color = config_id
  ) +
  geom_point() +
  geom_line() +
  theme_minimal()


run_host_data_fitness_cv_no_burnin_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_cv.png"), run_host_data_fitness_cv_no_burnin_plt, 3000, 1500)

# plt fitness cv
run_host_data_fitness_cv_box_plt <- ggplot(run_host_data_fitness_summary_no_burnin) +
  aes(
    x = config_id,
    y = fitness_cv
  ) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = config_id)) +
  theme_minimal()

run_host_data_fitness_cv_box_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_cv_box.png"), run_host_data_fitness_cv_box_plt, 3000, 1500)


#### Process pathogen fitness data ####

### pathogen fitness over time summary
run_pathogen_data_fitness_mean_over_time_burnin_infecting <- run_pathogen_data_combined_meta %>%
  filter(bBurninMode == FALSE & total_infections > 0) %>%
  group_by(config_id, generation, species) %>%
  summarise(fitness_mean = mean(fitness), fitness_sd = sd(fitness), fitness_cv = cv(fitness), fitness_median = median(fitness), .groups = "keep")


run_pathogen_data_fitness_mean_no_burnin_infecting_plt_autoscale <- ggplot(run_pathogen_data_fitness_mean_over_time_burnin_infecting) +
  aes(
    x = generation,
    y = fitness_mean,
    group = interaction(config_id, species),
    color = config_id
  ) +
  geom_hline(yintercept = as.numeric(run_config_common[1,"pathogens.fitness_minimum"])) +
  geom_point() +
  geom_line() +
  theme_minimal()

run_pathogen_data_fitness_mean_no_burnin_infecting_plt_autoscale
save_plot_defaults(file.path(output_dir, current_run_id, "pathogen_fitness_mean_no_burnin_infecting_autoscale.png"), run_pathogen_data_fitness_mean_no_burnin_infecting_plt_autoscale, 3000, 1500)


run_pathogen_data_fitness_mean_no_burnin_infecting_plt = run_pathogen_data_fitness_mean_no_burnin_infecting_plt_autoscale + ylim(0,1)
run_pathogen_data_fitness_mean_no_burnin_infecting_plt
save_plot_defaults(file.path(output_dir, current_run_id, "pathogen_fitness_mean_no_burnin_infecting.png"), run_pathogen_data_fitness_mean_no_burnin_infecting_plt, 3000, 1500)



#### Process host heterozygosity ####

# add allele frequency data 
# currently not used
# run_host_genome_data_combined_meta_count_freq = run_host_genome_data_combined_meta %>% 
#   left_join(run_host_allele_data_combined, 
#                  by = c("config_id", "generation", "species", "locus_id","allele_1_id" = "allele_id"), 
#                  suffix = c(".allele_1", ".allele_2")) %>%
#   left_join(run_host_allele_data_combined,
#             by = c("config_id", "generation", "species", "locus_id","allele_2_id" = "allele_id"), 
#             suffix = c(".allele_1", ".allele_2"))


# calculate Hexp 
# sum allele frequencies per generation/config/species/locus
run_host_allele_data_hexp = run_host_allele_data_combined_meta %>%
  mutate(frequency_sq = frequency * frequency) %>%
  group_by(config_id, generation, species, locus_id) %>%
  summarise(freq_sq_sum = sum(frequency_sq), .groups = "keep") %>%
  mutate(hexp = 1 - freq_sq_sum)


# add label to individuals
run_host_data_combined_meta_genome = run_host_data_combined_meta %>%
  left_join(run_host_genome_data_combined, by = c("config_id", "generation", "species", "id")) %>%
  mutate(zygosity = case_when(
    allele_1_id == allele_2_id ~ "homozygous",
    allele_1_id != allele_2_id ~ "heterozygous"
  ))

# filter by burnin
run_host_data_combined_meta_genome_no_burnin = run_host_data_combined_meta_genome %>%
  filter(bBurninMode == FALSE)

# summarize by zygosity
run_host_data_summary_by_zygosity = run_host_data_combined_meta_genome %>%
  group_by(config_id, generation, species, locus_id, zygosity, bBurninMode) %>%
  summarise(fitness_mean = mean(fitness), fitness_sd = sd(fitness), fitness_cv = cv(fitness), .groups = "keep")

run_host_data_summary_by_zygosity_no_burnin = run_host_data_combined_meta_genome %>%
  filter(bBurninMode == FALSE) %>%
  group_by(config_id, generation, species, locus_id, zygosity) %>%
  summarise(fitness_mean = mean(fitness), fitness_sd = sd(fitness), fitness_cv = cv(fitness), .groups = "keep")


# plt fitness over time by zygosity
run_host_data_summary_by_zygosity_no_burnin_plt = ggplot(run_host_data_summary_by_zygosity_no_burnin) + 
  aes(
    x = generation,
    y = fitness_mean,
    group = zygosity,
    color = zygosity,
    facet = config_id
  ) +
  geom_line() + 
  geom_vline(xintercept = as.numeric(run_config_common[1,"burnin_generations"])) + 
  theme_minimal() + 
  scale_color_discrete(name = "Zygosity", labels = c("Hetero","Homo")) +
  facet_wrap(
    vars(config_id),
  )

run_host_data_summary_by_zygosity_no_burnin_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_mean_by_zygosity_no_burnin.png"), run_host_data_summary_by_zygosity_no_burnin_plt, 5000, 3000)

# plt fitness over zygosity
run_host_data_combined_meta_genome_box_no_burnin_plt = ggplot(run_host_data_summary_by_zygosity_no_burnin) + 
  aes(
    x = zygosity,
    y = fitness_mean,
    color = zygosity,
    facet = config_id
  ) + 
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(
    vars(config_id)
  )

run_host_data_combined_meta_genome_box_no_burnin_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_mean_by_zygosity_box_no_burnin.png"), run_host_data_combined_meta_genome_box_no_burnin_plt, 5000, 3000)


## calculate hobserved vs hexpected
run_host_genome_zygosity_counts = run_host_genome_data_combined_meta %>%
  mutate(zygosity = case_when(
    allele_1_id == allele_2_id ~ "homozygous",
    allele_1_id != allele_2_id ~ "heterozygous"
  )) %>%
  group_by(config_id, generation, species, locus_id) %>%
  count(zygosity, name = "zygosity_count") %>%
  pivot_wider(names_from = zygosity, values_from = zygosity_count) %>%
  mutate(hobs = heterozygous/(homozygous+heterozygous))



run_host_genome_hobs_hexp = run_host_genome_zygosity_counts %>%
  left_join(run_host_allele_data_hexp, by = c("config_id", "generation", "species", "locus_id")) %>%
  mutate(hobsoverhexp = hobs/hexp) %>%
  mutate(hobsminushexp = hobs-hexp)



run_host_genome_hobs_hexp_plt = ggplot(run_host_genome_hobs_hexp) +
  aes(
    x = generation,
    facet = config_id
  ) +
  geom_line(aes(y = hobs, color = "blue")) + 
  geom_line(aes(y = hexp, color = "red")) + 
  scale_color_discrete(name = "Heterozygosity", labels = c("Hexp","Hobs")) +
  geom_vline(xintercept = as.numeric(run_config_common[1,"burnin_generations"])) + 
  theme_minimal() +
  facet_wrap(
    vars(config_id),
  )

run_host_genome_hobs_hexp_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_heterozygosity.png"), run_host_genome_hobs_hexp_plt, 5000, 3000)

run_host_genome_hobs_over_hexp_plt = ggplot(run_host_genome_hobs_hexp) +
  aes(
    x = generation,
    y = hobsoverhexp,
    group = config_id,
    color = config_id
  ) +
  geom_line() + 
  geom_vline(xintercept = as.numeric(run_config_common[1,"burnin_generations"])) + 
  theme_minimal()

run_host_genome_hobs_over_hexp_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_hobs_over_hexp.png"), run_host_genome_hobs_over_hexp_plt, 5000, 3000)


run_host_genome_hobs_minus_hexp_plt = ggplot(run_host_genome_hobs_hexp) +
  aes(
    x = generation,
    y = hobsminushexp,
    group = config_id,
    color = config_id,
    facet = config_id
  ) +
  geom_line() + 
  geom_vline(xintercept = as.numeric(run_config_common[1,"burnin_generations"])) + 
  theme_minimal() +
  facet_wrap(
    vars(config_id),
  )

run_host_genome_hobs_minus_hexp_plt
save_plot_defaults(file.path(output_dir, current_run_id, "host_hobs_minus_hexp.png"), run_host_genome_hobs_minus_hexp_plt, 5000, 3000)


#### calculate host marginal fitness of alleles ####
# WIP

run_host_data_combined_meta_genome_longer_by_allele = run_host_data_combined_meta_genome %>%
  pivot_longer(cols = c("allele_1_id", "allele_2_id"), names_to = "chromosome", values_to = "allele_id")

run_host_data_marginal_fitness_summary_no_burnin = run_host_data_combined_meta_genome_longer_by_allele %>%
  filter(bBurninMode == FALSE) %>%
  group_by(config_id, generation, species, allele_id) %>%
  summarise(marginal_fitness_mean = mean(fitness), marginal_fitness_sd = sd(fitness), marginal_fitness_cv = cv(fitness), .groups = "keep")

run_host_data_marginal_fitness_mean_no_burnin_plt = ggplot(run_host_data_marginal_fitness_summary_no_burnin) +
  aes(
    x = generation,
    y = marginal_fitness_mean,
    group = allele_id,
    color = allele_id,
    facet = config_id
  ) +
  geom_point() + 
  theme_minimal() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(
    vars(config_id)
  )

#run_host_data_marginal_fitness_mean_no_burnin_plt
# maybe not that interesting as its basically just the mean fitness(?)
#TODO(JAN): save plot

run_host_data_marginal_fitness_mean_box_no_burnin_plt = ggplot(run_host_data_marginal_fitness_summary_no_burnin) +
  aes(
    x = generation,
    y = marginal_fitness_mean,
    group = generation,
    color = config_id,
    facet = config_id
  ) +
  geom_violin() + 
  theme_minimal() +
  facet_wrap(
    vars(config_id)
  )

run_host_data_marginal_fitness_mean_box_no_burnin_plt


run_host_data_marginal_fitness_summary_frequencies = run_host_data_marginal_fitness_summary_no_burnin %>%
  left_join(run_host_allele_data_combined, by = c("config_id", "generation", "species", "allele_id"))




## Allele counts
# get allele counts by config, generation and species
run_allele_counts = run_host_allele_data_combined %>% count(config_id, generation, species)
run_allele_counts_grouped = run_allele_counts %>% group_by(config_id)

# add meta data to allele counts -> simulation mode most importantly
run_allele_counts_grouped_meta = add_data(run_allele_counts_grouped, run_meta_data_combined)


# select the last 100 burnin and last 100 post-burnin generations
run_allele_counts_last_100_by_burnin <- run_allele_counts_grouped_meta %>%
  group_by(config_id, bBurninMode) %>%
  filter(generation >= max(generation) - 100)


# build mean allele counts and SD by config and mode for the last 100 generations
run_allele_counts_last_100_mean <- run_allele_counts_last_100_by_burnin %>%
  summarise(mean_alleles = mean(n), sd = sd(n), .groups = "keep")



## Allele frequency distribution
#TODO(JAN): not working yet, fix
run_allele_data_first_last_1 <- run_host_allele_data_combined_meta %>%
  group_by(config_id, bBurninMode) %>%
  filter(generation == max(generation) | generation == 0)

summary(run_allele_data_first_last_1)


allele_dist_first_last_plt = ggplot(run_allele_data_first_last_1) +
  aes(
    x = count,
    facet = generation,
    group = bBurninMode
  ) +
  geom_histogram( color="#e9ecef", alpha=0.6, binwidth = 5) +
  facet_wrap(
    vars(generation),
  )

allele_dist_first_last_plt


# boxplot of the allele counts of the last 100 generations by simulation mode
run_allele_counts_last_100_box_plt = ggplot(run_allele_counts_last_100_by_burnin) +
  aes(
    x = bBurninMode,
    y = n,
    facet=config_id
  ) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = config_id)) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1000, by = 10)) +
  facet_wrap(
    vars(config_id),
    labeller=labeller(config_id = label_from_config_id)
  )

run_allele_counts_last_100_box_plt

save_plot_defaults(file.path(output_dir, current_run_id, "allele_counts_last_100_box.png"), run_allele_counts_last_100_box_plt, 3000, 3000)


# allele counts over time
# TODO(JAN): last 100 generations
run_allele_counts_plt = ggplot(run_allele_counts_grouped) +
  aes(
    x = generation,
    y = n,
    color = as.factor(species),
    group = species
  ) +
  geom_line(size = 0.2) +
  expand_limits(y=0) + 
  theme_minimal() +
  facet_wrap(
    vars(config_id),
  )


run_allele_counts_plt
save_plot_defaults(file.path(output_dir, current_run_id, "allele_counts.png"), run_allele_counts_plt, 20000, 5000)


run_allele_counts_last_n = run_allele_counts_grouped_meta %>%
  group_by(config_id,species) %>%
  filter(generation > max(generation) - 500)

run_allele_counts_last_100_plt = ggplot(run_allele_counts_last_n) +
  aes(
    x = generation,
    y = n,
    color = config_id,
    group = species,
  ) +
  geom_line(size = 0.5) +
  expand_limits(y=0) + 
  theme_minimal() +
  facet_wrap(
    vars(config_id),
  )


run_allele_counts_last_100_plt
save_plot_defaults(file.path(output_dir, current_run_id, "allele_counts_last_500.png"), run_allele_counts_last_100_plt, 20000, 5000)


# allele age
run_host_allele_data_age_summary = run_host_allele_data_combined_meta %>% 
  group_by(config_id, generation, species, locus_id) %>%
  summarise(age_mean = mean(age), .groups = "keep")


run_host_allele_data_age_plt = ggplot(run_host_allele_data_age_summary) +
  aes(
    x = generation,
    y = age_mean,
    color = config_id,
    facet = config_id
  ) + 
  geom_line() +
  theme_minimal() + 
  facet_wrap(
    vars(config_id)
  )
  
run_host_allele_data_age_plt
save_plot_defaults(file.path(output_dir, current_run_id, "allele_age_mean.png"), run_host_allele_data_age_plt, 5000, 3000)

# allele frequencies over time
# run_allele_freqs_plt = ggplot(run_host_allele_data_combined) +
#   aes(
#     x = generation,
#     y = frequency,
#     color = created_at,
#     group = allele_id,
#     facet = config_id
#   ) +
#   geom_line(size = 0.2) +
#   theme_minimal() +
#   scale_color_viridis_c(option = "turbo") +
#   facet_wrap(
#     vars(config_id)
#   )
# 
# run_allele_freqs_plt
# save_plot_defaults(file.path(output_dir, current_run_id, "allele_freqs.png"), run_allele_freqs_plt, 5000, 3000)


# only last n generations
run_host_allele_data_combined_last_n = run_host_allele_data_combined %>%
  filter(generation > max(generation) - 400)

run_allele_freqs_last_n_plt = ggplot(run_host_allele_data_combined_last_n) +
  aes(
    x = generation,
    y = frequency,
    color = created_at,
    group = allele_id,
    facet = config_id
  ) +
  geom_line(size = 0.2) +
  theme_minimal() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(
    vars(config_id)
  )

run_allele_freqs_last_n_plt
save_plot_defaults(file.path(output_dir, current_run_id, "allele_freqs_last_400.png"), run_allele_freqs_last_n_plt, 5000, 3000)


# run_allele_frequencies_plt = ggplot(run_host_allele_data_combined) +
#  aes(
#    x = generation,
#    y = count,
#    fill = allele_id,
#    color = config_id,
#    group = config_id
#  ) +
#  geom_bar(position = "stack",stat = "identity") +
#  geom_step(size = 0.5) +
#  #scale_color_distiller(palette = "Set1", direction = 1) +
#  expand_limits(y=0) + 
#  theme_minimal() +
#  facet_wrap(
#    vars(config_id),
#    labeller=labeller(config_id = label_from_config_id)
#  )



#ggsave(file.path(output_dir, current_run_id, "allele_counts_step.png"), plot = run_allele_frequencies_plt, width = 2000, height = 1000, units = "px")


#esquisse::esquisser(run_allele_counts_grouped)

## Reading data from config, config id 0 in this case
# config0_path = file.path(current_run_path, "config-0")
# 
# config0_config = read_config_config(config0_path)
# config0_data = read_config_data(config0_path)
# 
# # get allele counts
# config0_data$host_allele_counts = config0_data$host_allele_data %>% count(generation, species)
# 
# config0_data$meta_data
# View(config0_data$host_allele_data)
# esquisse::esquisser(config0_data$host_allele_counts)
