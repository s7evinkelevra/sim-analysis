library(tools)
library(rjson)
library(jsonlite)

library(tidyverse)
library(plotly)
library(Cairo)

library(esquisse)

# nice name lol
read_config_config <- function(config_path) {
  return(rjson::fromJSON(file=file.path(config_path, "config.json")))
}

read_config_data <- function(config_path) {
  csv_files = list.files(path=config_path, pattern="\\.csv$")
  csv_data = list()
  
  for(filename in csv_files){
    name = tools::file_path_sans_ext(filename)
    csv_data[[name]] = read_delim(file.path(config_path, filename), 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  }
  return(csv_data)
}

read_run_configs <- function(run_path) {
  configs = list()
  config_folders = list.dirs(run_path)[-1]
  
  for(config_folder in config_folders) {
    current_config = read_config_config(config_folder)
    configs[[basename(config_folder)]] = current_config
  }
  
  return(configs)
}

read_run_data <- function(run_path,csv_name) {
  data = list()
  data_folders = list.dirs(run_path)[-1]
  data_name = tools::file_path_sans_ext(csv_name)
  
  for(data_folder in data_folders) {
    current_data_file_path = file.path(data_folder, csv_name)
    current_data = read_delim(current_data_file_path, delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
    current_config = read_config_config(data_folder)
    data[[as.character(current_config$configId)]] = current_data
  }
  
  return(data)
}

generate_run_config_summary <- function(run_id, run_configs) {
  config_tibble = NULL
  
  for(config in run_configs){
    config["run_id"] = run_id
    if(is.null(config_tibble)) {
      config_tibble = config %>% unlist() %>% enframe() %>% pivot_wider() 
    }else{
      config_tibble = add_row(config_tibble, config %>% unlist() %>% enframe() %>% pivot_wider())
    }
  }
  config_tibble = config_tibble %>% column_to_rownames(var="configId")
  return(config_tibble)
}


label_from_config_id <- function(config_id){
  label = ""
  
  if(ncol(run_config_changed) == 0){
    return(paste0("No changes between configs\nconfig-",config_id))
  }
  
  for(col_i in 1:ncol(run_config_changed[config_id,])) {
    label = paste(label, colnames(run_config_changed[config_id,])[col_i], "=", run_config_changed[config_id, col_i], "\n")
  }
  
  return(label)
}

property_from_config_id <- function(config_id, property_name){
  return(run_config_summary[config_id, property_name])
}

add_data <- function(base, adding){
  return(left_join(base, adding, by = c("config_id", "generation")))
}

save_plot_defaults <- function(path, plt, width, height){
  ggsave(path, plot = plt, width = width, height = height, dpi = "retina", device = "png", type="cairo", units = "px", bg = "#ffffff", limitsize = FALSE)
}


data_dir <- "./data"
current_run_id <- "2022-07-06-3"
current_run_path <- file.path(data_dir, current_run_id)

output_dir <- "./output"


# Read run data
run_configs = read_run_configs(current_run_path)

# read allele frequency data
run_host_allele_data = read_run_data(current_run_path, "host_allele_data.csv")
run_host_allele_data_combined = bind_rows(run_host_allele_data, .id = "config_id")

# read host data (infection/fitness data + ancestry)
run_host_data = read_run_data(current_run_path, "host_data.csv")
run_host_data_combined = bind_rows(run_host_data, .id = "config_id")

# read host genome data (allele ids of all host individuals)
run_host_genome_data = read_run_data(current_run_path, "host_genome_data.csv")
run_host_genome_data_combined = bind_rows(run_host_genome_data, .id = "config_id")

# read run meta data
run_meta_data = read_run_data(current_run_path, "meta_data.csv")
run_meta_data_combined = bind_rows(run_meta_data, .id="config_id")

# build mode info from flags
run_meta_data_combined = run_meta_data_combined %>% 
  mutate(bInitialMode =
    bInfection == 0 & 
    bHostFitnessproportionalReproduction == 0 &
    bPathogenFitnessproportionalReproduction == 0 &
    bPathogenMutation == 0 &
    bHostMutation == 0) %>%
  mutate(bBurninMode =
    generation == 0 |
    bInfection == 0 & 
    bHostFitnessproportionalReproduction == 0 &
    bPathogenFitnessproportionalReproduction == 0 &
    bPathogenMutation == 1 &
    bHostMutation == 1) %>%
  mutate(bCoevolutionMode = 
    generation != 0 &
    bInfection == 1 & 
    bHostFitnessproportionalReproduction == 1 &
    bPathogenFitnessproportionalReproduction == 1 &
    bPathogenMutation == 1 &
    bHostMutation == 1
  ) %>%
  mutate(bNoCoevolutionMode = 
    generation != 0 &
    bInfection == 1 & 
    bHostFitnessproportionalReproduction == 1 &
    bPathogenFitnessproportionalReproduction == 0 &
    bPathogenMutation == 1 &
    bHostMutation == 1
)

## Housekeeping
# Create output folder for this run
dir.create(file.path(output_dir, current_run_id), showWarnings = FALSE)

## Configs
# Write summary of the configs to csv
run_config_summary = generate_run_config_summary(current_run_id, run_configs)
write.csv(run_config_summary, file.path(output_dir, current_run_id, "config_summary.csv"))

# Get all config vars that have been changed throughout the run
run_config_changed_mask = run_config_summary %>% sapply(function(x) !length(unique(x)) == 1)
run_config_changed = run_config_summary[, run_config_changed_mask]
run_config_common = run_config_summary[1, !run_config_changed_mask]

# write common config values in long format to file for better use later
run_config_key_value = run_config_common %>% pivot_longer(cols = everything(), names_to = "property", values_to = "values")
write.csv(run_config_key_value, file.path(output_dir, current_run_id, "config_common.csv"), row.names = FALSE)


## Fitness data
run_host_data_combined_meta = add_data(run_host_data_combined, run_meta_data_combined) 

# successful antigen presentations in the last generation
run_host_data_last_gen_successful_presentation <- run_host_data_combined_meta %>%
  group_by(config_id, species, successful_presentations) %>%
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

#plotly::ggplotly(run_host_data_last_gen_successful_presentation_dist_plt)
save_plot_defaults(file.path(output_dir, current_run_id, "host_presentation_counts.png"), run_host_data_last_gen_successful_presentation_dist_plt, 3000, 1500)


# average fitness over time
run_host_data_fitness_mean_over_time <- run_host_data_combined_meta %>%
  filter(bBurninMode == FALSE) %>%
  group_by(config_id, generation, species) %>%
  summarize(mean_fitness = mean(fitness), sd_fitness = sd(fitness))
 
run_host_data_fitness_mean_over_time_plt <- ggplot(run_host_data_fitness_mean_over_time) +
  aes(
    x = generation,
    y = mean_fitness,
    group = config_id,
    color = config_id
  ) +
  geom_hline(yintercept = as.numeric(run_config_common[1,"hosts.fitness_minimum"])) +
  geom_point() +
  geom_line() +
  theme_minimal()


run_host_data_fitness_mean_over_time_plt
plotly::ggplotly(run_host_data_fitness_mean_over_time_plt)
save_plot_defaults(file.path(output_dir, current_run_id, "host_fitness_mean.png"), run_host_data_fitness_mean_over_time_plt, 3000, 1500)


## Heterozygosity
run_host_genome_data_combined_meta = add_data(run_host_genome_data_combined, run_meta_data_combined)


# add allele frequency data 
run_host_genome_data_combined_meta_count_freq = run_host_genome_data_combined_meta %>% 
  left_join(run_host_allele_data_combined, 
                 by = c("config_id", "generation", "species", "locus_id","allele_1_id" = "allele_id"), 
                 suffix = c(".allele_1", ".allele_2")) %>%
  left_join(run_host_allele_data_combined,
            by = c("config_id", "generation", "species", "locus_id","allele_2_id" = "allele_id"), 
            suffix = c(".allele_1", ".allele_2"))



run_host_genome_zygosity = run_host_genome_data_combined_meta %>%
  mutate(zygosity = case_when(
    allele_1_id == allele_2_id ~ "homozygous",
    allele_1_id != allele_2_id ~ "heterozygous"
  )) %>%
  group_by(config_id, generation, species, locus_id) %>%
  count(zygosity, name = "zygosity_count") %>%
  pivot_wider(names_from = zygosity, values_from = zygosity_count) %>%
  mutate(hobs = heterozygous/(homozygous+heterozygous))

run_host_genome_zygosity_plt = ggplot(run_host_genome_zygosity) +
  aes(
    x = generation,
    y = hobs,
    group = config_id,
    color = config_id
  ) +
  geom_line() + 
  geom_vline(xintercept = as.numeric(run_config_common[1,"burnin_generations"])) + 
  theme_minimal()

run_host_genome_zygosity_plt


# marginal fitness of alleles


## Allele counts
# get allele counts by config, generation and species
run_allele_counts = run_host_allele_data_combined %>% count(config_id, generation, species)
run_allele_counts_grouped = run_allele_counts %>% group_by(config_id)

# add meta data to allele counts -> simulation mode most importantly
run_allele_counts_grouped_meta = add_data(run_allele_counts_grouped, run_meta_data_combined)


# select the last 100 burnin and last 100 post-burnin generations
run_allele_counts_last_100 <- run_allele_counts_grouped_meta %>%
  group_by(config_id, bBurninMode) %>%
  filter(generation >= max(generation) - 100)


# build mean allele counts and SD by config and mode for the last 100 generations
run_allele_counts_last_100_mean <- run_allele_counts_last_100 %>%
  summarise(mean_alleles = mean(n), sd = sd(n), .groups = "keep")



## Allele frequency distribution
run_host_allele_data_combined_meta <- add_data(run_host_allele_data_combined, run_meta_data_combined)


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
run_allele_counts_last_100_box_plt = ggplot(run_allele_counts_last_100) +
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
run_allele_counts_plt = ggplot(run_allele_counts_grouped) +
  aes(
    x = generation,
    y = n,
    color = species,
    group = species
  ) +
  geom_line(size = 0.2) +
  expand_limits(y=0) + 
  scale_color_distiller(palette = "Set1", direction = 1) +
  theme_minimal() +
  facet_wrap(
    vars(config_id),
    labeller=labeller(config_id = label_from_config_id)
  )


run_allele_counts_plt
save_plot_defaults(file.path(output_dir, current_run_id, "allele_counts.png"), run_allele_counts_plt, 20000, 5000)

# run_allele_frequencies_plt = ggplot(run_host_allele_data_combined) +
#   aes(
#     x = generation,
#     y = count,
#     #fill = allele_id,
#     color = created_at,
#     group = species
#   ) +
#   #geom_bar(position = "stack",stat = "identity") +
#   geom_step(size = 0.5) +
#   scale_color_distiller(palette = "Set1", direction = 1) +
#   expand_limits(y=0) + 
#   theme_minimal() +
#   facet_wrap(
#     vars(config_id),
#     labeller=labeller(config_id = label_from_config_id)
#   )

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
