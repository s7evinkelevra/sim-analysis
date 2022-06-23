library(tools)
library(rjson)
library(jsonlite)

library(tidyverse)

library(esquisse)

data_dir = "./data"
current_run_id = "run-2022-06-17-1"
current_run_path = file.path(data_dir, current_run_id)

output_dir = "./output"

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
  
  for(col_i in 1:ncol(run_config_changed[config_id,])) {
    label = paste(label, colnames(run_config_changed[config_id,])[col_i], "=", run_config_changed[config_id, col_i], "\n")
  }
  
  return(label)
}

property_from_config_id <- function(config_id, property_name){
  return(run_config_summary[config_id, property_name])
}


# Read run data
run_configs = read_run_configs(current_run_path)
run_host_allele_data = read_run_data(current_run_path, "host_allele_data.csv")
run_host_allele_data_combined = bind_rows(run_host_allele_data, .id = "config_id")

# Create output folder for this run
dir.create(file.path(output_dir, current_run_id), showWarnings = FALSE)

# Write summary of the configs to csv
run_config_summary = generate_run_config_summary(current_run_id, run_configs)
write.csv(run_config_summary, file.path(output_dir, current_run_id, "config_summary.csv"))


# Get all config vars that have been changed throughout the run
run_config_changed_mask = run_config_summary %>% sapply(function(x) !length(unique(x)) == 1)
run_config_changed = run_config_summary[, run_config_changed_mask]



run_allele_counts = run_host_allele_data_combined %>% count(config_id, generation, species)
run_allele_counts_grouped = run_allele_counts %>% group_by(config_id)



run_allele_counts_plt = ggplot(run_allele_counts_grouped) +
  aes(
    x = generation,
    y = n,
    color = species,
    group = species
  ) +
  geom_line(size = 0.5) +
  scale_color_distiller(palette = "Set1", direction = 1) +
  theme_minimal() +
  facet_wrap(
    vars(config_id),
    labeller=labeller(config_id = label_from_config_id)
  )


run_allele_counts_plt
ggsave(file.path(output_dir, current_run_id, "allele_counts.png"), plot = run_allele_counts_plt, width = 5000, height = 5000, units = "px")

run_allele_frequencies_plt = ggplot(run_host_allele_data_combined) +
  aes(
    x = generation,
    y = count,
    #fill = allele_id,
    color = created_at,
    group = species
  ) +
  #geom_bar(position = "stack",stat = "identity") +
  geom_step(size = 0.5) +
  scale_color_distiller(palette = "Set1", direction = 1) +
  theme_minimal() +
  facet_wrap(
    vars(config_id),
    labeller=labeller(config_id = label_from_config_id)
  )

ggsave(file.path(output_dir, current_run_id, "allele_counts_step.png"), plot = run_allele_frequencies_plt, width = 2000, height = 1000, units = "px")


esquisse::esquisser(run_allele_counts_grouped)

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
