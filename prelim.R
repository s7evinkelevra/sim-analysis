library(tools)
library(rjson)

library(tidyverse)

library(esquisse)

data_dir = "./data"
current_run_id = "run-2022-06-15-2"
current_run_path = file.path(data_dir, current_run_id)


# nice name lol
read_config_config <- function(config_path) {
  return(fromJSON(file=file.path(config_path, "config.json")))
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

read_run_data_data <- function(run_path,csv_name) {
  data = list()
  data_folders = list.dirs(run_path)[-1]
  data_name = tools::file_path_sans_ext(csv_name)
  
  for(data_folder in data_folders) {
    current_file_path = file.path(data_folder, csv_name)
    current_data = read_delim(current_file_path, delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
    data[[basename(data_folder)]] = current_data
  }
  
  return(data)
}


run_configs = read_run_configs(current_run_path)
run_host_allele_data = read_run_data_data(current_run_path, "host_allele_data.csv")
run_host_allele_data_combined = bind_rows(run_host_allele_data, .id = "config_id")

run_allele_counts = run_host_allele_data_combined %>% count(config_id, generation, species)
run_allele_counts_grouped = run_allele_counts %>% group_by(config_id)

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
