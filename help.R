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

add_sim_mode <- function(meta_data_tibble) {
  return(
    meta_data_tibble %>% 
      mutate(bInitialMode =
               bInfection == 0 & 
               bHostFitnessproportionalReproduction == 0 &
               bPathogenFitnessproportionalReproduction == 0 &
               bPathogenMutation == 0 &
               bHostMutation == 0) %>%
      mutate(bBurninMode =
               generation == 0 |
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
      ) %>%
      mutate(derived_sim_mode = case_when(
        bBurninMode ~ "Neutrality",
        bCoevolutionMode ~ "Coevolution",
        bNoCoevolutionMode ~ "No-coevolution"
      ))
  )
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
