# nice name lol
# prelim.R related helper
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



# thes.R related helper
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
    run_configs_tibble = generate_run_config_summary(id, configs) %>% rownames_to_column(var = "config_id")
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
    data_combined = bind_rows(data, .id = "config_id") %>% filter(generation >= min_generation & generation <= max_generation)
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

add_meta_sim_mode_analysis <- function(left, meta_data) {
  meta_data_sim_mode = meta_data %>% select(run_id, config_id, generation, derived_sim_mode)
  return(
    left_join(left, meta_data_sim_mode, by = c("run_id", "config_id", "generation"))
  )
}


add_run_config_analysis <- function(left, analysis_configs) {
  return(
    left_join(left, analysis_configs, by = c("run_id"))
  )
}

add_run_config_run_id_same_config_different_mode_analysis <- function(left, analysis_configs) {
  analysis_configs_reduced = analysis_configs %>% select(run_id, run_id_same_config_different_mode)
  return(
    left_join(left, analysis_configs_reduced, by = c("run_id"))
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

label_from_id <- function(id) {
  row_id_analysis_configs = analysis_configs_merged_unique_changed %>%
    column_to_rownames("id_same_config_different_mode") %>%
    select((!run_id_same_config_different_mode & !hash_no_sim_mode))
  
  label = ""
  
  for(col_i in 1:ncol(row_id_analysis_configs[id,])) {
    label = paste(label, colnames(row_id_analysis_configs[id,])[col_i], "=", row_id_analysis_configs[id, col_i], "\n")
  }
  
  return(label)
}


theme_Publication <- function(base_size=14, base_family="helvetica", legend_size_cm = 0.2, legend_direction = "horizontal") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = legend_direction,
            legend.key.size= unit(legend_size_cm, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

theme_Publication_Legend_Side <- function(base_size=14, base_family="helvetica", legend_size_cm = 0.2) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(legend_size_cm, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
