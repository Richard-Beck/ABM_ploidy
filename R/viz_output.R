

# Load necessary libraries
library(dplyr)
library(stringr)

main_dir <- "/Users/4477116/Documents/projects/ABM_ploidy/output/test_diffusion/data"
setwd("/Users/4477116/Documents/projects/ABM_ploidy/output/run_6")
#x <- read.csv("summary.csv")
o <- read.csv("oxygen.csv")
o$day <- round(o$time)

smm_o <- aggregate(list(nsteps=o$day),by=list(day=o$day),length)

p <- ggplot(o[498<o$time &o$time<499,],aes(x=time,y=O2))+
  geom_line()+
  scale_y_log10()+
  labs(x="time",y="Avg O2")+
  ylim(min(o[498<o$time &o$time<499,]$O2),min(o[498<o$time &o$time<499,]$O2)+.0005)
p
p <- ggplot(o,aes(x=time,y=O2))+
  geom_line()
p
p <- ggplot(smm_o,aes(x=day,y=nsteps))+
  geom_line()+
  scale_y_log10()
p

library(ggplot2)
#p <- ggplot(x,aes(x=Time,y=NCells_R))+
 # geom_line()
#p







#Changing Consumptin rates to see it's impact
# Get a list of all subdirectories in the main directory
subdirs <- list.dirs(main_dir, recursive = FALSE)

#Which variable what the seep run for?
config_file_path <-paste0(dirname(main_dir), "/test_homeostasis/config.txt") # Replace with the  path to config file
config_content <- readLines(config_file_path)
sweeptype_line_index <- which(str_detect(config_content, "sweeptype="))

parameter_line <- config_content[sweeptype_line_index + 1]
parameter_name <- str_extract(parameter_line, "^[^=]+")
#dynamic_value <- str_extract(parameter_line, "(?<=\\=).*")




# Initialize a data frame to store results
o_data <- data.frame()
summary<-data.frame()

# Loop through each subdirectory
for (subdir in subdirs) {
  # Construct the path to the params file and the oxygen.csv file
  params_file <- file.path(subdir, "params.txt")
  oxygen_file <- file.path(subdir, "oxygen.csv")
  summary_file<-file.path(subdir, "summary.csv")
  
  
 
  # Check if the params file exists
  if (file.exists(params_file)) {
    params_filename <- basename(params_file)  # Get just the file name from the path
    extraction_pattern <- parameter_name 
    
    # Read the params file as text
    params_content <- readLines(params_file)
    relevant_line <- params_content[str_detect(params_content, extraction_pattern)]
    
    # Extract the numeric value from the relevant line
    extracted_value <- as.numeric(str_extract(relevant_line, "\\d+\\.?\\d*"))
    
    
  }
  

  # Read the oxygen.csv file
  if (file.exists(oxygen_file)) {
    oxygen_data <- read.csv(oxygen_file)
    
    # Add a new column to identify the data
   oxygen_data <- mutate(oxygen_data, 
                         !!extraction_pattern:= extracted_value) # load dplyr library here
   
                          
    # Add the 'day' column by rounding the 'time'
    oxygen_data$day <- round(oxygen_data$time)
    
    # Summarize the data: count the number of steps (nsteps) for each day
    summarized_oxygen <- aggregate(list(nsteps = oxygen_data$day), by = list(day = oxygen_data$day), length)
    
    # Add the summarized data back as a separate column to the original dataset (if needed)
    oxygen_data <- left_join(oxygen_data, summarized_oxygen, by = "day")
    
    
    # Combine this data with the rest
    o_data <- bind_rows(o_data, oxygen_data)
  }
  
  if (file.exists(summary_file)) {
    summary_data <- read.csv(summary_file)
    summary_data <- mutate(summary_data, 
                           !!extraction_pattern:= extracted_value)
    summary_data$total_pop<- summary_data$nNormal+ summary_data$nCancer
    
    
   summary<-bind_rows(summary, summary_data)
  }
  
}
  

head(o_data)

o_data <- o_data %>%
  drop_na(time, maxDelta, all_of(parameter_name))


setLegendTitle=function(x){
  if(parameter_name=="normal_oxygen_consumption_rate"){
    legend_title=expression(O[2]~"Cons.")
  }else if(parameter_name=="initialTumorSize"){
    legend_title="init. Size"
  }
}
legend_title <- setLegendTitle()




# Create the ggplot
d<-ggplot(o_data[o_data$time>1&o_data$time<2 ,], aes(x = time, y = maxDelta, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1) +
  scale_y_log10() +  # Use a logarithmic scale for maxDelta
  scale_color_viridis_d(name =legend_title) +  # Use a continuous color gradient for clarity
  labs(x = "time (seconds)", y = "max delta", title = "Max delta") +
  theme(axis.text.x.top = element_text(hjust = 0.5),
        text=element_text(size=18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))
d

# Create the ggplot
o<-ggplot(o_data, aes(x = time, y = O2, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1.5) +
 #scale_y_log10() +  # Use a logarithmic scale for maxDelta
  scale_color_viridis_d(name = legend_title)+  # Use a continuous color gradient for clarity
  labs(x = "time (days)", y = "Average O2 (mmol/L)", title = "Average Oxygen") +
  theme(axis.text.x.top = element_text(hjust = 0.5),text=element_text(size=18))
o


plots_list <- list()  # Create an empty list to store individual plots

# Loop through time intervals from 1 to 10
for (i in 1:10:100) {
  # Create the ggplot for each interval (i < time <= i+1)
  p <- ggplot(o_data[o_data$time > i & o_data$time <= (i + 1), ], 
              aes(x = time, y = maxDelta, color = as.factor(.data[[parameter_name]]))) +
    geom_line(linewidth = 1) +
    scale_y_log10() +  # Use a logarithmic scale for maxDelta
    scale_color_viridis_d(name = legend_title) +  # Use a continuous color gradient for clarity
    labs(x = "time (seconds)", y = "max delta", title = paste("Max delta for", i, "< time <=", i + 1)) +
    theme(axis.text.x.top = element_text(hjust = 0.5))
  
  # Store the plot in the list
  plots_list[[i]] <- p
}

# Arrange all plots together using ggarrange from ggpubr
o_plot <- ggarrange(plotlist = plots_list, 
                           ncol = 2, nrow = 5)  # Arrange in a 2x5 grid
# Display the combined plot
print(o_plot)

   

#Count nsteps per timepoint

count_df <- o_data %>%
  mutate(time = floor(time)) %>%  # Group time into intervals (e.g., t = 0 for [0, 1), 1 for [1, 2), etc.)
  group_by(!!sym(parameter_name), time) %>%  # Dynamically group by the parameter_name
  summarize(nsteps = first(nsteps),  # Select first nsteps for each interval
            !!sym(parameter_name) := first(!!sym(parameter_name))) %>%  # Retain constant ConsumptionRate for each block
  ungroup() %>%
  select(time, nsteps, !!sym(parameter_name))  # Select relevant columns, dynamically include parameter_name


np<-ggplot(count_df %>% filter(time< 20),aes(x = time, y = nsteps, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1) +
  scale_y_log10() + 
  scale_color_viridis_d(name = legend_title) +  # Use a continuous color gradient for clarity
  labs(x = "time (days)", y = "# of diff. steps", title = "Diffusion steps ") 
np


  
n_plot<-ggplot(summary,aes(x = time, y = nNormal, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1.5) +
  #scale_y_log10() +  # Use a logarithmic scale for maxDelta
  scale_color_viridis_d(name = legend_title) +  # Use a continuous color gradient for clarity
  labs(x = "time (days)", y = "Population", title = "Normal population ") +
  theme(axis.text.x.top = element_text(hjust = 0.5),
        text=element_text(size=18))

n_plot

c_plot<-ggplot(summary,aes(x = time, y = nCancer, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1) +
  #scale_y_log10() +  # Use a logarithmic scale for maxDelta
  scale_color_viridis_d(name = legend_title) +  # Use a continuous color gradient for clarity
  labs(x = "time (days)", y = "Caner Size", title = "Cancer Population",color="con") +
  theme(axis.text.x.top = element_text(hjust = 0.5),
        text=element_text(size=18))

c_plot


p<-ggplot(summary,aes(x = time, y = total_pop, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1.5) +
  scale_y_log10() +  # Use a logarithmic scale for maxDelta
  scale_color_viridis_d(name = legend_title) +  # Use a continuous color gradient for clarity
  labs(x = "time (days)", y = "Size", title = "Cell Population") +
  theme(axis.text.x.top = element_text(hjust = 0.5),
        text=element_text(size=18))

p






#Ignore

#Overlay plot of homeostatic times
x_intercepts <- c(42, 123, 270, 375, 423, 450) #Times it too for the various setups to reach homeostasis
final_time_point <- max(summary$time)  # Get the final time point from the data

# Generate the base plot
p <- ggplot(summary, aes(x = time, y = total_pop, color = as.factor(.data[[parameter_name]]))) +
  geom_line(linewidth = 1.5) +
  scale_y_log10() +  # Use a logarithmic scale for total_pop
  scale_color_viridis_d(name = legend_title) +  # Color scale for the lines
  labs(x = "time (days)", y = "Size", title = "Cell Population") +
  theme(axis.text.x.top = element_text(hjust = 0.5), text = element_text(size = 18)) +
  
  # Set custom x-axis ticks (0, x_intercepts, final time point)
  scale_x_continuous(breaks = c(0, x_intercepts, final_time_point))

# Loop through the x_intercepts and add geom_vline for each with matching colors
for (i in seq_along(x_intercepts)) {
  p <- p + geom_vline(xintercept = x_intercepts[i], 
                      linetype = "dashed", size = 1, 
                      color = scales::viridis_pal()(length(unique(summary[[parameter_name]])))[i])
}

# Display the plot
p