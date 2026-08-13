# Manuscript section: Preprocessing utility for per-subject timecourse plotting
# Analysis family: preprocessing helper
# Original source path: scripts/plot_sub.R
# Primary input dataset(s): transformed per-session eye-tracking fixation data
# Primary output(s): per-subject timecourse plots
# Known TODOs: clarify whether these plots are exploratory only or manuscript-adjacent QC outputs
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(zoo)
library(generics)

plot_sub <- function(Sub_ID, session, result_fix) {
  #Marker start and end times
  start_end_times <- result_fix %>% 
    group_by(label) %>%
    summarize(start_time = first(pupil_timestamp), end_time = last(pupil_timestamp))
  
  #Time interval for bins in which rate in Hz is calculated
  time_interval_size <- 1
  
  #Dataframe for plotting
  result_plot <- result_fix %>% ungroup() %>% 
    pivot_longer(cols = fix_on_face:fix_on_background, 
                 names_to = "fix", values_to = "value") %>%
    mutate(fix = as.integer(factor(fix, levels = unique(fix)))) %>% 
    group_by(fix) %>% 
    mutate(time_interval = 
             floor(pupil_timestamp / time_interval_size) * time_interval_size) 
    #mutate(rolling_average = rollapply(value, width = 25, mean, na.rm = T, fill = NA, align = "right")) %>% 
    #mutate(change_rate = c(rep(NA, 2), diff(rolling_average, lag = 2)))
  
    plot_rate <- result_plot %>% ungroup() %>% 
      group_by(time_interval, fix) %>%
      distinct(id_fix, .keep_all = T) %>% 
      mutate(rate = sum(value, na.rm = T) / time_interval_size)
    
    plot_roll <- result_plot %>% ungroup() %>%
      group_by(time_interval, fix) %>% 
      distinct(id_fix, fix, .keep_all = T) %>% 
      group_by(fix) %>% 
      mutate(rolling_average = rollapply(value, width = 5, mean, na.rm = T, fill = NA, align = "right"))
      #mutate(change_rate = c(rep(NA, 2), diff(rolling_average, lag = 2)))
  
  
  #Plot rate (Hz)
   p <- ggplot(data = result_plot, 
         aes(x = pupil_timestamp, y = fix)) + 
    # geom_tile(data = result_plot %>% filter(value == T) %>% distinct(id_fix, .keep_all = T), aes(height = 0.25, width = 0.1)) +
    geom_rect(data = start_end_times, 
              aes(xmin = start_time, 
                  xmax = c(lead(start_time)[-length(start_time)], tail(end_time, 1)),
                  ymin = 0, ymax = Inf, 
                  fill = label, color = label),
              alpha = 0.25,
              inherit.aes = FALSE) +
    geom_step(data = plot_rate %>% 
                filter(fix == 1), aes(y = rate/10 + 1)) +
    geom_step(data = plot_rate %>% 
                filter(fix == 2), aes(y = rate/10 + 2)) +
    geom_step(data = plot_rate %>% 
                filter(fix == 3), aes(y = rate/10 + 3)) +
    geom_step(data = plot_rate %>% 
                filter(fix == 4), aes(y = rate/10 + 4)) +
    geom_hline(yintercept = c(1:4)) +
    scale_y_continuous(breaks = c(1:4), limits = c(0, 5),
                       labels = c("Face", "Eyes", "Mouth", "Background"),
                       sec.axis = sec_axis(~., breaks = c(1.5:4.5), 
                                           labels = rep(c(5), 4))) +
    labs(x = "Time (s)", y = "Fixation Rate (Hz)", 
         fill = "Questions", color = "Questions", 
         title = paste("Fixation Timecourse ", Sub_ID, " ", session)) +
    theme_classic() +
    theme(plot.title = element_text(face="bold", hjust = 0.5, size = 20),
          legend.position = "none")

   p
   
  return(p)
}
