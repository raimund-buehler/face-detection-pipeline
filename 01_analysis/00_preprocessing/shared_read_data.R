# Manuscript section: Preprocessing utility for reading raw eye-tracking exports
# Analysis family: preprocessing helper
# Original source path: scripts/read_data.R
# Primary input dataset(s): pupil_positions.csv; result.csv; annotations.csv
# Primary output(s): in-memory per-session raw tables and parsed IDs/session labels
# Known TODOs: parameterize path assumptions; document required export folder structure
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)

read_data <- function(time_path, result_path, marker_path) {
  #Extract Sub_ID from path (using regular expression)
  Sub_ID <- sub(".*/recordings_sorted/([^/]+)/.*", "\\1", result_path)
  session <- sub(".*/(Session \\d+)/.*", "\\1", result_path)
  
  print(paste("Sub_ID: ", Sub_ID, ", Session: ", session))
  
  #Read in Data
  time <- read_csv(time_path)
  result <- read_csv(result_path)
  marker <- read_csv(marker_path, col_types = cols("label" = col_character())) %>% select(index, label) %>% rename("frame_id" = index)
  
  return(list(time = time, result = result, marker = marker, Sub_ID = Sub_ID, session = session))
}
