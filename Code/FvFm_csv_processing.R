# R script to process csv files output from ImagingWin software reports from Walz Maxi-Imaging-PAM.
# Authors: Sabina Aitken, Alicia Cook, Pieter Arnold

# Aim: batch process to align csv files for one day of experiment at a time, assuming that each day has 
# samples named consistently: Temp_Time_Measurement, e.g. 40_15_i for 40˚C, 15-min, initial Fv/Fm.
# AOI = Area Of Interest (circle placed on a sample from which Fv/Fm is estimated)

library(tidyr)
library(dplyr)
library(stringr)
library(purrr)

# Below creates a function to read and clean each CSV file which can then be applied to the whole list of filenames 
clean_csv <- function(file) {
  data <- read.csv(file, sep = ";")
  cleaned_data <- as.data.frame(t(data)) %>%
    slice(-(1:4), -53) %>% # CHANGE last row value depending on number of areas of interest (samples)
    # Remove first (usually) four rows and last row that don't contain Fv/Fm values
    mutate(filename = rep(basename(file), each = nrow(data))) %>% 
    # add column for filename which has format: Temp_Time_Measurement
    separate(filename, into = c("Temp", "Time", "Measurement"), sep = "_") %>% 
    # Separate filename column in to 3 columns for each of time, temp and measurement (initial, post or final).
    mutate(Measurement = str_remove(Measurement, ".csv")) %>%  # clean measurement names
    mutate(AOI = rownames(.)) %>% # Make column of AOI number
    mutate(AOI = str_remove(AOI, "Y.II.")) %>% # clean up AOI number
    rename(FvFm = V1)
  
  return(cleaned_data)
}

setwd("your_directory") # CHANGE WORKING DIRECTORY OR SPECIFY FULL FILEPATH BELOW

# CHANGE FILEPATHS
# Read in 1 csv to check no. rows etc. 
csv <- read.csv("Exp 1 - CSVs Day 1/30_60_final.csv", sep = ";") 
file_path <- "Exp 1 - CSVs Day 1/30_60_final.csv"
filename <- tools::file_path_sans_ext(basename(file_path))

csv <- as.data.frame(t(csv)) %>%
  slice(-(1:4), -53) %>% # Remove first four rows (Date, Time, No., PAR) and last row (X) that don't contain Fv/Fm values
  mutate(filename = rep(filename, each = nrow(.))) %>% # Rename each row
  separate(filename, into = c("Temp", "Time", "Measurement"), sep="_") %>% # Separate filename column in to 3 columns for each of time, temp and measurement.
  mutate(Measurement = str_remove(Measurement, ".csv")) %>%
  mutate(AOI = rownames(.)) %>%
  mutate(AOI = str_remove(AOI, "Y.II."))

head(csv)

## READ ID FILE FOR THE EXPERIMENT THAT SHOWS WHICH SAMPLES ARE WHICH AOI
ID <- read.csv("Exp 1 ID.csv") # CHANGE FILEPATH
head(ID)

## READ AND CLEAN ALL CSV FILES IN FOLDER
folder_path <- "Exp 1 - CSVs Day 1" # CHANGE FILEPATH

file_names <- list.files(folder_path, full.names = TRUE) # Get a list of all CSV files in the folder
file_names

day1csv <- map_dfr(file_names, clean_csv) # apply above function called "clean_csv" to all file_names in list and build one large dataframe

day1csvID <- merge(day1csv, ID, by = "AOI") # merge function based on AOI column being common to both dataframes

day1csvID <- day1csvID %>%
  mutate(Day = rep(1, nrow(day1csvID)))

## IF NEEDED (MULTIPLE DAYS OF EXPERIMENT)
# Exp1 Day 2
# CHANGE FILEPATH
folder_path <- "Exp 1 - CSVs Day 2" 

file_names <- list.files(folder_path, full.names = TRUE) # Get a list of all CSV files in the folder

day2csv <- map_dfr(file_names, clean_csv)
day2csvID <- merge(day2csv, ID, by = "AOI")

day2csvID <- day2csvID %>%
  mutate(Day = rep(2, nrow(day2csvID))) # add column to denote two different days of experiment

Exp1_all_long <- rbind(day1csvID, day2csvID) # combine both days into single dataframe

# write.csv(Exp1_all_long, "Exp1_all_long.csv")

Exp1_all_wide <- Exp1_all_long %>%
  pivot_wider(names_from = Measurement, values_from = FvFm) # Transform data to wide format with separate columns for each of initial, post and final Fv/Fm
head(Exp1_all_wide)

# write.csv(Exp1_all_wide, "Exp1_all_wide.csv")

