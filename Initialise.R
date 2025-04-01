#!/bin/bash Rscript

################## Package Management ##################
# Check if required packages are installed and install if missing
required_packages <- c("stringr", "tidyr", "dplyr", "ggplot2", "ggpubr", "ez", "rstatix",
                       "psychReport")
to_be_installed <- required_packages[!(required_packages %in% rownames(installed.packages()))]
if(length(to_be_installed) > 0){
  for(package in to_be_installed){
    install.packages(package)
  }
}

################## Library Imports ##################
# Load all required libraries for data analysis and visualization
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ez)
library(rstatix)
library(psychReport)

################## Path Setup [UPDATE THESE FOR YOUR LOCAL MACHINE] ##################
# Define base paths for scripts and data
scripts_path <- "~/Downloads/CategorySpecificAssociativeInference/"

# Change base_path to the folder that contains the downloaded data from OpenNeuro
# All other paths are defined relative to base_path
base_path <- "~/Downloads/"
data_path <- paste0(base_path, "CategorySpecificAssociativeInference/")
derivatives_path <- paste0(data_path, "derivatives/")

################## Participant Management ##################
# Read and filter participant information based on task completion
participants <- read.table(paste0(data_path,
                                  "/participants.tsv"),
                           header = TRUE)
aim_participants <- participants %>% 
  filter(AIM == "Completed") %>% 
  pull(participant_id)
loc_participants <- participants %>% 
  filter(Localiser == "Completed") %>% 
  pull(participant_id)

################## Factor Level Definitions ##################
# Define factor levels and labels for different experimental conditions
# Each list element contains standardized levels, display labels, and (where applicable) color codes
factor_levels <- list(category = list(levels = c("scene", "face"),
                                      labels = c("Scene", "Face"),
                                      colours = c("#66CC00", "#0099FF")),
                      pairtype_study = list(levels = c("AB", "BC"),
                                             labels=c("AB", "BC")),
                      pairtype_test = list(levels = c("AB", "BC", "AC"),
                                            labels=c("AB", "BC", "AC")),
                      trialtype_test = list(levels = c("dir", "indir"),
                                            labels=c("Direct", "Indirect")),
                      category_localiser = list(levels = c("scene", "face", "obj"),
                                                labels = c("Scene", "Face", "Object"),
                                                colours = c("#66CC00", "#0099FF", "#ffd100")),
                      studyitem_pairtype = list(levels = c("Unique_AB", "Linking_AB", 
                                                           "Linking_BC", "Unique_BC"),
                                                labels = c("A\n(object)", "B\n(scene/face)", 
                                                           "B\n(scene/face)", "C\n(object)")),
                      rois = list(levels = c("Ant Hipp Head", "CA1", "CA23DG", "Sub", "Hipp Tail", 
                                             "alERC", "pmERC", "PRC", "PHC",
                                             "FFA"),
                                  roi_type = c(rep_len("HC", 5),
                                               rep_len("MTL", 4),
                                               "FFA")),
                      hemisphere = list(levels = c("Left", "Right", "Collapsed"),
                                        labels = c("Left", "Right", "Bilateral")),
                      memory_test = list(levels = c("Correct", "Incorrect"),
                                        labels = c("Correct", "Incorrect")))

################## Utility Functions ##################
# Function to correct p-values using different methods
correct_p_vals <- function(df, pvalcol, usemethod = "fdr"){
  df <- as.data.frame(df)
  #Is raw p-value significant?
  df[, "sig"] <- ifelse(df[, pvalcol] <= 0.05, "*", "")
  
  #Correct p's based on method entered
  df[, paste(pvalcol, "_", usemethod, "corrected", sep = "")] <- p.adjust(df[, pvalcol], method = usemethod)
  #Is corrected p significant?
  df[, paste("sig_", usemethod, "corrected", sep = "")] <- ifelse(df[, paste(pvalcol, "_", usemethod, "corrected", sep = "")] <= 0.05, "*", "")
  
  return(df)
}

# Function to calculate summary statistics for a given column
summarise_data <- function(df, col_name, rm_na = FALSE){
  df <- as.data.frame(df)
  M <- mean(df[,col_name], na.rm = rm_na)
  SD <- sd(df[,col_name], na.rm = rm_na)
  SE <- SD / sqrt(nrow(df))
  LCI <- M - 1.96*SE
  HCI <- M + 1.96*SE
  MeanPlusSE <- M+SE
  MeanMinusSE <- M-SE
  NumOfRows <- nrow(df)
  data.frame(Mean=M, SD=SD, SE=SE, LCI=LCI, HCI=HCI, MeanPlusSE=MeanPlusSE, MeanMinusSE, Rows=NumOfRows)
}

################## Plotting Themes ##################
# Define consistent themes for paper-quality visualizations
# Theme for facet labels and panel borders
paper_facet_theme <- theme(strip.text.x = element_text(size = 22, colour = "black"),
                           strip.text.y = element_text(size = 22, colour = "black"), 
                           strip.background = element_rect(color="white", fill="white", linewidth=1.5, linetype="solid"),
                           panel.border = element_rect(colour = "black", fill = NA, linewidth=1.5))
# Theme for x-axis formatting
x_axis_theme <- theme(axis.title.x = element_blank(), 
                      axis.text.x = element_text(colour = "#000000", size = 18)) 
# Theme for y-axis formatting
y_axis_theme <-   theme(axis.title.y = element_text(colour = "#000000", size = 22), 
                        axis.text.y = element_text(colour = "#000000", size = 18))
# Theme for clean background without grid lines
blank_bg_theme <- theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"))
# Theme for legend formatting
legend_theme <- theme(legend.text = element_text(face = "bold", size = 20), 
                      legend.title = element_text(face = "bold", size = 25))
