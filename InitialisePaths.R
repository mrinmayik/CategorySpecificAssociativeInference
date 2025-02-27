#!/bin/bash Rscript

# Check if packages need to be installed
required_packages <- c("stringr", "tidyr", "dplyr", "ggplot2", "ggpubr", "ez", "rstatix")
to_be_installed <- required_packages[!(required_packages %in% rownames(installed.packages()))]
if(length(to_be_installed) > 0){
  for(package in to_be_installed){
    install.packages(package)
  }
}

library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ez)
library(rstatix)

# Setup paths where scripts are
scripts_path <- "~/GitDir/CodeWithPapers/CategorySpecificAssociativeInference"

# Initialise paths for the local machine
base_path <- "~/Desktop/AssociativeInference/"
data_path <- paste0(base_path, "Defaced/")
derivatives_path <- paste0(data_path, "derivatives/")
deconvolve_path <- paste0(data_path, "DeconvolveOutput/")

# Get participants information
participants <- read.table(paste0(data_path,
                                  "/participants.tsv"),
                           header = TRUE)
aim_participants <- participants %>% 
  filter(AIM == "Completed") %>% 
  pull(subj_id)
loc_participants <- participants %>% 
  filter(Localiser == "Completed") %>% 
  pull(subj_id)

# Function to correct for p-values
correct_p_vals <- function(df, pvalcol, usemethod="fdr"){
  #Is raw p-value significant?
  df[, "sig"] <- ifelse(df[, pvalcol]<=0.05, "*", "")
  
  #Correct p's based on method entered
  df[, paste(pvalcol, "_", usemethod, "corrected", sep="")] <- p.adjust(df[, pvalcol], method=usemethod)
  #Is corrected p significant?
  df[, paste("sig_", usemethod, "corrected", sep="")] <- ifelse(df[, paste(pvalcol, "_", usemethod, "corrected", sep="")]<=0.05, "*", "")
  
  return(df)
}

