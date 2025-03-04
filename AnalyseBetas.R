# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/CategorySpecificAssociativeInference/")
source("Initialise.R")

mtl_rois <- factor_levels$rois$levels

########################### STEP 1: Localiser ###########################

################## Read in betas and organise them ##################

# Initialize empty dataframe to store all participants' data
locbetas_data <- c()

# Loop through all participants and read their task event files
for(part in loc_participants){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(deconvolve_path, part, "/", part, "_task-local_betas.tsv"),
                          header=TRUE)
  # Add participant ID column
  part_data$Participant <- part
  # Combine with existing data
  locbetas_data <- rbind(locbetas_data, part_data)
}

# Reorder columns to put Participant ID first
locbetas_data <- relocate(locbetas_data,
                          Participant)

# Verify all intended participants are included in the dataset
all(sort(unique(locbetas_data$Participant)) == sort(loc_participants))


################## Analyse betas ##################
beta_anova <- c()
beta_anova_detailed <- c()

#### FIRST: Run ANOVA ####
for(roi in unique(locbetas_data$Area)){
  beta_roi_data <- locbetas_data[locbetas_data$Area == roi, ]
  
  beta_anova_roi <- ezANOVA(
    data = beta_roi_data, 
    wid = Participant, 
    within = c(Category, Hemisphere),
    dv = Beta, 
    detailed = TRUE
  )
  beta_anova_df <- aovEffectSize(beta_anova_roi)$ANOVA
  
  beta_anova_detailed[[roi]] <- beta_anova_roi
  beta_anova <- rbind(beta_anova, 
                      cbind(data.frame(Area = roi, 
                                     Event = "Localiser"), 
                           beta_anova_df))
}

#Only keep handsegmented ROIs so we can correct for them
beta_anova_mtl <- beta_anova %>% 
  filter(Area %in% mtl_rois,
         Effect != "(Intercept)")
beta_anova_mtl <- beta_anova_mtl %>% 
  group_by(Effect) %>% 
  do(correct_p_vals(., "p")
  )

####SECOND run Posthocs: ####
beta_plots <- list()
beta_posthoc <- c()
for(roi in unique(locbetas_data$Area)){
  beta_data_roi <- locbetas_data %>%
    filter(Area == roi)
  beta_anova_roi <- beta_anova_mtl %>% 
    filter(Area == roi)
  
  #If there is any interaction with hemisphere, then keep hemisphere for posthocs. 
  #If not, collapse
  if(any(beta_anova_roi %>% 
         filter(grepl("Category:Hemisphere", Effect)) %>% 
         pull(p_fdrcorrected) 
         <= 0.05)){
    beta_data_roi <- beta_data_roi
    hemi_interac <- TRUE
  }else if(all(beta_anova_roi %>% 
               filter(grepl("Category:Hemisphere", Effect)) %>% 
               pull(p_fdrcorrected) 
               > 0.05)){
    beta_data_roi <- as.data.frame(
      beta_data_roi %>% 
        group_by(Participant, ExperimentPhase, 
                 Category, Area) %>% 
        summarise(Beta = mean(Beta)) %>% 
        ungroup() %>% 
        mutate(Hemisphere = "Collapsed")
    )
    hemi_interac <- FALSE
  }
  
  # Compute mean Betas per Area and Category
  mean_betas <- beta_data_roi %>%
    group_by(Area, Category) %>%
    summarise(Mean_Beta = mean(Beta, na.rm = TRUE), .groups = "drop")
  
  # Generate all pairwise combinations within each Area
  pairwise_diffs <- mean_betas %>%
    group_by(Area) %>%
    summarise(pairwise = list(combn(Category, 2, simplify = FALSE)), .groups = "drop") %>%
    unnest(pairwise) %>%
    mutate(Category1 = purrr::map_chr(pairwise, ~ as.character(.x[1])),
           Category2 = purrr::map_chr(pairwise, ~ as.character(.x[2]))) %>%
    select(-pairwise) %>%
    left_join(mean_betas, by = c("Area", "Category1" = "Category")) %>%
    rename(Mean_Beta1 = Mean_Beta) %>%
    left_join(mean_betas, by = c("Area", "Category2" = "Category")) %>%
    rename(Mean_Beta2 = Mean_Beta,
           group1 = Category1,
           group2 = Category2) %>%
    mutate(MeanDiff = Mean_Beta1 - Mean_Beta2)
  
  beta_posthoc_roi <- full_join(
    beta_data_roi %>% 
      group_by(Hemisphere) %>% 
      t_test(Beta ~ Category,
             paired = TRUE),
    beta_data_roi %>% 
      group_by(Hemisphere) %>% 
      cohens_d(Beta ~ Category,
               paired = TRUE),
    by = join_by(Hemisphere, .y., group1, group2, n1, n2)
  ) %>% 
    full_join(pairwise_diffs,
              by = c("group1", "group2")) %>% 
    relocate(Area) %>% 
    relocate(starts_with("Mean"), .after = n2) %>% 
    relocate(c(effsize, magnitude), .after = df) %>% 
    select(-starts_with("p.adj"))
  
  beta_posthoc <- bind_rows(
    beta_posthoc_roi,
    beta_posthoc
  )
  
  #### THIRD: Plot the Betas ####
  beta_data_roi <- beta_data_roi %>% 
    mutate(Category = factor(Category,
                             levels = factor_levels$category_localiser$levels,
                             labels = factor_levels$category_localiser$labels),
           Hemisphere = factor(Hemisphere,
                               levels = factor_levels$hemisphere$levels,
                               labels = factor_levels$hemisphere$labels),
           Hemi_Area = paste(beta_data_roi$Hemisphere, beta_data_roi$Area))
  
  roi_type <- factor_levels$rois$roi_type[factor_levels$rois$levels == roi]
  if(roi_type == "HC"){
    ylims <- c(-10, 20)
    ybreaks <- seq(ylims[1], ylims[2], by = 10)
    max_y_range <- diff(ylims)
  }else if(roi_type == "MTL"){
    ylims <- c(-10, 50)
    ybreaks <- c(-10, 0, seq(20, ylims[2], by = 20))
    max_y_range <- diff(ylims)
  }else if(roi_type == "LocCat"){
    ylims <- c(NA, (max(beta_data_roi$Beta)+10))
    ybreaks <- seq(-50, ylims[2], by = 50)
    max_y_range <- ylims[2] - min(beta_data_roi$Beta)
  }
  
  dotsize <- (max(beta_data_roi$Beta) - min(beta_data_roi$Beta))/max_y_range
  # Calculate a consistent binwidth based on the global y-range
  global_binwidth <- max_y_range/30  # Divide the max range into 30 bins, adjust as needed
  
  
  #Plot together
  beta_boxplot <- ggplot(data = beta_data_roi, 
                       aes(x = Category, 
                           y = Beta, 
                           fill = Category)) +
    geom_boxplot(position = position_dodge(1),
                 lwd = 1.5,
                 outlier.colour = NA) +
    geom_dotplot(binaxis = 'y',
                 stackdir = 'center',
                 position = position_dodge(1),
                 binwidth = global_binwidth,
                 dotsize = 0.75,
                 alpha = 0.4) +
    facet_wrap(~Hemi_Area) + 
    coord_cartesian(ylim = ylims) +
    scale_y_continuous(breaks = ybreaks) +
    labs(x = "Category",
         y = "Percent Signal Change",
         fill = "Category") +
    scale_fill_manual(values = factor_levels$category_localiser$colours) +
    geom_hline(yintercept=0, linewidth=0.5, linetype="dashed", alpha=0.5) + 
    x_axis_theme + y_axis_theme + paper_facet_theme + 
    blank_bg_theme + legend_theme
  beta_plots[[roi]] <- beta_boxplot
}

anova_tests_loc_cat <- beta_anova
anova_tests_detailed_loc_cat <- beta_anova_detailed
beta_plots_loc_cat <- beta_plots
beta_posthoc_loc_cat <- beta_posthoc

#Only keep handsegmented ROIs so we can correct for them
beta_posthoc_mtl <- beta_posthoc %>% 
  filter(Area %in% mtl_rois)
beta_posthoc_mtl <- correct_p_vals(beta_posthoc_mtl, "p") %>% 
  mutate(Greater = ifelse(Mean_Beta1 > Mean_Beta2,
                          group1,
                          ifelse(Mean_Beta2 > Mean_Beta1, 
                                 group2, 
                                 "Equal")))

category_selectivity <- beta_posthoc_mtl %>% 
  filter(p_fdrcorrected <= 0.05) %>% 
  group_by(Area, Hemisphere, Greater) %>% 
  summarise(ln = length(Greater))

scene_face_selective_loc <- c(
  category_selectivity %>% 
    filter(Greater == "scene",
           ln == 2) %>% 
    pull(Area) %>% 
    unique(),
  category_selectivity %>% 
    filter(Greater == "face",
           ln == 2) %>% 
    pull(Area) %>% 
    unique()
)



