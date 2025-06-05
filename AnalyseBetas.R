# Set working directory and source initialization file
setwd("~/Downloads/CategorySpecificAssociativeInference/")
source("Initialise.R")

# Make a list of ROIs available
mtl_rois <- factor_levels$rois$levels
# Exclude the FFA ROI
mtl_rois <- mtl_rois[!(mtl_rois == "FFA")]

####################### Read in betas and organise them #######################

############ Localiser ############
# Initialize empty dataframe to store all participants' data
locbetas_data <- c()
deconvolve_path <- paste0(derivatives_path, "LocalDeconvolveOutput/")

# Loop through all participants and read their task event files
for(part in loc_participants){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(deconvolve_path, part, "/", part, "_task-local_betas.tsv"),
                          header = TRUE,
                          sep = "\t")
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

############ AIM Task ############
# Initialize empty dataframe to store all participants' data
aimbetas_data <- c()
deconvolve_path <- paste0(derivatives_path, "AIMDeconvolveOutput/")

# Loop through all participants and read their task event files
for(part in aim_participants){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(deconvolve_path, part, "/", part, "_task-aim_betas.tsv"),
                          header = TRUE,
                          sep = "\t")
  # Add participant ID column
  part_data$Participant <- part
  # Combine with existing data
  aimbetas_data <- rbind(aimbetas_data, part_data)
}

# Reorder columns to put Participant ID first
aimbetas_data <- relocate(aimbetas_data,
                          Participant)

# Verify all intended participants are included in the dataset
all(sort(unique(aimbetas_data$Participant)) == sort(aim_participants))

####################### STEP 1: Find category-selective subregions #######################

############ Localiser ############
# Initialize empty vectors to store ANOVA results
beta_anova <- c()
beta_anova_detailed <- c()

######## FIRST: Run ANOVA ########
# Loop through each ROI in the dataset
for(roi in unique(locbetas_data$Area)){
  # Filter data for current ROI
  beta_roi_data <- locbetas_data %>% 
    filter(Area == roi)
  
  # Perform repeated measures ANOVA
  # Tests effects of Category and Hemisphere on Beta values
  beta_anova_test <- ezANOVA(
    data = beta_roi_data, 
    wid = Participant,              # Subject identifier
    within = c(Category, Hemisphere), # Within-subject factors
    dv = Beta,                      # Dependent variable
    detailed = TRUE                 # Request detailed output
  )
  # Calculate effect sizes for the ANOVA results
  beta_anova_df <- aovEffectSize(beta_anova_test)$ANOVA
  
  # Store results:
  # Detailed results stored by ROI name
  beta_anova_detailed[[roi]] <- beta_anova_test
  # Summary results combined into single dataframe
  beta_anova <- rbind(beta_anova, 
                      cbind(data.frame(Area = roi, 
                                       Event = "Localiser"), 
                            beta_anova_df))
}

# Filter to keep MTL ROIs and correct ANOVA results for multiple comparisons
beta_anova_mtl <- beta_anova %>% 
  filter(Area %in% mtl_rois,
         Effect != "(Intercept)")
beta_anova_mtl <- beta_anova_mtl %>% 
  group_by(Effect) %>% 
  do(correct_p_vals(., "p")
  )

# Recombine with FFA results to use for assessing the 
# Category x Hemisphere interaction in the posthoc analysis
beta_anova_mtl <- bind_rows(
  beta_anova_mtl,
  beta_anova %>% 
    filter(Area == "FFA",
           Effect != "(Intercept)") %>% 
    mutate(p_fdrcorrected = p) # Keep the uncorrected p for FFA
)


######## SECOND run Posthocs: ########
# Initialize empty lists/vectors for storing results
beta_plots <- list()
beta_posthoc <- c()

# Loop through each ROI
for(roi in unique(locbetas_data$Area)){
  # Filter data for current ROI
  beta_data_roi <- locbetas_data %>%
    filter(Area == roi)
  beta_anova_roi <- beta_anova_mtl %>% 
    filter(Area == roi)
  
  # Check for hemisphere interactions
  if(any(beta_anova_roi %>% 
         filter(grepl("Category:Hemisphere", Effect)) %>% 
         pull(p_fdrcorrected) 
         <= 0.05)){
    # If there's a significant interaction between Category and Hemisphere,
    # keep hemispheres separate for analysis
    beta_data_roi <- beta_data_roi
    hemi_interac <- TRUE
  }else if(all(beta_anova_roi %>% 
               filter(grepl("Category:Hemisphere", Effect)) %>% 
               pull(p_fdrcorrected) 
               > 0.05)){
    # If no significant hemisphere interaction, collapse across hemispheres
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
  
  # Calculate mean beta values for each Area and Category combination
  mean_betas <- beta_data_roi %>%
    group_by(Area, Category) %>%
    summarise(Mean_Beta = mean(Beta, na.rm = TRUE), .groups = "drop")
  
  # Create all possible pairwise category comparisons and calculate differences
  pairwise_diffs <- mean_betas %>%
    group_by(Area) %>%
    # Generate all possible category pairs
    summarise(pairwise = list(combn(Category, 2, simplify = FALSE)), 
              .groups = "drop") %>%
    unnest(pairwise) %>%
    # Extract category names from pairs
    mutate(Category1 = purrr::map_chr(pairwise, ~ as.character(.x[1])),
           Category2 = purrr::map_chr(pairwise, ~ as.character(.x[2]))) %>%
    select(-pairwise) %>%
    # Join mean beta values for each category
    left_join(mean_betas, by = c("Area", "Category1" = "Category")) %>%
    rename(Mean_Beta1 = Mean_Beta) %>%
    left_join(mean_betas, by = c("Area", "Category2" = "Category")) %>%
    rename(Mean_Beta2 = Mean_Beta,
           group1 = Category1,
           group2 = Category2) %>%
    # Calculate difference between categories
    mutate(MeanDiff = Mean_Beta1 - Mean_Beta2)
  
  # Perform post-hoc tests: t-tests and Cohen's d
  beta_posthoc_roi <- full_join(
    # Paired t-tests between categories
    beta_data_roi %>% 
      group_by(Hemisphere) %>% 
      t_test(Beta ~ Category,
             paired = TRUE),
    # Effect size calculations
    beta_data_roi %>% 
      group_by(Hemisphere) %>% 
      cohens_d(Beta ~ Category,
               paired = TRUE),
    by = join_by(Hemisphere, .y., group1, group2, n1, n2)
  ) %>% 
    # Add mean differences and reorganize columns
    full_join(pairwise_diffs,
              by = c("group1", "group2")) %>% 
    relocate(Area) %>% 
    relocate(starts_with("Mean"), .after = n2) %>% 
    relocate(c(effsize, magnitude), .after = df) %>% 
    select(-starts_with("p.adj"))
  
  # Combine results with previous ROIs
  beta_posthoc <- bind_rows(
    beta_posthoc_roi,
    beta_posthoc
  )
  
  ######## THIRD: Plot the Betas ########
  # Format factor levels for plotting
  beta_data_roi <- beta_data_roi %>% 
    mutate(Category = factor(Category,
                             levels = factor_levels$category_localiser$levels,
                             labels = factor_levels$category_localiser$labels),
           Hemisphere = factor(Hemisphere,
                               levels = factor_levels$hemisphere$levels,
                               labels = factor_levels$hemisphere$labels),
           Hemi_Area = paste(Hemisphere, Area))
  
  # Set plot parameters based on ROI type
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
  
  # Calculate plot parameters for consistent dot sizes
  dotsize <- (max(beta_data_roi$Beta) - min(beta_data_roi$Beta))/max_y_range
  global_binwidth <- max_y_range/30  # Divide the max range into 30 bins
  
  # Create boxplot with overlaid dot plot
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
         y = "Beta Estimates",
         fill = "") +
    scale_fill_manual(values = factor_levels$category_localiser$colours) +
    geom_hline(yintercept=0, linewidth=0.5, linetype="dashed", alpha=0.5) + 
    x_axis_theme + y_axis_theme + paper_facet_theme + 
    blank_bg_theme + legend_theme
  beta_plots[[roi]] <- beta_boxplot
}

# Store results in final variables
anova_tests_loc_cat <- beta_anova_mtl
anova_tests_detailed_loc_cat <- beta_anova_detailed
beta_plots_loc_cat <- beta_plots
beta_posthoc_loc_cat <- beta_posthoc

# Process MTL ROIs separately
beta_posthoc_loc_cat_mtl <- beta_posthoc %>% 
  filter(Area %in% mtl_rois)
# Apply FDR correction and determine which category shows greater activation
beta_posthoc_loc_cat_mtl <- correct_p_vals(beta_posthoc_loc_cat_mtl, "p") %>% 
  mutate(Greater = ifelse(Mean_Beta1 > Mean_Beta2,
                          group1,
                          ifelse(Mean_Beta2 > Mean_Beta1, 
                                 group2, 
                                 "Equal")))

# Summarize category selectivity
# Count significant differences where one category shows consistently greater activation
category_selectivity <- beta_posthoc_loc_cat_mtl %>% 
  filter(p_fdrcorrected <= 0.05) %>% 
  group_by(Area, Hemisphere, Greater) %>% 
  summarise(ln = length(Greater))

# Identify ROIs that are selective for scenes or faces
scene_face_selective_loc <- c(
  # Scene-selective ROIs (greater response to scenes than faces and objects)
  category_selectivity %>% 
    filter(Greater == "scene",
           ln == 2) %>% 
    pull(Area) %>% 
    unique(),
  # Face-selective ROIs (greater response to faces than scenes and objects)
  category_selectivity %>% 
    filter(Greater == "face",
           ln == 2) %>% 
    pull(Area) %>% 
    unique()
)

############ AIM Task: AB trials ############
studyABbetas_data <- aimbetas_data %>% 
  filter(ExperimentPhase == "study",
         PairType == "AB",
         Accuracy == "Remembered")

# Calculate mean beta values for each Area and Category combination
pairwise_diffs <- studyABbetas_data %>%
  group_by(Area, Category) %>%
  summarise(Mean_Beta = mean(Beta, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = Category,
              values_from = Mean_Beta) %>% 
  dplyr::rename(Mean_Beta1 = face,
                Mean_Beta2 = scene) %>% 
  mutate(group1 = "face",
         group2 = "scene",
         MeanDiff = Mean_Beta1 - Mean_Beta2)

######## FIRST: Run ANOVA ########
# Initialize empty vectors to store ANOVA results
beta_anova <- c()
beta_anova_detailed <- c()
# Loop through each ROI in the dataset
for(roi in unique(studyABbetas_data$Area)){
  # Filter data for current ROI
  beta_roi_data <- studyABbetas_data %>% 
    filter(Area == roi)
  
  # Perform repeated measures ANOVA
  # Tests effects of Category and Hemisphere on Beta values
  beta_anova_test <- ezANOVA(
    data = beta_roi_data, 
    wid = Participant,              # Subject identifier
    within = c(Category, Hemisphere), # Within-subject factors
    dv = Beta,                      # Dependent variable
    detailed = TRUE                 # Request detailed output
  )
  # Calculate effect sizes for the ANOVA results
  beta_anova_df <- aovEffectSize(beta_anova_test)$ANOVA
  
  # Store results:
  # Detailed results stored by ROI name
  beta_anova_detailed[[roi]] <- beta_anova_test
  # Summary results combined into single dataframe
  beta_anova <- rbind(beta_anova, 
                      cbind(data.frame(Area = roi, 
                                       Event = "StudyAB"), 
                            beta_anova_df))
}

# Filter to keep MTL ROIs and correct ANOVA results for multiple comparisons
beta_anova_mtl <- beta_anova %>% 
  filter(Area %in% mtl_rois,
         Effect != "(Intercept)")
beta_anova_mtl <- beta_anova_mtl %>% 
  group_by(Effect) %>% 
  do(correct_p_vals(., "p")
  )

######## SECOND run Posthocs: ########
# Initialize empty lists/vectors for storing results
beta_posthoc <- c()

# Loop through each ROI
for(roi in unique(studyABbetas_data$Area)){
  # Filter data for current ROI
  beta_data_roi <- studyABbetas_data %>%
    filter(Area == roi)
  beta_anova_roi <- beta_anova_mtl %>% 
    filter(Area == roi)
  
  # Check for hemisphere interactions
  if(any(beta_anova_roi %>% 
         filter(grepl("Category:Hemisphere", Effect)) %>% 
         pull(p_fdrcorrected) 
         <= 0.05)){
    # If there's a significant interaction between Category and Hemisphere,
    # keep hemispheres separate for analysis
    beta_data_roi <- beta_data_roi
    hemi_interac <- TRUE
  }else if(all(beta_anova_roi %>% 
               filter(grepl("Category:Hemisphere", Effect)) %>% 
               pull(p_fdrcorrected) 
               > 0.05)){
    # If no significant hemisphere interaction, collapse across hemispheres
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
  
  # Perform post-hoc tests: t-tests and Cohen's d
  beta_posthoc_roi <- full_join(
    # Paired t-tests between categories
    beta_data_roi %>% 
      group_by(Hemisphere) %>% 
      t_test(Beta ~ Category,
             paired = TRUE),
    # Effect size calculations
    beta_data_roi %>% 
      group_by(Hemisphere) %>% 
      cohens_d(Beta ~ Category,
               paired = TRUE),
    by = join_by(Hemisphere, .y., group1, group2, n1, n2)
  )
  
  # Combine results with previous ROIs
  beta_posthoc <- bind_rows(
    beta_posthoc_roi %>% 
      mutate(Area = roi) %>% 
      relocate(Area),
    beta_posthoc
  )
}

# Store results in final variables
anova_tests_studyAB_pairtype <- beta_anova
anova_tests_detailed_studyAB_pairtype <- beta_anova_detailed
beta_plots_studyAB_pairtype <- beta_plots
beta_posthoc_studyAB_pairtype <- beta_posthoc

# Process MTL ROIs separately
beta_posthoc_studyAB_pairtype_mtl <- beta_posthoc_studyAB_pairtype %>% 
  filter(Area %in% mtl_rois) %>% 
  # Apply FDR correction and determine which category shows greater activation
  do(correct_p_vals(., "p")) %>% 
  # Add mean differences and reorganize columns
  full_join(pairwise_diffs,
            by = c("Area", "group1", "group2")) %>% 
  relocate(Area) %>% 
  relocate(starts_with("Mean"), .after = n2) %>% 
  relocate(c(effsize, magnitude), .after = df) %>% 
  select(-starts_with("p.adj")) %>% 
  mutate(Greater = ifelse(Mean_Beta1 > Mean_Beta2,
                          group1,
                          ifelse(Mean_Beta2 > Mean_Beta1, 
                                 group2, 
                                 "Equal")))

# Identify ROIs that are selective for scenes or faces
scene_face_selective_aim <- c(
  # Scene-selective ROIs (greater response to scenes than faces and objects)
  beta_posthoc_studyAB_pairtype_mtl %>% 
    filter(p_fdrcorrected <= 0.05,
           Greater == "scene") %>% 
    pull(Area) %>% 
    unique(),
  # Face-selective ROIs (greater response to faces than scenes and objects)
  beta_posthoc_studyAB_pairtype_mtl %>% 
    filter(p_fdrcorrected <= 0.05,
           Greater == "face") %>% 
    pull(Area) %>% 
    unique()
)

# Get the areas that exhibit scene- and face-selectivity that 
# are common between the localiser and study AB trials
scene_face_selective <- intersect(scene_face_selective_loc,
                                  scene_face_selective_aim)

####################### STEP 2: Analyse AIM study phase #######################
# Filter study phase data to include only remembered trials
studybetas_data <- aimbetas_data %>% 
  filter(ExperimentPhase == "study",
         Accuracy == "Remembered")

######## FIRST: Run ANOVA ########
# Initialize empty vectors to store ANOVA results
beta_anova <- c()
beta_anova_detailed <- c()

# Loop through ROIs that showed category selectivity in both localizer and AIM task
for(roi in c(scene_face_selective, "FFA")){
  # Filter data for current ROI
  beta_roi_data <- studybetas_data %>% 
    filter(Area == roi)
  
  # Perform repeated measures ANOVA
  # Tests effects of Category and Hemisphere on Beta values
  beta_anova_test <- ezANOVA(
    data = beta_roi_data, 
    wid = Participant,              # Subject identifier
    within = c(Category, 
               PairType,
               Hemisphere), # Within-subject factors
    dv = Beta,                      # Dependent variable
    detailed = TRUE                 # Request detailed output
  )
  # Calculate effect sizes for the ANOVA results
  beta_anova_df <- aovEffectSize(beta_anova_test)$ANOVA
  
  # Store results:
  # Detailed results stored by ROI name
  beta_anova_detailed[[roi]] <- beta_anova_test
  # Summary results combined into single dataframe
  beta_anova <- rbind(beta_anova, 
                      cbind(data.frame(Area = roi, 
                                       Event = "Study"), 
                            beta_anova_df))
}

# Process ANOVA results for MTL ROIs
beta_anova_mtl <- beta_anova %>% 
  filter(Area %in% scene_face_selective,
         Effect != "(Intercept)")
beta_anova_mtl <- beta_anova_mtl %>% 
  group_by(Effect) %>% 
  do(correct_p_vals(., "p")
  )

# Recombine with FFA results to use for assessing the 
# Category x Hemisphere interaction in the posthoc analysis
beta_anova_mtl <- bind_rows(
  beta_anova_mtl,
  beta_anova %>% 
    filter(Area == "FFA",
           Effect != "(Intercept)") %>% 
    mutate(p_fdrcorrected = p) # Keep the uncorrected p for FFA
)

######## SECOND run Posthocs: ########
# Initialize empty lists/vectors for storing results
beta_plots <- list()
beta_posthoc <- c()

# Loop through category-selective ROIs and the FFA
for(roi in c(scene_face_selective, "FFA")){
  # Filter data for current ROI
  beta_data_roi <- studybetas_data %>%
    filter(Area == roi)
  beta_anova_roi <- beta_anova_mtl %>% 
    filter(Area == roi)
  
  # Check for hemisphere interactions
  if(any(beta_anova_roi %>% 
         filter(grepl("Category:PairType:Hemisphere", Effect)) %>% 
         pull(p_fdrcorrected) 
         <= 0.05)){
    # If there's a significant interaction between Category and Hemisphere,
    # keep hemispheres separate for analysis
    beta_data_roi <- beta_data_roi
    hemi_interac <- TRUE
  }else if(all(beta_anova_roi %>% 
               filter(grepl("Category:PairType:Hemisphere", Effect)) %>% 
               pull(p_fdrcorrected) 
               > 0.05)){
    # If no significant hemisphere interaction, collapse across hemispheres
    beta_data_roi <- as.data.frame(
      beta_data_roi %>% 
        group_by(Participant, ExperimentPhase, 
                 Category, PairType, Area) %>% 
        summarise(Beta = mean(Beta)) %>% 
        ungroup() %>% 
        mutate(Hemisphere = "Collapsed")
    )
    hemi_interac <- FALSE
  }
  
  # Perform post-hoc tests: t-tests and Cohen's d
  beta_posthoc_roi <- bind_rows(
    full_join(
      # Paired t-tests between categories
      beta_data_roi %>% 
        group_by(Hemisphere, Category) %>% 
        t_test(Beta ~ PairType,
               paired = TRUE),
      # Effect size calculations
      beta_data_roi %>% 
        group_by(Hemisphere, Category) %>% 
        cohens_d(Beta ~ PairType,
                 paired = TRUE),
      by = join_by(Hemisphere, Category, .y., group1, group2, n1, n2)
    ) %>% 
      dplyr::rename(within = Category),
    full_join(
      # Paired t-tests between categories
      beta_data_roi %>% 
        group_by(Hemisphere, PairType) %>% 
        t_test(Beta ~ Category,
               paired = TRUE),
      # Effect size calculations
      beta_data_roi %>% 
        group_by(Hemisphere, PairType) %>% 
        cohens_d(Beta ~ Category,
                 paired = TRUE),
      by = join_by(Hemisphere, PairType, .y., group1, group2, n1, n2)
    ) %>% 
      dplyr::rename(within = PairType))
  
  # Combine results with previous ROIs
  beta_posthoc <- bind_rows(
    beta_posthoc_roi %>% 
      mutate(Area = roi) %>% 
      relocate(Area),
    beta_posthoc
  )
  
  ######## THIRD: Plot the Betas ########
  # Format factor levels for plotting with appropriate labels
  beta_data_roi <- beta_data_roi %>%
    mutate(Category = factor(Category,
                             levels = factor_levels$category_localiser$levels,
                             labels = factor_levels$category_localiser$labels),
           PairType = factor(PairType,
                             levels = factor_levels$pairtype_study$levels,
                             labels = factor_levels$pairtype_study$labels),
           Hemisphere = factor(Hemisphere,
                               levels = factor_levels$hemisphere$levels,
                               labels = factor_levels$hemisphere$labels),
           Hemi_Area = paste(Hemisphere, Area))
  
  # Set y-axis limits and breaks based on ROI type
  roi_type <- factor_levels$rois$roi_type[factor_levels$rois$levels == roi]
  if(roi_type == "HC"){
    ylims <- c(-15, 30)
    ybreaks <- seq(ylims[1], ylims[2], by = 10)
    max_y_range <- diff(ylims)
  }else if(roi_type == "MTL"){
    ylims <- c(-10, 90)
    ybreaks <- c(-10, 0, seq(20, ylims[2], by = 20))
    max_y_range <- diff(ylims)
  }else if(roi_type == "LocCat"){
    ylims <- c(NA, (max(beta_data_roi$Beta)+10))
    ybreaks <- seq(-50, ylims[2], by = 50)
    max_y_range <- ylims[2] - min(beta_data_roi$Beta)
  }
  
  # Calculate plot parameters for consistent dot sizes
  dotsize <- (max(beta_data_roi$Beta) - min(beta_data_roi$Beta))/max_y_range
  global_binwidth <- max_y_range/30  # Divide the max range into 30 bins
  
  # Create visualization combining boxplots and individual data points
  beta_boxplot <- ggplot(data = beta_data_roi,
                         aes(x = PairType,
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
    labs(x = "Pair Type",
         y = "Beta Estimates",
         fill = "") +
    scale_fill_manual(values = factor_levels$category$colours) +
    geom_hline(yintercept=0, linewidth=0.5, linetype="dashed", alpha=0.5) +
    x_axis_theme + y_axis_theme + paper_facet_theme +
    blank_bg_theme + legend_theme
  beta_plots[[roi]] <- beta_boxplot
}

# Store results in final variables
anova_tests_study_pairtypecat <- beta_anova_mtl
anova_tests_detailed_study_pairtypecat <- beta_anova_detailed
beta_plots_study_pairtypecat <- beta_plots
beta_posthoc_study_pairtypecat <- beta_posthoc

# Process MTL ROIs separately
beta_posthoc_study_pairtypecat <- beta_posthoc_study_pairtypecat %>% 
  filter(Area %in% scene_face_selective) %>% 
  # Apply FDR correction and determine which category shows greater activation
  do(correct_p_vals(., "p"))

####################### STEP 3: Analyse AIM test phase #######################
# Filter test phase data to include only correct trials
testbetas_data <- aimbetas_data %>% 
  filter(ExperimentPhase == "test",
         TestType == "indir")

# Exclude participants who didn't have any incorrect trials
participants_to_exclude <- testbetas_data %>%
  mutate(Participant = factor(Participant),
         Category = factor(Category),
         Accuracy = factor(Accuracy)) %>% 
  group_by(Participant, Category, Accuracy, .drop = FALSE) %>% 
  summarise(len = length(Accuracy)) %>% 
  filter(len == 0) %>% 
  pull(Participant)

# Exclude the participants' data
testbetas_data <- testbetas_data %>% 
  filter(!(Participant %in% participants_to_exclude))

######## FIRST: Run ANOVA ########
# Initialize empty vectors to store ANOVA results
beta_anova <- c()
beta_anova_detailed <- c()
# Loop through each ROI in the dataset
for(roi in scene_face_selective){
  # Filter data for current ROI
  beta_roi_data <- testbetas_data %>% 
    filter(Area == roi)
  
  # Perform repeated measures ANOVA
  # Tests effects of Category and Hemisphere on Beta values
  beta_anova_test <- ezANOVA(
    data = beta_roi_data, 
    wid = Participant,              # Subject identifier
    within = c(Category, 
               Accuracy,
               Hemisphere), # Within-subject factors
    dv = Beta,                      # Dependent variable
    detailed = TRUE                 # Request detailed output
  )
  # Calculate effect sizes for the ANOVA results
  beta_anova_df <- aovEffectSize(beta_anova_test)$ANOVA
  
  # Store results:
  # Detailed results stored by ROI name
  beta_anova_detailed[[roi]] <- beta_anova_test
  # Summary results combined into single dataframe
  beta_anova <- rbind(beta_anova, 
                      cbind(data.frame(Area = roi, 
                                       Event = "Test"), 
                            beta_anova_df))
}

# Filter to keep MTL ROIs and correct ANOVA results for multiple comparisons
beta_anova_mtl <- beta_anova %>% 
  filter(Area %in% scene_face_selective,
         Effect != "(Intercept)")
beta_anova_mtl <- beta_anova_mtl %>% 
  group_by(Effect) %>% 
  do(correct_p_vals(., "p")
  )

######## SECOND run Posthocs: ########

# Initialize empty lists/vectors for storing results
beta_plots <- list()
beta_posthoc <- c()

# Loop through each ROI
for(roi in scene_face_selective){
  # Filter data for current ROI
  beta_data_roi <- testbetas_data %>%
    filter(Area == roi)
  beta_anova_roi <- beta_anova_mtl %>% 
    filter(Area == roi)
  
  # Check for hemisphere interactions
  if(any(beta_anova_roi %>% 
         filter(grepl("Category:Accuracy:Hemisphere", Effect)) %>% 
         pull(p_fdrcorrected) 
         <= 0.05)){
    # If there's a significant interaction between Category and Hemisphere,
    # keep hemispheres separate for analysis
    beta_data_roi <- beta_data_roi
    hemi_interac <- TRUE
  }else if(all(beta_anova_roi %>% 
               filter(grepl("Category:Accuracy:Hemisphere", Effect)) %>% 
               pull(p_fdrcorrected) 
               > 0.05)){
    # If no significant hemisphere interaction, collapse across hemispheres
    beta_data_roi <- as.data.frame(
      beta_data_roi %>% 
        group_by(Participant, ExperimentPhase, 
                 Category, Accuracy, Area) %>% 
        summarise(Beta = mean(Beta)) %>% 
        ungroup() %>% 
        mutate(Hemisphere = "Collapsed")
    )
    hemi_interac <- FALSE
  }
  
  # Perform post-hoc tests: t-tests and Cohen's d
  beta_posthoc_roi <- bind_rows(
    full_join(
      # Paired t-tests between categories
      beta_data_roi %>% 
        group_by(Hemisphere, Category) %>% 
        t_test(Beta ~ Accuracy,
               paired = TRUE),
      # Effect size calculations
      beta_data_roi %>% 
        group_by(Hemisphere, Category) %>% 
        cohens_d(Beta ~ Accuracy,
                 paired = TRUE),
      by = join_by(Hemisphere, Category, .y., group1, group2, n1, n2)
    ) %>% 
      dplyr::rename(within = Category),
    full_join(
      # Paired t-tests between categories
      beta_data_roi %>% 
        group_by(Hemisphere, Accuracy) %>% 
        t_test(Beta ~ Category,
               paired = TRUE),
      # Effect size calculations
      beta_data_roi %>% 
        group_by(Hemisphere, Accuracy) %>% 
        cohens_d(Beta ~ Category,
                 paired = TRUE),
      by = join_by(Hemisphere, Accuracy, .y., group1, group2, n1, n2)
    ) %>% 
      dplyr::rename(within = Accuracy))
  
  # Combine results with previous ROIs
  beta_posthoc <- bind_rows(
    beta_posthoc_roi %>% 
      mutate(Area = roi) %>% 
      relocate(Area),
    beta_posthoc
  )
  
  ######## THIRD: Plot the Betas ########
  # Format factor levels for plotting
  beta_data_roi <- beta_data_roi %>%
    mutate(Category = factor(Category,
                             levels = factor_levels$category_localiser$levels,
                             labels = factor_levels$category_localiser$labels),
           Accuracy = factor(Accuracy,
                             levels = factor_levels$acc_test$levels,
                             labels = factor_levels$acc_test$labels),
           Hemisphere = factor(Hemisphere,
                               levels = factor_levels$hemisphere$levels,
                               labels = factor_levels$hemisphere$labels),
           Hemi_Area = paste(Hemisphere, Area))
  
  # Set plot parameters based on ROI type
  roi_type <- factor_levels$rois$roi_type[factor_levels$rois$levels == roi]
  if(roi_type == "HC"){
    ylims <- c(-30, 30)
    ybreaks <- seq(ylims[1], ylims[2], by = 15)
    max_y_range <- diff(ylims)
  }else if(roi_type == "MTL"){
    ylims <- c(-25, 60)
    ybreaks <- seq(ylims[1], ylims[2], by = 15)
    max_y_range <- diff(ylims)
  }else if(roi_type == "LocCat"){
    ylims <- c(NA, (max(beta_data_roi$Beta)+10))
    ybreaks <- seq(-50, ylims[2], by = 50)
    max_y_range <- ylims[2] - min(beta_data_roi$Beta)
  }
  
  # Calculate plot parameters for consistent dot sizes
  dotsize <- (max(beta_data_roi$Beta) - min(beta_data_roi$Beta))/max_y_range
  global_binwidth <- max_y_range/30  # Divide the max range into 30 bins
  
  # Create boxplot with overlaid dot plot
  beta_boxplot <- ggplot(data = beta_data_roi,
                         aes(x = Accuracy,
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
    labs(x = "Test Type",
         y = "Beta Estimates",
         fill = "") +
    scale_fill_manual(values = factor_levels$category$colours) +
    geom_hline(yintercept=0, linewidth=0.5, linetype="dashed", alpha=0.5) +
    x_axis_theme + y_axis_theme + paper_facet_theme +
    blank_bg_theme + legend_theme
  beta_plots[[roi]] <- beta_boxplot
}

# Store results in final variables
anova_tests_test_acccat <- beta_anova_mtl
anova_tests_detailed_test_acccat <- beta_anova_detailed
beta_plots_test_acccat <- beta_plots
beta_posthoc_test_acccat <- beta_posthoc

# Process MTL ROIs separately
beta_posthoc_test_acccat <- beta_posthoc_test_acccat %>% 
  filter(Area %in% scene_face_selective) %>% 
  # Apply FDR correction and determine which category shows greater activation
  do(correct_p_vals(., "p"))

####################### Arrange Plots and View Results #######################
############ Localiser ############

# View ANOVAs from the localiser results
View(anova_tests_loc_cat %>%
       mutate(`F` = round(`F`, 3),
              p = round(p, 3),
              p_fdrcorrected = round(p_fdrcorrected, 3),
              pes = round(pes, 3)))
# View posthoc test from the localiser results
View(beta_posthoc_loc_cat_mtl %>%
       mutate(statistic = round(statistic, 3),
              p = round(p, 3),
              p_fdrcorrected = round(p_fdrcorrected, 3),
              effsize = round(effsize, 3)))
# Create combined figure for localizer results
(loc_plots <- ggarrange(beta_plots_loc_cat$`Ant Hipp Head`,
                        beta_plots_loc_cat$Sub +
                          theme(axis.title.y = element_blank()),
                        beta_plots_loc_cat$CA1 +
                          theme(axis.title.y = element_blank()),
                        beta_plots_loc_cat$CA23DG,
                        beta_plots_loc_cat$`Hipp Tail` +
                          theme(axis.title.y = element_blank()),
                        beta_plots_loc_cat$alERC +
                          theme(axis.title.y = element_blank()),
                        beta_plots_loc_cat$pmERC,
                        beta_plots_loc_cat$PRC +
                          theme(axis.title.y = element_blank()),
                        beta_plots_loc_cat$PHC +
                          theme(axis.title.y = element_blank()), 
                        nrow = 3, ncol = 3, 
                        common.legend = TRUE, 
                        legend = "bottom"))

############ AIM study phase ############

# View ANOVAs from the AIM study phase results
View(anova_tests_study_pairtypecat %>%
       mutate(`F` = round(`F`, 3),
              p = round(p, 3),
              p_fdrcorrected = round(p_fdrcorrected, 3),
              pes = round(pes, 3)))
# View posthoc test from the AIM study phase results
View(beta_posthoc_study_pairtypecat %>%
       mutate(statistic = round(statistic, 3),
              p = round(p, 3),
              p_fdrcorrected = round(p_fdrcorrected, 3),
              effsize = round(effsize, 3)))
# Create combined figure for study phase results
(aim_study_plots <- ggarrange(beta_plots_study_pairtypecat$`Ant Hipp Head`,
                              beta_plots_study_pairtypecat$alERC +
                                theme(axis.title.y = element_blank()),
                              beta_plots_study_pairtypecat$pmERC,
                              beta_plots_study_pairtypecat$PHC +
                                theme(axis.title.y = element_blank()), 
                              nrow = 2, ncol = 2, 
                              common.legend = TRUE, 
                              legend = "bottom"))

############ AIM test phase ############

# View ANOVAs from the AIM study phase results
View(anova_tests_test_acccat %>%
       mutate(`F` = round(`F`, 3),
              p = round(p, 3),
              p_fdrcorrected = round(p_fdrcorrected, 3),
              pes = round(pes, 3)))
# View posthoc test from the AIM study phase results
View(beta_posthoc_test_acccat %>%
       mutate(statistic = round(statistic, 3),
              p = round(p, 3),
              p_fdrcorrected = round(p_fdrcorrected, 3),
              effsize = round(effsize, 3)))
# Create combined figure for test phase results
(aim_test_plots <- ggarrange(beta_plots_test_acccat$`Ant Hipp Head`,
                             beta_plots_test_acccat$alERC +
                               theme(axis.title.y = element_blank()),
                             beta_plots_test_acccat$pmERC,
                             beta_plots_test_acccat$PHC +
                               theme(axis.title.y = element_blank()), 
                             nrow = 2, ncol = 2, 
                             common.legend = TRUE, 
                             legend = "bottom"))

