# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/CategorySpecificAssociativeInference/")
source("InitialisePaths.R")

################## Read in behavioural data and organise it ##################

# Initialize empty dataframe to store all participants' data
behavioural_data <- c()

# Loop through all participants and read their task event files
for(part in aim_participants){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(data_path, part, "/func/", part, "_task-aim_events.tsv"),
                          header=TRUE)
  # Add participant ID column
  part_data$Participant <- part
  # Combine with existing data
  behavioural_data <- rbind(behavioural_data, part_data)
}

# Reorder columns to put Participant ID first
behavioural_data <- relocate(behavioural_data,
                             Participant)

# Verify all intended participants are included in the dataset
all(sort(unique(behavioural_data$Participant)) == sort(aim_participants))


################## Analyse test phase accuracy data ##################

# Filter data to include only test phase trials
test_data <- behavioural_data %>% 
  filter(phase == "test") %>% 
  mutate(category = factor(category,
                           levels = factor_levels$category$labels),
         trialtype = factor(trialtype,
                            levels = factor_levels$trialtype_test$levels))


# Calculate accuracy metrics grouped by participant, trial type, and category
accuracy_trialtypecategory_data <- test_data %>% 
  group_by(Participant, trialtype, category) %>% 
  summarise(TotalTrials = length(Participant),
            CorrectTrials = sum(accuracy, na.rm = TRUE),
            PercentAccuracy = (CorrectTrials/TotalTrials)*100)

# Perform 2x2 repeated measures ANOVA on accuracy data
# Factors: category and trial type
accuracy_trialtypecategory_anova <- 
  ezANOVA(data = accuracy_trialtypecategory_data,
          dv = PercentAccuracy,
          within = c(category, trialtype),
          wid = Participant,
          detailed = TRUE)

# Calculate accuracy metrics grouped by participant and trial type
accuracy_trialtype_data <- test_data %>% 
  group_by(Participant, trialtype) %>% 
  summarise(TotalTrials = length(Participant),
            CorrectTrials = sum(accuracy, na.rm = TRUE),
            PercentAccuracy = (CorrectTrials/TotalTrials)*100)

# Perform post-hoc paired t-tests and calculate Cohen's d for trial type main effect
accuracy_trialtype_main_posthoc <- full_join(
  # Paired t-test between trial types
  accuracy_trialtype_data %>% 
    as.data.frame() %>% 
    t_test(PercentAccuracy ~ trialtype,
           paired = TRUE),
  # Calculate effect size (Cohen's d)
  accuracy_trialtype_data %>% 
    as.data.frame() %>% 
    cohens_d(PercentAccuracy ~ trialtype,
             paired = TRUE),
  by = join_by(.y., group1, group2, n1, n2)
)

# Perform post-hoc tests for category effects within each trial type
accuracy_category_interac_posthoc <- full_join(
  # Paired t-tests between categories within each trial type
  accuracy_trialtypecategory_data %>% 
    as.data.frame() %>% 
    group_by(trialtype) %>% 
    t_test(PercentAccuracy ~ category,
           paired = TRUE),
  # Calculate effect sizes (Cohen's d)
  accuracy_trialtypecategory_data %>% 
    as.data.frame() %>% 
    group_by(trialtype) %>% 
    cohens_d(PercentAccuracy ~ category,
             paired = TRUE),
  by = join_by(trialtype, .y., group1, group2, n1, n2)
) %>% 
  as.data.frame()

# Apply multiple comparison correction to p-values
accuracy_category_interac_posthoc <- 
  correct_p_vals(accuracy_category_interac_posthoc, 
                 "p")

################## Analyse test phase reaction time data ##################

# Calculate RT metrics grouped by participant, trial type, and category
rt_trialtypecategory_data <- test_data %>% 
  filter(accuracy == 1) %>% 
  group_by(Participant, trialtype, category) %>% 
  summarise(MeanRT = mean(test_RT, na.rm = TRUE))

# Perform 2x2 repeated measures ANOVA on RT data
# Factors: category and trial type
rt_trialtypecategory_anova <- 
  ezANOVA(data = rt_trialtypecategory_data,
          dv = MeanRT,
          within = c(category, trialtype),
          wid = Participant,
          detailed = TRUE)

# Calculate RT metrics grouped by participant and trial type
rt_trialtype_data <- test_data %>% 
  group_by(Participant, trialtype) %>% 
  summarise(MeanRT = mean(test_RT, na.rm = TRUE),
            SDRT = sd(test_RT, na.rm = TRUE))

# Perform post-hoc paired t-tests and calculate Cohen's d for trial type main effect
rt_trialtype_posthoc <- full_join(
  # Paired t-test between trial types
  rt_trialtype_data %>% 
    as.data.frame() %>% 
    t_test(MeanRT ~ trialtype,
           paired = TRUE),
  # Calculate effect size (Cohen's d)
  rt_trialtype_data %>% 
    as.data.frame() %>% 
    cohens_d(MeanRT ~ trialtype,
             paired = TRUE),
  by = join_by(.y., group1, group2, n1, n2)
) %>% 
  # Remove any existing adjusted p-values from the t_test operation
  select(-starts_with("p.adj")) %>% 
  # Need this for correct_p_vals to work
  as.data.frame()

# Apply Bonferroni correction to p-values for multiple comparisons
rt_trialtype_posthoc <- 
  correct_p_vals(rt_trialtype_posthoc, 
                 pvalcol = "p")


# Perform post-hoc tests for category effects within each trial type
rt_categorytrialtype_posthoc <- full_join(
  # Paired t-tests between categories within each trial type
  rt_trialtypecategory_data %>% 
    as.data.frame() %>% 
    group_by(trialtype) %>% 
    t_test(MeanRT ~ category,
           paired = TRUE),
  # Calculate effect sizes (Cohen's d)
  rt_trialtypecategory_data %>% 
    as.data.frame() %>% 
    group_by(trialtype) %>% 
    cohens_d(MeanRT ~ category,
             paired = TRUE),
  by = join_by(.y., trialtype, group1, group2, n1, n2)
) %>% 
  # Need this for correct_p_vals to work
  as.data.frame()

# Apply multiple comparison correction to p-values
rt_categorytrialtype_posthoc <- 
  correct_p_vals(rt_categorytrialtype_posthoc, 
                 pvalcol = "p")

################## Analyse eye-tracking data ##################

# Initialize empty dataframe to store all participants' data
et_data <- c()

# Loop through all participants and read their eye-tracking files
for(part in aim_participants){
  # Read TSV file containing eye-tracking data for each participant
  part_data <- read.table(paste0(data_path, part, "/func/", part, "_task-aim_eyetracking.tsv"),
                          header=TRUE)
  # Add participant ID column
  part_data$Participant <- part
  # Combine with existing data
  et_data <- rbind(et_data, part_data)
}

# Reorder columns to put Participant ID first
et_data <- relocate(et_data,
                    Participant)

# Verify all intended participants are included in the dataset
all(sort(unique(et_data$Participant)) == sort(aim_participants))

# Create new variables and convert categorical variables to factors
et_data <- et_data %>% 
  mutate(
    # Combine item type and trial type into a single identifier
    ia_itemtype_trialtype = paste(ia_item_type, trialtype, sep = "_"),
    # Convert category to factor with predefined levels
    category = factor(category,
                      levels = factor_levels$category$labels),
    # Convert trial type to factor with predefined study phase levels
    trialtype = factor(trialtype,
                       levels = factor_levels$trialtype_study$levels),
    # Convert combined identifier to factor with predefined levels and labels
    ia_itemtype_trialtype = factor(ia_itemtype_trialtype,
                                   levels = factor_levels$studyitem_pairtype$levels,
                                   labels = factor_levels$studyitem_pairtype$labels))

# Calculate number of fixations for each trial and interest area
fixnum_data <- et_data %>% 
  # Group by all relevant variables to maintain trial-level data
  group_by(Participant, block, trial, category, trialtype, ia_item_type, ia_itemtype_trialtype) %>% 
  # Count number of fixations (using length of ia_item_type as each row represents one fixation)
  summarise(num_fix = length(ia_item_type)) 

# Calculate mean number of fixations across trials for each condition
fixnum_categorytrialtypeitemtype_data <- fixnum_data %>% 
  # Group by participant and condition variables
  group_by(Participant, category, trialtype, ia_item_type, ia_itemtype_trialtype) %>% 
  # Calculate mean number of fixations
  summarise(mean_num_fix = mean(num_fix))

# Perform 3-way repeated measures ANOVA on fixation count data
fixnum_categorytrialtypeitemtype_anova <- 
  ezANOVA(data = fixnum_categorytrialtypeitemtype_data, 
          wid = Participant, 
          within = c(category, trialtype, ia_item_type),
          dv = mean_num_fix,
          detailed = TRUE)

# Posthocs for the category x item type interaction
fixnum_categoryitemtype_data <- fixnum_data %>% 
  group_by(Participant, category, ia_item_type) %>% 
  summarise(mean_num_fix = mean(num_fix))

fixnum_categoryitemtype_posthoc <- rbind(
  # First analysis: Compare categories within each item type
  full_join(
    # Conduct paired t-tests comparing categories for each item type
    fixnum_categoryitemtype_data %>% 
      group_by(ia_item_type) %>% 
      t_test(mean_num_fix ~ category,
             paired = TRUE),
    # Calculate effect sizes (Cohen's d) for the same comparisons
    fixnum_categoryitemtype_data %>% 
      group_by(ia_item_type) %>% 
      cohens_d(mean_num_fix ~ category,
               paired = TRUE),
    # Join t-test results with effect sizes
    by = join_by(.y., ia_item_type, group1, group2, n1, n2)
  ) %>% 
    # Rename ia_item_type column to 'within' for rbind
    rename(within = ia_item_type),
  # Second analysis: Compare item types within each category
  full_join(
    # Conduct paired t-tests comparing item types for each category
    fixnum_categoryitemtype_data %>% 
      group_by(category) %>% 
      t_test(mean_num_fix ~ ia_item_type,
             paired = TRUE),
    # Calculate effect sizes (Cohen's d) for the same comparisons
    fixnum_categoryitemtype_data %>% 
      group_by(category) %>% 
      cohens_d(mean_num_fix ~ ia_item_type,
               paired = TRUE),
    # Join t-test results with effect sizes
    by = join_by(.y., category, group1, group2, n1, n2)
  ) %>% 
    # Rename category column to 'within' for consistency rbind
    dplyr::rename(within = category)
) %>% 
  # Need this for correct_p_vals to work
  as.data.frame()

# Apply multiple comparison correction to p-values
fixnum_categoryitemtype_posthoc <- 
  correct_p_vals(fixnum_categoryitemtype_posthoc,
                 "p")

# Posthocs for the pair type x item type interaction
fixnum_trialtypeitemtype_data <- fixnum_data %>% 
  group_by(Participant, trialtype, ia_item_type) %>% 
  summarise(mean_num_fix = mean(num_fix))

# Perform post-hoc tests for the trial type x item type interaction
fixnum_trialtypeitemtype_posthoc <- rbind(
  # First analysis: Compare trial types within each item type
  full_join(
    # Conduct paired t-tests comparing trial types for each item type
    fixnum_trialtypeitemtype_data %>% 
      group_by(ia_item_type) %>% 
      t_test(mean_num_fix ~ trialtype,
             paired = TRUE),
    # Calculate effect sizes (Cohen's d) for the same comparisons
    fixnum_trialtypeitemtype_data %>% 
      group_by(ia_item_type) %>% 
      cohens_d(mean_num_fix ~ trialtype,
               paired = TRUE),
    # Join t-test results with effect sizes
    by = join_by(.y., ia_item_type, group1, group2, n1, n2)
  ) %>% 
    # Rename ia_item_type column to 'within' for rbind
    rename(within = ia_item_type),
  # Second analysis: Compare item types within each trial type
  full_join(
    # Conduct paired t-tests comparing item types for each trial type
    fixnum_trialtypeitemtype_data %>% 
      group_by(trialtype) %>% 
      t_test(mean_num_fix ~ ia_item_type,
             paired = TRUE),
    # Calculate effect sizes (Cohen's d) for the same comparisons
    fixnum_trialtypeitemtype_data %>% 
      group_by(trialtype) %>% 
      cohens_d(mean_num_fix ~ ia_item_type,
               paired = TRUE),
    # Join t-test results with effect sizes
    by = join_by(.y., trialtype, group1, group2, n1, n2)
  ) %>% 
    # Rename trialtype column to 'within' for rbind
    dplyr::rename(within = trialtype)
) %>% 
  # Need this for correct_p_vals to work
  as.data.frame()

# Apply multiple comparison correction to p-values
fixnum_trialtypeitemtype_posthoc <- 
  correct_p_vals(fixnum_trialtypeitemtype_posthoc,
                 "p")

################## Plot these ##################

# Calculate summary statistics for accuracy data
accuracy_trialtypecategory_summary <- accuracy_trialtypecategory_data %>% 
  group_by(trialtype, category) %>% 
  do(summarise_data(., "PercentAccuracy"))

# Create plot for accuracy
accuracy_categorytrialtype_plot <- 
  ggplot(data = accuracy_trialtypecategory_summary, 
         aes(x = trialtype, 
             y = Mean, 
             fill = category, 
             colour = category))  +
  # Add chance level line at 50%
  geom_hline(yintercept = 50, 
             size = 1, 
             linetype = "dashed") +
  geom_dotplot(data = accuracy_trialtypecategory_data, 
               position = position_dodge(0.7),
               mapping = aes(x = trialtype, 
                             y = PercentAccuracy, 
                             fill = category),
               binaxis = 'y', 
               stackdir = 'center', colour = "black", 
               dotsize = 0.7, alpha = 0.4) +
  geom_line(mapping = aes(group = category), 
            position = position_dodge(0.7),
            linewidth = 1.5) +
  geom_errorbar(mapping = aes(ymin = Mean-SE, 
                              ymax = Mean+SE), 
                width = 0.3, 
                size = 1.2, 
                colour = "black", 
                position = position_dodge(0.7)) +
  geom_point(size = 7, 
             shape = 21, 
             colour = "black", 
             position = position_dodge(0.7)) +
  labs(x="Test Type", y="Accuracy", fill="", colour="") + 
  coord_cartesian(ylim = c(35, 115)) +
  scale_y_continuous(breaks = c(40, 60, 80, 100)) +
  scale_color_manual(values = factor_levels$category$colours) +
  scale_fill_manual(values = factor_levels$category$colours) +
  x_axis_theme + y_axis_theme + paper_facet_theme + blank_bg_theme + legend_theme


# Calculate summary statistics for RT data
rt_trialtypecategory_summary <- rt_trialtypecategory_data %>% 
  group_by(trialtype, category) %>% 
  do(summarise_data(., "MeanRT"))

# Create plot for reaction time
rt_categorytrialtype_plot <- 
  ggplot(data = rt_trialtypecategory_summary, 
         aes(x = trialtype, 
             y = Mean, 
             fill = category, 
             colour = category))  +
  # Add individual participant data points
  geom_dotplot(data = rt_trialtypecategory_data, 
               position = position_dodge(0.7),
               mapping = aes(x = trialtype, 
                             y = MeanRT, 
                             fill = category),
               binaxis = 'y', 
               stackdir = 'center', 
               colour = "black", 
               dotsize = 0.7, 
               alpha = 0.4) +
  # Add lines connecting means
  geom_line(mapping = aes(group = category), 
            position = position_dodge(0.7),
            linewidth = 1.5) +
  # Add error bars
  geom_errorbar(mapping = aes(ymin = Mean-SE, 
                              ymax = Mean+SE), 
                width = 0.3, 
                size = 1.2, 
                colour = "black", 
                position = position_dodge(0.7)) +
  # Add mean points
  geom_point(size = 7, 
             shape = 21, 
             colour = "black", 
             position = position_dodge(0.7)) +
  # Add labels and customize appearance
  labs(x="Test Type", y="Reaction Time (ms)", fill="", colour="") + 
  scale_color_manual(values = factor_levels$category$colours) +
  scale_fill_manual(values = factor_levels$category$colours) +
  x_axis_theme + y_axis_theme + paper_facet_theme + blank_bg_theme + legend_theme



# Calculate summary statistics for fixation count data
fixnum_categorytrialtypeitemtype_summary <- fixnum_categorytrialtypeitemtype_data %>% 
  group_by(category, trialtype, 
           ia_item_type, ia_itemtype_trialtype) %>% 
  do(summarise_data(., "mean_num_fix"))

fixnum_categorytrialtypeitemtype_plot <- 
  ggplot(data = fixnum_categorytrialtypeitemtype_summary, 
         aes(x = ia_itemtype_trialtype, 
             y = Mean, 
             fill = category, 
             colour = category)) +
  facet_wrap(~trialtype, scales="free_x") +
  # Add individual participant data points
  geom_dotplot(data = fixnum_categorytrialtypeitemtype_data, 
               position = position_dodge(0.7),
               mapping = aes(x = ia_itemtype_trialtype, 
                             y = mean_num_fix, 
                             fill = category),
               binaxis = 'y', 
               stackdir = 'center', 
               colour = "black", 
               dotsize = 0.7, 
               alpha = 0.4) +
  # Add lines connecting means
  geom_line(mapping = aes(group = category), 
            position = position_dodge(0.7),
            linewidth = 1.5) +
  # Add error bars
  geom_errorbar(mapping = aes(ymin = Mean-SE, 
                              ymax = Mean+SE), 
                width = 0.3, 
                size = 1.2, 
                colour = "black", 
                position = position_dodge(0.7)) +
  # Add mean points
  geom_point(size = 7, 
             shape = 21, 
             colour = "black", 
             position = position_dodge(0.7)) +
  # Add labels and customize appearance
  labs(x="Test Type", y="Mean Fixation Count", fill="", colour="") + 
  scale_color_manual(values = factor_levels$category$colours) +
  scale_fill_manual(values = factor_levels$category$colours) +
  x_axis_theme + y_axis_theme + paper_facet_theme + blank_bg_theme + legend_theme


# Create final figure arrangement with all plots
blank <- ggplot() + theme_void()
ggarrange(
  # Top row: accuracy and RT plots with blank space for alignment
  ggarrange(accuracy_categorytrialtype_plot, rt_categorytrialtype_plot, blank, nrow = 1,
            widths = c(1, 1, 0), legend = "none"),
  # Bottom row: fixation count plot centered
  ggarrange(blank, fixnum_categorytrialtypeitemtype_plot, blank, nrow = 1,
            widths = c(0.5, 1, 0.5), legend = "none"),
  nrow = 2,
  legend = "none"
)
