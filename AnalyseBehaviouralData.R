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
accuracy_trialtypecategory_anova <- ezANOVA(data = accuracy_trialtypecategory_data,
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
accuracy_category_interac_posthoc <- correct_p_vals(accuracy_category_interac_posthoc, 
                                                    "p")

################## Analyse test phase reaction time data ##################

# Calculate RT metrics grouped by participant, trial type, and category
rt_trialtypecategory_data <- test_data %>% 
  filter(accuracy == 1) %>% 
  group_by(Participant, trialtype, category) %>% 
  summarise(MeanRT = mean(test_RT, na.rm = TRUE))

# Perform 2x2 repeated measures ANOVA on RT data
# Factors: category and trial type
rt_trialtypecategory_anova <- ezANOVA(data = rt_trialtypecategory_data,
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
rt_trialtype_posthoc <- correct_p_vals(rt_trialtype_posthoc, 
                                       pvalcol = "p")


# Perform post-hoc tests for category effects within each trial type
rt_category_posthoc <- full_join(
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
rt_category_posthoc <- correct_p_vals(rt_category_posthoc, 
                                      pvalcol = "p")

################## Analyse eye-tracking data ##################


################## Plot these ##################

# Calculate summary statistics for accuracy data
accuracy_trialtypecategory_summary <- accuracy_trialtypecategory_data %>% 
  group_by(trialtype, category) %>% 
  do(summarise_data(., "PercentAccuracy"))

# Create plot for accuracy
(accuracy_categorytrialtype_plot <- ggplot(data = accuracy_trialtypecategory_summary, 
                                           aes(x = trialtype, 
                                               y = Mean, 
                                               fill = category, 
                                               colour = category))  +
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
    x_axis_theme + y_axis_theme + paper_facet_theme + blank_bg_theme + legend_theme)


# Calculate summary statistics for RT data
rt_trialtypecategory_summary <- rt_trialtypecategory_data %>% 
  group_by(trialtype, category) %>% 
  do(summarise_data(., "MeanRT"))

# Create plot for reaction time
(rt_categorytrialtype_plot <- ggplot(data = rt_trialtypecategory_summary, 
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
    x_axis_theme + y_axis_theme + paper_facet_theme + blank_bg_theme + legend_theme)


ggarrange(accuracy_categorytrialtype_plot,
          rt_categorytrialtype_plot,
          common.legend = TRUE, 
          legend = "bottom") +
  legend_theme

