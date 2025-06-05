# Set working directory and source initialization file
setwd("~/Downloads/CategorySpecificAssociativeInference/")
source("Initialise.R")

################## Read in behavioural data and organise it ##################

# Initialize empty dataframe to store all participants' data
behavioural_data <- c()

# Loop through all participants and read their task event files
for(part in aim_participants){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(data_path, part, "/func/", part, "_task-aim_events.tsv"),
                          header=TRUE,
                          sep = "\t")
  # Add participant ID column
  part_data$Participant <- part
  # Combine with existing data
  behavioural_data <- bind_rows(behavioural_data, part_data)
}

# Reorder columns to put Participant ID first
behavioural_data <- relocate(behavioural_data,
                             Participant)

# Verify all intended participants are included in the dataset
all(sort(unique(behavioural_data$Participant)) == sort(aim_participants))


################## Analyse test phase trials without a response ##################

# Filter data to include only test phase trials and convert categorical variables to factors
test_data <- behavioural_data %>% 
  filter(phase == "test") %>% 
  mutate(category = factor(category,
                           levels = factor_levels$category$labels),
         trialtype = factor(trialtype,
                            levels = factor_levels$pairtype_test$levels))

# Count trials with no response (where RT = 0) for each participant
no_resp_data <- test_data %>% 
  mutate(Participant = factor(Participant)) %>% 
  filter(is.na(accuracy)) %>%          # Select only trials with no response
  group_by(Participant, .drop = F) %>%  # Group by participant, keeping all factor levels
  summarise(num_no_resp = length(Participant))  # Count no-response trials

# Calculate descriptive statistics for no-response trials across participants
no_resp_summary <- no_resp_data %>% 
  ungroup() %>% 
  do(summarise_data(., "num_no_resp"))

# Count no-response trials by participant, trial type, and category
noresp_trialtypecategory_data <- test_data %>% 
  mutate(Participant = factor(Participant)) %>% 
  filter(test_RT == 0) %>% 
  group_by(Participant, trialtype, category, .drop = F) %>% 
  summarise(num_no_resp = length(Participant)) 

# Perform 2x2 repeated measures ANOVA on no-response data
# Factors: trial type and category, within-subjects design
noresp_trialtypecategory_anova <-
  ezANOVA(data = noresp_trialtypecategory_data,
          dv = num_no_resp,           # Dependent variable: number of no-response trials
          within = c(trialtype, category),  # Within-subject factors
          wid = Participant,          # Participant identifier
          detailed = TRUE)            # Request detailed output

# Count no-response trials by participant and trial type (collapsed across categories)
# to examine the main effect of trialtype
noresp_trialtype_data <- test_data %>% 
  mutate(Participant = factor(Participant)) %>% 
  filter(test_RT == 0) %>% 
  group_by(Participant, trialtype, .drop = F) %>% 
  summarise(num_no_resp = length(Participant)) 

# Perform post-hoc tests comparing trial types for no-response trials
noresp_trialtype_main_posthoc <- full_join(
  noresp_trialtype_data %>% 
    as.data.frame() %>% 
    t_test(num_no_resp ~ trialtype,
           paired = TRUE),
  # Calculate effect sizes (Cohen's d)
  noresp_trialtype_data %>% 
    as.data.frame() %>% 
    cohens_d(num_no_resp ~ trialtype,
             paired = TRUE),
  by = join_by(.y., group1, group2, n1, n2)
) %>% 
  as.data.frame()

# Apply multiple comparison correction to p-values 
noresp_trialtype_main_posthoc <- 
  correct_p_vals(noresp_trialtype_main_posthoc, 
                 "p")


################## Analyse test phase accuracy data ##################

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
accuracy_trialtypecategory_interac_posthoc <- full_join(
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
accuracy_trialtypecategory_interac_posthoc <- 
  correct_p_vals(accuracy_trialtypecategory_interac_posthoc, 
                 "p")

# Calculate adjusted AC score from the appendix
# Filter test data to include only trials where both AB and BC were answered correctly
dir_correct_testdata <- test_data %>% 
  pivot_wider(id_cols = c(Participant, category, pairnum, block),  # Reshape data: one row per unique trial
              names_from = trialtype,                             # Columns for each trial type (AB, BC, AC)
              values_from = accuracy) %>% 
  filter(AB == 1 & BC == 1)                                       # Keep only trials where both AB and BC are correct

# Calculate adjusted AC accuracy for each participant and category
adjaccuracy_trialtypecategory_data <- dir_correct_testdata %>% 
  group_by(Participant, category) %>% 
  summarise(
    TotalTrials = length(Participant),                             # Total number of eligible trials
    CorrectTrials = sum(AC, na.rm = TRUE),                        # Number of correct AC responses
    PercentAccuracy = (CorrectTrials/TotalTrials)*100              # Percent accuracy for AC
  ) %>% 
  select(c(Participant, category, PercentAccuracy)) %>%            # Keep only relevant columns
  rename(AC = PercentAccuracy) %>%                                 # Rename for clarity
  # Join with raw AB and BC accuracy data for the same participants/categories
  full_join(
    accuracy_trialtypecategory_data %>% 
      filter(trialtype %in% c("AB", "BC")) %>% 
      pivot_wider(id_cols = c(Participant, category), 
                  names_from = trialtype,
                  values_from = PercentAccuracy),
    by = c("Participant", "category")
  ) %>% 
  # Reshape so AB, BC, and AC accuracy are all in one column, with a new column indicating trial type
  pivot_longer(
    cols = c(AB, BC, AC),
    names_to = "trialtype",
    values_to = "PercentAccuracy"
  ) %>% 
  # Convert trialtype to a factor with specified levels for plotting/analysis
  mutate(trialtype = factor(trialtype,
                            levels = factor_levels$pairtype_test$levels))

# Perform 2x2 repeated measures ANOVA on adjusted accuracy data
# Factors: category and trial type
adjaccuracy_trialtypecategory_anova <- 
  ezANOVA(data = adjaccuracy_trialtypecategory_data,
          dv = PercentAccuracy,
          within = c(category, trialtype),
          wid = Participant,
          detailed = TRUE)

# Calculate adjustedaccuracy metrics grouped by participant and trial type
adjaccuracy_trialtype_data <- dir_correct_testdata %>% 
  group_by(Participant) %>% 
  summarise(
    TotalTrials = length(Participant),           # Total number of eligible trials for this participant
    CorrectTrials = sum(AC, na.rm = TRUE),      # Number of correct AC responses
    PercentAccuracy = (CorrectTrials/TotalTrials)*100  # Percent accuracy for AC
  ) %>% 
  select(c(Participant, PercentAccuracy)) %>%    # Keep only relevant columns
  rename(AC = PercentAccuracy) %>%               # Rename for clarity
  # Join with raw AB and BC accuracy data for the same participants
  full_join(
    accuracy_trialtype_data %>% 
      filter(trialtype %in% c("AB", "BC")) %>% 
      pivot_wider(
        id_cols = c(Participant),                # One row per participant
        names_from = trialtype,                  # Columns for AB and BC
        values_from = PercentAccuracy
      ),
    by = c("Participant")
  ) %>% 
  # Reshape so AB, BC, and AC accuracy are all in one column, with a new column indicating trial type
  pivot_longer(
    cols = c(AB, BC, AC),
    names_to = "trialtype",
    values_to = "PercentAccuracy"
  ) %>% 
  # Convert trialtype to a factor with specified levels for plotting/analysis
  mutate(trialtype = factor(trialtype,
                            levels = factor_levels$pairtype_test$levels))

# Perform post-hoc paired t-tests and calculate Cohen's d for trial type main effect
adjaccuracy_trialtype_main_posthoc <- full_join(
  # Paired t-test between trial types
  adjaccuracy_trialtype_data %>%
    as.data.frame() %>%
    t_test(PercentAccuracy ~ trialtype,
           paired = TRUE),
  # Calculate effect size (Cohen's d)
  adjaccuracy_trialtype_data %>%
    as.data.frame() %>%
    cohens_d(PercentAccuracy ~ trialtype,
             paired = TRUE),
  by = join_by(.y., group1, group2, n1, n2)
)

# Perform post-hoc tests for trial type effects within each category
adjaccuracy_trialtypecategory_interac_posthoc <- full_join(
  # Paired t-tests between trial types within each categories
  adjaccuracy_trialtypecategory_data %>% 
    as.data.frame() %>% 
    group_by(category) %>% 
    t_test(PercentAccuracy ~ trialtype,
           paired = TRUE),
  # Calculate effect sizes (Cohen's d)
  adjaccuracy_trialtypecategory_data %>% 
    as.data.frame() %>% 
    group_by(category) %>% 
    cohens_d(PercentAccuracy ~ trialtype,
             paired = TRUE),
  by = join_by(category, .y., group1, group2, n1, n2)
) %>% 
  as.data.frame()

# Apply multiple comparison correction to p-values
adjaccuracy_trialtypecategory_interac_posthoc <- 
  correct_p_vals(adjaccuracy_trialtypecategory_interac_posthoc, 
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
rt_trialtype_main_posthoc <- full_join(
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
rt_trialtype_main_posthoc <- 
  correct_p_vals(rt_trialtype_main_posthoc, 
                 pvalcol = "p")


# Perform post-hoc tests for category effects within each trial type
rt_categorytrialtype_interac_posthoc <- full_join(
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
rt_categorytrialtype_interac_posthoc <- 
  correct_p_vals(rt_categorytrialtype_interac_posthoc, 
                 pvalcol = "p")

################## Analyse eye-tracking data ##################

# Initialize empty dataframe to store all participants' data
et_data <- c()

# Loop through all participants and read their eye-tracking files
for(part in aim_participants){
  # Read TSV file containing eye-tracking data for each participant
  part_data <- read.table(paste0(data_path, part, "/func/", part, "_task-aim_eyetracking.tsv"),
                          header = TRUE,
                          sep = "\t")
  # Add participant ID column
  part_data$Participant <- part
  # Combine with existing data
  et_data <- bind_rows(et_data, part_data)
}

# Reorder columns to put Participant ID first
et_data <- relocate(et_data,
                    Participant)

# Verify all intended participants are included in the dataset
all(sort(unique(et_data$Participant)) == sort(aim_participants))

# Create new variables and convert categorical variables to factors
et_data <- et_data %>% 
  mutate(
    # Combine item type and pair type into a single identifier
    ia_itemtype_pairtype = paste(ia_item_type, pairtype, sep = "_"),
    # Convert category to factor with predefined levels
    category = factor(category,
                      levels = factor_levels$category$labels),
    # Convert pair type to factor with predefined study phase levels
    pairtype = factor(pairtype,
                      levels = factor_levels$pairtype_study$levels),
    # Convert combined identifier to factor with predefined levels and labels
    ia_itemtype_pairtype = factor(ia_itemtype_pairtype,
                                  levels = factor_levels$studyitem_pairtype$levels,
                                  labels = factor_levels$studyitem_pairtype$labels))

# Calculate number of fixations for each trial and interest area
fixnum_data <- et_data %>% 
  # Group by all relevant variables to maintain trial-level data
  group_by(Participant, block, trial, category, pairtype, ia_item_type, ia_itemtype_pairtype) %>% 
  # Count number of fixations
  summarise(num_fix = length(ia_item_type)) 

# Calculate mean number of fixations across trials for each condition
fixnum_categorypairtypeitemtype_data <- fixnum_data %>% 
  # Group by participant and condition variables
  group_by(Participant, category, pairtype, ia_item_type, ia_itemtype_pairtype) %>% 
  # Calculate mean number of fixations
  summarise(mean_num_fix = mean(num_fix))

# Perform 3-way repeated measures ANOVA on fixation count data
fixnum_categorypairtypeitemtype_anova <- 
  ezANOVA(data = fixnum_categorypairtypeitemtype_data, 
          wid = Participant, 
          within = c(category, pairtype, ia_item_type),
          dv = mean_num_fix,
          detailed = TRUE)

# Posthocs for the category x item type interaction
fixnum_categoryitemtype_data <- fixnum_data %>% 
  group_by(Participant, category, ia_item_type) %>% 
  summarise(mean_num_fix = mean(num_fix))

fixnum_categoryitemtype_posthoc <- bind_rows(
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
    # Rename ia_item_type column to 'within' for bind_rows
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
    # Rename category column to 'within' for bind_rows
    dplyr::rename(within = category)
) %>% 
  # Need this for correct_p_vals to work
  as.data.frame()

# Apply multiple comparison correction to p-values
fixnum_categoryitemtype_posthoc <- 
  correct_p_vals(fixnum_categoryitemtype_posthoc,
                 "p")

# Posthocs for the pair type x item type interaction
fixnum_pairtypeitemtype_data <- fixnum_data %>% 
  group_by(Participant, pairtype, ia_item_type) %>% 
  summarise(mean_num_fix = mean(num_fix))

# Perform post-hoc tests for the pair type x item type interaction
fixnum_pairtypeitemtype_posthoc <- bind_rows(
  # First analysis: Compare pair types within each item type
  full_join(
    # Conduct paired t-tests comparing pair types for each item type
    fixnum_pairtypeitemtype_data %>% 
      group_by(ia_item_type) %>% 
      t_test(mean_num_fix ~ pairtype,
             paired = TRUE),
    # Calculate effect sizes (Cohen's d) for the same comparisons
    fixnum_pairtypeitemtype_data %>% 
      group_by(ia_item_type) %>% 
      cohens_d(mean_num_fix ~ pairtype,
               paired = TRUE),
    # Join t-test results with effect sizes
    by = join_by(.y., ia_item_type, group1, group2, n1, n2)
  ) %>% 
    # Rename ia_item_type column to 'within' for bind_rows
    rename(within = ia_item_type),
  # Second analysis: Compare item types within each pair type
  full_join(
    # Conduct paired t-tests comparing item types for each trial type
    # Conduct paired t-tests comparing item types for each pair type
    fixnum_pairtypeitemtype_data %>% 
      group_by(pairtype) %>% 
      t_test(mean_num_fix ~ ia_item_type,
             paired = TRUE),
    # Calculate effect sizes (Cohen's d) for the same comparisons
    fixnum_pairtypeitemtype_data %>% 
      group_by(pairtype) %>% 
      cohens_d(mean_num_fix ~ ia_item_type,
               paired = TRUE),
    # Join t-test results with effect sizes
    by = join_by(.y., pairtype, group1, group2, n1, n2)
  ) %>% 
    # Rename trialtype column to 'within' for bind_rows
    dplyr::rename(within = pairtype)
) %>% 
  # Need this for correct_p_vals to work
  as.data.frame()

# Apply multiple comparison correction to p-values
fixnum_pairtypeitemtype_posthoc <- 
  correct_p_vals(fixnum_pairtypeitemtype_posthoc,
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

# Calculate summary statistics for adjusted accuracy data
adjaccuracy_trialtypecategory_summary <- adjaccuracy_trialtypecategory_data %>% 
  group_by(trialtype, category) %>% 
  do(summarise_data(., "PercentAccuracy"))

# Create plot for accuracy
adjaccuracy_categorytrialtype_plot <- 
  ggplot(data = adjaccuracy_trialtypecategory_summary, 
         aes(x = trialtype, 
             y = Mean, 
             fill = category, 
             colour = category))  +
  # Add chance level line at 50%
  geom_hline(yintercept = 50, 
             size = 1, 
             linetype = "dashed") +
  geom_dotplot(data = adjaccuracy_trialtypecategory_data, 
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
fixnum_categorypairtypeitemtype_summary <- fixnum_categorypairtypeitemtype_data %>% 
  group_by(category, pairtype, 
           ia_item_type, ia_itemtype_pairtype) %>% 
  do(summarise_data(., "mean_num_fix"))

fixnum_categorytypeitemtype_plot <- 
  ggplot(data = fixnum_categorypairtypeitemtype_summary, 
         aes(x = ia_itemtype_pairtype, 
             y = Mean, 
             fill = category, 
             colour = category)) +
  facet_wrap(~pairtype, scales="free_x") +
  # Add individual participant data points
  geom_dotplot(data = fixnum_categorypairtypeitemtype_data, 
               position = position_dodge(0.7),
               mapping = aes(x = ia_itemtype_pairtype, 
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
  labs(x="Pair Type", y="Mean Fixation Count", fill="", colour="") + 
  scale_color_manual(values = factor_levels$category$colours) +
  scale_fill_manual(values = factor_levels$category$colours) +
  x_axis_theme + y_axis_theme + paper_facet_theme + blank_bg_theme + legend_theme


# Create final figure arrangement with all plots
blank <- ggplot() + theme_void()
(behavioural_plots <- ggarrange(
  # Top row: accuracy and RT plots with blank space for alignment
  ggarrange(accuracy_categorytrialtype_plot +
              labs(title = "A) Recognition Accuracy") + 
              title_theme, 
            rt_categorytrialtype_plot +
              labs(title = "B) Reaction Time") + 
              title_theme, 
            blank, 
            nrow = 1,
            widths = c(1, 1, 0), 
            legend = "none"),
  # Bottom row: fixation count plot centered
  ggarrange(blank, 
            fixnum_categorytypeitemtype_plot +
              labs(title = "C) Eye-movement\n Behaviour") + 
              title_theme,
            blank, nrow = 1,
            widths = c(0.5, 1, 0.5), 
            legend = "bottom"),
  nrow = 2,
  legend = "none"
))

# Create from appendix
(appendix_accuracy_plots <- ggarrange(accuracy_categorytrialtype_plot +
                                        labs(title = "Raw AC Accuracy") + 
                                        title_theme, 
                                      adjaccuracy_categorytrialtype_plot + 
                                        theme(axis.title.y = element_blank()) +
                                        labs(title = "Adjusted AC Accuracy") + 
                                        title_theme, 
                                      nrow = 1,
                                      legend = "bottom", common.legend = TRUE))
