# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/CategorySpecificAssociativeInference/")
source("../InitialisePaths.R")

################## Read in data and organise it ##################

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


################## Analyse Test Phase Data ##################

# Filter data to include only test phase trials
test_data <- behavioural_data %>% 
  filter(phase == "test")

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
             paired = TRUE)
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
             paired = TRUE)
) %>% 
  as.data.frame()

# Apply multiple comparison correction to p-values
accuracy_category_interac_posthoc <- correct_p_vals(accuracy_category_interac_posthoc, 
                                                    "p")

