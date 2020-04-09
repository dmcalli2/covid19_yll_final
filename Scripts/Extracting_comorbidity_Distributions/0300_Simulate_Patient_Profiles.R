library(tidybayes)
library(tidyverse)
library(MASS)

# Sample Posterior --------------------------------------------------------

# Read in runjags object
results <- readRDS("Data/JAGSresults")
# Convert runjags object to tidy dataframe
tidy_result <- tidy_draws(results)

# Useful parameters
n_dis <- 11 # Number of disease
sample_size <- 10000 # Number of patient profiles drawn

# Sample rows from posterior dataframe
posterior_samples <- sample_n(tidy_result, sample_size) 


# Define Profile Generating Function --------------------------------------

# This function takes a row from the dataframe above, extracts and reassembles 
# the correlation matrix and vector of means, passes it to a multivariate normal 
# and binarises the resulting vector.
get_patient_profile <- function(posterior_draw,n_dis){
  corr <- matrix(unlist(dplyr::select(posterior_draw, starts_with("Sigma"))), n_dis, n_dis)
  intercept <-as.vector(unlist(dplyr::select(posterior_draw, starts_with("intercept"))[1,]))
  cont_prof <- MASS::mvrnorm(n = 1, mu = intercept, Sigma = as.matrix(corr))
  bin_prof <- ifelse(cont_prof >= 0, 1, 0)
  bin_prof
}


# Generate Profiles -------------------------------------------------------

# Initialise data storage
patient_profiles <- matrix(NA, nrow = sample_size, ncol = n_dis)
# Loop over posterior samples and feed into patient profile generating function.
for(i in seq_along(1:sample_size)){
  patient_profiles[i,1:n_dis] <- get_patient_profile(posterior_samples[i,], n_dis)
}

# Fix column names
patient_profiles_df <- as.data.frame(patient_profiles)
names(patient_profiles_df) <- dis_names <- c("IHD", "Atr_Fib", "H_Fail", "Stroke", "Hypert", "Diab", "Demen", "COPD", "Cancer", "Liver", "Renal")

# Save Object
# saveRDS(patient_profiles_df, "Data/SimulatedCovidPatientProfiles")



