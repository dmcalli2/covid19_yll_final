# Aggregate Data Preparation ----------------------------------------------

# March 26th Italian data
full_dat <- read.table("Data/italy_update.txt",   sep="\t", header=TRUE)

## Extract disease counts
diseases <- full_dat[1:11,2] # Note, these values are sensitive to form in which data are imported
n_dis <- length(diseases)

## Extract comorbidity counts
# Number of patients in each comorbidity count bin
comorb_bins <- full_dat[13:16,2]
# Row numbers at which comorbidity counts change
markers  <- cumsum(comorb_bins)
# Extract population size
pop_size <- sum(comorb_bins)

## Generate vector of comorbidity counts
# Initialise data storage
comorb_count <- rep(NA, pop_size)
# Turn into comorbidity vector
comorb_count[1:markers[1]] <- 0
comorb_count[(markers[1]+1):markers[2]] <- 1
comorb_count[(markers[2]+1):markers[3]] <- 2
# comorb_count[markers[3]+1:markers[4]] <- NA

# Coding solution to make row sums and comorbidity counts match in JAGS 
zeroes <-  rep(0, pop_size)


# Individual Patient Data (IPD) --------------------------------------------
#' 
#' To protect patient confidentiality we are unable to share these data here,
#' however, the code below will allow the model that incorporates individual 
#' patient data to run for testing. We also used the IPD as presence only data,
#' i.e. where a disease was identified as present it was included in the data as
#' a 1, otherwise it was treated as an NA. This further protects patient ID and 
#' allows the model to account for undiagnosed diseases.
#'
data_ipd <- matrix(NA, nrow = 1, ncol = n_dis)

#' 
#' The following code would be identical for real IPD assuming the data are in 
#' the simulated form above
#' 
# Number of patients in individual patient dataset 
ipd_pop_size <- nrow(data_ipd)

# Coding solution to make row sums and comorbidity counts match in JAGS 
zeroes_ipd <- rep(0,ipd_pop_size)

#'
#'  In the IPD, comorbidity counts are unknown due to the treatment of the 
#'  data as presence-only described above. We therefore also treat the comorbidity
#'  counts as NAs to be imputed for the IPD.
comorb_count_ipd <- rep(NA,ipd_pop_size)


