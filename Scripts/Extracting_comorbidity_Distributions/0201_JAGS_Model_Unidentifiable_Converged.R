################## Run Model to Retrieve Patient Level Data #####################
library(runjags)
# Load libraries for post-hoc transforms
library(tidyverse)
library(posterior)
library(MASS)
source("Scripts/Extracting_comorbidity_Distributions/0100_Data_Preparation.R")

# Number of individual patients with 1 or more LTCs
ipd_pop_size <- 18
# Replace NAs from individual patient data with 0s
data_ipd <- ifelse(is.na(data_ipd), 0, 1)[1:ipd_pop_size,]
# Create row sums for IPD
comorb_count_ipd <- rowSums(data_ipd)

# Correct data to remove patients with zero LTCs
comorb_count <- comorb_count[(markers[1]+1):pop_size]
pop_size <- pop_size - markers[1]
markers <- markers - markers[1]
markers <- markers[2:4]

model_unid_scaled <- "model{
# Aggregate Data  ---------------------------------------------------------

# The following for-loop takes the marginal counts for comorbidity loads of 1-2
# and sets them equal to the row sums by stating that the difference between them 
# is zero (with some small error for flexibility that aids convergence) and states 
# that they are drawn from a Poisson distribution with shared rate lambda (as 
# multimorbidity load is known to follow a Poisson distribution). The former ties 
# the imputed diseases in each row to the known distribution of row sums and the 
# latter allows us to infer the distribution of comorbidity counts in the 3+ bin 
# below.

for(i in 1:markers[2]){
zeroes[i] ~ dnorm(sum(pres[i,]) - comorb_count[i], 4)
comorb_count[i] ~ dpois(lambda)
}

# For comorbidity counts greater than 3 the data are binned. We load these values
# as NAs and impute them by stating that they are drawn from a truncated Poisson 
# distribution (i.e. all the values are 3 or larger) with rate lambda shared with
# the known comorbidity count distribution below 3.

for(j in (markers[2]+1):(markers[3])){
zeroes[j] ~ dnorm(sum(pres[j,]) - comorb_count[j], 4)
comorb_count[j] ~ dpois(lambda)T(3,11)
}

# For the column counts, there is no need to tie the values to two distributions
# (i.e. there is no functional form for the disease presences, we do not need to
# link to the column sums and a Poisson distribution, for instance). Thus, we say
# that the disease counts (data) are drawn from a peaked Normal distribution
# with a mean equal to the column sums. By adding a small amount of stochasticity
# convergence is faster.

for(D in 1:n_dis){ 
diseases[D] ~ dnorm(sum(pres[,D]), 0.5)
}

# In this for-loop, we draw a vector of latent disease values from a multivariate
# normal distribution (MVN). The MVN is a computationally efficient distribution 
# for drawing correlated values. As these values are continuous and unbounded, 
# we need to convert them to 1s or 0s to reflect the data (disease presence or 
# absence). We do this using the step function which takes values of 0 or higher
# in the latent space and converts them to 1, and values below 0 to 0. This
# is equivalent to a Probit link. Note discussion of MVN priors below. 

for(p in 1:pop_size){ 
Z[p, 1:n_dis] ~ dmnorm(mu[p, 1:n_dis], Sigma[,]) # Precision formulation of MVN
for(d in 1:n_dis){ 
pres[p,d] <- step(Z[p,d]) # If latent state is 0 or higher, disease is present
mu[p,d] <- intercept[d] # Disease specific, intercept only linear predictor
}
} 


# Individual Patient Data -------------------------------------------------

# These data are structured differently, with 1s and 0s rather than all NAs. To connect
# connect to data we need to introduce a stochastic step - a Bernoulli trial (see 
# innermost) loop. Beyond this adjustment, the process is identical above. The
# shared model attributes between the IPD and aggregated data are the correlation
# matrix (Sigma), the disease specific intercepts, and the comorbidity count rate 
# (lambda).

for(ipd in 1:ipd_pop_size){ 
Z_ipd[ipd, 1:n_dis] ~ dmnorm(mu_ipd[ipd, 1:n_dis], Sigma[,]) # Precision formulation of MVN
for(d_ipd in 1:n_dis){ 
data_ipd[ipd,d_ipd] ~ dbern(ifelse(step(Z_ipd[ipd,d_ipd]), 0.999, 0.001))
mu_ipd[ipd,d_ipd] <- intercept[d_ipd] # Disease specific, intercept only linear predictor
}

# Comorbidity counts are both equal to row sums and Poisson distributed with an 
# unknown rate - lambda
comorb_count_ipd[ipd] ~ dpois(lambda)
}


# Priors ------------------------------------------------------------------

# Disease specific intercept - minimally informative prior for variance around 1
for (D in 1:n_dis){
intercept[D] ~ dnorm(0, 0.2) 
}

# A vague prior for the rate of comorbidity counts which can be a maximum of 11
lambda ~ dgamma(0.1, 0.1)

# Rather than using the covariance parameterisation of the MV normal we are
# going to use the precision parameterisation with a Wishart distribution
# prior. We will *NOT* normalise the variances (i.e. set them to 1). Even 
# though the variances and means are non-identifiable in the probit model,
# the model is able to converge better with the added flexibility of not being 
# normalised. The variances and means are then  rescaled post hoc to make the means
# identifiable). The idea behind this is to improve mixing by making the means more
# flexible. This is very counter intuitive!!!

Sigma ~ dscaled.wishart(rep(1,n_dis), 50)

#monitor# Sigma, intercept, lambda
#data# pop_size, n_dis, diseases, comorb_count,zeroes, markers,ipd_pop_size,data_ipd, comorb_count_ipd
}"

# Run JAGS code
results_unid_scaledwish <- run.jags(model_unid_scaled, 
                         n.chains = 10,
                         # Adapatation in this model does not seem possible and makes model output less predictable
                         adapt = 0, 
                         burnin = 3000, 
                         sample = 50000,
                         modules = 'glm', 
                         method = "parallel")

## RESCALE TO MAKE IDENTIFIABLE
# Convert from JAGS to dataframe
results_unid_df <- as_draws_df(coda::as.mcmc.list(results_unid_scaledwish))

# Extract samples for sigmas and intercepts
sigmas <- results_unid_df %>% dplyr::select(starts_with("Sigma")) %>% as.data.frame()
means <- results_unid_df %>% dplyr::select(starts_with("intercept")) %>% as.data.frame()

# Slow function which takes samples and rescales them to satisfy var = 1 identifiability constraint
rescale_mvprobit <- function(index, prec_vec, mean_vec){
  prec_vec <- slice(prec_vec, index) %>% unlist() # Remove extraneous info
  mean_vec <- slice(mean_vec, index) %>% unlist() # Remove extraneous info
  dimensions <- 11 # Number of diseases
  # Assemble precisions into matrix
  prec_mat <- matrix(c(prec_vec), nrow = dimensions, ncol = dimensions, byrow = FALSE)
  # Invert precision to get covariance
  cov_mat <- ginv(prec_mat)
  # Rescale covariance to get correlation
  cor_mat <- cov2cor(cov_mat)
  # Transform to vector of correlations
  cor_vec <- c(cor_mat)
  # Rescale means based on square root of covariance diagonal
  mean_vec_resc <- (mean_vec/sqrt(diag(cov_mat)))
  attributes(mean_vec_resc) <- NULL  # Remove extraneous info
  # Combine results into single vector
  combi_vec <- c(cor_vec, mean_vec_resc) 
  combi_vec
}

# Iterate over samples and rescale one row at a time
rescaled_samples_wide <- 1:nrow(results_unid_df) %>% 
  map_dfc(~rescale_mvprobit(.,(sigmas), (means)))

# Convert from one column per iteration to one row per iteration
rescaled_samples <- t(rescaled_samples_wide)

# Place in dataframe
rescaled_samples_df <- as.data.frame(rescaled_samples)
# Add column names
names(rescaled_samples_df) <- c(names(sigmas), names(means))

# Add chain and iteration info
rescaled_samples_df$.chain <-results_unid_df$.chain
rescaled_samples_df$.iteration <-results_unid_df$.iteration

# Produce rescaled summary
summ_resc <- summarise_draws(rescaled_samples_df)

# Extract only the upper triangle of correlations
sigm <- summ_resc %>% filter(str_detect(variable, "Sigma"))
sigm_tri <- summ_resc %>% 
  filter(variable %in% matrix(sigm$variable, 11,11)[upper.tri(matrix(sigm$variable, 11,11))])
# Extract intercepts
intc <- summ_resc %>% filter(str_detect(variable, "inter"))
# Extract lambda
lmd <- results_unid_df %>% dplyr::select(starts_with("lambda")) %>% as.data.frame() %>% summarise_draws()
# Glue back together
summ_resc <- rbind(sigm_tri, intc, lmd)
