## Fit model to age data to obtain an age distribution, and to allow us to
## assume an association between age and comorbidity

library(tidyverse)
library(runjags)

## Run separately for men and women ----
sex <- "male"
## run again with larger SD around effect size of comorbidity on age
eff_size_sd <- 0.5
# eff_size_sd <- 1

## Read data ----
italy <- readRDS("Data/SimulatedProfiles.Rds")
names(italy) <- str_to_lower(names(italy))
italy <- as.matrix(italy)
italy_count <- rowSums(italy)
rm(italy)
age_sex <- read_tsv("Data/age_distribution_italian_deaths.txt")

if(sex == "male") {
  age_obs <- age_sex %>% 
    filter(sex == "men")
  sex_eff_size <- 5
} else {
  age_obs <- age_sex %>% 
    filter(sex == "women")
  sex_eff_size <- 3
}

age_obs <- rep(age_obs$age_point, age_obs$n)

## Run model with estimate for association between comorbidity count and mean age.
## Estimate age distribution in people without comorbidity (alpha) and 
## use this to predict the age distribution for each combination of comorbidities
## Derive variance for each comorbidity group
## Note uses normal distribution. Similar results for the mean and SD were obtained where the likelihood was a truncated normal (below)
n_como <- table(italy_count) %>% 
  as.vector()
p_como <- n_como/sum(n_como)
como_df <- tibble(como = (0:11), p_como = p_como, n_como) %>% 
  mutate(como_trnc = if_else(como <=6, como, 6L))  %>% 
  group_by(como_trnc) %>% 
  summarise(p_como = sum(p_como),
            n_como = sum(n_como)) %>% 
  ungroup() %>% 
  mutate(como_final = como_trnc* p_como) 
rm(n_como, p_como)

## Calculate component of linear predictor based on comorbidity count
## and proportion of patients with that comorbidity count
weighted_sum <- como_df %>% 
  pull(como_final) %>% 
  sum()
## Need number of groups to calcualte variance within each comorbidity group
n_grps <- nrow(como_df)

## Estimate age distribution, ignoring comorbidity ----
modelstring <- 
  "model{
  ## Loop through patients
 for(i in 1:n){
  age_obs[i] ~ dnorm(mean_age[i], prec)
  mean_age[i] <- alpha 
 }
## Priors on prec and alpha
alpha  ~ dnorm(0, 0.01)
prec <- 1/tau
tau <- s^2
s ~ dnorm(0, 0.01)T(0,)
## Obtain prediction not conditional on comorbidity
pred ~ dnorm(alpha, prec)T(50, 110)
}"

b <- run.jags(modelstring, data = list(age_obs = age_obs,
                                       n = length(age_obs),
                                       weighted_sum = weighted_sum),
              monitor = c("pred" ), n.chains = 2)
b
saveRDS(b, paste0("Data/jags_model_age_flu1", sex,".Rds"))
age_select_no_comorbid <- b$mcmc %>% as.matrix()
age_select_no_comorbid <- age_select_no_comorbid [1:10000,1]

saveRDS(age_select_no_comorbid, paste0("Data/age_selection_no_comorbid_", sex,".Rds"))

## Estimate age distribution, ignoring comorbidity, but this time with truncation ----
modelstring <- 
  "model{
  ## Loop through prediction matrix
  ## 
 for(i in 1:n){
  age_obs[i] ~ dnorm(mean_age[i], prec)T(20, 95)
  mean_age[i] <- alpha 
 }
## Priors on sd and alpha
alpha  ~ dnorm(0, 0.01)
prec <- 1/tau
tau <- s^2
s ~ dnorm(0, 0.01)T(0,)

pred ~ dnorm(alpha, prec)T(50, 110)

}"

b_trunc <- run.jags(modelstring, data = list(age_obs = age_obs,
                                       n = length(age_obs),
                                       n_grps = n_grps,
                                       weighted_sum = weighted_sum),
              monitor = c("pred" ), n.chains = 2)
b_trunc
saveRDS(b_trunc, paste0("Data/jags_model_age_flu1b_trunc", sex,".Rds"))
        
## Estimate age-comorbidity distribution, using knowledge of association between age and comorbidity count for influenza deaths ----
modelstring <- 
  "model{
  ## Loop through ISS patients
 for(i in 1:n){
  age_obs[i] ~ dnorm(mean_age[i], prec)
  mean_age[i] <- alpha + como_effect * weighted_sum
 }
## Priors on sd and alpha
alpha  ~ dnorm(0, 0.01)
prec <- 1/tau
tau <- s^2
s ~ dnorm(0, 0.01)T(0,)
# effect of comorbidity on age
como_effect ~ dnorm(eff_size, 1/eff_size_sd^2)

## Mean for each comorbidity group
m[1] <- alpha
for(g in 2:n_grps){
 m[g] <- alpha + (g-1)*como_effect 
}

## Calculate SD in each group, assumes variance is the same in each comorbidity group
global_mean <- alpha + como_effect * weighted_sum
total_ss <- tau * (n-1)
group_ss <- total_ss - sum(n_grps*(m - global_mean)^2)
group_var <- group_ss/ (n-1)
group_prec <- 1/group_var

## sample to obtain 10000 patients in each comorbidity group
for(g in 1:n_grps){
age[g] ~ dnorm(m[g], group_prec)T(50,110)
}

}"

d <- run.jags(modelstring, data = list(age_obs = age_obs,
                                       n = length(age_obs),
                                       n_grps = n_grps,
                                       eff_size = sex_eff_size,
                                       eff_size_sd = eff_size_sd,
                                       weighted_sum = weighted_sum),
              monitor = c("tau", "group_var", "age", "global_mean" ), n.chains = 2)
d
saveRDS(d, paste0("Data/jags_model_age_flu2", sex,".Rds"))
# plot(d, plot.type = "trace")

## select patients by comorbidity  ----
res <- as.matrix(d$mcmc)
res <- res [ , paste0("age[", 1:7, "]")]

italy_count_indx <- if_else(italy_count >=6, 6, italy_count)
italy_count_indx <- italy_count_indx +1
## Select using matrix indexing
italy_count_indx <- cbind(1:10000, italy_count_indx)
age_slct <- res[italy_count_indx]

saveRDS(age_slct, paste0("Data/age_selection_", sex, "_", str_replace(eff_size_sd, "\\.", "_"),".Rds"))
