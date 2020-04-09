################## Run Model to Retrieve Patient Level Data #####################
library(runjags)
source("Scripts/Extracting_Joint_Distributions/0100_Data_Preparation.R")

model <- "model{
# Aggregate Data  ---------------------------------------------------------

# The following for-loop takes the marginal counts for comorbidity loads of 0-2
# and sets them equal to the row sums by stating that the difference between them 
# is zero (with some small error for flexibility that aids convergence) and states 
# that they are drawn from a Poisson distribution with shared rate lambda (as 
# multimorbidity load is known to follow a Poisson distribution). The former ties 
# the imputed diseases in each row to the known distribution of row sums and the 
# latter allows us to infer the distribution of comorbidity counts in the 3+ bin 
# below.

for(i in 1:markers[3]){
zeroes[i] ~ dnorm(sum(pres[i,]) - comorb_count[i], 100)
comorb_count[i] ~ dpois(lambda)
}

# For comorbidity counts greater than 3 the data are binned. We load these values
# as NAs and impute them by stating that they are drawn from a truncated Poisson 
# distribution (i.e. all the values are 3 or larger) with rate lambda shared with
# the known comorbidity count distribution below 3.

for(j in (markers[3]+1):(markers[4])){
zeroes[j] ~ dnorm(sum(pres[j,]) - comorb_count[j], 100)
comorb_count[j] ~ dpois(lambda)T(3,11)
}

# For the column counts, there is no need to tie the values to two distributions
# (i.e. there is no functional form for the disease presences, we do not need to
# link to the column sums and a Poisson distribution, for instance). Thus, we say
# that the disease counts (data) are drawn from a highly peaked Normal distribution
# with a mean equal to the column sums. By adding a small amount of stochasticity
# convergence is faster.

for(D in 1:n_dis){ 
diseases[D] ~ dnorm(sum(pres[,D]),100)
}

# In this for-loop, we draw a vector of latent disease values from a multivariate
# normal distribution (MVN). The MVN is a computationally efficient distribution 
# for drawing correlated values. As these values are continuous and unbounded, 
# we need to convert them to 1s or 0s to reflect the data (disease presence or 
# absence). We do this using the step function which takes values of 0 or higher
# in the latent space and converts them to 1, and values below 0 to 0. This
# is equivalent to a Probit link. Note discussion of MVN priors below. 

for(p in 1:pop_size){ 
Z[p, 1:n_dis] ~ dmnorm.vcov(mu[p, 1:n_dis], Sigma[,]) # Variance formulation of MVN
for(d in 1:n_dis){ 
pres[p,d] <- step(Z[p,d]) # If latent state is 0 or higher, disease is present
mu[p,d] <- intercept[d] # Disease specific, intercept only linear predictor
}
} 


# Individual Patient Data -------------------------------------------------

# These data are structured differently, with 1s rather than all NAs. To connect
# the 1s as data we need to introduce a stochastic step - a Bernoulli trial (see 
# innermost) loop. Beyond this adjustment, the process is identical above. The
# shared model attributes between the IPD and aggregated data are the correlation
# matrix (Sigma), the disease specific intercepts, and the comorbidity count rate 
# (lambda).

for(ipd in 1:ipd_pop_size){ 
Z_ipd[ipd, 1:n_dis] ~ dmnorm.vcov(mu_ipd[ipd, 1:n_dis], Sigma[,]) # Variance formulation of MVN
for(d_ipd in 1:n_dis){ 
data_ipd[ipd,d_ipd] ~ dbern(pres_ipd[ipd,d_ipd]) # Individual level data  brought in here
pres_ipd[ipd,d_ipd] <- step(Z_ipd[ipd,d_ipd]) # If latent state is 0 or higher, disease is present
mu_ipd[ipd,d_ipd] <- intercept[d_ipd] # Disease specific, intercept only linear predictor

}

# Comorbidity counts are both equal to row sums and Poisson distributed.
zeroes_ipd[ipd] ~ dnorm(sum(pres_ipd[ipd,]) - comorb_count_ipd[ipd], 100)
comorb_count_ipd[ipd] ~ dpois(lambda)
} 


# Priors ------------------------------------------------------------------

# Disease specific intercept - vague prior
for (D in 1:n_dis){
intercept[D] ~ dnorm(0, 10^-6) 
}

# A vague prior for the rate of comorbidity counts 
lambda ~ dgamma(1.0E-3, 1.0E-3)

# We have used to the variance formulation of the MVN (rather than the precision
# formulation) as this allows us to normalise the variances (i.e. parameterise 
# with a correlation matrix rather than a variance-covariance matrix). Without
# fixing the variance to 1, the model is unidentifiable with binary data as it
# can achieve the same fit to the data by adjusting either the variance or the 
# mean. Unfortunately, this means we can't use the conjugate prior, the inverse 
# Wishart as it is not implemented in JAGS. Instead, we fix the diagonal values
# of the correlation matrix to 1 and draw the off diagonal values from a scaled
# Beta distribution. The Beta distribution has two parameters, alpha and beta,
# when alpha = beta, the distribution is symmetric. under these conditions, the
# closer alpha is to 1, the flatter (and thus vaguaer) the prior. This allows us
# to use alpha to induce shrinkage if wanted. The current value of 10 is quite flat.

alpha <- 10 # Shrinkage factor
smin <- -1 # Lower rescaling bound
smax <- 1 # Upper rescaling bound

# Diagonal elements
for(i in 1:n_dis) {Sigma[i,i]<-1 } # Fix variance to 1
# Above-diagonal elements
for(i in 1:(n_dis-1)) 
{
  for(j in (i+1):n_dis)
  {
  Sigmad[i,j]~ dbeta(alpha,alpha)
  Sigma[i,j]<-smin+Sigmad[i,j]*(smax-smin)
  }
}

# Below-diagonal elements - will be equal to corresponding upper-diagonal elements.
for(i in 2:n_dis) 
{
  for(j in 1:(i-1))
  {
  Sigma[i,j]<-Sigma[j,i]
  }
}

#monitor# Sigma, intercept, lambda
#data# pop_size, n_dis, diseases, comorb_count,zeroes, markers,ipd_pop_size,data_ipd, zeroes_ipd, comorb_count_ipd
}"

results<-run.jags(model, n.chains=4, adapt=1000, burnin=2000, sample=20000,
                  modules = 'glm', method = "parallel")

summary(results)
plot(results)
