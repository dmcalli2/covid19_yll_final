Combine component models to estimate Years of Life Lost
================
David McAllister
17 December, 2020

# Years of Life Lost for women assuming age and comorbidity are independent. The correlation between comorbidities was modelled and the model did/was converged .

This document combines the comorbidity, age and survival model outputs
to estimate years of life lost. It can be rendered for men and women,
and assuming that age and comorbidity are associated or independent.
This particular report is for women and assumes that age and comorbidity
are independent.

It also allows for sensitivity analyses on the correlation between
comorbidities. This analysis is

.

``` r
library(tidyverse)
```

# Functions

``` r
SurvGomp <- function(t, start, ...) {
  # This is a modified version of the summary function (survival part)
  # from flexsurv::summary.flexsurvreg, only difference is triple colon as 
  # workhorse function undereath (written in c++) is hidden from users
  ret <- (1 - flexsurv:::pgompertz_work(t, ...))/(1-flexsurv:::pgompertz_work(start, ...))
  ret[t<start] <- 1
  ret
}
 # example of this function
a <- SurvGomp(t = 41:130, start = 50, shape = 0.10345, rate = exp(-11.0345),
              lower_tail = TRUE, give_log = FALSE)
```

# Inputs to synthesis

Read in simulated patients from comorbidity model, age-sex distribution
from the Istituto Superiore di Sanità (ISS report), the WHO global
burden of disease life tables and the age distributions from the age
models.

``` r
if(params$correlation_como == "converged")   italy <- readRDS("Data/SimulatedProfilesConverged.Rds")
if(params$correlation_como == "unconverged") italy <- readRDS("Data/SimulatedProfilesUnconverged.Rds")
if(params$correlation_como == "maximal") italy <- readRDS("Data/highest_cor_como.Rds")
if(params$correlation_como == "independent") italy <- readRDS("Data/independent_como.Rds")

names(italy) <- str_to_lower(names(italy))
italy <- as.matrix(italy)

age_sex <- read_tsv("Data/age_distribution_italian_deaths.txt")

if(params$sex == "men") {
  my_vcov <- read_csv("Data/sail_outputs/male_vcov.csv")[,-1]
  my_coef <- read_csv("Data/sail_outputs/male_coef.csv")
  if(params$correlation == "associated") {
    if(params$sd_como_age == "low") {age_smpls <- readRDS("Data/age_selection_male_0_5.Rds")}
    if(params$sd_como_age == "high"){age_smpls <- readRDS("Data/age_selection_male_1.Rds")}
    }
  if(params$correlation == "independent") age_smpls <- readRDS("Data/age_selection_no_comorbid_male.Rds")
  who <- read_csv("Data/sail_outputs/who_compare_men.csv")

}

if(params$sex == "women") {
  my_vcov <- read_csv("Data/sail_outputs/female_vcov.csv")[,-1]
  my_coef <- read_csv("Data/sail_outputs/female_coef.csv")
  if(params$correlation == "associated") {
    if(params$sd_como_age == "low") {age_smpls <- readRDS("Data/age_selection_female_0_5.Rds")}
    if(params$sd_como_age == "high"){age_smpls <- readRDS("Data/age_selection_female_1.Rds")}
  }
  if(params$correlation == "independent") age_smpls <- readRDS("Data/age_selection_no_comorbid_female.Rds")
  who <- read_csv("Data/sail_outputs/who_compare_women.csv")

}

my_vcov <- as.matrix(my_vcov)
rownames(my_vcov) <- colnames(my_vcov)
my_coef_vect <- my_coef[,2, drop = TRUE]
names(my_coef_vect) <- my_coef[,1, drop = TRUE]
my_coef <- my_coef_vect
rm(my_coef_vect)


who_reference <- read_tsv("Data/who_mortality_table.txt")
```

# WHO life tables based YLL

The following uses the proportions in each age who died and the WHO GBD
2010 life tables to estimate the YLL.

``` r
age_sex %>% select(-age_point) %>% 
  spread(sex, n) %>% 
  knitr::kable()
```

| age\_bands |  men | women |
|:-----------|-----:|------:|
| 30-39      |   14 |     3 |
| 40-49      |   49 |    18 |
| 50-59      |  190 |    53 |
| 60-69      |  607 |   154 |
| 70-79      | 1847 |   556 |
| 80-89      | 1808 |   894 |
| 90+        |  274 |   334 |

``` r
age_sex_whole <- age_sex %>% 
    filter(sex == params$sex) %>% 
    mutate(prp = n/sum(n), 
           age_mean = prp*age_point)

age_sex_restr <- age_sex %>% 
    filter(sex == params$sex, age_point >=50) %>% 
    mutate(prp = n/sum(n), 
           age_mean = prp*age_point)

## Mean age
mean_age_whole <- age_sex_whole %>% 
  pull(age_mean) %>% 
  sum() %>% 
  round(1)

mean_age_restricted <- age_sex_restr %>% 
  pull(age_mean) %>% 
  sum() %>% 
  round(1)

## estimate yll using WHO life expectancy and age distribution alone ----
who_reference <- who_reference %>% 
  group_by(map2italy) %>% 
  summarise(expected_remain = mean(expected_remain)) %>% 
  ungroup() 
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
yll_mean_who_whole <- who_reference %>% 
  inner_join(age_sex_whole %>% select(map2italy = age_point, prp)) %>% 
  mutate(res = expected_remain * prp) %>% 
  pull(res) %>% 
  sum() %>% 
  round(1)
```

    ## Joining, by = "map2italy"

``` r
yll_mean_who_restricted <- who_reference %>% 
  inner_join(age_sex_restr %>% select(map2italy = age_point, prp)) %>% 
  mutate(res = expected_remain * prp) %>% 
  pull(res) %>% 
  sum() %>% 
  round(1)
```

    ## Joining, by = "map2italy"

The mean age for the ISS deaths was 81.1 overall and 81.5 when people
who died aged &lt;50 were excluded from the calculation.

The estimated years of life lost obtained by applying the age
distributions from the ISS report to these tables was 12.2 per person
using the whole cohort and 11.8 after excluding those aged under 50.

# Combining samples from survival models, age models and comorbidity models

The following combines 10,000 samples from survival models, age models
and comorbidity models to estimate the age, sex and comorbidity specific
life expectancy and subsequently years of life lost.

## Age distributions

The following shows the modelled age-distribution of the ISS data (based
on the published counts of patients in each age band). The grey bars are
for the simulated data and the coloured bars are the original data
(assuming a uniform distribution within each age-band).

``` r
## Sum conditions to Calculate comorbidity count
italy_count <- rowSums(italy)
age_tibble <- tibble(comorbidity_count = italy_count, age = age_smpls)

## spread ages along range of possible values to allow comparison
iss <- age_sex %>% 
  filter(sex == params$sex) %>% 
  select(age_bands, age_point, n) %>% 
  group_by(age_bands) %>% 
  nest()
iss$age <- map(iss$data, ~ rep(.x$age_point, .x$n))
iss$age <- map(iss$age, ~ .x + seq(-5, 4, length.out = length(.x)))

iss <- iss %>% 
  select(age) %>% 
  unnest(cols = c(age)) 
```

    ## Adding missing grouping variables: `age_bands`

``` r
plot_ages1 <- ggplot(age_tibble, aes(x = age)) +
  geom_histogram() +
  geom_histogram(data = iss, mapping = aes(fill = age_bands), alpha = 0.9)
plot_ages1
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Report_women_independent_converged_high_files/figure-gfm/agedistr-1.png)<!-- -->

Here the above histogram is broken down by the comorbidity count.

``` r
plot_ages2 <- ggplot(age_tibble, aes(x = age)) +
  geom_histogram(mapping = aes(fill = factor(comorbidity_count, ordered = TRUE))) 
plot_ages2
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Report_women_independent_converged_high_files/figure-gfm/repeategraphclearer-1.png)<!-- -->

## Comorbidity counts

``` r
como_tbl <- tibble(`Comorbidity count` = italy_count) %>% 
  group_by(`Comorbidity count`) %>% 
  summarise(N = length(`Comorbidity count`)) %>% 
  mutate(Percentage = round(100* N /sum(N),2))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
knitr::kable(como_tbl)
```

| Comorbidity count |    N | Percentage |
|------------------:|-----:|-----------:|
|                 0 |  110 |       1.10 |
|                 1 | 1387 |      13.87 |
|                 2 | 2267 |      22.67 |
|                 3 | 2222 |      22.22 |
|                 4 | 1747 |      17.47 |
|                 5 | 1165 |      11.65 |
|                 6 |  694 |       6.94 |
|                 7 |  294 |       2.94 |
|                 8 |   91 |       0.91 |
|                 9 |   22 |       0.22 |
|                10 |    1 |       0.01 |

Here the above histogram is shown as small multiples having stratified
by comorbidity count.

``` r
plot_ages3 <- ggplot(age_tibble, aes(x = age)) +
  geom_histogram() +
  facet_wrap( ~ comorbidity_count) +
  scale_x_continuous("Age (years)") +
  scale_y_continuous(name = "Number of simulated patients (of 20,000)")
plot_ages3
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Report_women_independent_converged_high_files/figure-gfm/agestratcomo-1.png)<!-- -->

The following is a check that the names of the comorbidities are the
same in the samples from the survival models and comorbidity estimates.
The comorbidity samples should have no differences (character 0) and the
coefficients and variance-covariance matrix from the survival model
shoudl be identical apart from the shape and rate parameters.

``` r
name_cols <- c("atr_fib", "cancer", "copd", "demen", "diab", 
               "h_fail", "hypert", "ihd", "renal", "liver", "stroke")
name_cols_age <- c("atr_fib:age", 
                   "cancer:age", "copd:age", "demen:age", "diab:age", "h_fail:age", 
                   "hypert:age", "ihd:age", "renal:age", "liver:age", "stroke:age")

setdiff(colnames(italy), name_cols)
```

    ## character(0)

``` r
setdiff(names(my_coef), c(name_cols, name_cols_age))
```

    ## [1] "shape" "rate"

``` r
setdiff(colnames(my_vcov), c(name_cols, name_cols_age))
```

    ## [1] "shape" "rate"

## Design matrix

The following code multiplies the comorbidity samples, which is a
prediction matrix by age to get age:comorbidity interactions

``` r
italy_ages <- age_smpls * italy
head(italy_ages)
```

    ##          ihd atr_fib  h_fail  stroke  hypert    diab demen    copd  cancer   liver   renal
    ## [1,]  0.0000 72.7021  0.0000  0.0000 72.7021  0.0000     0  0.0000  0.0000 72.7021  0.0000
    ## [2,]  0.0000  0.0000 86.4479 86.4479 86.4479  0.0000     0  0.0000  0.0000  0.0000  0.0000
    ## [3,] 68.0758  0.0000 68.0758 68.0758  0.0000 68.0758     0 68.0758 68.0758  0.0000 68.0758
    ## [4,] 86.0029 86.0029  0.0000  0.0000  0.0000 86.0029     0 86.0029  0.0000  0.0000 86.0029
    ## [5,]  0.0000  0.0000  0.0000  0.0000 75.2458  0.0000     0  0.0000  0.0000  0.0000  0.0000
    ## [6,]  0.0000  0.0000  0.0000  0.0000 73.7530  0.0000     0  0.0000  0.0000 73.7530  0.0000

## Sampling from survival model

THe following code samples from the survival model to obtain 10,000
samples of the coefficients.

Check data are correctly aligned then obtain 10,000 samples.

``` r
all(names(my_coef) == colnames(my_vcov))
```

    ## [1] TRUE

``` r
surv_smpls <- mvtnorm::rmvnorm(10000, my_coef, my_vcov)

a <- (colMeans(surv_smpls) %>% round(3))
b <- apply(surv_smpls, 2, sd) %>% round(3)
knitr::kable(tibble("_" = names(a), "est" = a, se = b))
```

| \_           |     est |    se |
|:-------------|--------:|------:|
| shape        |   0.109 | 0.001 |
| rate         | -12.263 | 0.049 |
| atr\_fib     |  -0.108 | 0.108 |
| cancer       |   2.345 | 0.073 |
| copd         |   2.733 | 0.078 |
| demen        |   2.719 | 0.123 |
| diab         |   0.781 | 0.074 |
| h\_fail      |   1.296 | 0.115 |
| hypert       |  -0.412 | 0.064 |
| ihd          |   0.861 | 0.133 |
| renal        |   0.364 | 0.079 |
| liver        |   2.606 | 0.194 |
| stroke       |   0.704 | 0.097 |
| atr\_fib:age |   0.005 | 0.001 |
| cancer:age   |  -0.026 | 0.001 |
| copd:age     |  -0.028 | 0.001 |
| demen:age    |  -0.022 | 0.001 |
| diab:age     |  -0.007 | 0.001 |
| h\_fail:age  |  -0.011 | 0.001 |
| hypert:age   |   0.006 | 0.001 |
| ihd:age      |  -0.008 | 0.002 |
| renal:age    |  -0.003 | 0.001 |
| liver:age    |  -0.028 | 0.003 |
| stroke:age   |  -0.005 | 0.001 |

## Estimating the linear predictor for each of the 10,000 iterations

Note that these multiplications are not matrix multiplications. They are
element-wise multiplications of identically positioned elements in each
matrix. Each part of the linear predictor is calculated separately,
these are then summed to obtain the linear predictor based on the age
and comorbidity distributions of the simulated patients.

``` r
## shape parameter
shpe       <- surv_smpls[ , "shape"]

## Intercept (rate parameter)
cept_lp    <- surv_smpls[ , "rate"]
```

The following shows the linear predictor for covariates without
interactions by age. These are mostly positive, reflecting the fact that
almost all of the comorbidites are associated with increased hazard
rates for death.

``` r
## This multiplies the covariate level (binary) with the simulated coefficeints
non_age_lp <- italy[,name_cols] * (surv_smpls[, name_cols])
## This adds these together to get the non-age dependent part of the linear predictor
non_age_lp <- rowSums(non_age_lp)
hist(non_age_lp, breaks = 100)
```

![](Report_women_independent_converged_high_files/figure-gfm/nonagelp-1.png)<!-- -->

The following shows the linear predictor for covariates with
interactions by age. All are negative, reflecting the fact that there is
attenuation in the hazard ratios for the comorbid diseases with
increasing age, as is common in chronic disease epidemiology.

``` r
## Part of lp dependent on age
# ensure_matrices_match
all(paste0(name_cols, ":age") == name_cols_age )
```

    ## [1] TRUE

``` r
## Loop through the different age groups to get the part of the LP with age interactions
ages_lp <-  rowSums(italy_ages[, name_cols]  * surv_smpls[, name_cols_age])
hist(ages_lp, breaks = 100)
```

![](Report_women_independent_converged_high_files/figure-gfm/agelp-1.png)<!-- -->

The following sums the components of the linear predictor.

``` r
## total lp
tot_lp <- cept_lp + non_age_lp + ages_lp
hist(tot_lp, breaks = 50, xlim = c(-12, -5))
```

![](Report_women_independent_converged_high_files/figure-gfm/sumlp-1.png)<!-- -->

The following transforms the linear predictor onto the rate scale.

``` r
tot_rat <- exp(tot_lp)
hist(10000*tot_rat, breaks = 500, xlim = c(0,5), xlab = "Rate per 10,000 person-years")
```

![](Report_women_independent_converged_high_files/figure-gfm/rate_scale-1.png)<!-- -->

## Estimate the survival at each time-point (3-monthly)

The following uses a function from the flexsurv package called by
summary.flexsurv to obtain survival probabilities at each follow-up
time. Where the time is before the patients age at baseline (in this
case the age at death in the ISS data) the survival probability is by
definition one.

In this function all notional patients have the same shape parameter,
but have patient-specific rate parameters. This function is not
vectorised with respect to the rate and shape parameters, but it is
nonetheless very fast.

For one simulated patient in the “low” SD age group the rate and age
variables are such that the function returns a missing value, this is
the only sample that does so. Presumably there is a numerical reason in
the internals of the function. Give a warning and drop from the summary
function at the end. Return an error if there is more than one missing.
Only occured in most recent set of simulated patients which is a fresh
simulation.

``` r
times_vect <- seq(50, 150, 0.25)
tot_surv <- matrix(NA, nrow = 10000, ncol = length(times_vect))
for(i in 1:10000){
  a <- SurvGomp(t = times_vect, start = age_smpls[i], shape = shpe[i], rate = tot_rat[i],
                lower_tail = TRUE, give_log = FALSE)
  tot_surv[i,] <- t(a)
}

tot_surv_na <- matrix(is.na(tot_surv), nrow = 10000, ncol = length(times_vect))
tot_surv_na <- rowSums(tot_surv_na)
tot_surv_na <- any(tot_surv_na >=1)
tot_surv_na <- sum(tot_surv_na)
if(tot_surv_na >=1) warning(paste0("The SurvGomp function returned NaN for ", tot_surv_na, " simulated patients"))
if(tot_surv_na >=5) stop()
```

### The following plots show the survival distribution for a random sample of 3 notional patients.

``` r
par(mfrow=c(1,3))
for(k in 1:3){
  i <- sample(1:10000, 1)
  plot(times_vect, tot_surv[i,], xlab = "age", ylab = "survival", main = paste0("Age ", age_smpls[i] %>%  round(), " years"))
}
```

![](Report_women_independent_converged_high_files/figure-gfm/example_plots-1.png)<!-- -->

``` r
par(mfrow=c(1,1))
```

## Survival distribution for all patients stratified by age and comorbidity count

The following figures show this for all patients according to age, and
comorbidity.

``` r
surv_all <- as.data.frame(tot_surv)
names(surv_all) <- times_vect
surv_all <- surv_all[ , times_vect %in% seq(50, 120, 2.5)]
surv_all <- surv_all %>% 
  mutate(pt = 1:10000)
surv_all <-  surv_all %>% 
  as_tibble() %>% 
  mutate(comorbidity_count = italy_count,
         age = age_smpls) %>% 
  gather("times", "survival", -pt, -comorbidity_count, -age) %>% 
  mutate(times = as.integer(times),
         comorbidity_count = factor(comorbidity_count, ordered = TRUE)) %>% 
  arrange(pt, times)
surv_all <- surv_all %>% 
  filter(age <= times, times <= 120) %>% 
  mutate(age_bands = Hmisc::cut2(age, cuts = seq(50, 80, 10))) 

plot_surv <- ggplot(surv_all, aes(x = times, y = survival, colour = comorbidity_count, group = pt)) +
  geom_line() +
  facet_wrap(~ age_bands)
plot_surv
```

![](Report_women_independent_converged_high_files/figure-gfm/survival_all-1.png)<!-- -->

## Calculate life expectancy (mean survival) for each patient

The following sums over the distribution of survival on time to obtain
the mean survival, which is the life expectancy.

``` r
le <- apply(tot_surv, 1, function(x) {
    x <- 1-x
    # calculate proportion who died in the interval
    x <- x - lag(x, default = 0)
    # Calculate the expectation of death
    sum(x * times_vect)
})

le_check <- apply(tot_surv, 1, function(x) {
    x <- 1-x
    # calculate proportion who died in the interval
    x <- x - lag(x, default = 0)
    # Calculate the expectation of death
    sum(x)
})
all(le_check ==1)
```

    ## [1] TRUE

``` r
saveRDS(le,  paste0("Data/", params$sex, "_", params$correlation, "_", params$correlation_como, "_le_pre_offset.Rds"))
```

## Offset life expectancy to allow for the fact that the LE in Wales is lower than the WHO levels

The following corrects for the fact that the intercept-only life
expectancy from the models in Wales are slightly lower than those in the
WHO GBD life-expectancy tables. The Wales life expectancies are
multiplied by between 1.07 and 1.02 depending on the age.

``` r
who <- who %>%
  mutate(who_ratio = who_remain/expect_remain_years,
         who_le_ratio = (who_remain+age)/(expect_remain_years + age))
mean(who$who_le_ratio)
```

    ## [1] 1.026913

``` r
if(params$sex == "men"){
  who_offset <- case_when(
  age_smpls <= 55 ~ 1.07,
  age_smpls <= 65 ~ 1.06,
  age_smpls <= 75 ~ 1.04,
  age_smpls <= 85 ~ 1.02,
  TRUE ~ 1
)
}
if(params$sex == "women"){
  who_offset <- case_when(
  age_smpls <= 55 ~ 1.04,
  age_smpls <= 65 ~ 1.03,
  age_smpls <= 75 ~ 1.02,
  age_smpls <= 85 ~ 1.01,
  TRUE ~ 1
)
}
le <- le * who_offset
```

## Years of life lost

The following calculates the years of life lost from the age, sex and
comorbidity defined survival and the age and comorbidities at death from
the ISS data.

``` r
yll <- le - age_smpls
par(mfrow = c(1,2))
plot(age_smpls, yll, xlab = "Age (years)", ylab = "Years of life lost")
hist(yll, breaks = 50, xlab = "Years of life lost", ylab = "Number of simulated patients")
```

![](Report_women_independent_converged_high_files/figure-gfm/yll_wales-1.png)<!-- -->

## Compare YLL using WHO tables and comorbidity-based assesment

Among those over 50 for which both methods were used, the WHO tables
based estimates of YLL, and our estimates which additionally accomodate
comorbidity, yielded similar per capita YLL; the values were 11.8 and
9.16 respectively.

## YLL by age and comorbidity

The following shows the distribution of YLL by age, sex and comorbidity
by age categories.

### Table YLL per capita by comorbidity count, and age

``` r
yll_df <- tibble(yll = yll) %>% 
  mutate(comorbidity_count = factor(italy_count, ordered = TRUE),
         age = round(age_smpls),
         age_cat = Hmisc::cut2(age, cuts = seq(50, 80, 10)))

yll_smry <- yll_df %>% 
  group_by(comorbidity_count, age_cat) %>% 
  summarise(yll = mean(yll) %>% round(2)) %>% 
  rename(`Comorbidity count` = comorbidity_count) %>% 
  spread(age_cat, yll)
```

    ## `summarise()` regrouping output by 'comorbidity_count' (override with `.groups` argument)

``` r
knitr::kable(yll_smry)
```

| Comorbidity count | \[ 50, 60) | \[ 60, 70) | \[ 70, 80) | \[ 80,110\] |
|:------------------|-----------:|-----------:|-----------:|------------:|
| 0                 |         \- |      25.23 |      16.60 |        8.38 |
| 1                 |      33.94 |      24.83 |      16.39 |        7.68 |
| 2                 |      29.73 |      21.33 |      14.00 |        6.67 |
| 3                 |      25.38 |      17.69 |      11.87 |        5.76 |
| 4                 |      21.11 |      14.68 |       9.88 |        4.92 |
| 5                 |      16.55 |      12.74 |       8.21 |        4.17 |
| 6                 |      11.89 |       9.43 |       6.67 |        3.41 |
| 7                 |      11.91 |       7.10 |       5.16 |        2.81 |
| 8                 |       8.80 |       5.47 |       3.90 |        2.40 |
| 9                 |         \- |         \- |       2.95 |        1.56 |
| 10                |         \- |         \- |         \- |        0.99 |

## Plot distribution of years of life lost by age and comrobidity count

Plot distributions of YLL overall, then stratified by age, with colour
indicating the comorbidity count.

``` r
plot1_yll <- ggplot(yll_df, aes(x = yll,  fill = comorbidity_count)) + 
  geom_histogram(bins = 50) +
  scale_y_continuous("Number of simulated patients") +
  scale_x_continuous("Years of life lost")
plot1_yll
```

![](Report_women_independent_converged_high_files/figure-gfm/plot_yll_overall-1.png)<!-- -->

``` r
plot2_yll <- plot1_yll +
  facet_wrap(~ age_cat) 
plot2_yll
```

![](Report_women_independent_converged_high_files/figure-gfm/plot_yll_age_como-1.png)<!-- -->

## Crude approach to obtain 95% credible interval

The following approach uses empirical bootstrapping to calculate an
uncertainty interval around the mean YLL. This

``` r
res <- map_dbl(1:1000, ~  mean(sample(yll, 404, replace = TRUE)))
hist(res, xlab = "Mean Years of Life Lost")
```

![](Report_women_independent_converged_high_files/figure-gfm/credible%20intervals-1.png)<!-- -->

``` r
est <- mean(res, na.rm = TRUE) %>% round(1)
lci_uci <- quantile(res, probs = c(0.025, 0.975), na.rm = TRUE) %>% round(1) %>% paste(collapse = "-")
est
```

    ## [1] 9.2

``` r
lci_uci
```

    ## [1] "8.6-9.7"

The mean and (95% credible interval) for the years of life lost using
this approach was 9.2 (95% CI 8.6-9.7).

``` r
mylist <- list(est, lci_uci, mean_age_restricted, mean_age_whole,
               plot_ages1, plot_ages2, plot_ages3, plot_surv, plot1_yll, plot2_yll,
               who, yll_df, yll_mean_who_restricted, yll_mean_who_whole, yll_smry)
saveRDS(mylist, paste0("Output_rmd/", 
                       params$sex, "_", 
                       params$correlation, "_", 
                       params$correlation_como, "_", 
                       params$sd_como_age,
                       ".Rds"))
```
