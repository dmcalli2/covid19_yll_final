## Analysis to support addendum

## Packages ----
library(tidyverse)

## Functions ----
SurvGomp <- function(t, start, ...) {
  # This is a modified version of the summary function (survival part)
  # from flexsurv::summary.flexsurvreg, only difference is triple colon as 
  # workhorse function undereath (written in c++) is hidden from users
  ret <- (1 - flexsurv:::pgompertz_work(t, ...))/(1-flexsurv:::pgompertz_work(start, ...))
  ret[t<start] <- 1
  ret
}

MakeONSCompare <- function(myshape, myrate, interval_yrs = 1){
  res <- map(seq(50, 100, interval_yrs), ~ SurvGomp(t = seq(.x, 130),
                                                   start = .x,
                                                   shape = myshape,
                                                   rate = myrate,
                                                   lower_tail = TRUE,
                                                   give_log = FALSE))
  names(res) <- paste0("aged", seq(50, 100, interval_yrs))
  res <- map2(res, seq(50, 100, interval_yrs), ~
                cbind(srval = .x, times = seq(.y, 130), age= .y))
  res <- map(res, as_tibble)
  res <- bind_rows(res, .id = "age_start")
  
  res <- res %>% 
    group_by(age_start) %>% 
    mutate(est = 1 - srval,
           prb = est - lag(est, default = 0),
           prb_t = times * prb) %>% 
    ungroup()
  res %>% 
    group_by(age_start, age) %>% 
    summarise(mean_dist = sum(prb_t)) %>% 
    ungroup() %>% 
    mutate(expect_remain_years = mean_dist - age) %>% 
    select(age, age_start, expect_remain_years, expect_death = mean_dist)
  
}


## Obtain and plot life tables ----
## Wales coefficeints data, use to create life table
wales_coef_m <- read_csv("Data/sail_outputs/cept_coef_gomp_m.csv")
wales_coef_f <- read_csv("Data/sail_outputs/cept_coef_gomp_f.csv")

shape_m <- wales_coef_m %>% filter(param == "shape") %>% pull(value)
rate_m  <- wales_coef_m %>% filter(param == "rate") %>% pull(value)
wales_m <- MakeONSCompare(shape_m, exp(rate_m))
shape_f <- wales_coef_f %>% filter(param == "shape") %>% pull(value)
rate_f  <- wales_coef_f %>% filter(param == "rate") %>% pull(value)
wales_f <- MakeONSCompare(shape_f, exp(rate_f))
rm(rate_m, rate_f, shape_m, shape_f, wales_coef_m, wales_coef_f)
wales <- bind_cols(wales_m %>% select(age, male = expect_remain_years),
                   wales_f %>% select(female = expect_remain_years)) %>% 
  mutate(lifetable = "wales_gomp_cept")
rm(wales_f, wales_m)

## READ in ONS and USA data, can share these on github
ons_usa <- read_tsv("Data/sensitivity_lifetables.txt") %>% 
  select(lifetable, age, male, female)
## Read in italy data, cannot share this on github as per www.mortality.org user agreement
italy_lt <- read_tsv("../covid19_yll/Data/sensitivity_lifetables_italy.txt") %>% 
  select(lifetable, age, male, female)

lts <- bind_rows(ons_usa, italy_lt, wales)
rm(ons_usa, italy_lt, wales)

lts_lng <- lts %>% 
  gather("sex", "expected_remain", -lifetable, -age) %>% 
  filter(age >=50)
plot_le <- ggplot(lts_lng, aes(x = age, y = expected_remain, colour = lifetable)) +
  geom_point() +
  geom_line() +
  facet_grid(~sex) +
  scale_y_continuous("Expected years remaining", limits = c(0,NA))
plot_le

## Read in GBD data
gbd <- read_tsv("Data/who_mortality_table.txt")
gbd_plot <- gbd %>% 
  mutate(age = as.integer(str_sub(age, 1, 2)) + 2.5,
         age = if_else(age == 85, 95, age)) %>% 
  select(age, expected_remain) %>% 
  mutate(lifetable = "gbd2010")
gbd_plot <- bind_rows(gbd_plot %>% mutate(sex = "male"),
                      gbd_plot %>% mutate(sex = "female"))

lts_lng_gbd <- bind_rows(lts_lng, gbd_plot)

lts_lng_gbd <- lts_lng_gbd %>% 
  mutate(lifetable2 = factor(lifetable, 
                             levels = c("gbd2010", "italy_2017", "uk2016_2018", "us_2017","wales_gomp_cept"),
                             labels = c("GBD 2010", "Italy 2017", "UK 2016-2018", "US 2017", "Wales (see paper)"), ordered = TRUE),
         sex2 = if_else(sex == "female", "Female", "Male")
  )

## Plot life expectancy based on different life tables
plot_le <- ggplot(lts_lng_gbd %>% filter(age >=50), aes(x = age, y = expected_remain, colour = lifetable2)) +
  geom_point() +
  geom_line() +
  facet_grid(~sex2) +
  scale_y_continuous("Expected years remaining", limits = c(0,NA)) +
  scale_x_continuous("Age in years") +
  scale_colour_brewer("Life table", palette = "Paired")
plot_le
saveRDS(plot_le, "Outputs/compare_lifetables.Rds")

## Estimate age-sex specific ratio between LE in Wales data and LE in life tables ----
## read in modelled Italy data
italy <- readRDS("Data/SimulatedProfiles.Rds")
names(italy) <- str_to_lower(names(italy))
italy <- as.matrix(italy)
italy_count <- rowSums(italy)
rm(italy)

## add values for 101 to 110 for countries other than italy which is already available up to 110
lts_lng_101_pls <- lts_lng %>% 
  filter(age == 100 & lifetable != "italy_2017") %>% 
  slice(rep(1:nrow(.), each = 10))
lts_lng_101_pls$age <- 1:10 + lts_lng_101_pls$age
lts_lng <- bind_rows(lts_lng,
                     lts_lng_101_pls) %>% 
  arrange(lifetable, sex, age)

## Calculate ratios
ons_usa_italy <- lts_lng %>% 
  filter(lifetable != "wales_gomp_cept") 
wales <- lts_lng %>% 
  filter(lifetable == "wales_gomp_cept") %>% 
  rename(wales = expected_remain) %>% 
  select(-lifetable)
lts_wide <- ons_usa_italy %>% 
  inner_join(wales) %>% 
  mutate(expected_remain_ratio = expected_remain/wales,
         le_ratio = (expected_remain + age) /(wales + age))
lts_ratio <- lts_wide %>% 
  select(lifetable, age, sex, le_ratio)
quantile(lts_ratio$le_ratio)
mean(lts_ratio$le_ratio)

## Calcualte YLL based on age alone for different life tables ----
men_age_smpls   <- readRDS("Data/age_selection_male.Rds")
women_age_smpls <- readRDS("Data/age_selection_female.Rds")
age_both <- bind_rows(tibble(sex = "male", age = men_age_smpls, comorbidity_count = italy_count),
                      tibble(sex = "female", age = women_age_smpls, comorbidity_count = italy_count)) %>% 
  mutate(age = round(age))
age_both_cmpr <- lts_lng %>% 
  inner_join(age_both) 

## Calcualte YLL based on age and comorbidities using different life tables ----
le_multi_m <- readRDS("Data/men_associated_modelled_le_pre_offset.Rds")
le_multi_f <- readRDS("Data/women_associated_modelled_le_pre_offset.Rds")

le_multi <- bind_rows(tibble(sex = "male",   age = men_age_smpls,   le_pre = le_multi_m, comorbidity_count = italy_count),
                      tibble(sex = "female", age = women_age_smpls, le_pre = le_multi_f, comorbidity_count = italy_count)) %>% 
  mutate(age = round(age))
rm(le_multi_f, le_multi_m)

multi_both_cmpr <- lts_ratio %>% 
  inner_join(le_multi) %>% 
  mutate(le_post = le_pre * le_ratio,
         expected_remain = le_post -age) 

multi_both_cmpr <- bind_rows(multi_both_cmpr,
                             le_multi %>% mutate(lifetable = "wales_gomp_cept",
                                                 expected_remain = le_pre - age)) %>% 
  select(lifetable, age, sex, comorbidity_count, expected_remain)

## Summarise results for multiple life tables ----
multi_both_cmpr_smry <- multi_both_cmpr %>% 
  group_by(lifetable, sex) %>% 
  summarise(yll = mean(expected_remain)) %>% 
  ungroup() 
age_both_cmpr_smry <- age_both_cmpr %>% 
  group_by(lifetable, sex) %>% 
  summarise(yll = mean(expected_remain)) %>% 
  ungroup()
both_cmp_smry <- bind_rows(`Not LTC adjusted` = age_both_cmpr_smry, `LTC adjusted` = multi_both_cmpr_smry, .id = "LCT_adjusted") %>% 
  spread(LCT_adjusted, yll) %>% 
  arrange(sex) %>% 
  mutate(`Difference with adjustment (months)` = 12*(`Not LTC adjusted` - `LTC adjusted`),
         lifetable = factor(lifetable, 
                levels = c("gbd2010", "italy_2017", "uk2016_2018", "us_2017","wales_gomp_cept"),
                labels = c("GBD 2010", "Italy 2017", "UK 1016-2018", "US 2017", "Wales (see paper)"), ordered = TRUE),
         sex = if_else(sex == "female", "Female", "Male")) %>% 
  select(Sex = sex, `Life Table` = lifetable, everything()) 

men_gbd <- readRDS("Output_rmd/men_associated_modelled.Rds")
names(men_gbd) <- c("est", "lci_uci", "mean_age_restricted", "mean_age_whole",
  "plot_ages1", "plot_ages2", "plot_ages3", "plot_surv", "plot1_yll", "plot2_yll",
  "who", "yll_df", "yll_mean_who_restricted", "yll_mean_who_whole", "yll_smry")
men_gbd <- men_gbd[c("est", "yll_mean_who_restricted")] %>% 
  as.data.frame() 
women_gbd <- readRDS("Output_rmd/women_associated_modelled.Rds")
names(women_gbd) <- c("est", "lci_uci", "mean_age_restricted", "mean_age_whole",
                    "plot_ages1", "plot_ages2", "plot_ages3", "plot_surv", "plot1_yll", "plot2_yll",
                    "who", "yll_df", "yll_mean_who_restricted", "yll_mean_who_whole", "yll_smry")
women_gbd <- women_gbd[c("est", "yll_mean_who_restricted")] %>% 
  as.data.frame() 

gbd <- bind_rows(Male = men_gbd, Female = women_gbd, .id = "Sex") %>% 
  mutate(`Life Table` = "GBD 2010",
         `Difference with adjustment (months)` = 12*(yll_mean_who_restricted -est)) %>% 
  select(Sex, `Life Table`, `Not LTC adjusted` = yll_mean_who_restricted, `LTC adjusted` = est, everything())
both_cmp_smry2 <- bind_rows(gbd, both_cmp_smry) %>% 
  arrange(desc(Sex)) 
sex_n <- read_tsv("Data/age_distribution_italian_deaths.txt") %>% 
  mutate(Sex = if_else(sex == "men", "Male", "Female")) %>% 
  group_by(Sex) %>% 
  summarise(n = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/sum(n)) %>% 
  select(-n)
both_cmp_smry2b <- both_cmp_smry2 %>% 
  inner_join(sex_n) %>% 
  group_by(`Life Table`) %>% 
  summarise(`Not LTC adjusted` = weighted.mean(`Not LTC adjusted`, prop),
            `LTC adjusted` = weighted.mean(`LTC adjusted`, prop)) %>% 
  ungroup() %>% 
  mutate(Sex = "Total",
         `Difference with adjustment (months)` = round(12 *(`Not LTC adjusted` - `LTC adjusted`)))
both_cmp_smry3 <- bind_rows(both_cmp_smry2,
                            both_cmp_smry2b) %>%
  mutate_at(vars(`LTC adjusted`, `Not LTC adjusted`), function(x) round(x, 1)) %>% 
  mutate(`Difference with adjustment (months)` = round(`Difference with adjustment (months)`),
         `Not LTC adjusted` = if_else(`Life Table` == "GBD 2010", paste0(`Not LTC adjusted`, "*"), as.character(`Not LTC adjusted` ))) %>% 
  arrange(desc(Sex)) %>% 
  rename(`LTC number and type unadjusted`  =  `Not LTC adjusted`, `LTC number and type adjusted` = `LTC adjusted`)
saveRDS(both_cmp_smry3, "Outputs/Compare_lifetables_yll.Rds")

## Replicate Table 2 of manuscript (version 1) with UK life tables ----
yll_df <- multi_both_cmpr %>% 
  filter(lifetable == "uk2016_2018", sex == "male") %>%
  rename(yll = expected_remain) %>% 
  mutate(age_cat = Hmisc::cut2(age, cuts = seq(50, 80, 10)))
yll_smry <- yll_df %>% 
  group_by(comorbidity_count, age_cat) %>% 
  summarise(yll = mean(yll) %>% round(2)) %>% 
  rename(`Comorbidity count` = comorbidity_count) %>% 
  spread(age_cat, yll)

yll_df_f <- multi_both_cmpr %>% 
  filter(lifetable == "uk2016_2018", sex == "female") %>%
  rename(yll = expected_remain) %>% 
  mutate(age_cat = Hmisc::cut2(age, cuts = seq(50, 80, 10)))
yll_smry_f <- yll_df_f %>% 
  group_by(comorbidity_count, age_cat) %>% 
  summarise(yll = mean(yll) %>% round(2)) %>% 
  rename(`Comorbidity count` = comorbidity_count) %>% 
  spread(age_cat, yll)


age_lkp <- c(`[ 50, 60)` = "50-59",
`[ 60, 70)` = "60-69",
`[ 70, 80)` = "70-79",
`[ 80,110]` = "80 plus")
names(yll_smry)[-1] <- age_lkp[names(yll_smry)[-1]]
names(yll_smry_f)[-1] <- age_lkp[names(yll_smry_f)[-1]]
saveRDS(yll_smry, "Outputs/yll_uklt.Rds")
saveRDS(yll_smry_f, "Outputs/yll_uklt_f.Rds")

yll_df_nocomo <- age_both_cmpr %>% 
  filter(lifetable == "uk2016_2018", sex == "male") %>%
  rename(yll = expected_remain) %>% 
  mutate(age_cat = Hmisc::cut2(age, cuts = seq(50, 80, 10)))
yll_smry_nocomo <- yll_df_nocomo %>% 
  group_by(comorbidity_count, age_cat) %>% 
  summarise(yll = mean(yll) %>% round(2)) %>% 
  rename(`Comorbidity count` = comorbidity_count) %>% 
  spread(age_cat, yll)

## Summarise the comorbidity distribution by age ----
como_age <- multi_both_cmpr %>% 
  filter(lifetable == "uk2016_2018") %>%
  mutate(age_cut = Hmisc::cut2(age, cuts = seq(50,90, 10))) %>% 
  group_by(age_cut, sex) %>% 
  count(comorbidity_count) %>% 
  ungroup() %>% 
  group_by(age_cut, sex) %>% 
  mutate(prcnt = round(100*n/sum(n), 1)) %>% 
  ungroup()

como_age_smry <- como_age %>% 
  select(age_cut, sex, comorbidity_count, prcnt) %>% 
  spread(comorbidity_count, prcnt, fill = 0)
como_age_smry

sail_como <- read_csv("Data/sail_outputs/comorbidity_count_summary_sail.csv")
sail_como <- sail_como %>% 
  select(-n)
como_cmpr <- bind_rows(sail = sail_como, italy = como_age_smry, .id = "data_source") 
como_cmpr <- como_cmpr %>% 
  gather("como_count", "percentage", -age_cut, -sex, -data_source) %>% 
  mutate(como_count = as.integer(como_count))
## Add n for each age-sex combination in Italian deaths
italy_counts <- read_tsv("Data/age_distribution_italian_deaths.txt")
italy_counts <- italy_counts %>% 
  mutate(age_cut = paste("Age = " ,age_point-5, "-" ,age_point + 4, sep = ""),
         sex = if_else(sex == "women", "Women", "Men")) %>% 
  select(age_cut, sex, n) %>% 
  mutate(age_cut = if_else(age_cut == "Age = 90-99", "Age = 90 plus", age_cut))

age_cat_lkp <- c(`[ 50, 60)` = "50-59", `[ 60, 70)` = "60-69", `[ 70, 80)` = "70-79", `[ 80, 90)` = "80-89", 
                 `[ 90,110]` = "90 plus")
age_cat_lkp[] <- paste0("Age = ", age_cat_lkp)

como_cmpr <- como_cmpr %>% 
  mutate(age_cut = as.character(age_cat_lkp[age_cut]),
         sex = if_else(sex == "female", "Women", "Men"))
como_cmpr2 <- como_cmpr %>% 
  inner_join(italy_counts) %>% 
  mutate(n = paste0("N. deaths = ", n ),
         new_facet = paste(sex, age_cut, n, sep = "\n"),
         data_source = if_else(data_source == "sail", "SAIL databank, Wales", "Italian (modelled)"))

plot_como <- ggplot(como_cmpr2, aes(x = como_count, y = percentage, colour = data_source)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ new_facet, nrow = 2) +
  scale_color_discrete("") +
  scale_x_continuous("Count of long term condition") +
  scale_y_continuous("Percentage with count (%)")
plot_como
saveRDS(plot_como, "Outputs/comparison_comorbidity_sail_italy.Rds")

## Calculate marginal distribution of comorbidity counts ----
age_both_weight <- age_both %>% 
  select(-comorbidity_count) %>% 
  mutate(age_cut = Hmisc::cut2(age, cuts = seq(50, 90, 10))) %>% 
  group_by(age_cut, sex) %>% 
  count()

sail_como_marg <- sail_como %>% 
  inner_join(age_both_weight)
sail_como_marg <- map(sail_como_marg %>% 
      select(`0`:`11`),
    ~ weighted.mean(.x, sail_como_marg$n))  %>% 
  stack() 
names(sail_como_marg) <- c("prcnt", "como_count")
sail_como_marg <- sail_como_marg %>% 
  mutate(como_count = as.integer(como_count))
sail_como_marg$como_count <- sail_como_marg$como_count -1

italy_como_marg <- tibble(como_count = italy_count) %>% 
  count(como_count) %>% 
  mutate(prcnt = 100*n/sum(n))
both_marg <- bind_rows(`ISS Report` = italy_como_marg,
                       `SAIL data`  = sail_como_marg,
                       .id = "data_source") 
both_marg2 <- both_marg %>%
  select(-n) %>% 
  spread(como_count, prcnt, fill = 0L)

## Tabulate effect estimates ----
male <- read_csv("Data/sail_outputs/male_coef.csv")
female <- read_csv("Data/sail_outputs/female_coef.csv")
coefs <- bind_rows(male = male, female = female, .id = "sex")
names(coefs) <- c("sex", "parameter", "value")
coefs <- coefs %>% 
  filter(!parameter %in% c("shape", "rate"))
coefs_age <- coefs %>% 
  filter(str_detect(parameter, "age")) %>% 
  mutate(parameter = str_replace(parameter, fixed(":age"), ""))
coefs <- coefs %>% 
  filter(!str_detect(parameter, "age"))

ages <- tibble(age = seq(50, 90, 10))

coefs_age <- merge(coefs_age, ages) 

coefs_age <- coefs_age %>% 
  mutate(est_inter = age*value) %>% 
  select(sex, parameter, age, est_inter) %>% 
  inner_join(coefs) 

coefs_age <- coefs_age %>%
  mutate(est = value + est_inter,
         hr = exp(est))
  
coefs_age_wide <- coefs_age %>% 
  mutate(hr = round(hr, 2)) %>% 
  select(sex, parameter, age, hr) %>% 
  spread(age, hr)

condition_name_lkp <- c(atr_fib = "Atrial fibrillation",
                        cancer = "Cancer",
                        copd = "Chronic obstructive pulmonary disease",
                        demen = "Dementia",
                        diab = "Diabetes",
                        h_fail = "Heart failure",
                        hypert = "Hypertension",
                        ihd = "Ischaemic heart disease",
                        liver = "Liver failure",
                        renal = "Renal failure",
                        stroke = "Stroke")
coefs_age_wide <- coefs_age_wide %>% 
  mutate(parameter = condition_name_lkp[parameter],
         sex = if_else(sex == "male", "Male", "Female")) %>% 
  rename(Sex = sex, Condition = parameter)
names(coefs_age_wide)[-(1:2)] <- paste0("Age ", names(coefs_age_wide)[-(1:2)])
saveRDS(coefs_age_wide, "Outputs/Effect_estimates_by_age.Rds")

