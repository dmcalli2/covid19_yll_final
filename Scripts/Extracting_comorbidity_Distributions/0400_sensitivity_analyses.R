## Fully independent sampling then sample in proportion to comorbidity count
library(tidyverse)
library(matrixStats)

a <- read_tsv("data/italy_update.txt")

como <- a %>% 
  slice(1:11) 

como <- como %>% 
  mutate(prcnt = `%` /100,
         x = 710*prcnt,
         no = N-x)

cnts <- a %>% 
  slice(13:16) %>% 
  mutate( prcnt = `%`/100) %>% 
  select(cnt = Diseases, N, prcnt)

## Create most highly correlated disease structure ----
## Except for hypertension which is associated with lower mortalit, so treat hypertension as no hypertension
como <- como %>% 
  mutate(prcnt_bad = if_else(Diseases == "Hypertension", 1 - prcnt, prcnt))

## Select distribution that maximises correlation between conditions/ not having hypertension
highest <- map(como$prcnt, function(x) {
  got <- rep(1, as.integer(710*x))
  not <- rep(0, 710 - length(got))
  c(got, not)
})
highest <- bind_cols(highest)

## plot pairwise correlations
cor(highest)

# sum count distribution (this will be wrong )
rowSums(highest) %>% table()

## label to match modelled data
names(highest) <- c("ihd", "atr_fib", "h_fail", "stroke",
                  "hypert",  "diab",  "demen",
                  "copd", "cancer", "renal", "liver" )
highest <- highest[ , c("atr_fib", "cancer", "copd", "demen", "diab", 
                    "h_fail", "hypert", "ihd", "renal", "liver", "stroke")]


## resample from 710 patients to get 10,000 samples
highest_smpls <- highest[sample(1:710, 10000, replace = TRUE),]
saveRDS(highest_smpls, "Data/highest_cor_como.Rds")

## Next sample for independence
indepen_smpls <- map(highest, ~ .x[sample(1:710, 10000, replace = TRUE)])
indepen_smpls <- bind_cols(indepen_smpls) %>% 
  as.matrix()
cor(indepen_smpls)
saveRDS(indepen_smpls, "Data/independent_como.Rds")
