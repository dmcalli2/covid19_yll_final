# collect_plots


library(ggplot2)
library(tidyverse)
library(cowplot)

a <- list.files("Output_rmd/")
b <- map(a, ~ readRDS(paste0("Output_rmd/", .x)))
names(b) <- a %>% str_sub(1,-5)
b <- transpose(b)
names(b) <- c("est", "lci_uci", "mean_age_restricted", "mean_age_whole",
              "plot_ages1", "plot_ages2", "plot_ages3", "plot_surv", "plot1_yll", "plot2_yll",
              "who", "yll_df", "yll_mean_who_restricted", "yll_mean_who_whole", "yll_smry")



## Age plots ----
plot_age2 <- b$plot_ages2
plot_age2 <- plot_age2[c("men_associated_modelled", "men_independent_modelled",
"women_associated_modelled", "women_independent_modelled")]
plot_age2 <- map(plot_age2, ~ .x + scale_fill_ordinal(guide = FALSE) + scale_y_continuous("Number of simulated patients"))
names(plot_age2) <- paste(c("Men", "Men", "Women", "Women"), rep(c("Associated", "Independent"),2), sep = ", ")
plot_age2 <- map2(plot_age2, names(plot_age2), ~ .x + ggtitle(.y))
plot_overall <- plot_grid(plot_age2[[1]], plot_age2[[2]], plot_age2[[3]], plot_age2[[4]])

tiff("Outputs/Age_plots2.tiff")
plot_overall
dev.off()

## Survival plots ----
srvs <- b$plot_surv
names(srvs)
srvs <- srvs[c("men_associated_modelled", "women_associated_modelled")]

srvs <- map2(srvs, c("Men", "Women"), ~ .x +
               scale_colour_ordinal(guide = FALSE) + 
               scale_y_continuous("Cumulative survival probability") +
               scale_x_continuous("Age (years)") +
               ggtitle(.y))
plot_overall <- plot_grid(Men = srvs$men_associated_modelled, Women = srvs$women_associated_modelled, align = "h", labels = "auto")
plot_overall
win.metafile("Outputs/surv_plots.emf", width = 8)
plot_overall
dev.off()

## YLL plots overall
ylls <- b$plot1_yll
names(ylls)
ylls <- ylls[c("men_associated_modelled", "women_associated_modelled")]

ylls <- map2(ylls, c("Men", "Women"), ~ .x +
               scale_fill_ordinal(guide = FALSE) + 
               ggtitle(.y))
plot_overall <- plot_grid(Men = ylls$men_associated_modelled, Women = ylls$women_associated_modelled, align = "h", labels = "auto")

tiff("Outputs/yll_plots_overall.tiff")
plot_overall
dev.off()

## YLL stratified
ylls <- b$plot2_yll
names(ylls)
ylls <- ylls[c("men_associated_modelled", "women_associated_modelled")]

age_labs <- list("50-59", "60-69", "70-79", "80+")
names(age_labs) <- list("[50,60]", "[60,70]", "[70,80]", "[80,110]")

age_labeller <- function(variable,value){
  return(age_labs[value])
}

ylls <- map2(ylls, c("Men", "Women"), ~ .x +
               scale_fill_ordinal(guide = FALSE) + 
               ggtitle(.y) +
               facet_wrap(~ age_cat, 
                          labeller = age_labeller))
plot_overall <- plot_grid(Men = ylls$men_associated_modelled, Women = ylls$women_associated_modelled, align = "h", labels = "auto")

tiff("Outputs/yll_plots_stratified.tiff")
plot_overall
dev.off()
