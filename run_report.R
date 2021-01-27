library(tidyverse)
## Make sure that your un this on a fresh restart of R or else
## you may get a conflict if there is a conflict wiht dplyr::select (because the code does not use the double-colon notation)
render_report = function(sex, correlation, correlation_como, sd_como_age) {
  rmarkdown::render(
    "Scripts/simulation_combine.Rmd", params = list(
      sex = sex,
      correlation = correlation,
      correlation_como = correlation_como,
      sd_como_age = sd_como_age
    ),
    output_file = paste0("Report_", 
                         sex, "_", 
                         correlation, "_", 
                         correlation_como, "_",
                         sd_como_age, ".md"))
}

for(sexes in c("men", "women")){
  for(correlation in c("associated", "independent")){
    for(correlation_como in c("independent", "maximal", "converged", "unconverged")){
      for(sd_como_age in c("low","high")){
        print(paste(c(sexes, correlation, correlation_como, sd_como_age), sep = " "))
        render_report(sexes, correlation, correlation_como, sd_como_age)
   }
  }
 }
}

if(!dir.exists("Outputs/reports")) dir.create("Outputs/reports")
a <- (list.files("Scripts/", patt = "^Rep") )
b <- paste0("Outputs/reports/", a)
a <- paste0("Scripts/", a)
file.copy(from=a, to=b, 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)
file.remove(a)
