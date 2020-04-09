library(tidyverse)

render_report = function(sex, correlation, correlation_como) {
  rmarkdown::render(
    "Scripts/simulation_combine.Rmd", params = list(
      sex = sex,
      correlation = correlation,
      correlation_como = correlation_como
    ),
    output_file = paste0("Report_", sex, "_", correlation, "_", correlation_como,".md"))
}

for(sexes in c("men", "women")){
  for(correlation in c("associated", "independent")){
    for(correlation_como in c("independent", "maximal", "modelled")){
      print(paste(c(sexes, correlation, correlation_como), sep = " "))
      render_report(sexes, correlation, correlation_como)
    }
  }
}

a <- (list.files("Scripts/", patt = "^Rep") )
b <- paste0("Outputs/reports/", a)
a <- paste0("Scripts/", a)
file.copy(from = a, to = b)

file.copy(from=a, to=b, 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)

file.remove(a)
