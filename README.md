# covid_yll
Estimation of years of life lost according to age, sex and comorbidity.
Code for pre-peer reviewed article posted at https://wellcomeopenresearch.org/articles/5-75.

Please see [this link for the recent addendum](Scripts/Addendum.md)

Please only include comments related to the code - eg errors, difficulty running the code etc on the issues here. **Please discuss any issues related to the modelling assumptions, interpretation etc to the comments section of the paper which can be found at at https://wellcomeopenresearch.org/articles/5-75**.


The following files are included in this repository:-

- [Data/sail_outputs](Data/sail_outputs) - model outputs from the SAIL repository alongside model diagnostics
- [Data/Case_series_and_comorbidity.xlsx](Data/Case_series_and_comorbidity.xlsx) - original exploration of data sources indicating comorbidity distribution among people with COVID-19 who died
- [Data/](Data/) - SimulatedProfilesConverged.Rds and SimulatedProfilesUnconverged.Rds - simulated patients based on earlier unconverged model and new converged model
- [Data/age_distribution_italian_deaths.txt](Data/age_distribution_italian_deaths.txt) - age distribution of Italian peopel dying with COVID-19 from ISS report
- [Data/](Data/) - files starting with "age_selection" show age distributions under different assumptions about relation of age and comorbidity
- [Data/highest_cor_como.Rds](Data/highest_cor_como.Rds) - comorbidity dataframe assuming highly correlated comorbidities
- [Data/independent_como.Rds](Data/independent_como.Rds) - comorbidity dataframe assuming independent comorbidities
- [Data/italy_update](Data/italy_update.txt) - aggregate comorbidity data from Italy
- [Data/](Data/) - starting with "jags_model_age" - MCMC from model fitting age to comorbidity
- [Data/](Data/) - starting with "men_" or "women_" - life expectancy estimates for men and women under different modelling assumptions
- [Data/results_converged_cutdown.rds](Data/results_converged_cutdown.rds) - MCMC output from comrobidity model
- [Data/sensitivity_lifetables.txt](Data/sensitivity_lifetables.txt) - US and UK lifetables
- [Data/who_mortality_table.txt](Data/who_mortality_table.txt) - WHO GBD lifetable

- [Output_rmd/](Output_rmd/) - Contains R objects created when run [Scripts/simulation_combine.Rmd](Scripts/simulation_combine.Rmd) for different parameters
- [Outputs/reports/](Outputs/reports/) - Results in html format created by running [Scripts/simulation_combine.Rmd](Scripts/simulation_combine.Rmd) for different parameters
- [Outputs/](Outputs) - Other files are YLL and other results under different assumptions

- [Scripts/Estimating_age](Scripts/Estimating_age) - Fit model to age data to obtain a simulated age distribution for calculations under assumptions about association between age and comorbidity
- [Scripts/Extracting_comorbidity_Distributions/](Scripts/Extracting_comorbidity_Distributions/) - Bayesian model for estimating comorbidity distributions as well as sensitivity analysis

- [Supporting/COVID_icd10_codes_for_sharing.csv](Supporting/COVID_icd10_codes_for_sharing.csv) - ICD-10 codes used to define comorbidities in SAIL dataset
- [Supporting/COVID_readcodes_for_sharing.csv](Supporting/COVID_readcodes_for_sharing.csv) - Read codes used to define comorbidities in SAIL dataset
- [Supporting/age_distribution_italian_deaths.xlsx](Supporting/age_distribution_italian_deaths.xlsx) - age distribution of deaths from ISS report

- [run_report.R](run_report.R) - meta-script to run analyses
