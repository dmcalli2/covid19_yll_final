Addendum to COVID-19 – exploring the implications of long-term condition
type and extent of multimorbidity on years of life lost: a modelling
study
================
Peter Hanlon, Fergus Chadwick, Anoop Shah, Rachael Wood, Jon Minton,
Gerry McCartney, Colin Fischbacher, Frances S Mair, Dirk Husmeier, Jason
Matthiopoulos and David A McAllister

# Executive summary

  - Following additional analyses, using a range of national life tables
    in addition to the Global Burden of Disease (GBD) 2010 life tables
    which we used previously, the average years of life lost (YLL) due
    to COVID-19 remained above 10, even after adjusting for the number
    of long-term conditions

  - The number and type of long-term conditions has a large impact on
    YLL for individual patients, but a minimal impact on the overall
    average YLL

  - The comparatively small average overall impact is largely due to the
    fact that long-term conditions in general and multimorbidity in
    particular (the presence of two or more long-term conditions) are
    common in the older general adult population, not just among people
    who died with COVID-19

# Introduction

We received several interesting and useful suggestions, as well as some
press reports, concerning publication of
[version](https://wellcomeopenresearch.org/articles/5-75) one of our
recent manuscript. In response to these, we have undertaken further
analyses and produced additional tables and plots.

We will shortly incorporate these into another version of the manuscript
which will be available via the Wellcome Open Research journal, but in
view of the attention that the work received we have opted to post this
additional work rapidly via placing it immediately on our [project
github repository](https://github.com/dmcalli2/covid19_yll_final).

Specifically this addendum addresses:-

  - The impact of whether using different life tables would change the
    results

  - Why accounting for the number and type of long-term conditions does
    not have a large impact on the average years of life lost

  - The limited generalisability of these findings for special
    populations such as care homes.

# The impact of using different life tables

In our original analysis we used the Global Burden of Disease 2010
(GBD-2010) life tables rather than UK or Italian life tables. We did so
to allow comparison of the burden of COVID-19 against other causes of
death within an international framework. This was true both for the
simple “standard” comparison, and in the more complex analysis
additionally accounting for the number and type of long-term conditions
(as a post-modelling correction).

Some commentators interpreted our estimated YLL, however, as being the
YLL for the UK. However, our intention was for our estimates to be
comparable with other causes of death where the YLL was also benchmarked
against the GBD 2010 life tables. Although the age-distribution for
deaths from COVID-19 that we used originated from the Italian Istituto
Superiore di Sanita (ISS) report, this was understandable given that we
(i) used data from Wales to estimate the impact of number of long-term
conditions in the general population, (ii) included some (albeit a small
amount) of data from Scotland to inform our modelling and (iii) are
UK-based researchers.

We have therefore added a series of sensitivity analyses using national
lifetables for the UK and other countries, to examine whether greater
consistency in the populations and time period across our data sources,
and the closer alignment with the current UK context, makes a
substantial difference to the overall YLL estimates.

Figure A1 shows the remaining expected years of life by age and sex
using standard life tables from GBD-2010 (which we used in our original
analysis), [Italy (2017)](https://www.mortality.org/), the [UK
(2016-2018)](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables)
and the [US
(2017)](https://www.cdc.gov/nchs/data/nvsr/nvsr68/nvsr68_07-508.pdf). It
also shows the remaining expected years of life directly obtained from a
sex-specific intercept-only Gompertz model fitted to the Wales data
within the SAIL repository (the same data and type of model we used in
[our recent manuscript](https://wellcomeopenresearch.org/articles/5-75),
but without covariates for the long-term conditions). The Italian life
tables are more similar to the GBD estimates than to either the UK or US
for women (as expected given the high life expectancy in Italy), but are
as expected higher for men, since the GBD uses the same life tables for
men and women. The Wales estimates have lower expected years of life
remaining at younger ages, but higher expected years remaining at older
ages. The UK and US life tables are available in our github repository,
and the Italian life tables are available at
<https://www.mortality.org/>.

### Figure A1 Life tables

![](Addendum_files/figure-gfm/le-1.png)<!-- -->

Table A1 shows the recalculated YLL estimates for COVID-19 deaths using
each of these lifetables (Table A1). While lower than the GBD-2010 based
estimates, particularly for men, the interpretation of our findings does
not change with the use of national life tables. Even after adjusting
for the number and type of long-term conditions, the YLL remained above
10; this is consistent with both our press release and the majority of
the recent press-coverage

### Table A1 Years of life lost by sex and life table, unadjusted and adjusted for number and type of long-term conditions

| Sex    | Life Table        | LTC number and type unadjusted | LTC number and type adjusted | Difference with adjustment (months) |
| :----- | :---------------- | :----------------------------- | ---------------------------: | ----------------------------------: |
| Total  | GBD 2010          | 13.3\*                         |                         12.4 |                                  11 |
| Total  | Italy 2017        | 11.3                           |                         10.8 |                                   6 |
| Total  | UK 1016-2018      | 11.1                           |                         10.6 |                                   6 |
| Total  | US 2017           | 10.9                           |                         10.4 |                                   6 |
| Total  | Wales (see paper) | 11                             |                         10.5 |                                   6 |
| Male   | GBD 2010          | 14\*                           |                         13.2 |                                  10 |
| Male   | Italy 2017        | 11.5                           |                         11.1 |                                   4 |
| Male   | UK 1016-2018      | 11.4                           |                         11.0 |                                   4 |
| Male   | US 2017           | 11.1                           |                         10.7 |                                   4 |
| Male   | Wales (see paper) | 11.2                           |                         10.9 |                                   4 |
| Female | GBD 2010          | 11.8\*                         |                         10.5 |                                  16 |
| Female | Italy 2017        | 10.9                           |                         10.2 |                                   9 |
| Female | UK 1016-2018      | 10.5                           |                          9.8 |                                   9 |
| Female | US 2017           | 10.5                           |                          9.7 |                                   9 |
| Female | Wales (see paper) | 10.3                           |                          9.6 |                                   9 |

*LTC: long-term conditions. \* The GBD-2010 life tables were only
available in 5-year bands, hence we estimated the YLL in the age-bands
as presented in the ISS report. For the remaining tables, since these
were available for single-years, we used single years of age derived
from our age models.*

We would continue to argue that national agencies should estimate the
YLL using their local age, sex and comorbidity data, where available. We
recommend that the choice of life tables should reflect the specific
question – for example national life tables should be used for
within-country comparisons and GBD or similar life tables for
international comparisons. However these results indicate that the
choice of life table does not materially change the overall finding,
that, on average, people dying from COVID-19 are losing around 10 years
of life.

# Why accounting for the number and type of long-term conditions does not have a large impact on the total years of life lost

Some commentators have been surprised that adjusting for the number and
type of long-term conditions does not have a larger impact on the YLL.
Here, we present some additional results to help explain this finding.
First, as has previously been described
(<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2424085/>), the effect of
most of the long-term conditions on mortality are attenuated with age
(Table A2). This may partly explain why the estimated YLL was not
greatly affected by adjusting for the number and type of long-term
conditions.

### Table A2 Hazard ratios at selected ages estimated separately for men and women having mutually adjusted for all other variables included in the model

| Sex    | Condition                             | Age 50 | Age 60 | Age 70 | Age 80 | Age 90 |
| :----- | :------------------------------------ | -----: | -----: | -----: | -----: | -----: |
| Female | Atrial fibrillation                   |   1.18 |   1.24 |   1.31 |   1.38 |   1.46 |
| Female | Cancer                                |   2.80 |   2.15 |   1.66 |   1.27 |   0.98 |
| Female | Chronic obstructive pulmonary disease |   3.89 |   2.95 |   2.24 |   1.70 |   1.29 |
| Female | Dementia                              |   4.98 |   3.99 |   3.19 |   2.56 |   2.05 |
| Female | Diabetes                              |   1.52 |   1.41 |   1.31 |   1.22 |   1.13 |
| Female | Heart failure                         |   2.09 |   1.87 |   1.67 |   1.50 |   1.34 |
| Female | Hypertension                          |   0.89 |   0.95 |   1.00 |   1.07 |   1.13 |
| Female | Ischaemic heart disease               |   1.55 |   1.42 |   1.31 |   1.20 |   1.10 |
| Female | Liver failure                         |   3.29 |   2.48 |   1.87 |   1.41 |   1.06 |
| Female | Renal failure                         |   1.27 |   1.24 |   1.20 |   1.17 |   1.14 |
| Female | Stroke                                |   1.58 |   1.51 |   1.44 |   1.37 |   1.30 |
| Male   | Atrial fibrillation                   |   1.13 |   1.18 |   1.24 |   1.30 |   1.37 |
| Male   | Cancer                                |   2.10 |   1.80 |   1.54 |   1.31 |   1.12 |
| Male   | Chronic obstructive pulmonary disease |   3.03 |   2.47 |   2.01 |   1.64 |   1.34 |
| Male   | Dementia                              |   3.49 |   3.17 |   2.89 |   2.63 |   2.39 |
| Male   | Diabetes                              |   1.54 |   1.39 |   1.26 |   1.14 |   1.03 |
| Male   | Heart failure                         |   2.09 |   1.88 |   1.68 |   1.51 |   1.35 |
| Male   | Hypertension                          |   0.98 |   1.01 |   1.04 |   1.07 |   1.11 |
| Male   | Ischaemic heart disease               |   1.28 |   1.24 |   1.20 |   1.17 |   1.13 |
| Male   | Liver failure                         |   5.65 |   3.55 |   2.23 |   1.40 |   0.88 |
| Male   | Renal failure                         |   1.21 |   1.22 |   1.23 |   1.24 |   1.25 |
| Male   | Stroke                                |   1.67 |   1.57 |   1.48 |   1.39 |   1.31 |

*Note that we also fit a model within the SAIL data repository for the
age-covariate interactions, treating age as a categorical variable
rather than as a continuous variable and this yielded almost identical
estimates for YLL in men and women as the model shown in Table A1
(GBD-2010 - YLL 13.2 for men and 10.5 for women).*

However, we do not think that this “attenuation effect” is the most
important factor for two reasons. Firstly, *within* the set of simulated
patients, the number of long-term conditions (multimorbidity count)
strongly impacted the YLL, even among the older groups. This can be seen
in Table 2 in [version one of our recent
manuscript](https://wellcomeopenresearch.org/articles/5-75), and in the
similar table below (Table A3) modified to use the UK life table rather
than the GBD-2010 life table.

### Table A3 Mean years of life lost among people dying from COVID-19, accounting for type of long-term conditions, by age-band, and multimorbidity count - with UK life table reference

#### Men

| Comorbidity count | 50-59 | 60-69 | 70-79 | 80 plus |
| ----------------: | ----: | ----: | ----: | ------: |
|                 0 | 30.98 | 23.06 | 15.79 |    9.35 |
|                 1 | 30.29 | 22.39 | 14.97 |    8.56 |
|                 2 | 25.28 | 18.59 | 12.29 |    6.81 |
|                 3 | 20.88 | 15.78 | 10.21 |    5.42 |
|                 4 | 19.56 | 13.19 |  8.42 |    4.02 |
|                 5 | 14.43 | 10.21 |  6.58 |    2.81 |
|                 6 |    \- |  4.05 |  5.09 |    1.90 |
|                 7 |    \- |  5.31 |  4.38 |    1.52 |
|                 8 |    \- |  4.24 |  2.81 |    1.06 |
|                 9 |    \- |  2.78 |  1.99 |    0.86 |
|                10 |    \- |    \- |  1.36 |    0.62 |
|                11 |    \- |    \- |    \- |    0.58 |

#### Women

| Comorbidity count | 50-59 | 60-69 | 70-79 | 80 plus |
| ----------------: | ----: | ----: | ----: | ------: |
|                 0 | 33.09 | 24.28 | 16.71 |    9.75 |
|                 1 | 32.86 | 24.32 | 16.22 |    8.26 |
|                 2 | 27.33 | 20.17 | 13.36 |    6.60 |
|                 3 | 24.36 | 16.93 | 11.19 |    5.29 |
|                 4 | 19.02 | 14.49 |  9.21 |    4.05 |
|                 5 | 15.20 | 10.65 |  7.48 |    3.09 |
|                 6 | 16.66 |  9.18 |  5.79 |    2.32 |
|                 7 |    \- |  7.18 |  4.14 |    1.95 |
|                 8 |    \- |  5.32 |  3.26 |    1.50 |
|                 9 |    \- |    \- |  2.28 |    1.17 |
|                10 |    \- |  2.27 |  1.83 |    0.77 |
|                11 |    \- |    \- |  0.91 |    0.86 |

Instead, we suspect that the main reason for the comparatively modest
impact of adjusting for multimorbidity count on the overall YLL is that
multimorbidity is common in older people among the general population,
not solely in people who die from COVID-19. Figure A2 shows the observed
(empirical) age-sex specific distribution of multimorbidity count in the
Welsh population alongside the modelled age-sex distribution (see
manuscript, this was derived from the marginal distributions from the
ISS report as well as from models of the association between
multimorbidity count and age among a small set of people coded as having
died from influenza in Wales) for the ISS deaths. For most age-groups,
multimorbidity counts were lower in the population than in people dying
from COVID-19, but among older people were nonetheless high for both.

### Figure A2 Distribution of multimorbidity by age and sex

![](Addendum_files/figure-gfm/compare_comorbidity-1.png)<!-- -->

If we assume that age and long-term condition count are independent
rather than associated, this does lead to a larger attenuation in YLL
(see Table 1 of the main manuscript). However, even then the attenuation
in YLL was less than two years, despite this being a somewhat extreme
and implausible assumption.

In conclusion, a major driver of the fact that accounting for long-term
condition count has only a modest effect on overall YLL is that, as well
as being common in people who die from COVID-19, long-term conditions
and multimorbidity are common in older people in the general population.

# The limitations of these findings for special populations such as care homes

Our work was completed before the recent concerns over COVID-19 deaths
in care homes became prominent and was not designed to address this
issue. It was designed instead to determine whether, since the presence
of long-term conditions and multimorbidity are common among people who
died from COVID-19, the years of life lost among people dying from
COVID-19 can be assumed to be low. Our findings demonstrate that this is
not the case.

Since inclusion in the Welsh dataset simply required having been being
registered with a participating general practice (GP), and since all
care home residents in the UK are registered with a local GP, people
resident in care homes will have been included in our survival models.
Moreover, although the comorbidity prevalence data from the ISS reports
was based solely on hospitalised patients, the age-sex distribution from
the ISS reports included anyone who died from COVID-19 provided they had
tested positive for SARS CoV2. Therefore, care home residents will have
contributed to the average YLL in the estimates we produced.

Nonetheless, we strongly agree that care home residents, are a special
population, in whom more severe disease, multimorbidity and frailty are
likely to be commoner. We also agree that there are good biological
reasons for suspecting that care home residents may be over-represented
among COVID-19 deaths compared to other causes of death (because it is
an infectious disease and people in care homes are in a communal
residence), and that the inclusion of more care home residents would
likely have lowered the YLL in our analyses. Therefore, we would argue
that the best approach for determining life expectancy in this group
would be to estimate it using data which includes care home residency
status.