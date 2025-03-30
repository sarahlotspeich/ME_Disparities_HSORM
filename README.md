Examining the Impacts of Measurement Error in Quantifying Health
Disparities: A Case Study on Diabetes and the Food Environment
================

This repository contains `R` code to reproduce results from the
manuscript by Hung, Green, and Lotspeich (2025+).

## Packages

Our implementation of the health disparities methods in the paper rely
on the following `R` packages, which can be installed and loaded as
follows.

``` r
# Run once: install.packages(c("dplyr", "ggplot2", "scales", "epitools", "rineq"))
library(dplyr) ## for data wrangling
library(ggplot2) ## for concentration curves
library(scales) ## for % labels on concentration curves
library(epitools) ## for rate ratios
library(rineq) ## for concentration indices
```

## Data

With permission from the Institutional Review Board at Wake Forest
University School of Medicine, a de-identified version of the study
dataset have been made available with the paper as Supplementary
Material. Once downloaded, these data can be loaded as follows.

``` r
# Read in the data from CSV
ehr_food = read.csv("https://raw.githubusercontent.com/sarahlotspeich/ME_Disparities_HSORM/refs/heads/main/me_disparities_analysis_2020.csv")

## Inspect data summary
ehr_food |> 
  summary()
```

    ##        ID           DIABETES       MAP_PROXIMITY       STRAIGHT_PROXIMITY 
    ##  Min.   :  1.0   Min.   :0.00000   Min.   : 0.000029   Min.   : 0.000029  
    ##  1st Qu.:231.8   1st Qu.:0.00000   1st Qu.: 0.950276   1st Qu.: 0.607415  
    ##  Median :462.5   Median :0.00000   Median : 1.701083   Median : 1.140446  
    ##  Mean   :462.5   Mean   :0.06818   Mean   : 2.511577   Mean   : 1.737758  
    ##  3rd Qu.:693.2   3rd Qu.:0.00000   3rd Qu.: 3.572429   3rd Qu.: 2.290866  
    ##  Max.   :924.0   Max.   :1.00000   Max.   :15.474724   Max.   :10.193788  
    ##  MAP_PROXIMITY_BINNED STRAIGHT_PROXIMITY_BINNED EST_FOOD_INSECURITY
    ##  Length:924           Length:924                Min.   : 5.30      
    ##  Class :character     Class :character          1st Qu.:11.90      
    ##  Mode  :character     Mode  :character          Median :16.10      
    ##                                                 Mean   :16.08      
    ##                                                 3rd Qu.:20.90      
    ##                                                 Max.   :29.20      
    ##  LB_FOOD_INSECURITY UB_FOOD_INSECURITY EST_FOOD_INSECURITY_BINNED
    ##  Min.   : 4.60      Min.   : 6.00      Length:924                
    ##  1st Qu.:10.60      1st Qu.:13.20      Class :character          
    ##  Median :14.50      Median :18.00      Mode  :character          
    ##  Mean   :14.46      Mean   :17.76                                
    ##  3rd Qu.:19.00      3rd Qu.:22.80                                
    ##  Max.   :26.70      Max.   :31.60                                
    ##  LB_FOOD_INSECURITY_BINNED UB_FOOD_INSECURITY_BINNED
    ##  Length:924                Length:924               
    ##  Class :character          Class :character         
    ##  Mode  :character          Mode  :character         
    ##                                                     
    ##                                                     
    ## 

``` r
## Define levels for binned factors
prox_levels = c("More than 5 Miles", "Between 4-5 Miles", "Between 3-4 Miles", 
                "Between 2-3 Miles", "Between 1-2 Miles", "Up to 1 Mile")
insec_levels = c("More than 20%", "Between 17.5-20%", "Between 15-17.5%", 
                 "Between 12.5-15%", "Between 10-12.5%", "Up to 10%")
### Order from most disadvantaged --> least disadvantaged
ehr_food = ehr_food |> 
  mutate(MAP_PROXIMITY_BINNED = factor(x = MAP_PROXIMITY_BINNED, 
                                       levels = prox_levels), 
         STRAIGHT_PROXIMITY_BINNED = factor(x = STRAIGHT_PROXIMITY_BINNED, 
                                            levels = prox_levels), 
         EST_FOOD_INSECURITY_BINNED = factor(x = EST_FOOD_INSECURITY_BINNED, 
                                             levels = insec_levels), 
         LB_FOOD_INSECURITY_BINNED = factor(x = LB_FOOD_INSECURITY_BINNED, 
                                            levels = insec_levels),
         UB_FOOD_INSECURITY_BINNED = factor(x = UB_FOOD_INSECURITY_BINNED, 
                                            levels = insec_levels))
```

The dataset contains the following columns.

- `ID`: a deidentified patient identifier
- `DIABETES`: a binary indicator of being diagnosed with type-2 diabetes
  (`=0` if no, `=1` if yes)
- `MAP_PROXIMITY`: unbinned map-based proximity to healthy foods (in
  miles)
- `STRAIGHT_PROXIMITY`: unbinned straight-line proximity to healthy
  foods (in miles)
- `MAP_PROXIMITY_BINNED`: binned map-based proximity to healthy foods
  (categorical)
- `STRAIGHT_PROXIMITY_BINNED`: binned straight-line proximity to healthy
  foods (categorical)
- `EST_FOOD_INSECURITY`: unbinned survey estimate of neighborhood rate
  of food insecurity (percent, between `0` and `100`)
- `LB_FOOD_INSECURITY`: unbinned lower bound on neighborhood rate of
  food insecurity, taken from the 95% confidence interval (percent,
  between `0` and `100`)
- `UB_FOOD_INSECURITY`: unbinned upper bound on neighborhood rate of
  food insecurity, taken from the 95% confidence interval (percent,
  between `0` and `100`)
- `EST_FOOD_INSECURITY_BINNED`: binned survey estimate of neighborhood
  rate of food insecurity (categorical)
- `LB_FOOD_INSECURITY_BINNED`: binned lower bound on neighborhood rate
  of food insecurity, taken from the 95% confidence interval
  (categorical)
- `UB_FOOD_INSECURITY_BINNED`: binned upper bound on neighborhood rate
  of food insecurity, taken from the 95% confidence interval
  (categorical)

Each row represents data for one patient.

## Useful Functions

The `logbin_rr` function, defined below, calculates the rate ratio (RR)
and its 95% confidence interval from a log-binomial regression model.

``` r
# Write a function that calculates RR (95% CI) from a log-binomial model
### formula: Model formula, like Y ~ X 
### data: Dataset where variables in the model formula can be found
logbin_rr = function(formula, data) {
  fit = glm(formula = formula, 
            data = data, 
            family = binomial(link = "log"))
  to_return = suppressMessages(
    c(exp(coefficients(fit)[2]), 
      exp(confint(object = fit))[2, ])
  )
  names(to_return) = c("RR", "LB", "UB")
  return(to_return)
}
```

The `mod_rii` function, defined below, calculates the relative index of
inequality (RII) and its 95% confidence interval using the Poisson
regression approach by [Moreno-Betancur et
al. (2015)](https://journals.lww.com/epidem/abstract/2015/07000/relative_index_of_inequality_and_slope_index_of.12.aspx).

``` r
# Write a function that calculates RII (95% CI) from a Poisson model
## Arguments: 
### health_var: name of health outcome 
### group_by_var: name of the socioeconomic exposure to group by (ordered from most- to least-disadvantaged)
### data: dataset where health_var and group_by_var can be found
mod_rii = function(health_var, group_by_var, data) {
  prev_by_group = data |> 
    group_by({{group_by_var}}) |> 
    summarize(sum_health = sum({{health_var}}), 
              num = n()) |> 
    ungroup() |> 
    arrange({{group_by_var}}) |> 
    mutate(est_prev = sum_health / num, 
           prop_num = num / sum(num), 
           cum_prop_num = cumsum(prop_num))
  cdf = data.frame(
      sum_health = 0,
      num = 0,
      est_prev = 0, 
      prop_num = 0, 
      cum_prop_num = 0) |> 
      bind_rows(prev_by_group) |> 
      mutate(rank = (1 - cum_prop_num) + 1 / 2 * prop_num) |> 
      slice(-1)
  rii = glm(formula = sum_health ~ rank, 
            offset = log(num), 
            family = "poisson", 
            data = cdf)
  to_return = suppressMessages(exp(c(rii$coefficients[2], confint(rii)[2, ])))
  names(to_return) = c("RII", "LB", "UB")
  return(to_return)
}
```

The `make_cc_dat` function, defined below, constructs the ranked
dataset, which can be used to build a concentration curve.

``` r
# Write a function that create the dataset for a concentration curve
## Arguments: 
### health_var: name of health outcome 
### group_by_var: name of the socioeconomic exposure to group by (ordered from most- to least-disadvantaged)
### data: dataset where health_var and group_by_var can be found
make_cc_dat = function(health_var, group_by_var, data) {
  data |> 
    group_by({{group_by_var}}) |> 
    summarize(num = n(), 
              sum_health = sum({{health_var}})) |> 
    arrange({{group_by_var}}) |> 
    mutate(cumsum_num = cumsum(num), 
           cumprop_num = cumsum_num / sum(num),
           cumsum_health = cumsum(sum_health), 
           cumprop_health = cumsum_health / sum(sum_health)) |>
    bind_rows(
      data.frame(num = 0, 
                 sum_health = 0, 
                 cumsum_num = 0, 
                 cumprop_num = 0, 
                 cumsum_health = 0, 
                 cumprop_health = 0)
      ) |> 
    mutate(CI = cumprop_num * lead(x = cumprop_health, n = 1, default = 0) -
                    cumprop_health * lead(x = cumprop_num, n = 1, default = 0))
}
```

The `calc_grouped_ci` function, defined below, calculates the
concentration index with its 95% confidence interval from grouped data.

``` r
# Write a function that calculate the concentration index (C) from grouped data
## Arguments: 
### health_var: name of health outcome 
### group_by_var: name of the socioeconomic exposure to group by (ordered from most- to least-disadvantaged)
### data: dataset where health_var and group_by_var can be found
calc_grouped_ci = function(group_by_var, health_var, data) {
  # Save useful constants
  n = nrow(data) ## sample size
  
  # Build summary dataset by group
  cc_dat = data |> 
    group_by({{group_by_var}}) |> 
    summarize(num = dplyr::n(), 
              prop = num / n,
              cases = sum({{health_var}}), 
              mean_health = mean({{health_var}}), 
              var_health = var({{health_var}})) |> 
    arrange({{group_by_var}}) |> 
    mutate(cumsum_num = cumsum(num), 
           cumprop_num = cumsum_num / sum(num),
           cumsum_cases = cumsum(cases), 
           cumprop_cases = cumsum_cases / sum(cases)) |> 
    mutate(CI = cumprop_num * lead(x = cumprop_cases, n = 1, default = 0) -
                    cumprop_cases * lead(x = cumprop_num, n = 1, default = 0), 
           R = lag(x = cumprop_num, n = 1, default = 0) + 1 / 2 * prop)
  
  # Compute concentration index
  C = sum(cc_dat$CI)
  
  # Compute variance for C 
  Ybar = data |> 
    pull({{health_var}}) |> 
    mean()
  q = c(0, 1 / Ybar * cumsum(cc_dat$mean_health * cc_dat$prop)) ## augment with q0 = 0 
  a = (cc_dat$mean_health / Ybar) * (2 * cc_dat$R - 1 - C) + 2 - q[-length(q)] - q[-1] 
  varC = 1 / n * (sum(cc_dat$prop * a ^ 2) - (1 + C) ^ 2) + 
    1 / (n * Ybar ^ 2) * sum(cc_dat$prop * cc_dat$var_health * (2 * cc_dat$R - 1 - C) ^ 2)
  
  # Return C with its 95% CI 
  to_return = c(C, C + c(-1.96, 1.96) * sqrt(varC))
  names(to_return) = c("CI", "LB", "UB")
  return(to_return)
}
```

## Analysis of Household Proximity to Healthy Foods

### Rate Ratios

#### Binned Exposure

``` r
## Prevalence rate ratios based on straight-line proximity (binned)
tab_straight = with(ehr_food, table(STRAIGHT_PROXIMITY_BINNED, DIABETES))
rr_straight = data.frame(riskratio(tab_straight)$measure)
rr_straight
```

    ##                    estimate     lower    upper
    ## More than 5 Miles 1.0000000        NA       NA
    ## Between 4-5 Miles 0.0000000 0.0000000      NaN
    ## Between 3-4 Miles 1.0305556 0.3459918 3.069567
    ## Between 2-3 Miles 0.6934579 0.2309697 2.082022
    ## Between 1-2 Miles 0.7571429 0.2923353 1.960986
    ## Up to 1 Mile      0.6929782 0.2788354 1.722231

``` r
## Prevalence rate ratios based on map-based proximity (binned)
tab_map = with(ehr_food, table(MAP_PROXIMITY_BINNED, DIABETES))
rr_map = data.frame(riskratio(tab_map)$measure)
rr_map
```

    ##                    estimate     lower    upper
    ## More than 5 Miles 1.0000000        NA       NA
    ## Between 4-5 Miles 1.9154930 0.7255041 5.057329
    ## Between 3-4 Miles 1.1086957 0.3856808 3.187107
    ## Between 2-3 Miles 0.6746032 0.2201190 2.067471
    ## Between 1-2 Miles 1.4060150 0.6176463 3.200664
    ## Up to 1 Mile      1.0200000 0.4272495 2.435111

#### Unbinned Exposure

``` r
## Prevalence rate ratio based on straight-line proximity (unbinned)
logbin_rr(formula = DIABETES ~ STRAIGHT_PROXIMITY, 
          data = ehr_food)
```

    ##        RR        LB        UB 
    ## 1.0012044 0.8569016 1.1427091

``` r
## Prevalence rate ratio based on map-based proximity (unbinned)
logbin_rr(formula = DIABETES ~ MAP_PROXIMITY,
          data = ehr_food)
```

    ##        RR        LB        UB 
    ## 1.0136741 0.9055229 1.1152475

### Relative Index of Inequality

``` r
## Relative index of inequality based on straight-line proximity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = STRAIGHT_PROXIMITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 1.1196108 0.4485444 2.7463805

``` r
## Relative index of inequality based on straight-line proximity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = MAP_PROXIMITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 1.0966057 0.4541682 2.6378360

### Concentration Curve

#### Binned Exposure

``` r
## Define plot colors
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6", "#ff914d")

## Make dataset for concentration curves (both distance types)
cc_dat = make_cc_dat(health_var = DIABETES, 
                     group_by_var = STRAIGHT_PROXIMITY_BINNED, 
                     data = ehr_food) |> 
  mutate(PROXIMITY = "Straight-Line (Binned)") |> 
  bind_rows(
    make_cc_dat(health_var = DIABETES, 
                group_by_var = MAP_PROXIMITY_BINNED, 
                data = ehr_food) |> 
      mutate(PROXIMITY = "Map-Based (Binned)")
  ) |> 
  mutate(PROXIMITY = factor(x = PROXIMITY, 
                            levels = c("Straight-Line (Binned)", "Map-Based (Binned)")))

## Use ggplot() to build concentration curves (both distance types)
cc_dat|> 
  ggplot(aes(x = cumprop_num, y = cumprop_health, color = PROXIMITY)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black") + 
  geom_line(linewidth = 1.2, alpha = 1) + 
  coord_equal() + 
  theme_minimal(base_size = 14) + 
  labs(x = "Cumulative Proportion of Patients\n(Ranked by Proximity)", 
       y = "Cumulative Proportion of Diabetes Cases") + 
  theme(legend.text = element_text(color = "black"), 
        legend.title = element_text(color = "black", face = "bold"),
        axis.title = element_text(color = "black", face = "bold"),
        legend.position = "top") + 
  scale_x_continuous(labels = percent_format()) + 
  scale_y_continuous(labels = percent_format()) + 
  scale_color_manual(values = cols, name = "Distance")
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

#### Unbinned Exposure

``` r
## Create negative proximity to rank from farthest (most disadvantaged) --> nearest (least disadvantaged)
ehr_food = ehr_food |> 
  dplyr::mutate(NEG_STRAIGHT_PROXIMITY = - STRAIGHT_PROXIMITY, 
                NEG_MAP_PROXIMITY = - MAP_PROXIMITY)

## Make dataset for concentration curves (both distance types)
cc_dat = make_cc_dat(health_var = DIABETES, 
                     group_by_var = NEG_STRAIGHT_PROXIMITY, 
                     data = ehr_food) |> 
  mutate(PROXIMITY = "Straight-Line") |> 
  bind_rows(
    make_cc_dat(health_var = DIABETES, 
                group_by_var = NEG_MAP_PROXIMITY, 
                data = ehr_food) |> 
      mutate(PROXIMITY = "Map-Based")
  ) |> 
  mutate(PROXIMITY = factor(x = PROXIMITY, 
                            levels = c("Straight-Line", "Map-Based")))

## Use ggplot() to build concentration curves (both distance types)
cc_dat |> 
  ggplot(aes(x = cumprop_num, y = cumprop_health, color = PROXIMITY)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black") + 
  geom_line(linewidth = 1.2, alpha = 1) + 
  coord_equal() + 
  theme_minimal(base_size = 14) + 
  labs(x = "Cumulative Proportion of Patients\n(Ranked by Proximity)", 
       y = "Cumulative Proportion of Diabetes Cases") + 
  theme(axis.title = element_text(color = "black", face = "bold"), 
        legend.title = element_text(color = "black", face = "bold"), 
        legend.position = "top") + 
  scale_x_continuous(labels = percent_format()) + 
  scale_y_continuous(labels = percent_format()) + 
  scale_color_manual(values = cols, name = "Distance")
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Concentration Index

#### Binned Exposure

``` r
## Concentration index based on straight-line proximity (binned)
calc_grouped_ci(group_by_var = STRAIGHT_PROXIMITY_BINNED, 
                health_var = DIABETES, 
                data = ehr_food)
```

    ##         CI         LB         UB 
    ## -0.0168522 -0.1470394  0.1133350

``` r
## Concentration index based on map-based proximity (binned)
calc_grouped_ci(group_by_var = MAP_PROXIMITY_BINNED, 
                health_var = DIABETES, 
                data = ehr_food)
```

    ##          CI          LB          UB 
    ## -0.01461898 -0.14701374  0.11777578

#### Unbinned Exposure

``` r
## Concentration index based on straight-line proximity (unbinned)
ci_straight = ci(ineqvar = ehr_food$NEG_STRAIGHT_PROXIMITY, ### Food environment exposure
                 outcome = ehr_food$DIABETES, ### Health outcome 
                 type = "CIw", ### Wagstaff's correction for binary health outcome
                 robust_se = TRUE) ### Use sandwich standard errors for confidence interval 
summary(ci_straight)
```

    ## Call:
    ## [1] "ci(ineqvar = ehr_food$NEG_STRAIGHT_PROXIMITY, outcome = ehr_food$DIABETES, "
    ## [2] "    type = \"CIw\", robust_se = TRUE)"                                      
    ## 
    ## Type of Concentration Index:
    ##  CIw 
    ## 
    ## Health Concentration Index:
    ##  -0.003078738 
    ## 
    ## Variance:
    ##  0.005611874 
    ## 
    ## 95% Confidence Interval:
    ##  -0.1499044 0.1437469

``` r
## Concentration index based on map-based proximity (unbinned)
ci_map = ci(ineqvar = ehr_food$NEG_MAP_PROXIMITY, ### Food environment exposure
                 outcome = ehr_food$DIABETES, ### Health outcome 
                 type = "CIw", ### Wagstaff's correction for binary health outcome
                 robust_se = TRUE) ### Use sandwich standard errors for confidence interval 
summary(ci_map)
```

    ## Call:
    ## [1] "ci(ineqvar = ehr_food$NEG_MAP_PROXIMITY, outcome = ehr_food$DIABETES, "
    ## [2] "    type = \"CIw\", robust_se = TRUE)"                                 
    ## 
    ## Type of Concentration Index:
    ##  CIw 
    ## 
    ## Health Concentration Index:
    ##  -0.02026068 
    ## 
    ## Variance:
    ##  0.00553609 
    ## 
    ## 95% Confidence Interval:
    ##  -0.1660916 0.1255703

## Analysis of Neighborhood Food Insecurity

### Rate Ratios

#### Binned Exposure

``` r
## Prevalence rate ratios based on survey estimate food insecurity (binned)
tab_est = with(ehr_food, table(EST_FOOD_INSECURITY_BINNED, DIABETES))
rr_est = data.frame(riskratio(tab_est)$measure)
rr_est
```

    ##                   estimate      lower    upper
    ## More than 20%    1.0000000         NA       NA
    ## Between 17.5-20% 1.0133445 0.48205805 2.130173
    ## Between 15-17.5% 0.7280431 0.31983259 1.657263
    ## Between 12.5-15% 0.7988166 0.39756781 1.605029
    ## Between 10-12.5% 0.3030303 0.07280309 1.261311
    ## Up to 10%        0.8320493 0.42264900 1.638017

``` r
## Prevalence rate ratios based on lower bound food insecurity (binned)
tab_lb = with(ehr_food, table(LB_FOOD_INSECURITY_BINNED, DIABETES))
rr_lb = data.frame(riskratio(tab_lb)$measure)
rr_lb
```

    ##                   estimate      lower     upper
    ## More than 20%    1.0000000         NA        NA
    ## Between 17.5-20% 0.5929468 0.29929127 1.1747283
    ## Between 15-17.5% 0.4920635 0.11724109 2.0652015
    ## Between 12.5-15% 0.7635468 0.38260805 1.5237623
    ## Between 10-12.5% 0.2242315 0.07569018 0.6642837
    ## Up to 10%        0.5507032 0.26345554 1.1511392

``` r
## Prevalence rate ratios based on upper bound food insecurity (binned)
tab_ub = with(ehr_food, table(UB_FOOD_INSECURITY_BINNED, DIABETES))
rr_ub = data.frame(riskratio(tab_ub)$measure)
rr_ub
```

    ##                   estimate      lower    upper
    ## More than 20%    1.0000000         NA       NA
    ## Between 17.5-20% 1.0444444 0.49484619 2.204451
    ## Between 15-17.5% 0.8837607 0.45442442 1.718730
    ## Between 12.5-15% 0.2433657 0.05913981 1.001472
    ## Between 10-12.5% 0.5632959 0.20365731 1.558020
    ## Up to 10%        0.9641026 0.45584225 2.039069

#### Unbinned Exposure

``` r
## Prevalence rate ratio based on survey estimate food insecurity (unbinned)
logbin_rr(formula = DIABETES ~ EST_FOOD_INSECURITY, 
          data = ehr_food)
```

    ##        RR        LB        UB 
    ## 1.0325328 0.9857904 1.0819680

``` r
## Prevalence rate ratio based on lower bound food insecurity (unbinned)
logbin_rr(formula = DIABETES ~ LB_FOOD_INSECURITY, 
          data = ehr_food)
```

    ##        RR        LB        UB 
    ## 1.0356549 0.9852063 1.0891817

``` r
## Prevalence rate ratio based on upper bound food insecurity (unbinned)
logbin_rr(formula = DIABETES ~ UB_FOOD_INSECURITY, 
          data = ehr_food)
```

    ##        RR        LB        UB 
    ## 1.0302874 0.9867722 1.0761693

### Relative Index of Inequality

``` r
## Relative index of inequality based on survey estimate food insecurity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = EST_FOOD_INSECURITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 1.6530293 0.6871525 4.0523816

``` r
## Relative index of inequality based on lower bound food insecurity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = LB_FOOD_INSECURITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 2.2605885 0.9371285 5.5732896

``` r
## Relative index of inequality based on upper bound food insecurity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = UB_FOOD_INSECURITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 1.8074586 0.7338411 4.6124180

### Concentration Curve

#### Binned Exposure

``` r
## Make dataset for concentration curves (all values)
cc_dat = make_cc_dat(group_by_var = LB_FOOD_INSECURITY_BINNED, 
                     health_var = DIABETES, 
                     data = ehr_food)  |> 
  mutate(FOOD_INSECURITY = "Lower Bound (Binned)") |> 
  bind_rows(
    make_cc_dat(group_by_var = EST_FOOD_INSECURITY_BINNED, 
                     health_var = DIABETES, 
                     data = ehr_food) |> 
      mutate(FOOD_INSECURITY = "Survey Estimate (Binned)")
  ) |> 
  bind_rows(
    make_cc_dat(group_by_var = UB_FOOD_INSECURITY_BINNED, 
                     health_var = DIABETES, 
                     data = ehr_food) |> 
      mutate(FOOD_INSECURITY = "Upper Bound (Binned)")
  ) |> 
  mutate(FOOD_INSECURITY = factor(x = FOOD_INSECURITY, 
                                  levels = c("Lower Bound (Binned)", 
                                             "Survey Estimate (Binned)", 
                                             "Upper Bound (Binned)")))

## Use ggplot() to build concentration curves (both distance types)
cc_dat |> 
  ggplot(aes(x = cumprop_num, y = cumprop_health, color = FOOD_INSECURITY)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black") + 
  geom_line(linewidth = 1.2, alpha = 1) + 
  coord_equal() + 
  theme_minimal(base_size = 14) + 
  labs(x = "Cumulative Proportion of Patients\n(Ranked by Food Insecurity)", 
       y = "Cumulative Proportion of Diabetes Cases") + 
  theme(axis.title = element_text(color = "black", face = "bold"), 
        legend.title = element_text(color = "black", face = "bold"), 
        legend.position = "top") + 
  scale_x_continuous(labels = percent_format()) + 
  scale_y_continuous(labels = percent_format()) + 
  scale_color_manual(values = rev(cols[-length(cols)]), name = "Distance")
```

![](README_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

#### Unbinned Exposure

``` r
## Create negative neighborhood food insecurity to rank from highest (most disadvantaged) --> lowest (least disadvantaged)
ehr_food = ehr_food |> 
  dplyr::mutate(NEG_EST_FOOD_INSECURITY = - EST_FOOD_INSECURITY, 
                NEG_LB_FOOD_INSECURITY = - LB_FOOD_INSECURITY, 
                NEG_UB_FOOD_INSECURITY = - UB_FOOD_INSECURITY)

## Make dataset for concentration curves (all values)
cc_dat = make_cc_dat(group_by_var = NEG_LB_FOOD_INSECURITY, 
                     health_var = DIABETES, 
                     data = ehr_food) |> 
  mutate(FOOD_INSECURITY = "Lower Bound") |> 
  bind_rows(
    make_cc_dat(group_by_var = NEG_EST_FOOD_INSECURITY, 
                     health_var = DIABETES, 
                     data = ehr_food) |> 
      mutate(FOOD_INSECURITY = "Survey Estimate")
  ) |> 
  bind_rows(
    make_cc_dat(group_by_var = NEG_UB_FOOD_INSECURITY, 
                     health_var = DIABETES, 
                     data = ehr_food) |> 
      mutate(FOOD_INSECURITY = "Upper Bound")
  ) |> 
  mutate(FOOD_INSECURITY = factor(x = FOOD_INSECURITY, 
                                  levels = c("Lower Bound", 
                                             "Survey Estimate", 
                                             "Upper Bound")))

## Use ggplot() to build concentration curves (both distance types)
cc_dat |> 
  ggplot(aes(x = cumprop_num, y = cumprop_health, color = FOOD_INSECURITY)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black") + 
  geom_line(linewidth = 1.2, alpha = 1) + 
  coord_equal() + 
  theme_minimal(base_size = 14) + 
  labs(x = "Cumulative Proportion of Patients\n(Ranked by Food Insecurity)", 
       y = "Cumulative Proportion of Diabetes Cases") + 
  theme(axis.title = element_text(color = "black", face = "bold"), 
        legend.title = element_text(color = "black", face = "bold"), 
        legend.position = "top") + 
  scale_x_continuous(labels = percent_format()) + 
  scale_y_continuous(labels = percent_format()) + 
  scale_color_manual(values = rev(cols[-length(cols)]), name = "Distance")
```

![](README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Concentration Index

#### Binned Exposure

``` r
## Concentration index based on survey estimate food insecurity (binned)
calc_grouped_ci(group_by_var = EST_FOOD_INSECURITY_BINNED, 
                health_var = DIABETES, 
                data = ehr_food)
```

    ##          CI          LB          UB 
    ## -0.07950251 -0.21598215  0.05697713

``` r
## Concentration index based on lower bound food insecurity (binned)
calc_grouped_ci(group_by_var = LB_FOOD_INSECURITY_BINNED, 
                health_var = DIABETES, 
                data = ehr_food)
```

    ##           CI           LB           UB 
    ## -0.128908129 -0.267642306  0.009826048

``` r
## Concentration index based on upper bound food insecurity (binned)
calc_grouped_ci(group_by_var = UB_FOOD_INSECURITY_BINNED, 
                health_var = DIABETES, 
                data = ehr_food)
```

    ##          CI          LB          UB 
    ## -0.08888202 -0.22069648  0.04293245

#### Unbinned Exposure

``` r
## Concentration index based on survey estimate food insecurity (unbinned)
ci_est = ci(ineqvar = ehr_food$NEG_EST_FOOD_INSECURITY, ### Food environment exposure
            outcome = ehr_food$DIABETES, ### Health outcome 
            type = "CIw", ### Wagstaff's correction for binary health outcome
            robust_se = TRUE) ### Use sandwich standard errors for confidence interval 
summary(ci_est)
```

    ## Call:
    ## [1] "ci(ineqvar = ehr_food$NEG_EST_FOOD_INSECURITY, outcome = ehr_food$DIABETES, "
    ## [2] "    type = \"CIw\", robust_se = TRUE)"                                       
    ## 
    ## Type of Concentration Index:
    ##  CIw 
    ## 
    ## Health Concentration Index:
    ##  -0.1603709 
    ## 
    ## Variance:
    ##  0.006078574 
    ## 
    ## 95% Confidence Interval:
    ##  -0.3131799 -0.007561918

``` r
## Concentration index based on lower bound food insecurity (unbinned)
ci_lb = ci(ineqvar = ehr_food$NEG_LB_FOOD_INSECURITY, ### Food environment exposure
           outcome = ehr_food$DIABETES, ### Health outcome 
           type = "CIw", ### Wagstaff's correction for binary health outcome
           robust_se = TRUE) ### Use sandwich standard errors for confidence interval 
summary(ci_lb)
```

    ## Call:
    ## [1] "ci(ineqvar = ehr_food$NEG_LB_FOOD_INSECURITY, outcome = ehr_food$DIABETES, "
    ## [2] "    type = \"CIw\", robust_se = TRUE)"                                      
    ## 
    ## Type of Concentration Index:
    ##  CIw 
    ## 
    ## Health Concentration Index:
    ##  -0.1602603 
    ## 
    ## Variance:
    ##  0.006087245 
    ## 
    ## 95% Confidence Interval:
    ##  -0.3131783 -0.007342358

``` r
## Concentration index based on upper bound food insecurity (unbinned)
ci_ub = ci(ineqvar = ehr_food$NEG_UB_FOOD_INSECURITY, ### Food environment exposure
           outcome = ehr_food$DIABETES, ### Health outcome 
           type = "CIw", ### Wagstaff's correction for binary health outcome
           robust_se = TRUE) ### Use sandwich standard errors for confidence interval 
summary(ci_ub)
```

    ## Call:
    ## [1] "ci(ineqvar = ehr_food$NEG_UB_FOOD_INSECURITY, outcome = ehr_food$DIABETES, "
    ## [2] "    type = \"CIw\", robust_se = TRUE)"                                      
    ## 
    ## Type of Concentration Index:
    ##  CIw 
    ## 
    ## Health Concentration Index:
    ##  -0.1601128 
    ## 
    ## Variance:
    ##  0.006077284 
    ## 
    ## 95% Confidence Interval:
    ##  -0.3129056 -0.007320041
