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
# Run once: install.packages(c("dplyr", "ggplot2", "epitools", "rineq"))
library(dplyr) ## for data wrangling
library(ggplot2) ## for concentration curves
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
and its 95% from a log-binomial regression model.

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
    ## Between 1-2 Miles 1.0000000        NA       NA
    ## Between 2-3 Miles 0.9158879 0.3914167 2.143114
    ## Between 3-4 Miles 1.3611111 0.5877271 3.152183
    ## Between 4-5 Miles 0.0000000 0.0000000      NaN
    ## More than 5 Miles 1.3207547 0.5099476 3.420730
    ## Up to 1 Mile      0.9152542 0.5096420 1.643684

``` r
## Prevalence rate ratios based on map-based proximity (binned)
tab_map = with(ehr_food, table(MAP_PROXIMITY_BINNED, DIABETES))
rr_map = data.frame(riskratio(tab_map)$measure)
rr_map
```

    ##                    estimate     lower    upper
    ## Between 1-2 Miles 1.0000000        NA       NA
    ## Between 2-3 Miles 0.4797980 0.1860006 1.237663
    ## Between 3-4 Miles 0.7885375 0.3300255 1.884071
    ## Between 4-5 Miles 1.3623560 0.6335271 2.929652
    ## More than 5 Miles 0.7112299 0.3124352 1.619050
    ## Up to 1 Mile      0.7254545 0.3851491 1.366443

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

    ##      RII       LB       UB 
    ## 1.143123 0.458853 2.813689

``` r
## Relative index of inequality based on straight-line proximity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = MAP_PROXIMITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 1.3751221 0.5711289 3.3464406

### Concentration Curve

### Concentration Index

## Analysis of Neighborhood Food Insecurity

### Rate Ratios

#### Binned Exposure

``` r
## Prevalence rate ratios based on survey estimate food insecurity (binned)
tab_est = with(ehr_food, table(EST_FOOD_INSECURITY_BINNED, DIABETES))
rr_est = data.frame(riskratio(tab_est)$measure)
rr_est
```

    ##                  estimate     lower    upper
    ## Between 10-12.5% 1.000000        NA       NA
    ## Between 12.5-15% 2.636095 0.5981745 11.61700
    ## Between 15-17.5% 2.402542 0.5120758 11.27218
    ## Between 17.5-20% 3.344037 0.7424673 15.06138
    ## More than 20%    3.300000 0.7928257 13.73568
    ## Up to 10%        2.745763 0.6290040 11.98595

``` r
## Prevalence rate ratios based on lower bound food insecurity (binned)
tab_lb = with(ehr_food, table(LB_FOOD_INSECURITY_BINNED, DIABETES))
rr_lb = data.frame(riskratio(tab_lb)$measure)
rr_lb
```

    ##                  estimate     lower     upper
    ## Between 10-12.5% 1.000000        NA        NA
    ## Between 12.5-15% 3.405172 1.1544189 10.044187
    ## Between 15-17.5% 2.194444 0.4179449 11.522062
    ## Between 17.5-20% 2.644351 0.9006491  7.763950
    ## More than 20%    4.459677 1.5053810 13.211754
    ## Up to 10%        2.455959 0.8078658  7.466256

``` r
## Prevalence rate ratios based on upper bound food insecurity (binned)
tab_ub = with(ehr_food, table(UB_FOOD_INSECURITY_BINNED, DIABETES))
rr_ub = data.frame(riskratio(tab_ub)$measure)
rr_ub
```

    ##                   estimate      lower    upper
    ## Between 10-12.5% 1.0000000         NA       NA
    ## Between 12.5-15% 0.4320388 0.08104507 2.303133
    ## Between 15-17.5% 1.5689103 0.51479836 4.781444
    ## Between 17.5-20% 1.8541667 0.57832076 5.944684
    ## More than 20%    1.7752660 0.64184010 4.910209
    ## Up to 10%        1.7115385 0.53313022 5.494650

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
    ## 0.6137938 0.2513117 1.4754482

``` r
## Relative index of inequality based on lower bound food insecurity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = LB_FOOD_INSECURITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 0.5600026 0.2294239 1.3453199

``` r
## Relative index of inequality based on upper bound food insecurity (binned)
mod_rii(health_var = DIABETES, 
        group_by_var = UB_FOOD_INSECURITY_BINNED, 
        data = ehr_food)
```

    ##       RII        LB        UB 
    ## 0.4569237 0.1798752 1.1279066

### Concentration Curve

### Concentration Index
