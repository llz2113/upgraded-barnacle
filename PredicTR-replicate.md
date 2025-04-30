PredicTR replicate
================
2025-03-18

## Install/load packages and import data

## Table S4

``` r
summary(predictr_data[biomarkers])
```

    ##       BCL2             COX2          CyclinD1         EGFRext      
    ##  Min.   :  0.00   Min.   :  0.0   Min.   :  0.00   Min.   :  0.00  
    ##  1st Qu.:  0.00   1st Qu.:100.0   1st Qu.: 26.67   1st Qu.: 76.67  
    ##  Median :  8.33   Median :146.1   Median : 99.58   Median :135.00  
    ##  Mean   : 62.64   Mean   :147.2   Mean   :126.54   Mean   :137.80  
    ##  3rd Qu.:100.00   3rd Qu.:190.8   3rd Qu.:225.00   3rd Qu.:195.83  
    ##  Max.   :300.00   Max.   :300.0   Max.   :300.00   Max.   :300.00  
    ##  NA's   :102      NA's   :99      NA's   :105      NA's   :106     
    ##     HPVISHHR           p16             PLK1           Survivin     
    ##  Min.   :0.0000   Min.   :  0.0   Min.   :  0.00   Min.   :  0.00  
    ##  1st Qu.:0.0000   1st Qu.:  0.0   1st Qu.: 15.00   1st Qu.: 40.00  
    ##  Median :0.0000   Median :220.0   Median : 27.50   Median : 62.50  
    ##  Mean   :0.3081   Mean   :158.7   Mean   : 35.57   Mean   : 66.52  
    ##  3rd Qu.:1.0000   3rd Qu.:280.0   3rd Qu.: 47.50   3rd Qu.: 90.00  
    ##  Max.   :1.0000   Max.   :300.0   Max.   :180.00   Max.   :192.50  
    ##  NA's   :99       NA's   :111     NA's   :143      NA's   :109     
    ##       TILS      
    ##  Min.   :1.000  
    ##  1st Qu.:2.000  
    ##  Median :2.000  
    ##  Mean   :2.077  
    ##  3rd Qu.:3.000  
    ##  Max.   :3.000  
    ##  NA's   :247

## Table S5

n are different - many more observations than in supplemental materials.

``` r
# Function for univariable Cox
univariable_cox <- function(data, vars) {
  results <- lapply(vars, function(var) {
    
    # Subset to non-NA for this variable, t.os, and Outcome
    valid_data <- data[!is.na(data[[var]]) & !is.na(data$t.os) & !is.na(data$Outcome), ]
    n_obs <- nrow(valid_data)
    
    formula <- as.formula(paste("Surv(t.os, Outcome) ~", var))
    model <- coxph(formula, data = data)
    coef_summary <- summary(model)$coefficients
    return(data.frame(Variable = var, 
                      Coefficient = coef_summary[, "coef"],
                      Lower = coef_summary[, "coef"] - 1.96 * coef_summary[, "se(coef)"],
                      Upper = coef_summary[, "coef"] + 1.96 * coef_summary[, "se(coef)"],
                      P = coef_summary[, "Pr(>|z|)"],
                      n = n_obs))
  })
  results_df <- bind_rows(results)
  # Adjust p-values using Benjamini & Hochberg method
  results_df$Q <- p.adjust(results_df$P, method = "BH")
  return(bind_rows(results_df))
}

colnames(data_processed)[colnames(data_processed) == 'Age at Diagnosis'] <- 'Age.at.Diagnosis'
colnames(data_processed)[colnames(data_processed) == 'Smoking status'] <- 'Smoking.status'

all_vars1 <- c("BCL2", "COX2", "CyclinD1", "EGFRext", "HIF1alpha", "PLK1","Survivin", 
              "p16_binary", "HPVISHHR","TILS", "Gender","Age.at.Diagnosis", "Surgery", 
              "RT", "CT", "T", "N", "Smoking.status")

# Summary of key columns
str(data_processed[, c("t.os", "Outcome", all_vars1)])
```

    ## tibble [985 × 20] (S3: tbl_df/tbl/data.frame)
    ##  $ t.os            : num [1:985] 6.84 7.62 9.79 1.39 5.08 ...
    ##  $ Outcome         : num [1:985] 0 1 0 1 0 1 1 1 1 1 ...
    ##  $ BCL2            : num [1:985] -0.209 2.268 -0.677 -0.587 -0.677 ...
    ##  $ COX2            : num [1:985] -0.493 0.656 0.388 -0.187 -0.723 ...
    ##  $ CyclinD1        : num [1:985] -0.554 -0.84 1.112 -0.126 -1.093 ...
    ##  $ EGFRext         : num [1:985] -0.899 0.192 1.942 -0.498 -0.765 ...
    ##  $ HIF1alpha       : num [1:985] 0.391 -0.611 -0.338 1.138 -0.192 ...
    ##  $ PLK1            : num [1:985] 0.311 -1.089 -0.842 -0.431 -0.842 ...
    ##  $ Survivin        : num [1:985] 0.1811 -0.3511 -0.3732 0.0925 1.4676 ...
    ##  $ p16_binary      : Factor w/ 2 levels "0","1": 2 2 1 2 2 2 1 1 1 2 ...
    ##  $ HPVISHHR        : Factor w/ 2 levels "0","1": 1 1 1 2 2 2 1 1 1 1 ...
    ##  $ TILS            : Factor w/ 3 levels "1","2","3": NA 1 2 2 3 2 2 3 3 2 ...
    ##  $ Gender          : Factor w/ 2 levels "F","M": 2 2 1 2 2 2 2 2 2 2 ...
    ##  $ Age.at.Diagnosis: Factor w/ 2 levels "0","1": 1 2 1 2 1 2 2 1 2 1 ...
    ##  $ Surgery         : Factor w/ 2 levels "0","1": 2 2 1 1 1 1 1 1 1 1 ...
    ##  $ RT              : Factor w/ 2 levels "0","1": 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ CT              : Factor w/ 2 levels "0","1": 1 2 2 2 2 2 2 2 2 2 ...
    ##  $ T               : Factor w/ 2 levels "1","2": 1 1 2 2 2 2 2 2 2 2 ...
    ##  $ N               : Factor w/ 2 levels "1","2": 2 2 2 2 1 2 1 1 2 2 ...
    ##  $ Smoking.status  : Factor w/ 3 levels "0","1","2": 2 NA 1 NA NA NA NA 3 2 3 ...

``` r
sapply(data_processed[, c("t.os", "Outcome", all_vars1)], function(x) sum(is.na(x)))
```

    ##             t.os          Outcome             BCL2             COX2 
    ##              150              107              102               99 
    ##         CyclinD1          EGFRext        HIF1alpha             PLK1 
    ##              105              106              116              143 
    ##         Survivin       p16_binary         HPVISHHR             TILS 
    ##              109              111               99              247 
    ##           Gender Age.at.Diagnosis          Surgery               RT 
    ##               64               66               79               83 
    ##               CT                T                N   Smoking.status 
    ##              138               99              100              182

``` r
complete_only <- data_processed[!is.na(data_processed$t.os) & !is.na(data_processed$Outcome),]

print(univariable_cox(data_processed, all_vars1))
```

    ##                         Variable  Coefficient       Lower       Upper
    ## ...1                        BCL2 -0.332222659 -0.47800722 -0.18643810
    ## ...2                        COX2 -0.087911935 -0.20804547  0.03222160
    ## ...3                    CyclinD1  0.557246735  0.43832239  0.67617108
    ## ...4                     EGFRext  0.255814916  0.13684608  0.37478375
    ## ...5                   HIF1alpha -0.016257808 -0.13416537  0.10164975
    ## ...6                        PLK1  0.090386867 -0.02338653  0.20416027
    ## ...7                    Survivin  0.040107761 -0.07923392  0.15944944
    ## ...8                  p16_binary -1.296072552 -1.54497634 -1.04716877
    ## ...9                    HPVISHHR -1.229607272 -1.57308297 -0.88613157
    ## TILS2                       TILS -0.653187599 -0.93596335 -0.37041185
    ## TILS3                       TILS -1.641254039 -2.05831809 -1.22418999
    ## ...12                     Gender  0.002812274 -0.25789742  0.26352197
    ## ...13           Age.at.Diagnosis  0.674595488  0.35075632  0.99843466
    ## ...14                    Surgery -0.634840553 -0.87040771 -0.39927340
    ## ...15                         RT -0.950833523 -1.27650917 -0.62515788
    ## ...16                         CT -0.156156763 -0.39513775  0.08282422
    ## ...17                          T  0.895390553  0.65677607  1.13400504
    ## ...18                          N  0.152093236 -0.08307214  0.38725861
    ## Smoking.status1   Smoking.status  0.452369976  0.05617262  0.84856733
    ## Smoking.status2   Smoking.status  1.189567793  0.83292367  1.54621192
    ##                            P   n            Q
    ## ...1            7.948500e-06 734 1.589700e-05
    ## ...2            1.514870e-01 737 2.019826e-01
    ## ...3            4.152899e-20 732 4.152899e-19
    ## ...4            2.503033e-05 730 4.550970e-05
    ## ...5            7.869628e-01 722 8.283819e-01
    ## ...6            1.194432e-01 701 1.706331e-01
    ## ...7            5.100839e-01 728 5.667599e-01
    ## ...8            1.864669e-24 729 3.729339e-23
    ## ...9            2.273385e-12 739 9.093542e-12
    ## TILS2           5.970510e-06 640 1.326780e-05
    ## TILS3           1.227950e-14 640 8.186335e-14
    ## ...12           9.831320e-01 809 9.831320e-01
    ## ...13           4.447468e-05 807 7.412446e-05
    ## ...14           1.277168e-07 805 3.192919e-07
    ## ...15           1.050536e-08 803 3.001532e-08
    ## ...16           2.002927e-01 755 2.410941e-01
    ## ...17           1.911913e-13 787 9.559565e-13
    ## ...18           2.049300e-01 790 2.410941e-01
    ## Smoking.status1 2.522839e-02 725 3.881290e-02
    ## Smoking.status2 6.256514e-11 725 2.085505e-10
