# script for changes for the revision 1 for the paper

# get the libraries


library(easypackages)
libraries(c("tidyverse","survival","rms","riskRegression","Hmisc","janitor",
            "haven","readxl","cmprsk","etm","survminer","mstate","geepack","prodlim",
            "boot","rsample","webshot"))



#=======================================================
# loading some functions 



etm_to_df <- function(object, ci.fun = "cloglog", level = 0.95, ...) {
  l.X <- ncol(object$X)
  l.trans <- nrow(object[[1]]$trans)
  res <- list()
  for (i in seq_len(l.X)) {
    temp <- summary(object[[i]], ci.fun = ci.fun, level = level)
    res[[i]] <- data.table::rbindlist(
      temp[object$failcode + 1], idcol = "CIF"
    )[, CIF := paste0("CIF", CIF, "; ", names(object)[i])]
  }
  do.call(rbind, res)
}


win = function(a, na.rm = T){
  low = quantile(a, 0.01, na.rm = T)
  high = quantile(a, 0.99, na.rm = T)
  
  b = ifelse(a <= low, low, 
             ifelse(a >= high, high, a))
  
  return(b)
}



cumu_inc_summary = function(sum_obj, deci = 3,time.points = 1){
  
  class = attr(sum_obj,"class")
  
  
  #  stopifnot(class == "summary.survfit")
  
  if(class != "summary.survfit") 
    stop("Error: Not a surfit object")
  
  cum_inc = format((1 - sum_obj$surv)*100, nsmall = deci+2) 
  upper = format((1 - sum_obj$lower)*100, nsmall = deci+2) 
  lower = format((1 - sum_obj$upper)*100, nsmall = deci+2) 
  
  cum_inc2 = cum_inc
  u2 = upper
  l2 = lower
  
  #  cum_inc2 = as.numeric(format(round(cum_inc,4), nsmall = deci))
  #  u2 = as.numeric(format(round(upper,4), nsmall = deci))
  #  l2 = as.numeric(format(round(lower,4), nsmall = deci))
  
  x = time.points
  
  
  
  result = tibble(
    group = sum_obj$strata,
    time = sum_obj$time,
    total_patients = rep(sum_obj$n,time.points),
    nrisk = sum_obj$n.risk,
    nevent = sum_obj$n.event,
    ncensor = sum_obj$n.censor,
    `cumulative_incidence(%)` = cum_inc2,
    `lower95(%)` = l2,
    `upper95(%)` = u2)
  
  
  result
}


# t_col ##############################

t_col <- function(color, percent = 80, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
#======================================end


# now to get the data 


setwd("P:\\ORD_Deo_202205002D\\DM\\DM_2010\\")

# get many functions now.

source('ckbplotr_theme.r')
source('cuminc_plot_function_new.r')
source('survplot_fn.r')
source('cif_plot.r')
source('code_ci.R')

# now to get the main dataset and then to fit a model in the main dataset.

setwd("P:\\ORD_Deo_202205002D\\DM\\watch_dm_validate\\")

df = read_csv("data/df_use.csv")


# changes outlined & done for revision 1.

# 1. age range for the data -

glimpse(df)


summary(df$age)

summary(df$age_w)

# 2. fit model using only sdi as covariate & then report the R2 and goodness of fit.

tabyl(df$sdi_cat)
summary(df$sdi_score)


df2 = df %>% filter(!is.na(sdi_cat))



model = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points,
              data = df2, x = T)

summary(model)

model2 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points + sdi_cat ,
               data = df2, x = T)

summary(model2)

anova(model, model2)

AIC(model); AIC(model2)

# comparing models 
# > model = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points,
#                 +               data = df2, x = T)
# > 
#   > summary(model)
# Call:
#   coxph(formula = Surv(hf_admit_comp_years_5y, event_cox_5y) ~ 
#           watch_dm_points, data = df2, x = T)
# 
# n= 1012315, number of events= 54356 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# watch_dm_points 0.117839  1.125063 0.001128 104.4   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# watch_dm_points     1.125     0.8888     1.123     1.128
# 
# Concordance= 0.621  (se = 0.001 )
# Likelihood ratio test= 10568  on 1 df,   p=<2e-16
# Wald test            = 10907  on 1 df,   p=<2e-16
# Score (logrank) test = 10884  on 1 df,   p=<2e-16
# 
# > 
#   > model2 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points + sdi_cat ,
#                    +                data = df2, x = T)
# > 
#   > summary(model2)
# Call:
#   coxph(formula = Surv(hf_admit_comp_years_5y, event_cox_5y) ~ 
#           watch_dm_points + sdi_cat, data = df2, x = T)
# 
# n= 1012315, number of events= 54356 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)    
# watch_dm_points 0.120231  1.127758 0.001127 106.68   <2e-16 ***
#   sdi_cat         0.131415  1.140441 0.003309  39.71   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# watch_dm_points     1.128     0.8867     1.125     1.130
# sdi_cat             1.140     0.8769     1.133     1.148
# 
# Concordance= 0.63  (se = 0.001 )
# Likelihood ratio test= 12166  on 2 df,   p=<2e-16
# Wald test            = 12541  on 2 df,   p=<2e-16
# Score (logrank) test = 12523  on 2 df,   p=<2e-16
# 
# > 
#   > anova(model, model2)
# Analysis of Deviance Table
# Cox model: response is  Surv(hf_admit_comp_years_5y, event_cox_5y)
# Model 1: ~ watch_dm_points
# Model 2: ~ watch_dm_points + sdi_cat
# loglik  Chisq Df Pr(>|Chi|)    
# 1 -740140                         
# 2 -739342 1597.4  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > 
#   > AIC(model); AIC(model2)
# [1] 1480283
# [1] 1478687

BIC(model)
BIC(model2)
# 
# > BIC(model)
# [1] 1480292
# > BIC(model2)
# [1] 1478705



#- get the various social  & zip code level information to report the general characteristics

library(zipcodeR)

glimpse(df)

glimpse(zip_info)

# get the distribution of median income for the VA pop



zipcodes = unique(df$ZIP)

zip_info = zipcodeR::reverse_zipcode(zip_code = zipcodes)

summary(zip_info$median_household_income)

summary(zip_info$median_home_value)

# get the sdi components 

sdic = read_csv("data/sdi_components.csv")

glimpse(sdic)

# limit to only my zip codes 

sdic2 = sdic %>% filter(ZCTA5_FIPS %in% zipcodes)

summary(sdic2$pct_Poverty_LT100)
summary(sdic2)


# SHOW THAT CALIBRATION IS OFF EVEN WITH DECILES 

summary(df$sdi_score)

# remove the missing 

df2 = df %>% drop_na(sdi_score)

# now to see the decile 

df2$sdi.2 = df2$sdi_score/10

df2$sdi10 = with(df2, 
                 ifelse(sdi_score <= 1 , 1,
                        ifelse(
                          sdi_score >1 & sdi_score <= 2, 2, 
                          ifelse(
                            sdi_score > 2 & sdi_score <= 3, 3, 
                            ifelse(
                              sdi_score > 3 & sdi_score <= 4, 4, 
                              ifelse(
                                sdi_score > 4 & sdi_score <= 5, 5, 
                                ifelse(
                                  sdi_score > 5 & sdi_score <= 6, 6, 
                                  ifelse(
                                    sdi_score > 6 & sdi_score <= 7, 7, 
                                    ifelse(
                                      sdi_score > 7 & sdi_score <= 8, 8, 
                                      ifelse(
                                        sdi_score > 8 & sdi_score <= 9, 9, 10
                                      )
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )))


# now to get the predicted risk with original watch dm & then add sdi groups 

model10 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points,
              data = df2, x = T)

model3 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points + sdi10 ,
               data = df2, x = T)

summary(model3)

anova(model10, model3)

# REPORT AGE ADJUSTED MORTALITY RATES



m = coxph(Surv(fupyears_5y, died_5y) ~ age, data = df, x = T, y = T)

predictCox(object = m,times = c(1,3,5),
           newdata = tibble(age = 67),type = "survival",se = T)


# for each SDI - not going to report 


tapply(INDEX = df$sdi_cat,X = df$age_w,FUN = summary)



m1 = coxph(Surv(fupyears_5y, died_5y) ~ age, data = df[df$sdi_cat == 1, ], x = T, y = T)

predictCox(object = m1,times = 5,newdata = tibble(age = 69.58),type = "survival",se = T)




m2 = coxph(Surv(fupyears_5y, died_5y) ~ age, data = df[df$sdi_cat == 2, ], x = T, y = T)

predictCox(object = m2,times = 5,newdata = tibble(age = 68.48),type = "survival",se = T)




m3 = coxph(Surv(fupyears_5y, died_5y) ~ age, data = df[df$sdi_cat == 3, ], x = T, y = T)

predictCox(object = m3,times = 5,newdata = tibble(age = 67.59),type = "survival",se = T)




m4 = coxph(Surv(fupyears_5y, died_5y) ~ age, data = df[df$sdi_cat == 4, ], x = T, y = T)

predictCox(object = m4,times = 5,newdata = tibble(age = 66.82),type = "survival",se = T)



m5 = coxph(Surv(fupyears_5y, died_5y) ~ age, data = df[df$sdi_cat == 5, ], x = T, y = T)

predictCox(object = m5,times = 5,newdata = tibble(age = 65.71),type = "survival",se = T)




