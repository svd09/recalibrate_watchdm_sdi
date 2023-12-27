################################################################################
#                                                                              #
#                          WATCH-DM validation study                           #
#                                                                              #
################################################################################

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



setwd("P:\\ORD_Deo_202205002D\\DM\\DM_2010\\")

# get many functions now.

source('ckbplotr_theme.r')
source('cuminc_plot_function_new.r')
source('survplot_fn.r')
source('cif_plot.r')
source('code_ci.R')

setwd("P:\\ORD_Deo_202205002D\\DM\\watch_dm_validate\\")



# get the data 

df = read_csv("data/df_use.csv")

# now to get the variables to make the original 5 year predicted risk for HF
# first impute qrs duration, change hba1c to average blood glucose

df$bg = (28.7*df$hba1c_wi) - 46.7

summary(df$bg)

df$qrs = rnorm(n = length(df$scrssn),mean = 94,sd = 16.6)

summary(df$qrs)

df$qrs_long = with(df, ifelse(qrs > 120, 1, 0))

tabyl(df$qrs_long)

# df = tibble(
#   age = rnorm(n = 100, mean = 60, sd = 5),
#   fpg = rnorm(n = 100, mean = 75, sd = 4),
#   bmi = rnorm(n = 100, mean = 26, sd = 2),
#   sbp = rnorm(n = 100, mean = 130, sd = 4),
#   dbp = rnorm(n = 100, mean = 75, sd = 3),
#   creat = rnorm(n = 100, mean = 1.2, sd = 0.4),
#   hdl = rnorm(n = 100, mean = 30, sd = 4),
#   qrs = rbinom(n = 100,size = 1, prob = 0.05),
#   mi = rbinom(n = 100, size = 1, prob = 0.1),
#   cabg = rbinom(n = 100, size = 1, prob = 0.05))


summary(df)


S0 = 0.97358

# get baseline hazard 

s = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
            data = df)

b = basehaz(s,centered = F)

b[b$time == 5, ]

hr = c(1.04, 1.00, 1.06, 1.02, 0.98, 2.03, 
       0.98, 1.79, 1.82, 1.96)

coef = log(hr)

coef[2] <- 0.002
coef[5]<- -0.025
coef[7]<- -0.023

coef

d = df %>% select(age, bg, bmi_w, systolic_wi, diastolic_wi, creat_wi, hdl_wi, qrs_long,p_mi,p_cabg)


des_matr2 <- as.data.frame(model.matrix(~ age + bg + bmi_w + systolic_wi + 
                                          diastolic_wi + creat_wi + hdl_wi + 
                                          qrs_long + p_mi + p_cabg,data = d))


des_matr2$`(Intercept)` <- NULL


d$PI <- as.vector(as.matrix(des_matr2) %*% cbind(coef))

summary(d$PI)

bs = 0.97358

d$predict = as.vector(1 - bs**exp(d$PI))

summary(1 - d$predict)



df$orig_prediction = d$predict

mean(df$orig_prediction)

# fit cox model 

i = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_prediction,
          data = df, x = T)


df$y5rpred = 1 - pec::predictSurvProb(object = i,newdata = df,times = 5)


summary(df$y5rpred)



# see the group according to this 

df$q_original = cut2(df$orig_prediction,g = 5)

# see the cum incidence according to quintiles 

df$q_original = as.factor(df$q_original)

s = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ q_original,
            data = df)

s2 = summary(s, times = 5)

cumu_inc_summary(summary(s, times = 5), deci = 5,time.points = 1)

mean(df$orig_prediction)

sd(df$orig_prediction)

# calculate o/e


h_c = concordance(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_prediction, 
                  data = df,reverse = T)

h_c




# see actual predictd risk 

p = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df)

df$predictedrisk = 1 - pec::predictSurvProb(object = p,times = 5,newdata = df)

summary(df$predictedrisk)


