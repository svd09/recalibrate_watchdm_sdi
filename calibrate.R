####### code to create a recalibration equation for the WATCHDM score.

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


# first fit the cox model on the main dataaset and include only watchdm points
# in the model

model = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points,
data = df, x = T)

# now to create a model object from this model 


model_info = list(
	"coef" = stats::coef(model),
	"baseline_haz" = survival::basehaz(model, centered = F),
	"model_terms" = model[["terms"]])



str(model_info)

# obtain the baseline hazard at 5 years from this original model
# convert the baseline hazard to the S0t.5y 


t = model_info$baseline_haz

t$St.5y = exp(-t$hazard)

# this code can help to get the results that are very close to the value

t %>% filter(near(time, 2, tol = 1e-2)|
               near(time, 3, tol = 1e-2)|
               near(time, 4, tol = 1e-2)|
               near(time, 5, tol = 1e-2))

# so the S0.5y = 0.9904424

S0.5y = 0.9904424

# now that we have this information, we can obtain the predicted risk 
# in the cohort with least deprivation


df_sdi1 = df %>% filter(sdi_cat == 1)

# now to get the model matrix 

model_matrix = as.data.frame(model.matrix( ~ watch_dm_points, data = df_sdi1))

colnames(model_matrix)

model_matrix$`(Intercept)` <- NULL

# get the coefficients for the covariates in the model 

coef <- unname(model_info$coef)

coef

df_sdi1$linp = as.vector(as.matrix(model_matrix) %*% cbind(coef))

summary(df_sdi1$linp)

# get the 5-year predicted risk from the model 

df_sdi1$pred_5y = as.vector(1 - S0.5y**exp(df_sdi1$linp))

summary(df_sdi1$pred_5y)

# now to get the EO ration for recalibrating this estimate 
# 1. get the 5-year observed risk for HFH from df_sdi1


km = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df_sdi1)

# get the observed incidence of HF at 5 years

summary(km, times = 5)


obs = 1 - summary(km, times = 5)$surv

obs

# now to get the EO ratio

eo_ratio = mean(df_sdi1$pred_5y)/obs

eo_ratio

log(eo_ratio)


# use jf to get the new predicted estimates

df_sdi1$eo_pred = as.vector(1 - S0.5y**exp(df_sdi1$linp - log(eo_ratio)))

obs

# now as the observed is lower than the predicted the jf adjusted should be lower 
# than the raw model predicted

summary(df_sdi1$pred_5y)
summary(df_sdi1$eo_pred)

length(df_sdi1$pred_5y)
length(df_sdi1$eo_pred)

tiff("figures/q1hist.tiff", height = 6, width = 6, units = "in",res = 600)
hist(df_sdi1$pred_5y, col = t_col("red"), xlim = c(0, 0.25),
     main = "Quintile 1",xlab = "5-year HFH risk", ylab = "Count")
hist(df_sdi1$eo_pred, add = T, col = t_col("blue"))
dev.off()

# see the new calibration plot and do the GND test 
# do the GND test on the new data 

# for the GND test create groups of the predicted risk

source("scripts/gnd.R")

df_sdi1$eo_groups = as.numeric(cut2(df_sdi1$eo_pred, g = 10))



sdi1_gnd = GND.calib(pred = df_sdi1$eo_pred,
          tvar = df_sdi1$hf_admit_comp_years_5y,
          out = df_sdi1$event_cox_5y,
          cens.t = 5,
          groups = df_sdi1$eo_groups,
          adm.cens = 5)


sdi1_gnd$chi2gw
sdi1_gnd$pvalgw

# to evaluate re-calibration we will plot the calibration plot
# without adjustment

#===== create calibration plot with confidence intervals #====

df_sdi1$pred_5y_g = as.numeric(cut2(df_sdi1$pred_5y, g = 10))



sdi1_gnd_orig = GND.calib(pred = df_sdi1$pred_5y,
                     tvar = df_sdi1$hf_admit_comp_years_5y,
                     out = df_sdi1$event_cox_5y,
                     cens.t = 5,
                     groups = df_sdi1$pred_5y_g,
                     adm.cens = 5)



tiff("figures/recalib1.tiff", height = 5, width = 5, 
     units = "in",res = 600)
plot(x = sdi1_gnd$expectedperc,
     y = sdi1_gnd$kmperc, cex = 2, frame = F,
     col = "blue",
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH",
     ylab = "Observed 5-year HFH", main = "Quintile 1")
abline(0,1,col = "red", lty = 2)
points(x = sdi1_gnd_orig$expectedperc,
       y = sdi1_gnd_orig$kmperc, cex = 2, pch = 3, 
       col = "black")

legend("topleft",pch = c(1,3), col = c("blue","black"),
       legend = c("re-calibrated WATCH-DM score",
                  "original WATCH-DM score"))
dev.off()

#==== quintile 2 


df_sdi2 = df %>% filter(sdi_cat == 2)

# now to get the model matrix 

model_matrix2 = as.data.frame(model.matrix( ~ watch_dm_points, data = df_sdi2))

colnames(model_matrix2)

model_matrix2$`(Intercept)` <- NULL

#== the model_info the same for all the quintiles of the SDI.
# get the coefficients for the covariates in the model 

coef <- unname(model_info$coef)

coef

df_sdi2$linp = as.vector(as.matrix(model_matrix2) %*% cbind(coef))

summary(df_sdi2$linp)

# get the 5-year predicted risk from the model 

df_sdi2$pred_5y = as.vector(1 - S0.5y**exp(df_sdi2$linp))

summary(df_sdi2$pred_5y)

# now to get the EO ration for recalibrating this estimate 
# 1. get the 5-year observed risk for HFH from df_sdi1


km2 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df_sdi2)

# get the observed incidence of HF at 5 years

summary(km2, times = 5)


obs2 = 1 - summary(km2, times = 5)$surv

obs2




#== using the eoo method

eo_ratio2 = mean(df_sdi2$pred_5y)/obs2

eo_ratio

log(eo_ratio2)

df_sdi2$eo_pred = as.vector(1 - S0.5y**exp(df_sdi2$linp -log(eo_ratio2)))

obs2

# now as the observed is lower than the predicted the jf adjusted should be lower 
# than the raw model predicted

summary(df_sdi2$pred_5y)
summary(df_sdi2$eo_pred)
obs2


tiff("figures/q2hist.tiff", height = 6, width = 6, units = "in",res = 600)
hist(df_sdi2$pred_5y, col = t_col("red"), xlim = c(0, 0.25),
     main = "Quintile 2",xlab = "5-year HFH risk", ylab = "Count")
hist(df_sdi2$eo_pred, add = T, col = t_col("blue"))
dev.off()

# see the new calibration plot and do the GND test 
# do the GND test on the new data 

# for the GND test create groups of the predicted risk

source("scripts/gnd.R")

df_sdi2$eo_groups = as.numeric(cut2(df_sdi2$eo_pred, g = 10))



sdi2_gnd = GND.calib(pred = df_sdi2$eo_pred,
                     tvar = df_sdi2$hf_admit_comp_years_5y,
                     out = df_sdi2$event_cox_5y,
                     cens.t = 5,
                     groups = df_sdi2$eo_groups,
                     adm.cens = 5)


sdi2_gnd$chi2gw
sdi2_gnd$pvalgw

# to evaluate re-calibration we will plot the calibration plot
# without adjustment

#===== create calibration plot with confidence intervals #====

df_sdi2$pred_5y_g = as.numeric(cut2(df_sdi2$pred_5y, g = 10))



sdi2_gnd_orig = GND.calib(pred = df_sdi2$pred_5y,
                          tvar = df_sdi2$hf_admit_comp_years_5y,
                          out = df_sdi2$event_cox_5y,
                          cens.t = 5,
                          groups = df_sdi2$pred_5y_g,
                          adm.cens = 5)



tiff("figures/recalib2.tiff", height = 5, width = 5, 
     units = "in",res = 600)
plot(x = sdi2_gnd$expectedperc,
     y = sdi2_gnd$kmperc, cex = 2, frame = F,
     col = "blue",
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH",
     ylab = "Observed 5-year HFH", main = "Quintile 2")
abline(0,1,col = "red", lty = 2)
points(x = sdi2_gnd_orig$expectedperc,
       y = sdi2_gnd_orig$kmperc, cex = 2, pch = 3, 
       col = "black")

legend("topleft",pch = c(1,3), col = c("blue","black"),
       legend = c("re-calibrated WATCH-DM score",
                  "original WATCH-DM score"))
dev.off()

#==== quintile 3 


df_sdi3 = df %>% filter(sdi_cat == 3)

# now to get the model matrix 

model_matrix3 = as.data.frame(model.matrix( ~ watch_dm_points, data = df_sdi3))

colnames(model_matrix3)

model_matrix3$`(Intercept)` <- NULL

#== the model_info the same for all the quintiles of the SDI.
# get the coefficients for the covariates in the model 

coef <- unname(model_info$coef)

coef

df_sdi3$linp = as.vector(as.matrix(model_matrix3) %*% cbind(coef))

summary(df_sdi3$linp)

# get the 5-year predicted risk from the model 

df_sdi3$pred_5y = as.vector(1 - S0.5y**exp(df_sdi3$linp))

summary(df_sdi3$pred_5y)

# now to get the EO ration for recalibrating this estimate 
# 1. get the 5-year observed risk for HFH from df_sdi1


km3 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df_sdi3)

# get the observed incidence of HF at 5 years

summary(km3, times = 5)


obs3 = 1 - summary(km3, times = 5)$surv

obs3




#== using the eoo method

eo_ratio3 = mean(df_sdi3$pred_5y)/obs3

eo_ratio3

log(eo_ratio3)

df_sdi3$eo_pred = as.vector(1 - S0.5y**exp(df_sdi3$linp -log(eo_ratio3)))

obs3

# now as the observed is lower than the predicted the jf adjusted should be lower 
# than the raw model predicted

summary(df_sdi3$pred_5y)
summary(df_sdi3$eo_pred)
obs3


tiff("figures/q3hist.tiff", height = 6, width = 6, units = "in",res = 600)
hist(df_sdi3$pred_5y, col = t_col("red"), xlim = c(0, 0.25),
     main = "Quintile 3",xlab = "5-year HFH risk", ylab = "Count")
hist(df_sdi3$eo_pred, add = T, col = t_col("blue"))
dev.off()

# see the new calibration plot and do the GND test 
# do the GND test on the new data 

# for the GND test create groups of the predicted risk

source("scripts/gnd.R")

df_sdi3$eo_groups = as.numeric(cut2(df_sdi3$eo_pred, g = 10))



sdi3_gnd = GND.calib(pred = df_sdi3$eo_pred,
                     tvar = df_sdi3$hf_admit_comp_years_5y,
                     out = df_sdi3$event_cox_5y,
                     cens.t = 5,
                     groups = df_sdi3$eo_groups,
                     adm.cens = 5)


sdi3_gnd$chi2gw
sdi3_gnd$pvalgw

# to evaluate re-calibration we will plot the calibration plot
# without adjustment

#===== create calibration plot with confidence intervals #====

df_sdi3$pred_5y_g = as.numeric(cut2(df_sdi3$pred_5y, g = 10))



sdi3_gnd_orig = GND.calib(pred = df_sdi3$pred_5y,
                          tvar = df_sdi3$hf_admit_comp_years_5y,
                          out = df_sdi3$event_cox_5y,
                          cens.t = 5,
                          groups = df_sdi3$pred_5y_g,
                          adm.cens = 5)



tiff("figures/recalib3.tiff", height = 5, width = 5, 
     units = "in",res = 600)
plot(x = sdi3_gnd$expectedperc,
     y = sdi3_gnd$kmperc, cex = 2, frame = F,
     col = "blue",
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH",
     ylab = "Observed 5-year HFH", main = "Quintile 3")
abline(0,1,col = "red", lty = 2)
points(x = sdi3_gnd_orig$expectedperc,
       y = sdi3_gnd_orig$kmperc, cex = 2, pch = 3, 
       col = "black")

legend("topleft",pch = c(1,3), col = c("blue","black"),
       legend = c("re-calibrated WATCH-DM score",
                  "original WATCH-DM score"))
dev.off()
#== end q3

#=== quintile 4 


df_sdi4 = df %>% filter(sdi_cat == 4)

# now to get the model matrix 

model_matrix4 = as.data.frame(model.matrix( ~ watch_dm_points, data = df_sdi4))

colnames(model_matrix4)

model_matrix4$`(Intercept)` <- NULL

#== the model_info the same for all the quintiles of the SDI.
# get the coefficients for the covariates in the model 

coef <- unname(model_info$coef)

coef

df_sdi4$linp = as.vector(as.matrix(model_matrix4) %*% cbind(coef))

summary(df_sdi4$linp)

# get the 5-year predicted risk from the model 

df_sdi4$pred_5y = as.vector(1 - S0.5y**exp(df_sdi4$linp))

summary(df_sdi4$pred_5y)

# now to get the EO ration for recalibrating this estimate 
# 1. get the 5-year observed risk for HFH from the df quintile


km4 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df_sdi4)

# get the observed incidence of HF at 5 years

summary(km4, times = 5)


obs4 = 1 - summary(km4, times = 5)$surv

obs4




#== using the eoo method

eo_ratio4 = mean(df_sdi4$pred_5y)/obs4

eo_ratio4

log(eo_ratio4)

df_sdi4$eo_pred = as.vector(1 - S0.5y**exp(df_sdi4$linp -log(eo_ratio4)))

# now as the observed is lower than the predicted the jf adjusted should be lower 
# than the raw model predicted

summary(df_sdi4$pred_5y)
summary(df_sdi4$eo_pred)
obs4


tiff("figures/q4hist.tiff", height = 6, width = 6, units = "in",res = 600)
hist(df_sdi4$pred_5y, col = t_col("red"), xlim = c(0, 0.40),
     main = "Quintile 4",xlab = "5-year HFH risk", ylab = "Count")
hist(df_sdi4$eo_pred, add = T, col = t_col("blue"))
dev.off()

# see the new calibration plot and do the GND test 
# do the GND test on the new data 

# for the GND test create groups of the predicted risk

source("scripts/gnd.R")

df_sdi4$eo_groups = as.numeric(cut2(df_sdi4$eo_pred, g = 10))



sdi4_gnd = GND.calib(pred = df_sdi4$eo_pred,
                     tvar = df_sdi4$hf_admit_comp_years_5y,
                     out = df_sdi4$event_cox_5y,
                     cens.t = 5,
                     groups = df_sdi4$eo_groups,
                     adm.cens = 5)


sdi4_gnd$chi2gw
sdi4_gnd$pvalgw

# to evaluate re-calibration we will plot the calibration plot
# without adjustment

#===== create calibration plot with confidence intervals #====

df_sdi4$pred_5y_g = as.numeric(cut2(df_sdi4$pred_5y, g = 10))



sdi4_gnd_orig = GND.calib(pred = df_sdi4$pred_5y,
                          tvar = df_sdi4$hf_admit_comp_years_5y,
                          out = df_sdi4$event_cox_5y,
                          cens.t = 5,
                          groups = df_sdi4$pred_5y_g,
                          adm.cens = 5)



tiff("figures/recalib4.tiff", height = 5, width = 5, 
     units = "in",res = 600)

plot(x = sdi4_gnd$expectedperc,
     y = sdi4_gnd$kmperc, cex = 2, frame = F,
     col = "blue",
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH",
     ylab = "Observed 5-year HFH", main = "Quintile 4")
abline(0,1,col = "red", lty = 2)
points(x = sdi4_gnd_orig$expectedperc,
       y = sdi4_gnd_orig$kmperc, cex = 2, pch = 3, 
       col = "black")

legend("topleft",pch = c(1,3), col = c("blue","black"),
       legend = c("re-calibrated WATCH-DM score",
                  "original WATCH-DM score"))


dev.off()

#= end for quintile 4 

# now to recalibrate for the most deprived patients.

# first we need to get the data

df_sdi5 = df %>% filter(sdi_cat == 5)

# now to get the model matrix 

model_matrix5 = as.data.frame(model.matrix(~ watch_dm_points, data = df_sdi5))

model_matrix5

model_matrix5$`(Intercept)` <- NULL

# get the linear predictor for this data
# the coef is going to be same as before from the main model


df_sdi5$linp = as.vector(as.matrix(model_matrix5) %*% cbind(coef))

# now to get the predicted risk from the linp 

S0.5y = 0.9904424

df_sdi5$pred_5y = as.vector(1 - S0.5y**exp(df_sdi5$linp))

summary(df_sdi5$pred_5y)

# now to see the observed 


km5 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, 
data = df_sdi5)


summary(km5, times = 5)


obs5 = 1 - summary(km5, times = 5)$surv

obs5

# to use the eo correction 


eo_ratio5 = mean(df_sdi5$pred_5y)/obs5

eo_ratio5

log(eo_ratio5)


df_sdi5$eo_pred = as.vector(1 - S0.5y**exp(df_sdi5$linp -log(eo_ratio5)))

summary(df_sdi5$pred_5y)
summary(df_sdi5$eo_pred)
obs5

library(ggdist)

a = tibble(value = df_sdi5$pred_5y, group = 1)
b = tibble(value = df_sdi5$eo_pred, group = 2)

c = rbind(a,b)

ggplot(data = c, aes(x = value, y = factor(group))) + stat_dots() + xlim(0, 0.3)
