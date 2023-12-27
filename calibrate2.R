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

#= we need to the figure for WATCHDM points overall 

df$orig_cat = factor(df$orig_cat)



h = df %>% ggplot(aes(x = watch_dm_points, fill = orig_cat)) + geom_bar() + 
     theme_ckb() + theme(legend.position = "none")  

ggsave(plot = h, filename = "figures/new_bar.tiff", height = 8, width = 6,
       units = "in", device = "tiff",dpi = 600)


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


# use the eo ratio to get the new predicted estimates

df_sdi1$eo_pred = as.vector(1 - S0.5y**exp(df_sdi1$linp - log(eo_ratio)))

obs

# now as the observed is lower than the predicted the jf adjusted should be lower 
# than the raw model predicted

summary(df_sdi1$pred_5y)
summary(df_sdi1$eo_pred)
obs

length(df_sdi1$pred_5y)
length(df_sdi1$eo_pred)

# tiff("figures/q1hist.tiff", height = 6, width = 6, units = "in",res = 600)
# hist(df_sdi1$pred_5y, col = t_col("red"), xlim = c(0, 0.25),
#      main = "Quintile 1",xlab = "5-year HFH risk", ylab = "Count")
# hist(df_sdi1$eo_pred, add = T, col = t_col("blue"))
# dev.off()

#- make a mirror histogram for the predicted values
#- using this in the paper...


a = df_sdi1$pred_5y
b = df_sdi1$eo_pred


tiff("figures/q1hist_new.tiff", height = 6, width = 6, units = "in",res = 600)

par(mfrow = c(2,1))

par(mar = c(0,5,3,3))
hist(a*100,  xlim = c(0, 25), ylim = c(0,40000), xaxt = "n",
     las = 1, breaks = 25, xlab = " ", ylab= " ", col = "gray40",
     main = "Quintile 1")
abline(v = mean(a)*100, lty = 2, col = "red")

par(mar = c(5,5,0,3))

hist(b*100,  
     xlim = c(0,25),ylim = c(40000,0),
     xlab = "5-year HFH risk (%)", las = 1, breaks = 25,
 col = "gray60", main = " ", ylab = " ")
abline(v = mean(b)*100, lty = 2, col = "red")

legend("bottomright",
       title = "5-Year HFH Predicted Risk (%)", fill = c("gray40","gray60"), 
       legend = c("Original WATCH-DM Score", "Re-calibrated WATCH-DM"))


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

sdi1_gnd_orig$chi2gw


tiff("figures/recalib1.tiff", height = 5, width = 5, 
     units = "in",res = 600)

plot(x = sdi1_gnd$expectedperc,
     y = sdi1_gnd$kmperc, cex = 1, frame = F,
     col = "red",pch = 16,
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH (From the WATCH-DM Score)",
     ylab = "Observed 5-year HFH", main = "Quintile 1")
lines(x = sdi1_gnd$expectedperc,
      y = sdi1_gnd$kmperc, lty = 2, col = "red")

abline(0,1,col = "blue", lty = 2)

points(x = sdi1_gnd_orig$expectedperc,
       y = sdi1_gnd_orig$kmperc, cex = 1, pch = 15, 
       col = "black")
lines(x = sdi1_gnd_orig$expectedperc,
      y = sdi1_gnd_orig$kmperc, lty = 2, col = "black")



legend("topleft",pch = c(16,15), col = c("red","black"),
       legend = c("WATCH-DM score\n re-calibrated for the SDI",
                  "WATCH-DM score"), bty = "n")

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


# tiff("figures/q2hist.tiff", height = 6, width = 6, units = "in",res = 600)
# hist(df_sdi2$pred_5y, col = t_col("red"), xlim = c(0, 0.25),
#      main = "Quintile 2",xlab = "5-year HFH risk", ylab = "Count")
# hist(df_sdi2$eo_pred, add = T, col = t_col("blue"))
# dev.off()


#= new histogram for the paper


tiff("figures/q2hist_new.tiff", height = 6, width = 6, units = "in",res = 600)

par(mfrow = c(2,1))

par(mar = c(0,5,3,3))
hist(df_sdi2$pred_5y*100,  xlim = c(0, 25), ylim = c(0,60000), xaxt = "n",
     las = 1, breaks = 25, xlab = " ", ylab= " ", col = "gray40",
     main = "Quintile 2")
abline(v = mean(df_sdi2$pred_5y)*100, lty = 2, col = "red")

par(mar = c(5,5,0,3))

hist(df_sdi2$eo_pred*100,  
     xlim = c(0,25),ylim = c(60000,0),
     xlab = "5-year HFH risk (%)", las = 1, breaks = 25,
     col = "gray60", main = " ", ylab = " ")
abline(v = mean(df_sdi2$eo_pred)*100, lty = 2, col = "red")

# legend("bottomright",
#        title = "5-Year HFH Predicted Risk (%)", fill = c("gray40","gray60"), 
#        legend = c("Original WATCH-DM Score", "Re-calibrated WATCH-DM"))


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

sdi2_gnd_orig$chi2gw


tiff("figures/recalib2.tiff", height = 5, width = 5, 
     units = "in",res = 600)


plot(x = sdi2_gnd$expectedperc,
     y = sdi2_gnd$kmperc, cex = 1, frame = F,
     col = "red",pch = 15,
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH (From the WATCH-DM Score)",
     ylab = "Observed 5-year HFH", main = "Quintile 2")

lines(x = sdi2_gnd$expectedperc,
      y = sdi2_gnd$kmperc, col = "red", lty = 2)

abline(0,1,col = "blue", lty = 2)


points(x = sdi2_gnd_orig$expectedperc,
       y = sdi2_gnd_orig$kmperc, cex = 1, pch = 16, 
       col = "black")
lines(x = sdi2_gnd_orig$expectedperc,
      y = sdi2_gnd_orig$kmperc, col = "black", lty = 2)


# legend("topleft",pch = c(16,15), col = c("blue","black"),
#        legend = c("WATCH-DM score\n re-calibrated for the SDI",
#                   "WATCH-DM score"), bty = "n")

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


# tiff("figures/q3hist.tiff", height = 6, width = 6, units = "in",res = 600)
# hist(df_sdi3$pred_5y, col = t_col("red"), xlim = c(0, 0.25),
#      main = "Quintile 3",xlab = "5-year HFH risk", ylab = "Count")
# hist(df_sdi3$eo_pred, add = T, col = t_col("blue"))
# dev.off()


#= mirror histogram for the paper 


tiff("figures/q3hist_new.tiff", height = 6, width = 6, units = "in",res = 600)

par(mfrow = c(2,1))

par(mar = c(0,5,3,3))
hist(df_sdi3$pred_5y*100,  xlim = c(0, 25), ylim = c(0,80000), xaxt = "n",
     las = 1, breaks = 25, xlab = " ", ylab= " ", col = "gray40",
     main = "Quintile 3")
abline(v = mean(df_sdi3$pred_5y)*100, lty = 2, col = "red")

par(mar = c(5,5,0,3))

hist(df_sdi3$eo_pred*100,  
     xlim = c(0,25),ylim = c(80000,0),
     xlab = "5-year HFH risk (%)", las = 1, breaks = 25,
     col = "gray60", main = " ", ylab = " ")
abline(v = mean(df_sdi3$eo_pred)*100, lty = 2, col = "red")

# legend("bottomright",
#        title = "5-Year HFH Predicted Risk (%)", fill = c("gray40","gray60"), 
#        legend = c("Original WATCH-DM Score", "Re-calibrated WATCH-DM"))


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

sdi3_gnd_orig$chi2gw


tiff("figures/recalib3.tiff", height = 5, width = 5, 
     units = "in",res = 600)


plot(x = sdi3_gnd$expectedperc,
     y = sdi3_gnd$kmperc, cex = 1, frame = F,
     col = "red",pch = 16,
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH (From the WATCH-DM Score)",
     ylab = "Observed 5-year HFH", main = "Quintile 3")

lines(x = sdi3_gnd$expectedperc,
      y = sdi3_gnd$kmperc, col = "red", lty = 2)

abline(0,1,col = "blue", lty = 2)


points(x = sdi3_gnd_orig$expectedperc,
       y = sdi3_gnd_orig$kmperc, cex = 1, pch = 15, 
       col = "black")

lines(x = sdi3_gnd_orig$expectedperc,
      y = sdi3_gnd_orig$kmperc, col = "black", lty = 2)

# legend("topleft",pch = c(16,15), col = c("blue","black"),
#        legend = c("WATCH-DM score\n re-calibrated for the SDI",
#                   "WATCH-DM score"), bty = "n")


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


# tiff("figures/q4hist.tiff", height = 6, width = 6, units = "in",res = 600)
# hist(df_sdi4$pred_5y, col = t_col("red"), xlim = c(0, 0.40),
#      main = "Quintile 4",xlab = "5-year HFH risk", ylab = "Count")
# hist(df_sdi4$eo_pred, add = T, col = t_col("blue"))
# dev.off()


#= mirror histogram 


tiff("figures/q4hist_new.tiff", height = 6, width = 6, units = "in",res = 600)

par(mfrow = c(2,1))

par(mar = c(0,5,3,3))
hist(df_sdi4$pred_5y*100,  xlim = c(0, 25), ylim = c(0,90000), xaxt = "n",
     las = 1, breaks = 25, xlab = " ", ylab= " ", col = "gray40",
     main = "Quintile 4")
abline(v = mean(df_sdi4$pred_5y)*100, lty = 2, col = "red")

par(mar = c(5,5,0,3))

hist(df_sdi4$eo_pred*100,  
     xlim = c(0,25),ylim = c(90000,0),
     xlab = "5-year HFH risk (%)", las = 1, breaks = 25,
     col = "gray60", main = " ", ylab = " ")
abline(v = mean(df_sdi4$eo_pred)*100, lty = 2, col = "red")

# legend("bottomright",
#        title = "5-Year HFH Predicted Risk (%)", fill = c("gray40","gray60"), 
#        legend = c("Original WATCH-DM Score", "Re-calibrated WATCH-DM"))


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


sdi4_gnd_orig$chi2gw


tiff("figures/recalib4.tiff", height = 5, width = 5, 
     units = "in",res = 600)

plot(x = sdi4_gnd$expectedperc,
     y = sdi4_gnd$kmperc, cex = 1, frame = F,
     col = "red",pch = 16,
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH (From the WATCH-DM Score)",
     ylab = "Observed 5-year HFH", main = "Quintile 4")

lines(x = sdi4_gnd$expectedperc,
      y = sdi4_gnd$kmperc, 
      col = "red", lty = 2)

abline(0,1,col = "blue", lty = 2)


points(x = sdi4_gnd_orig$expectedperc,
       y = sdi4_gnd_orig$kmperc, cex = 1, pch = 15, 
       col = "black")

lines(x = sdi4_gnd_orig$expectedperc,
      y = sdi4_gnd_orig$kmperc, lty = 2,
      col = "black")




dev.off()

#= end for quintile 4 

#==== Quintile 5 +++++++++++++++++++++++++++++++++++++++++

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





# tiff("figures/q5hist.tiff", height = 6, width = 6, units = "in",res = 600)
# 
# hist(df_sdi5$pred_5y, col = t_col("red"), xlim = c(0, 0.40),
#      main = "Quintile 5",xlab = "5-year HFH risk", ylab = "Count")
# hist(df_sdi5$eo_pred, add = T, col = t_col("blue"))
# 
# dev.off()

#= mirror histogram 


tiff("figures/q5hist_new.tiff", height = 6, width = 6, units = "in",res = 600)

par(mfrow = c(2,1))

par(mar = c(0,5,3,3))
hist(df_sdi5$pred_5y*100,  xlim = c(0, 25), ylim = c(0,80000), xaxt = "n",
     las = 1, breaks = 25, xlab = " ", ylab= " ", col = "gray40",
     main = "Quintile 5")
abline(v = mean(df_sdi5$pred_5y)*100, lty = 2, col = "red")

par(mar = c(5,5,0,3))

hist(df_sdi5$eo_pred*100,  
     xlim = c(0,25),ylim = c(80000,0),
     xlab = "5-year HFH risk (%)", las = 1, breaks = 25,
     col = "gray60", main = " ", ylab = " ")
abline(v = mean(df_sdi5$eo_pred)*100, lty = 2, col = "red")

# legend("bottomright",
#        title = "5-Year HFH Predicted Risk (%)", fill = c("gray40","gray60"), 
#        legend = c("Original WATCH-DM Score", "Re-calibrated WATCH-DM"))


dev.off()


# see the new calibration plot and do the GND test 
# do the GND test on the new data 

# for the GND test create groups of the predicted risk

source("scripts/gnd.R")

df_sdi5$eo_groups = as.numeric(cut2(df_sdi5$eo_pred, g = 10))


sdi5_gnd = GND.calib(pred = df_sdi5$eo_pred,
                          tvar = df_sdi5$hf_admit_comp_years_5y,
                          out = df_sdi5$event_cox_5y,
                          cens.t = 5,
                          groups = df_sdi5$eo_groups,
                          adm.cens = 5)




sdi5_gnd$chi2gw
sdi5_gnd$pvalgw

# to evaluate re-calibration we will plot the calibration plot
# without adjustment

#===== create calibration plot with confidence intervals #====

df_sdi5$pred_5y_g = as.numeric(cut2(df_sdi5$pred_5y, g = 10))



sdi5_gnd_orig = GND.calib(pred = df_sdi5$pred_5y,
                          tvar = df_sdi5$hf_admit_comp_years_5y,
                          out = df_sdi5$event_cox_5y,
                          cens.t = 5,
                          groups = df_sdi5$pred_5y_g,
                          adm.cens = 5)

sdi5_gnd_orig$chi2gw

tiff("figures/recalib5.tiff", height = 5, width = 5, 
     units = "in",res = 600)

plot(x = sdi5_gnd$expectedperc,
     y = sdi5_gnd$kmperc, cex = 1, frame = F,
     col = "red",pch = 16,
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH (From the WATCH-DM Score)",
     ylab = "Observed 5-year HFH", main = "Quintile 5")


lines(x = sdi5_gnd$expectedperc,
      y = sdi5_gnd$kmperc, lty = 2,
      col = "red")

abline(0,1,col = "blue", lty = 2)


points(x = sdi5_gnd_orig$expectedperc,
       y = sdi5_gnd_orig$kmperc, cex = 1, pch = 15, 
       col = "black")

lines(x = sdi5_gnd_orig$expectedperc,
      y = sdi5_gnd_orig$kmperc, lty = 2,
      col = "black")


# legend("topleft",pch = c(1,3), col = c("blue","black"),
#        legend = c("re-calibrated WATCH-DM score",
#                   "original WATCH-DM score"))


dev.off()

#=== check watch dm score distribution according to quintiles of the SDI

df %>% drop_na(sdi_cat) %>% group_by(sdi_cat) %>% summarise(
  mean = mean(watch_dm_points),
  median = median(watch_dm_points), 
  sd = sd(watch_dm_points),
  q1 = quantile(watch_dm_points, 0.25),
  q3 = quantile(watch_dm_points, 0.75)
)

# plot watch dm according to quintiles of the SDI 



d = df %>% drop_na(sdi_cat) %>% ggplot(aes(x = watch_dm_points, fill = factor(sdi_cat))) + 
  geom_bar() + facet_wrap(~sdi_cat) + xlab("WATCH-DM Score") + ylab( "") + xlim(5,30) + 
  theme_ckb() + scale_fill_grey() + labs(fill = "SDI Quintile")

ggsave(plot = d,filename = "figures/sdi_dist.tiff", height = 5, width = 8, units = "in",device = "tiff",dpi = 600)

#= get max and min watchdm

min(df$watch_dm_points); max(df$watch_dm_points)

# we need to get the watch dm values for each sdi with 2 decimals 


tapply(df$watch_dm_points,df$sdi_cat, mean)




#== check whether the mean decrease in watch dm is statistically signficant 

df$sdi_cat = factor(df$sdi_cat)

anova = aov(watch_dm_points ~ sdi_cat,data = df)

anova

TukeyHSD(anova,ordered = T)

#=== get the calibration plot and GND test for the whole cohort 

# model and then predicted scores 

glimpse(df)

df$pred5y = 1 - pec::predictSurvProb(object = model,
                                     times = 5, newdata = df)


summary(df$pred5y)

df$pred5y_g = as.numeric(cut2(df$pred5y, g = 10))

source("scripts/gnd.R")

df_gnd = GND.calib(pred = df$pred5y,
                          tvar = df$hf_admit_comp_years_5y,
                          out = df$event_cox_5y,
                          cens.t = 5,
                          groups = df$pred5y_g,
                          adm.cens = 5)


tiff("figures/calib_plot_w.tiff", height = 6, width = 6,units = "in",res = 600)

plot(x = df_gnd$expectedperc,
     y = df_gnd$kmperc,pch = 15, cex = 1, frame = F,
     col = "blue",
     xlim = c(0, 0.12),
     ylim = c(0, 0.12),
     xlab = "Predicted 5-year risk of HFH",
     ylab = "Observed 5-year HFH", main = "Whole Cohort")
abline(0,1,col = "red", lty = 2)
lines(x = df_gnd$expectedperc,
      y = df_gnd$kmperc, lty = 2, col = "black")

dev.off()



#== get the plot for the cumulative incidence of HFH as a cox model 

# surv object 

surv_obj = Surv(df$hf_admit_comp_years_5y,
                         df$event_cox_5y)

inc = survfit(surv_obj ~ 1, data = df)


hf_km = cuminc_plot_group_comb(surv.obj = inc,.bgcolor = "white",.ylab = "Cumulative HFH",.xlab = "Years",
                       .plotsize = 4,.risktablesize = 1.5,.legend_labs = F,.xlim = 5,.ylim = 0.08,
                       .break.x.by = 1,
                       .break.y.by = 0.02,.legend = F,.group = F)


ggsave(plot = hf_km,filename = "figures/rev_hf_km.tiff",height = 8, width = 8, 
       units = "in",device = "tiff",dpi = 600)


# HFH by groups 

inc_g = survfit(surv_obj ~ orig_cat, data = df)


hf_km_groups = cuminc_plot_group_comb(surv.obj = inc_g,
            .bgcolor = "white",.ylab = "Cumulative HFH",.xlab = "Years",
    .plotsize = 4,.risktablesize = 2,.legend_labs = c("I","II","III","IV","V"),.xlim = 5,.ylim = 0.12,
                               .break.x.by = 1,
                               .break.y.by = 0.02,.legend = T,
    .group = T)



ggsave(plot = hf_km_groups,
       filename = "figures/rev_hf_km_groups.tiff",height = 8, width = 8, 
       units = "in",device = "tiff",dpi = 600)



#== we need to create a table of predicted risk using the re-calibrated model
#== for each WATCH-DM group in each quintile of SDI .

# Q1

tapply(df_sdi1$eo_pred, df_sdi1$orig_cat, mean)

tapply(df_sdi2$eo_pred, df_sdi2$orig_cat, mean)

tapply(df_sdi3$eo_pred, df_sdi3$orig_cat, mean)

tapply(df_sdi4$eo_pred, df_sdi4$orig_cat, mean)

tapply(df_sdi5$eo_pred, df_sdi5$orig_cat, mean)


#= complete by doing DCA

library(dcurves)

# dca for the sdi 1

length(df_sdi1$eo_pred)
length(df_sdi1$pred_5y)

library(dcurves)

dca1 = dcurves::dca(formula = Surv(hf_admit_comp_years_5y,
                                   event_cox_5y) ~ pred_5y + eo_pred,
                    data = df_sdi1,thresholds = seq(0,0.15,0.1),time = 5)


plot(dca1)

df_sdi1$event_cox_5y = as.numeric(df_sdi1$event_cox_5y) - 1

df_sdi1$event_cox_5y

df_sdi1$e = as.numeric(df_sdi1$event_cox_5y)

df_sdi1 = as.data.frame(df_sdi1)

a = stdca(data = df_sdi1,
      outcome = 'e',ttoutcome = 'hf_admit_comp_years_5y',
      timepoint = 5,
      predictors = c("pred_5y","eo_pred"),
      xstart = 0,xstop = 0.15,
      xby = 0.01,
      ymin = -0.01)

b = a$net.benefit

tiff("figures/dca1.tiff", height = 6, width = 8, units = "in",res = 600)


plot(x = b$threshold, y = b$all, type = "l", frame = F, 
     ylim = c(-0.01, 0.05),
     ylab = "Net Benefit", xlab = "Risk Threshold" )
lines(x = b$threshold, y = b$pred_5y, col = "red")
lines(x = b$threshold, y = b$eo_pred, col = "blue")
abline(h = 0, col = "black")

legend("topright", fill = c("red","blue"), 
       legend = c("WATCH-DM Score", "Re-calibrated WATCH-DM Score"))

dev.off()

#= dca for sdi 2

df_sdi2$e = as.numeric(df_sdi2$event_cox_5y) 

df_sdi2 = as.data.frame(df_sdi2)

a2 = stdca(data = df_sdi2,
          outcome = 'e',
          ttoutcome = 'hf_admit_comp_years_5y',
          timepoint = 5,
          predictors = c("pred_5y","eo_pred"),
          xstart = 0,xstop = 0.15,
          xby = 0.01,
          ymin = -0.01)

b2 = a2$net.benefit


tiff("figures/dca2.tiff", height = 6, width = 8, units = "in",res = 600)


plot(x = b2$threshold, y = b2$all, type = "l", frame = F, 
     ylim = c(-0.01, 0.05),
     ylab = "Net Benefit", xlab = "Risk Threshold" )
lines(x = b2$threshold, y = b2$pred_5y, col = "red")
lines(x = b2$threshold, y = b2$eo_pred, col = "blue")
abline(h = 0, col = "black")

legend("topright", fill = c("red","blue"), 
       legend = c("WATCH-DM Score", "Re-calibrated WATCH-DM Score"))

dev.off()

#== sdi 3


df_sdi3 = as.data.frame(df_sdi3)

df_sdi3$pred_5y

a3 = stdca(data = df_sdi3,
          outcome = 'event_cox_5y',
          ttoutcome = 'hf_admit_comp_years_5y',
          timepoint = 5,
          predictors = c("pred_5y","eo_pred"),
          xstart = 0,xstop = 0.15,
          xby = 0.01,
          ymin = -0.01)

b3 = a3$net.benefit

tiff("figures/dca3.tiff", height = 6, width = 8, units = "in",res = 600)


plot(x = b3$threshold, y = b3$all, type = "l", frame = F, 
     ylim = c(-0.01, 0.05),
     ylab = "Net Benefit", xlab = "Risk Threshold" )
lines(x = b3$threshold, y = b3$pred_5y, col = "red")
lines(x = b3$threshold, y = b3$eo_pred, col = "blue")
abline(h = 0, col = "black")

legend("topright", fill = c("red","blue"), 
       legend = c("WATCH-DM Score", "Re-calibrated WATCH-DM Score"))

dev.off()


#= dca quintile 4 


df_sdi4 = as.data.frame(df_sdi4)

df_sdi4$pred_5y

a4 = stdca(data = df_sdi4,
           outcome = 'event_cox_5y',
           ttoutcome = 'hf_admit_comp_years_5y',
           timepoint = 5,
           predictors = c("pred_5y","eo_pred"),
           xstart = 0,xstop = 0.15,
           xby = 0.01,
           ymin = -0.01)

b4 = a4$net.benefit

tiff("figures/dca4.tiff", height = 6, width = 8, units = "in",res = 600)


plot(x = b4$threshold, y = b4$all, type = "l", frame = F, 
     ylim = c(-0.01, 0.05),
     ylab = "Net Benefit", xlab = "Risk Threshold" )
lines(x = b4$threshold, y = b4$pred_5y, col = "red")
lines(x = b4$threshold, y = b4$eo_pred, col = "blue")
abline(h = 0, col = "black")

legend("topright", fill = c("red","blue"), 
       legend = c("WATCH-DM Score", "Re-calibrated WATCH-DM Score"))

dev.off()

#= sdi quintile 5 



df_sdi5 = as.data.frame(df_sdi5)

df_sdi5$pred_5y

a5 = stdca(data = df_sdi5,
           outcome = 'event_cox_5y',
           ttoutcome = 'hf_admit_comp_years_5y',
           timepoint = 5,
           predictors = c("pred_5y","eo_pred"),
           xstart = 0,xstop = 0.20,
           xby = 0.01,
           ymin = -0.01)

b5 = a5$net.benefit

tiff("figures/dca5.tiff", height = 6, width = 8, units = "in",res = 600)


plot(x = b5$threshold, y = b5$all, type = "l", frame = F, 
     ylim = c(-0.01, 0.05),
     ylab = "Net Benefit", xlab = "Risk Threshold" )
lines(x = b5$threshold, y = b5$pred_5y, col = "red")
lines(x = b5$threshold, y = b5$eo_pred, col = "blue")
abline(h = 0, col = "black")

legend("topright", fill = c("red","blue"), 
       legend = c("WATCH-DM Score", "Re-calibrated WATCH-DM Score"))

dev.off()


#= making the for table for the paper for baseline characteristics
#= according to the SDI

df$female = with(df, ifelse( sex == "F", 1, 0))

fvars = c("p_anemia","p_cad","p_dialysis","p_pvd",'p_liver',
         "p_htn","p_cancer","p_lung","p_cva","race.mod","hispanic","female")


cvars = c("age_w", 
          
          "systolic_wi",  
          "hba1c_wi",  "diastolic_wi", 
          "creat_wi",  "hdl_wi",  "watch_dm_points")



all_vars = c(fvars, cvars)

t1_sdi = tableone::CreateCatTable(vars = vars, 
                              strata = "sdi_cat", data = df,
                              addOverall = F)

write.csv(print(t1_sdi),"results/t1_sdi.csv")


cvars = c("age_w", 
          
          "systolic_wi",  
          "hba1c_wi",  "diastolic_wi", 
          "creat_wi",  "hdl_wi",  "watch_dm_points")


t2_sdi = tableone::CreateContTable(vars = cvars, data = df,
                               strata = "sdi_cat", 
                               addOverall = F)

write.csv(print(t2_sdi),"results/t2_sdi.csv")

# get race, sex, ethnicity for groups



vars = c("race.mod","hispanic","female")


t3_sdi = tableone::CreateCatTable(vars = vars, strata = "sdi_cat",
                              data = df, addOverall = T)

write.csv(print(t3_sdi),"results/t3_sdi.csv")

# chronic kidney disease 

tabyl(df$ckd)

t_ckd_sdi = tableone::CreateCatTable(vars = "ckd", strata = "sdi_cat",
                                 data = df, addOverall = T)

write.csv(print(t_ckd_sdi),"results/ckd_t_sdi.csv")


#== also to get the revised table 1 using the original cat for WATCH-DM score

df$female = with(df, ifelse( sex == "F", 1, 0))


fvars = c("p_anemia","p_cad","p_dialysis","p_pvd",'p_liver',
         "p_htn","p_cancer","p_lung","p_cva","race.mod","hispanic","female")


cvars = c("age_w", 
          
          "systolic_wi",  
          "hba1c_wi",  "diastolic_wi", 
          "creat_wi",  "hdl_wi",  "watch_dm_points")


all_vars = c(fvars, cvars)

all_vars

t1_mod = tableone::CreateTableOne(vars = all_vars,factorVars = fvars, 
                                  strata = "orig_cat", data = df,
                                  addOverall = T)



write.csv(print(t1_mod), "results/mod_table1.csv")


#= going to see the CIF for 5 year predicted HFH in SDI1 according to race

tabyl(df_sdi1$race.mod)

#== see the HFH according to race 

with(df_sdi1, tapply(pred_5y, race.mod, mean))

#= see the observed 5 year HFH according to the race 

sdi1_race = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ race.mod, data = df_sdi1)

race1 = summary(sdi1_race, times = 5)

t1 = cumu_inc_summary(sum_obj = race1,2,time.points = 1)

# this shows that after adjusting for the SDI, race does not matter much

#= for the SDI 5

sdi5_race = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ race.mod, 
              data = df_sdi5)

race5 = summary(sdi5_race, times = 5)

t2 = cumu_inc_summary(sum_obj = race5, deci = 2,time.points = 1)

r1 = t1 %>% select(group, `cumulative_incidence(%)`, `lower95(%)`, `upper95(%)`) %>% 
  mutate(SDI = 1)
r2 = t2 %>% select(group, `cumulative_incidence(%)`, `lower95(%)`, `upper95(%)`) %>% 
  mutate(SDI = 5)

r = rbind(r1, r2)

r$`cumulative_incidence(%)` = as.numeric(r$`cumulative_incidence(%)`)
r$`lower95(%)` = as.numeric(r$`lower95(%)`)
r$`upper95(%)` = as.numeric(r$`upper95(%)`)

#= save this to make a bar plot

write_csv(r, "results/race_inc_sdi.csv")