# script for rev 2.
# do the cross validation using the rms library for the bootstrapped validation.


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

# fit model using rms package and then crossvalidate.

df2 = df %>% select(hf_admit_comp_years_5y,event_cox_5y,watch_dm_points,sdi_cat)

dd <- datadist(df2)
options(datadist = "dd")


# fit the model using both score + sdi 


model = cph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points + sdi_cat ,
 x = T, y = T, data = df2, surv = T)


v <- validate(model,B = 100)


v

# Brier score to demonstrate improvement with sdi 

df2 = data.frame(df2)

df3 = df2 %>% drop_na(sdi_cat)

model = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points,
              data = df3, x = T)


model2 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points + sdi_cat ,
               data = df3, x = T)

# get predicted risk for each model.

df3$m1 = predictRisk(object = model,times = 5,newdata = df3)

df3$m2 = predictRisk(object = model2, times = 5,newdata = df3)


summary(df2$m1);summary(df2$m2)

p = riskRegression::Score(list("model1" = df3$m1, "model2" = df3$m2), 
          formula = Surv(hf_admit_comp_years_5y,event_cox_5y)~1,
          data = df3,times = 5,metrics = "brier")

p


