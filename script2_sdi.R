####### External validation of the WATCH-DM score. 
####### The social determinants of health.

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
#=== to get the descriptive analysis for the SDI groups.
#=== then to do the rest of the analysis for the SDI groups.

#= get the source files and then get the dataset.
#=


setwd("P:\\ORD_Deo_202205002D\\DM\\DM_2010\\")

# get many functions now.

source('ckbplotr_theme.r')
source('cuminc_plot_function_new.r')
source('survplot_fn.r')
source('cif_plot.r')
source('code_ci.R')

#= get back to wd

setwd("P:\\ORD_Deo_202205002D\\DM\\watch_dm_validate\\")

df = read_csv("data/df_use.csv")
dim(df)
#= 1065691

df = df[!is.na(df$sdi_cat), ]; dim(df)

#= 1012315 patients for SDI 

#= see the dist of the categories of SDI

tabyl(df$sdi_cat)

df %>% tabyl(sdi_cat, orig_cat)

map_dbl(df$sdi_cat,.x = df$watch_dm_points,.f = mean)

df %>% select(watch_dm_points, sdi_cat) %>% 
  split(.$sdi_cat) %>% 
    map( summary)

# `1`
# watch_dm_points    sdi_cat 
# Min.   : 2.00   Min.   :1  
# 1st Qu.:13.00   1st Qu.:1  
# Median :15.00   Median :1  
# Mean   :15.17   Mean   :1  
# 3rd Qu.:17.00   3rd Qu.:1  
# Max.   :32.00   Max.   :1  
# 
# $`2`
# watch_dm_points    sdi_cat 
# Min.   : 2.00   Min.   :2  
# 1st Qu.:13.00   1st Qu.:2  
# Median :15.00   Median :2  
# Mean   :15.05   Mean   :2  
# 3rd Qu.:17.00   3rd Qu.:2  
# Max.   :33.00   Max.   :2  
# 
# $`3`
# watch_dm_points    sdi_cat 
# Min.   : 2.00   Min.   :3  
# 1st Qu.:12.00   1st Qu.:3  
# Median :15.00   Median :3  
# Mean   :14.94   Mean   :3  
# 3rd Qu.:17.00   3rd Qu.:3  
# Max.   :33.00   Max.   :3  
# 
# $`4`
# watch_dm_points    sdi_cat 
# Min.   : 2.00   Min.   :4  
# 1st Qu.:12.00   1st Qu.:4  
# Median :15.00   Median :4  
# Mean   :14.77   Mean   :4  
# 3rd Qu.:17.00   3rd Qu.:4  
# Max.   :33.00   Max.   :4  
# 
# $`5`
# watch_dm_points    sdi_cat 
# Min.   : 2.00   Min.   :5  
# 1st Qu.:12.00   1st Qu.:5  
# Median :14.00   Median :5  
# Mean   :14.45   Mean   :5  
# 3rd Qu.:17.00   3rd Qu.:5  
# Max.   :32.00   Max.   :5  


