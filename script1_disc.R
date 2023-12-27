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


# now to get the data 


setwd("P:\\ORD_Deo_202205002D\\DM\\DM_2010\\")

# get many functions now.

source('ckbplotr_theme.r')
source('cuminc_plot_function_new.r')
source('survplot_fn.r')
source('cif_plot.r')
source('code_ci.R')


df = read_csv("data/df_watchdm_5y.csv")
glimpse(df)

# get the zip code and join 

z = read_sas("SAS_DATA/DM2010_ZIP.sas7bdat")

# remove duplicates of scrssn in z and then join 

names(z) = tolower(names(z))


z2 = z %>% select(scrssn,zip_1)

z2 = distinct(z2, scrssn, .keep_all = T)

df2 = left_join(df, z2, by = "scrssn")

summary(df2$zip_1)

df2$zip_1 = as.numeric(df2$zip_1)



# bring back the wd to this project again.

setwd("P:\\ORD_Deo_202205002D\\DM\\watch_dm_validate\\")

# now to get the CDI 

cdi = read_excel("data\\SVI.xlsx")

glimpse(cdi)

df2 = df2 %>% rename(ZIP = zip_1)

df3 = left_join(df2, cdi, by = "ZIP")

summary(df3$sdi_score)

df3$sdi_cat = with(df3, 
                   ifelse(
                     sdi_score < 20, 1, 
                     ifelse(
                       sdi_score >= 20 & sdi_score < 40, 2, 
                       ifelse(
                         sdi_score >= 40 & sdi_score < 60, 3, 
                         ifelse(
                           sdi_score >= 60 & sdi_score < 80, 4, 5
                         )
                       )
                     )
                   ))


# convert df3 to df to do the rest of the analysis

df = df3

# going to get the original categories are reported in the WATCH-DM paper
# < 11, 12 - 13, 14 - 15, 16 - 18, >= 19 



df$orig_cat = with(df, 
                   ifelse(
                     watch_dm_points <= 11, 1, 
                     ifelse(watch_dm_points %in% c(12,13), 2,
                            ifelse(watch_dm_points %in% c(14,15), 3, 
                                   ifelse(watch_dm_points %in% c(16, 17, 18), 4, 5)))
                   ))


df %>% tabyl(orig_cat)


# going to make race and ethnicity into proper var.

tabyl(df$race)

df$race.mod = with(df, ifelse(race  %in% c("WHITE"), "white",
                              ifelse(race %in% c("BLACK OR AFRICAN AMERICAN"), "black", "others")))


df %>% tabyl(race.mod)

# declined to answer 3% , moved to white as mode category


tabyl(df$ethnicity)

df$hispanic = with(df, ifelse(
  ethnicity == "HISPANIC OR LATINO", 1,0
))

tabyl(df$hispanic)



# those NA or not reported, imputed as not hispanic 3.3%

dim(df)

summary(df$age)

# > summary(df$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.00   61.00   66.00   67.45   75.00  108.00

sd(df$age)

# > sd(df$age)
# [1] 9.938348

tabyl(df$race.mod)

# > tabyl(df$race.mod)
# df$race.mod      n    percent
# black 173491 0.16279672
# others  87010 0.08164656
# white 805190 0.75555672

tabyl(df$hispanic)

# df$hispanic      n   percent
# 0 998944 0.9373674
# 1  66747 0.0626326

tabyl(df$sex)


# df$sex       n    percent
# F   35432 0.03324791
# M 1030259 0.96675209


# making table 1 in parts because of the large N 
# making variables to get the table 1

df = df %>% mutate(
  
  p_cancer = ifelse(
    is.na(cancer), 0, 1),
  
  p_anemia = if_else(
    is.na(anemia), 0,1),
  
  p_cad = if_else(
    is.na(cad), 0, 1),
  
  p_htn = if_else(
    is.na(htn),0,1
  ),
  
  p_liver = if_else(
    is.na(liver), 0, 1
  ),
  
  p_pvd = if_else(
    is.na(pvd),0,1
  ),
  
  p_dialysis = if_else(
    is.na(kiddial), 0, 1
  ),
  
  p_lung = if_else(
    is.na(lung), 0,1
  ),
  
  p_cva = if_else(
    is.na(cva), 0, 1
  )
  
)

# 


vars = c("p_anemia","p_cad","p_dialysis","p_pvd",'p_liver',
         "p_htn","p_cancer","p_lung","p_cva")



t1 = tableone::CreateCatTable(vars = vars, 
                              strata = "orig_cat", data = df,
                              addOverall = T)

write.csv(print(t1),"results/t1.csv")


cvars = c("age_w", 
          
          "systolic_wi",  
          "hba1c_wi",  "diastolic_wi", 
          "creat_wi",  "hdl_wi",  "watch_dm_points")


t2 = tableone::CreateContTable(vars = cvars, data = df,
                               strata = "orig_cat", 
                               addOverall = T)

write.csv(print(t2),"results/t2.csv")

# get race, sex, ethnicity for groups


df$female = with(df, ifelse( sex == "F", 1, 0))

vars = c("race.mod","hispanic","female")


t3 = tableone::CreateCatTable(vars = vars, strata = "orig_cat",
                              data = df, addOverall = T)

write.csv(print(t3),"results/t3.csv")

# chronic kidney disease 

tabyl(df$ckd)

t_ckd = tableone::CreateCatTable(vars = "ckd", strata = "orig_cat",
                                 data = df, addOverall = T)

write.csv(print(t_ckd),"results/ckd_t.csv")

# plotting the WATCH-DM point distribution according to the categories


df$orig_cat = factor(df$orig_cat)

h = df %>% ggplot(aes(x = watch_dm_points, fill = orig_cat)) + 
  geom_bar() + ggtitle("Distribution of the WATCH-DM score") + 
  xlab("WATCH-DM score")+ ylab("Number of Patients") + theme_ckb()


h


ggsave(h, filename = "figures/hist.tiff", 
       height = 6, width = 5, units = "in", device = "tiff",dpi = 600)



########################### Without Competing risks ############################

# If we consider that although mortality is a competing outcome, covariates that influence the incidence of HF also influence 
# the death in the same direction.
# so we willl do the analysis without considering competing risks situation for now.


# get the HF incidence at 5 years 

glimpse(df)

# get event indicator as only 0/1

df$event_cox_5y = with(df, ifelse(hf_admit_comp_5y == 1, 1, 0))

# save this dataset so that we do not need to repeat these steps 

write_csv(df, "data/df_use.csv")


# now to get the 5 yr incidence for each group

inc = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_cat, data = df)

cumu_inc_summary(summary(inc, times = 5), deci = 5,time.points = 1)


# # A tibble: 5 x 9
# group       time total_patients  nrisk nevent ncensor `cumulative_incidence(%)` `lower95(%)` `upper95(%)`
# <fct>      <dbl>          <int>  <dbl>  <dbl>   <dbl> <chr>                     <chr>        <chr>       
#   1 orig_cat=1     5         192113 168204   5820  186293 " 3.1893017"              " 3.1085907" " 3.2699455"
# 2 orig_cat=2     5         194214 165017   7436  186778 " 4.0650995"              " 3.9744604" " 4.1556530"
# 3 orig_cat=3     5         236129 188827  11075  225054 " 5.0926038"              " 4.9999736" " 5.1851437"
# 4 orig_cat=4     5         277903 206019  16995  260908 " 6.8011659"              " 6.7020393" " 6.9001872"
# 5 orig_cat=5     5         165332 108386  15881  149451 "11.0030524"              "10.8403046" "11.1655031"



hf_plot_cat = cuminc_plot_group_comb(surv.obj = inc,.bgcolor = "gray90",.ylab = "Cumulative Mortality",.xlab = "Years",.plotsize = 4,
                                  .risktablesize = 1,.legend_labs = c("I","II","III","IV","V"),
                                  .xlim = 5,.ylim = 0.15,
                                  .break.x.by = 1,.break.y.by = 0.05,
                                  .legend = T,.group = T)

hf_plot_cat


ggsave(plot = hf_plot_cat,
               filename = "figures/m_plot_cat.tiff",
               height = 10, width = 16, units = "in", device = "tiff", dpi = 600)


# validation metrics 
# c-statistics Harrell
# uno c-statistic

# calibration 
# calibration plot
# ICI/E50/E90
# GND test 

library(timeROC)



h_c = concordance(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points, 
                  data = df,reverse = T)

h_c

h_c$concordance + 1.96*0.001882

h_c$concordance - 1.96*0.001882

# > h_c
# Call:
#   concordance.formula(object = Surv(hf_admit_comp_years_5y, event_cox_5y) ~ 
#                         watch_dm_points, data = df, reverse = T)
# 
# n= 1065691 
# Concordance= 0.6221 se= 0.001182
# concordant  discordant      tied.x      tied.y     tied.xy 
# 31971939511 18642154766  3949282541      841387       66989 
# > 
#   > h_c$concordance + 1.96*0.001882
# [1] 0.6258383
# > 
#   > h_c$concordance - 1.96*0.001882
# [1] 0.6184608


# get the predicted value for HF incidence at 5 years
# use that to calculate the O/E ratio.

m = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points, data = df, x = T, y = T)

df$pred_5y = 1 - pec::predictSurvProb(m, newdata = df, times = 5)

summary(df$pred_5y)

# get the mean predicted value for each category

df %>% group_by(orig_cat) %>% 
  summarise(
    mean = mean(pred_5y),
    median = median(pred_5y),
    sd = sd(pred_5y),
    lq = quantile(pred_5y, 0.25),
    uq = quantile(pred_5y, 0.75)
  )

# 

# # A tibble: 5 x 6
# orig_cat   mean   median      sd     lq     uq
# <dbl>       <dbl>  <dbl>   <dbl>  <dbl>  <dbl>
# 1        1 0.0301 0.0310 0.00480 0.0276 0.0348
# 2        2 0.0417 0.0439 0.00240 0.0391 0.0439
# 3        3 0.0523 0.0554 0.00302 0.0493 0.0554
# 4        4 0.0689 0.0697 0.00640 0.0621 0.0781
# 5        5 0.108  0.0979 0.0261  0.0875 0.123 

# now to calculate the O/E ratio

horizon = 5

obj = summary(survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df), data = df, times = horizon)

obj

OE <- (1 - obj$surv)/mean(df$pred_5y)

alpha <- 0.05

oe_ratio_res = list("OE" = OE,
"2.5%" = OE*exp(-qnorm(1 - alpha/2)*sqrt(1/obj$n.event)), 
"97.5%" = OE*exp(qnorm(1 - alpha/2)*sqrt(1/obj$n.event)))

oe_ratio_res

# 
# $OE
# [1] 0.9854872
# 
# $`2.5%`
# [1] 0.9774446
# 
# $`97.5%`
# [1] 0.993596
# 

# Calibration according to the heirarchy - 

# 1. calibration in the large (weak) = O/E ratio 

# 2. moderate calibration - calibration in deciles, GND test,
# calibration plot 
# GND test 

glimpse(df)

# split the pred prob into deciles

df$dec = as.numeric(cut2(df$pred_5y, g = 10))

table(df$dec, df$event_cox_5y)


# fit the test 

GND = GND.calib(pred = df$pred_5y, tvar = df$hf_admit_comp_years_5y,out = df$event_cox_5y,
                cens.t = 5, groups = df$dec, adm.cens = 5)


# GND test result

# > GND = GND.calib(pred = df$pred_5y, tvar = df$hf_admit_comp_years_5y,out = df$event_cox_5y,
#                   +                 cens.t = 5, groups = df$dec, adm.cens = 5)
# group totaln  censn numevents kmnum expected  kmperc expectedperc     kmvar
# 1      1 120397 601985      3465  3644     3287 0.03026      0.02730 2.567e-07
# 2      2 161193 805965      5534  5842     5995 0.03624      0.03719 2.292e-07
# 3      3 104737 523685      4257  4539     4600 0.04334      0.04392 4.232e-07
# 4      4 117900 589500      5338  5765     5815 0.04890      0.04932 4.275e-07
# 5      5 118229 591145      5737  6262     6545 0.05297      0.05536 4.653e-07
# 6      6 110156 550780      6062  6700     6842 0.06082      0.06211 5.766e-07
# 7      7  93634 468170      5767  6421     6522 0.06858      0.06966 7.653e-07
# 8      8  74113 370565      5166  5792     5787 0.07815      0.07808 1.100e-06
# 9      9  95346 476730      7978  9056     8761 0.09498      0.09188 1.035e-06
# 10    10  69986 349930      7903  9175     9106 0.13110      0.13012 1.926e-06
# GND_component
# 1      34.276974
# 2       3.936062
# 3       0.801840
# 4       0.413577
# 5      12.269086
# 6       2.896092
# 7       1.520052
# 8       0.003944
# 9       9.254433
# 10      0.498050
# > GND
# df       chi2gw       pvalgw 
# 9.000000e+00 6.587011e+01 9.768075e-11 

# Therneau et al 

# df$p = log(predict(m, newdata = df, type = "expected"))
# 
# df$lp = predict(m, newdata = df, type = "lp")
# 
# summary(df$lp)
# 
# df$logbase = df$p - df$lp
# 
# fit1 = glm(event_cox_5y ~ offset(p), 
#            family = "poisson", data = df)
# 
# 
# fit1
# 
# exp(fit1$coefficients)


# Calibration plot 

sc = riskRegression::Score(object = list("score" = df$pred_5y),
     formula = Hist(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                           data = df,
                           summary = "IPA",
                           plots = "calibration",
                           cause = 1,
                           times = 5,
                           cens.method = "jackknife",
                           cens.model = "km")

tiff("figures/calib_whole.tiff", height = 6, width = 8, 
     units = "in",res = 600)
plotCalibration(sc, bar = T, ylim = c(0, 0.14),
                xlab = " ", yaxt = "n")
dev.off()

################ SDI #################
# look at the HF incidence per category in SDI 1

df = read_csv("data/df_use.csv")

tabyl(df$sdi_cat)

# fit the model to the whole data 

m = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points, 
          data = df, x = T, y = T)


# for the models of SDI, remove those with missing information 

df_n = df[!is.na(df$sdi_cat), ]; dim(df_n)



# see the incidence of HF 
# SD1
df_sdi1 = df_n %>% filter(sdi_cat == 1)

# get the linear predictor 

df_sdi1$lp = predict(m, newdata = df_sdi1)

summary(df_sdi1$lp)

dim(df_sdi1) # sdi1 is the least deprived.

inc_sdi1 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_cat,
              data = df_sdi1)

cumu_inc_summary(summary(inc_sdi1, times = 5), deci = 5,time.points = 1)


# SD2
df_sdi2 = df_n %>% filter(sdi_cat == 2)

dim(df_sdi2) 

inc_sdi2 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_cat,
                   data = df_sdi2)

cumu_inc_summary(summary(inc_sdi2, times = 5), 
                 deci = 5,time.points = 1)


#SD3
df_sdi3 = df_n %>% filter(sdi_cat == 3)

dim(df_sdi3) 

inc_sdi3 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_cat,
                   data = df_sdi3)

cumu_inc_summary(summary(inc_sdi3, times = 5), 
                 deci = 5,time.points = 1)

#SD4
df_sdi4 = df_n %>% filter(sdi_cat == 4)

dim(df_sdi4) 

inc_sdi4 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_cat,
                   data = df_sdi4)

cumu_inc_summary(summary(inc_sdi4, times = 5), 
                 deci = 5,time.points = 1)


# see the incidence of HF 

df_sdi5 = df_n %>% filter(sdi_cat == 5)

inc_sdi5 = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ orig_cat,
                   data = df_sdi5)

cumu_inc_summary(summary(inc_sdi5, times = 5), deci = 5,time.points = 1)

# see c-statistic for each quintile of the SDI

library(timeROC)

c_sdi1 = concordance(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points, 
                  data = df_sdi1,
                  reverse = T)

c_sdi1

c_sdi1$concordance + 1.96*0.003704

c_sdi1$concordance - 1.96*0.001882






c_sdi5 = concordance(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points, 
                     data = df_sdi5,
                     reverse = T)

c_sdi5

c_sdi5$concordance + 1.96*0.002

c_sdi1$concordance - 1.96*0.002



df_w = df %>% filter(race.mod  == 'black')

glimpse(df_w)

# look at the event rate per categories in whites


inc_w = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y)
                ~ orig_cat, data = df_w)

cumu_inc_summary(summary(inc_w, times = 5), deci = 5,time.points = 1)




library(timeROC)



h_c = concordance(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ watch_dm_points, 
                  data = df_,
                  reverse = T)

h_c

h_c$concordance + 1.96*0.001882

h_c$concordance - 1.96*0.001882

# calibration // this is the more important part of the study.

# do calibration according to the hierarchy developed by van Calster etal 
# calibration in the large 
# fit model to the whole cohort 

m = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 
                watch_dm_points, data = df, x = T, y = T)

# use this model object to predict 5-year outcome in the SDI1 group.

df_sdi1$pred_5y = 1 - pec::predictSurvProb(m,
                                           newdata = df_sdi1, times = 5)


horizon = 5

obj = summary(survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                      data = df_sdi1), data = df_sdi1, times = horizon)

obj

OE <- (1 - obj$surv)/mean(df_sdi1$pred_5y)

alpha <- 0.05

oe_ratio_res = list("OE" = OE,
                    "2.5%" = OE*exp(-qnorm(1 - alpha/2)*sqrt(1/obj$n.event)), 
                    "97.5%" = OE*exp(qnorm(1 - alpha/2)*sqrt(1/obj$n.event)))

oe_ratio_res

# weak calibration 
# calibration slope 

# we have the linear predictor 

hist(df_sdi1$lp)

cal_slope1 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ lp, data = df_sdi1)


alpha = 0.05

calib_slope_sdi1 = c(
'slope' = cal_slope1$coefficients,
 "2.5%" = cal_slope1$coefficients - qnorm(1 - alpha/2)*sqrt(cal_slope1$var),
"97.5%" = cal_slope1$coefficients + qnorm(1 - alpha/2)*sqrt(cal_slope1$var))


calib_slope_sdi1



# moderate calibration = GND test 
# split the prediction score into deciles 

df_sdi1$dec = as.numeric(cut2(df_sdi1$pred_5y, g = 10))

table(df_sdi1$dec, df_sdi1$event_cox_5y)

# fit the test 

source("scripts/gnd.R")


GND_sdi1 = GND.calib(pred = df_sdi1$pred_5y, 
                tvar = df_sdi1$hf_admit_comp_years_5y,
                out = df_sdi1$event_cox_5y,
                cens.t = 5, 
                groups = df_sdi1$dec, 
                adm.cens = 5)

GND_sdi1$chi2gw
GND_sdi1$pvalgw

# moderate calibration and calibration plot 

# get the observed HF incidence at 5 years from the 

# fit the model and then get the predicted risk at 5 years for HF

m = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 
            watch_dm_points, data = df, x = T, y = T)


df_sdi1$pred5y = 1 - pec::predictSurvProb(object = m, newdata = df_sdi1,times = 5)


df_sdi1_score = riskRegression::Score(object = list("score" = df_sdi1$pred5y),
      formula = Hist(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
      data = df_sdi1,
      summary = "IPA",
      plots = "calibration",
      cause = 1,
      times = 5,
      cens.method = "jackknife",
      cens.model = "km")

tiff("figures/calplot1.tiff", height = 5, 
     width = 5, units = "in",res = 600)
plotCalibration(df_sdi1_score, xlim = c(0, 0.20), 
                ylim = c(0, 0.20),
                auc.in.legend = F,
                brier.in.legend = F, col = 1)

dev.off()

# 
# # now to get the plot for calibration 
# 
# 
# # calibration plot 
# 
# 
# m_cph = rms::cph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 
#             watch_dm_points, data = df, x = T, y = T)
# 
# 
# dat_cal1 = cbind.data.frame(
#   
#   
#   "obs" = 1 - rms::survest(m_cph, 
#                            newdata = df_sdi1,
#                            times = 5)$surv,
#   "lower" = 1 - rms::survest(m_cph,
#                              newdata = df_sdi1,
#                              times = 5)$lower,
#   "upper" = 1 - rms::survest(m_cph,
#                              newdata = df_sdi1,
#                              times = 5)$upper,
#   "pred" = df_sdi1$pred_5y,
#   
#   "cll" = log(-log(1 - df_sdi1$pred_5y)),
#   "id" = 1:nrow(df_sdi1))
#   
#                   
# dat_cal1 = dat_cal1[order(dat_cal1$pred), ] 
# 
# plot(
#   x = dat_cal1$pred,
#   y = dat_cal1$obs,
#   type = "l",
#   lty = 1,
#   xlim = c(0, 0.25),
#   ylim = c(0,0.25)
# )
# abline(0,1,col = "red")
# 

# for SDI group 2

df_sdi2 = df_n %>% filter(sdi_cat == 2)

df_sdi2$pred_5y = 1 - pec::predictSurvProb(m,
                                           newdata = df_sdi2, times = 5)


horizon = 5

obj = summary(survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                      data = df_sdi1), data = df_sdi2, times = horizon)

obj

OE <- (1 - obj$surv)/mean(df_sdi2$pred_5y)

alpha <- 0.05

oe_ratio_res = list("OE" = OE,
                    "2.5%" = OE*exp(-qnorm(1 - alpha/2)*sqrt(1/obj$n.event)), 
                    "97.5%" = OE*exp(qnorm(1 - alpha/2)*sqrt(1/obj$n.event)))

oe_ratio_res

# weak calibration - calibration slope 

df_sdi2$lp = predict(m, newdata = df_sdi2)

cal_slope2 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ lp, data = df_sdi2)


alpha = 0.05

calib_slope_sdi2 = c(
  'slope' = cal_slope2$coefficients,
  "2.5%" = cal_slope2$coefficients - qnorm(1 - alpha/2)*sqrt(cal_slope2$var),
  "97.5%" = cal_slope2$coefficients + qnorm(1 - alpha/2)*sqrt(cal_slope2$var))


calib_slope_sdi2


# GND test 

df_sdi2$dec = as.numeric(cut2(df_sdi2$pred_5y, g = 10))

table(df_sdi2$dec, df_sdi2$event_cox_5y)

# fit the test 

source("scripts/gnd.R")


GND_sdi2 = GND.calib(pred = df_sdi2$pred_5y, 
                     tvar = df_sdi2$hf_admit_comp_years_5y,
                     out = df_sdi2$event_cox_5y,
                     cens.t = 5, 
                     groups = df_sdi2$dec, 
                     adm.cens = 5)

GND_sdi2$chi2gw
GND_sdi2$pvalgw

# calibration plot 


m = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 
            watch_dm_points, data = df, x = T, y = T)


df_sdi2$pred5y = 1 - pec::predictSurvProb(object = m, newdata = df_sdi2,times = 5)


df_sdi2_score = riskRegression::Score(object = list("score" = df_sdi2$pred5y),
                                      formula = Hist(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                                      data = df_sdi2,
                                      summary = "IPA",
                                      plots = "calibration",
                                      cause = 1,
                                      times = 5,
                                      cens.method = "jackknife",
                                      cens.model = "km")



tiff("figures/calplot2.tiff", height = 5, width = 5,
     units = "in", res = 600)
plotCalibration(df_sdi2_score, xlim = c(0, 0.2),
                ylim = c(0, 0.2),
                auc.in.legend = F,
                brier.in.legend = F, col = 2)
dev.off()




# for group 3

sd_sdi3 = df_n %>% filter(sdi_cat == 3)

df_sdi3$pred_5y = 1 - pec::predictSurvProb(m,
                                           newdata = df_sdi3, times = 5)


horizon = 5

obj = summary(survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                      data = df_sdi3), data = df_sdi3, times = horizon)

obj

OE <- (1 - obj$surv)/mean(df_sdi3$pred_5y)

alpha <- 0.05

oe_ratio_res = list("OE" = OE,
                    "2.5%" = OE*exp(-qnorm(1 - alpha/2)*sqrt(1/obj$n.event)), 
                    "97.5%" = OE*exp(qnorm(1 - alpha/2)*sqrt(1/obj$n.event)))

oe_ratio_res

# weak calibration - calibration slope 


df_sdi3$lp = predict(m, newdata = df_sdi3)

cal_slope3 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ lp, data = df_sdi3)

cal_slope3$coefficients

alpha = 0.05

calib_slope_sdi3 = c(
  'slope' = cal_slope3$coefficients,
  "2.5%" = cal_slope3$coefficients - qnorm(1 - alpha/2)*sqrt(cal_slope3$var),
  "97.5%" = cal_slope3$coefficients + qnorm(1 - alpha/2)*sqrt(cal_slope3$var))


calib_slope_sdi3

# moderate calibration and calibration plot 



df_sdi3$pred5y = 1 - pec::predictSurvProb(object = m, newdata = df_sdi3,times = 5)


df_sdi3_score = riskRegression::Score(object = list("score" = df_sdi3$pred5y),
                                      formula = Hist(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                                      data = df_sdi3,
                                      summary = "IPA",
                                      plots = "calibration",
                                      cause = 1,
                                      times = 5,
                                      cens.method = "jackknife",
                                      cens.model = "km")




plotCalibration(df_sdi3_score, xlim = c(0, 0.3), ylim = c(0, 0.3),
                auc.in.legend = F,brier.in.legend = F, col = 3, lwd = 1)



# GND test 

df_sdi3$dec = as.numeric(cut2(df_sdi3$pred_5y, g = 10))

table(df_sdi3$dec, df_sdi3$event_cox_5y)

# fit the test 

source("scripts/gnd.R")


GND_sdi3 = GND.calib(pred = df_sdi3$pred_5y, 
                     tvar = df_sdi3$hf_admit_comp_years_5y,
                     out = df_sdi3$event_cox_5y,
                     cens.t = 5, 
                     groups = df_sdi3$dec, 
                     adm.cens = 5)

GND_sdi3$chi2gw

GND_sdi3$pvalgw

# calibration plot 

# SD4


df_sdi4$pred_5y = 1 - pec::predictSurvProb(m,
                                           newdata = df_sdi4, times = 5)


horizon = 5

obj = summary(survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                      data = df_sdi4), data = df_sdi4, times = horizon)

obj

OE <- (1 - obj$surv)/mean(df_sdi4$pred_5y)

alpha <- 0.05

oe_ratio_res = list("OE" = OE,
                    "2.5%" = OE*exp(-qnorm(1 - alpha/2)*sqrt(1/obj$n.event)), 
                    "97.5%" = OE*exp(qnorm(1 - alpha/2)*sqrt(1/obj$n.event)))

oe_ratio_res

# weak calibration - calibration slope 


df_sdi4$lp = predict(m, newdata = df_sdi4)

cal_slope4 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ lp, data = df_sdi4)

cal_slope4$coefficients

alpha = 0.05

calib_slope_sdi4 = c(
  'slope' = cal_slope4$coefficients,
  "2.5%" = cal_slope4$coefficients - qnorm(1 - alpha/2)*sqrt(cal_slope4$var),
  "97.5%" = cal_slope4$coefficients + qnorm(1 - alpha/2)*sqrt(cal_slope4$var))


calib_slope_sdi4


# GND test 

df_sdi4$dec = as.numeric(cut2(df_sdi4$pred_5y, g = 10))

table(df_sdi4$dec, df_sdi4$event_cox_5y)

# fit the test 

source("scripts/gnd.R")


GND_sdi4 = GND.calib(pred = df_sdi4$pred_5y, 
                     tvar = df_sdi4$hf_admit_comp_years_5y,
                     out = df_sdi4$event_cox_5y,
                     cens.t = 5, 
                     groups = df_sdi4$dec, 
                     adm.cens = 5)

GND_sdi4$chi2gw

GND_sdi4$pvalgw

# calibration plot for sdi4



df_sdi4$pred5y = 1 - pec::predictSurvProb(object = m, newdata = df_sdi4,times = 5)


df_sdi4_score = riskRegression::Score(object = list("score" = df_sdi4$pred5y),
                                      formula = Hist(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                                      data = df_sdi4,
                                      summary = "IPA",
                                      plots = "calibration",
                                      cause = 1,
                                      times = 5,
                                      cens.method = "jackknife",
                                      cens.model = "km")




plotCalibration(df_sdi4_score, 
                xlim = c(0, 0.3), ylim = c(0, 0.3),
                auc.in.legend = F,brier.in.legend = F, col = 4, lwd = 1)






# now for sdi group5


df_sdi5$pred_5y = 1 - pec::predictSurvProb(m,
                                           newdata = df_sdi5, times = 5)


horizon = 5

obj = summary(survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                      data = df_sdi5), data = df_sdi5, times = horizon)

obj

OE <- (1 - obj$surv)/mean(df_sdi5$pred_5y)

alpha <- 0.05

oe_ratio_res = list("OE" = OE,
                    "2.5%" = OE*exp(-qnorm(1 - alpha/2)*sqrt(1/obj$n.event)), 
                    "97.5%" = OE*exp(qnorm(1 - alpha/2)*sqrt(1/obj$n.event)))

oe_ratio_res

# weak calibration - calibration slope 


df_sdi5$lp = predict(m, newdata = df_sdi5)

cal_slope5 = coxph(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ lp, data = df_sdi5)

cal_slope5$coefficients

alpha = 0.05

calib_slope_sdi5 = c(
  'slope' = cal_slope5$coefficients,
  "2.5%" = cal_slope5$coefficients - qnorm(1 - alpha/2)*sqrt(cal_slope5$var),
  "97.5%" = cal_slope5$coefficients + qnorm(1 - alpha/2)*sqrt(cal_slope5$var))


calib_slope_sdi5


# GND test 

df_sdi5$dec = as.numeric(cut2(df_sdi5$pred_5y, g = 10))

table(df_sdi5$dec, df_sdi5$event_cox_5y)

# fit the test 

source("scripts/gnd.R")


GND_sdi5 = GND.calib(pred = df_sdi5$pred_5y, 
                     tvar = df_sdi5$hf_admit_comp_years_5y,
                     out = df_sdi5$event_cox_5y,
                     cens.t = 5, 
                     groups = df_sdi5$dec, 
                     adm.cens = 5)

GND_sdi5$chi2gw

GND_sdi5$pvalgw

# calibration plot 

df_sdi5$pred5y = 1 - pec::predictSurvProb(object = m, newdata = df_sdi5,times = 5)


df_sdi5_score = riskRegression::Score(object = list("score" = df_sdi5$pred5y),
                                      formula = Hist(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                                      data = df_sdi5,
                                      summary = "IPA",
                                      plots = "calibration",
                                      cause = 1,
                                      times = 5,
                                      cens.method = "jackknife",
                                      cens.model = "km")




plotCalibration(df_sdi5_score, 
                xlim = c(0, 0.3), ylim = c(0, 0.3),
                auc.in.legend = F,brier.in.legend = F, col = 5)




# calibration intercept

pseudos = data.frame(df_sdi5_score$Calibration$plotframe)

pseudos$cll_pred = log(-log(1 - pseudos$risk))

fit_cal_int5 = geese(
  pseudovalue ~ offset(cll_pred),
  data = pseudos,
  id = ID,
  scale.fix = T,
  family = gaussian,
  mean.link = "cloglog",
  corstr = "independence",
  jack = T
)

fit_cal_int5

summary(fit_cal_int5)

# DCA for the whole group 

df = as.data.frame(df)

df$pred_5y = 1 - pec::predictSurvProb(object = m, newdata = df, times = 5)

df$event_cox_5y

df$event_cox_5y = factor(df$event_cox_5y)

dca = stdca_new(data = df,outcome = 'event_cox_5y',
              ttoutcome = 'hf_admit_comp_years_5y',timepoint = 5,
              predictors = "pred_5y",xstop = 0.30,cmprsk = F,graph = T,)

library(dcurves)

dcurves::dca(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ pred_5y,
  data = df,
  thresholds = c(0, 0.5),
time = 5)

# dca for the whole cohort 

threshold = seq(0, 0.3, by = 0.01)

h = 5

survfit_a = survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1,
                      data = df)

b = summary(survfit_a, times = 5)
                      
f_all = 1 - b$surv                      
                      
# calculate NB accross all thresholds 

list_nb = lapply(threshold, function(ps){
  
  # treat all 
  
  NB_all = f_all - ((1 - f_all)*(ps/1- ps))

  p_exceed <- mean(df$pred_5y > ps)
  survfit_among_exceed = try(
    summary(
      survfit(Surv(hf_admit_comp_years_5y,event_cox_5y) ~ 1, data = df[df$pred_5y> ps, ]),
      times = h
    ), silent = T
  )

# if a)no more observations above threshold, or b)among subset exceeding..
  #-- no indiv as event time >= horizon, the NB = 0
  if(class(survfit_among_exceed) == "try-error"){
    NB = 0
  }else{
    f_given_exceed <- 1 - survfit_among_exceed$surv
    TP = f_given_exceed*p_exceed
    FP <- ( 1- f_given_exceed)*p_exceed
    NB <- TP - FP*(ps/(1- ps))
  }
  
  # retugn together
  
  df_res <- data.frame(
    
    "NB" = NB,
    "treat_all" = NB_all)
  
    return(df_res)
  })


df_nb <- do.call(rbind.data.frame, list_nb)

smooth_nb = smooth(df_nb$NB, twiceit = T)



plot(
  x = df_nb$NB
)

###### END SCRIPT ##############################################################
# # now to evaluate the 5-year mean observed risk for each category.
# 
# 
# hf_cif_group = cuminc_ci_got(ftime = df$hf_admit_comp_years_5y,fstatus = df$hf_admit_comp_5y,
#               group = df$orig_cat,times = c(1,3,5),
#               print_event = T,print_var = T,cencode = 0)
# 
# 
# 
# hf_cif_group5 = cuminc_ci_got(ftime = df$hf_admit_comp_years_5y,fstatus = df$hf_admit_comp_5y,
#                              group = df$orig_cat,times = 5,
#                              print_event = T,print_var = T,cencode = 0)
# 
# 
# 
# # get the cumulative mortality at 5 years for each WATCHDM category
# 
# 
# 
# m_q = survfit(Surv(fupyears_5y, died_5y)~ orig_cat, data = df)
# 
# m_q2 = summary(m_q, times = 5)
# 
# r_q = cumu_inc_summary(m_q2, deci = 5,time.points = 1)
# 
# r_q
# 
# # plot mortality 
# 
# 
# m_plot_q = cuminc_plot_group_comb(surv.obj = m_q,.bgcolor = "gray90",.ylab = "Cumulative Mortality",.xlab = "Years",.plotsize = 4,
#                                   .risktablesize = 1,.legend_labs = c("I","II","III","IV","V"),
#                                   .xlim = 5,.ylim = 0.3,
#                                   .break.x.by = 1,.break.y.by = 0.05,
#                                   .legend = T,.group = T)
# 
# 
# m_plot_q
# 
# 
# 
# 
# ggsave(plot = m_plot_q,
#        filename = "figures/m_plot_group.tiff",
#        height = 10, width = 16, units = "in", device = "tiff", dpi = 600)
# 
# 
# # plot CIF for HFH 
# 
# 
# hf_q = mstate::Cuminc(time = df$hf_admit_comp_years_5y,
#                       status = df$hf_admit_comp_5y,group = df$orig_cat,
#                       data = df,failcodes = 1)
# 
# 
# hf2 = cif_plot(.cuminc_obj = hf_q,.ylab = "Cumulative Incident Heart Failure Hospitalization",
#                .xlab = "Years",.ylim = 0.1,.xlim = 5,.bgcolor = "gray90",.legend = T,
#                .legend_labs = c("I","II","III","IV","V"),.break.x.by = 1,.breaky.by = 0.02,.group = T)
# 
# hf2
# 
# 
# ggsave(plot = hf2, 
#        filename = "figures/hfh_cif.tiff", height = 10, width = 16,
#        
#        units = "in",device = "tiff",dpi = 600)
# 
# # calculate the HFH CSC and then the crr subHR 
# 
# 
# fit_hr = CSC(Hist(hf_admit_comp_years_5y, hf_admit_comp_5y) ~ orig_cat,
#              data = df, fitter = "cph")
# 
# fit_hr
# 
# 
# 
# 
# calc_res = function(group, coef, se){
#   hr = exp(coef)
#   lower = exp(coef - 1.96*se)
#   upper = exp(coef + 1.96*se)
#   
#   tb = data.frame(group, hr,  lower, upper)
#   return(tb)
# }
# 
# 
# calc_res(group = c(2:5), coef = c(0.2472, 0.4783, 0.7770, 1.2840), 
#          se = c(0.0175, 0.0162, 0.0152, 0.0153))
# 
# # getting the sub HR now 
# 
# crr = FGR(Hist(hf_admit_comp_years_5y, hf_admit_comp_5y) ~ orig_cat, 
#           data = df,cause = 1,y = T)
# 
# library(tidycmprsk)
# 
# # use this package 
# 
# df$hf_admit_comp_5y = factor(df$hf_admit_comp_5y)
# 
# crr_model = crr(Surv(hf_admit_comp_years_5y, hf_admit_comp_5y) ~ orig_cat, 
#                 data = df)
# 
# # am able to obtain validation metrics using CIF models with riskRegression
# # leave calibration for now and do validation. 




