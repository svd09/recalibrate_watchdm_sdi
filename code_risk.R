library(survival)
library(tidyverse)


fit = coxph(Surv(time, status) ~ sex, data = lung)

p = log(predict(fit, newdata = lung, type = "expected"))

p

age = 0.039
fpg = 0.002
bmi = 0.057
sbp = 0.015
dbp = -0.025
scr = 0.710
hdl = -0.023
qrs = 0.583
mi = 0.600
cabg = 0.673

 age  = 70*0.039 
 fpg = 100*0.002 
 bmi = 40*0.057 
 sbp = 120*0.015 
 dbp = 80*0.025 
 cr = 1.4*0.710 
 qrs = 1*0.583 
 hdl = 40*0.023 
 mi = 1*0.600 
 cabg = 1*0.673
 
 )70*0.039 + 80*0.002 + 30*0.057 + 120*0.015 - 80*0.025 


lp = 70*0.039 + 100*0.002 + 40*0.057 + 120*0.015 - 80*0.025 
    + 1.4*0.710 + 1*0.583 - 40*0.023 + 1*0.600 + 1*0.673


bh = 1 - 0.97358


S0 = 0.97358

S0

pred5 <-  as.vector(1 - S0**exp(lp))

# create a df 
# 

df = tibble(
    age = rnorm(n = 100, mean = 60, sd = 5),
    fpg = rnorm(n = 100, mean = 75, sd = 4),
    bmi = rnorm(n = 100, mean = 26, sd = 2),
    sbp = rnorm(n = 100, mean = 130, sd = 4),
    dbp = rnorm(n = 100, mean = 75, sd = 3),
    creat = rnorm(n = 100, mean = 1.2, sd = 0.4),
    hdl = rnorm(n = 100, mean = 30, sd = 4),
    qrs = rbinom(n = 100,size = 1, prob = 0.05),
    mi = rbinom(n = 100, size = 1, prob = 0.1),
    cabg = rbinom(n = 100, size = 1, prob = 0.05))
    

summary(df)


S0 = 0.97358*5

coef <- c(0.342, 0.574, 0.304, -0.811,  0.362)

hr = c(1.04, 1.00, 1.06, 1.02, 0.98, 2.03, 
       0.98, 1.79, 1.82, 1.96)

coef = log(hr)

coef

des_matr2 <- as.data.frame(model.matrix(~ age + fpg + bmi + sbp + dbp + creat + hdl + 
                        qrs + mi + cabg,data = df))


des_matr2$`(Intercept)` <- NULL

df$PI <- as.vector(as.matrix(des_matr2) %*% cbind(coef))

df$pred5 <- 1 - as.vector(1 - S0**exp(df$PI))

summary(df$pred5)

hist(df$pred5)


