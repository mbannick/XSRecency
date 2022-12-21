rm(list=ls())
library(data.table)
library(reshape2)
source("./R/phi-functions.R")
library(tidyverse)
library(magrittr)
library(dplyr)
library(tidyr)

# Get the input and output directories
ISINFILE <- TRUE

if(ISINFILE){
  # in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/02-11-2022-11-40-45"
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/14-11-2022-14-35-34"
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/29-11-2022-20-58-03"
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/02-12-2022-17-13-19" # run with the fixes, but 5000 sims
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/05-12-2022-11-19-38" # run with the fixes, but 50,000 sims
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/08-12-2022-22-09-35" # run with bugfixes, 5000 sims (W1, W4)
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/15-12-2022-17-11-12"

  # Read in files
  f <- list.files(in.dir, full.names=T)
  f <- f[!grepl("summary", f)]
  f <- f[!grepl("detail", f)]
  f <- f[!grepl("README.md", f)]
  df <- lapply(f, fread) %>% rbindlist(fill=T)
  df[, V1 := NULL]
}

id.vars <- c("truth", "n_sims", "n", "p", "inc", "tau", "bigT", "itype",
             "window", "shadow")

for(var in c("rho", "phi_frr", "phi_tfrr", "phi_norm_mu",
             "phi_norm_sd", "phi_norm_div",
             "phi_pnorm_mu", "phi_pnorm_sd", "phi_pnorm_div",
             "frr_mix_start", "frr_mix_end",
             "ext_FRR", "duong_scale", "max_FRR", "last_point",
             "pt", "t_min", "t_max", "q", "gamma", "eta", "nu",
             "xi", "mech2", "exclude_pt_bigT")){
  if(var %in% colnames(df)){
    id.vars <- c(id.vars, var)
  }
}

wbar <- df[, lapply(.SD, mean), .SDcols=c("W1", "W2", "W3", "W4", "W5"), by=id.vars]
west <- df[, lapply(.SD, mean), .SDcols=c("EW1", "EW2", "EW3", "EW4", "EW5"), by=id.vars]

# Empirical means of W's
setnames(wbar, c("W1", "W2", "W3", "W4", "W5"),
         c("EW1", "EW2", "EW3", "EW4", "EW5"))
# Mean of estimated means of W's
setnames(west, c("EW1", "EW2", "EW3", "EW4", "EW5"),
         c("EW1E", "EW2E", "EW3E", "EW4E", "EW5E"))

df[, VW42 := n_p * p_A * (omega_TA**2 + omega_TA_var + var_TA + mu_TA**2 -
                            p * p_A * (mu_TA - omega_TA)**2 +
                            r_TA + p_A * r_TAprime * (n_p - (n_p)/n) - 2 * omega_TAstar)]

vbar <- df[, lapply(.SD, var), .SDcols=c("W1", "W2", "W3", "W4", "W5"), by=id.vars]
vest <- df[, lapply(.SD, mean), .SDcols=c("VW1", "VW2", "VW3", "VW4", "VW42", "VW5"), by=id.vars]

# Empirical variance of W's
setnames(vbar, c("W1", "W2", "W3", "W4", "W5"),
         c("VW1", "VW2", "VW3", "VW4", "VW5"))
# Mean of estimated variance of W's
setnames(vest, c("VW1", "VW2", "VW3", "VW4", "VW42", "VW5"),
         c("VW1E", "VW2E", "VW3E", "VW4E", "VW42E", "VW5E"))

# Empirical covariance of W's -- no need for name change
cbar <- df %>%
  group_by_at(id.vars) %>%
  summarise(
    C12 = cov(W1, W2),
    C13 = cov(W1, W3),
    C14 = cov(W1, W4),
    C15 = cov(W1, W5),
    C23 = cov(W2, W3),
    C24 = cov(W2, W4),
    C25 = cov(W2, W5),
    C34 = cov(W3, W4),
    C35 = cov(W3, W5),
    C45 = cov(W4, W5),
    # C14_A = cov(C14_1, C14_2),
    # C14_B = cov(C14_2, C14_3),
    .groups="drop"
  )

# Mean of the estimated covariance of W's
cest <- df[, lapply(.SD, mean), .SDcols=c("C12", "C13", "C14", "C15", "C23", "C24", "C25",
                                          "C34", "C35", "C45"), by=id.vars]

df[, p_est := n_p/n]
df[, pr := n_r/n_p]
df[, c14_1 := n * p_est * (
  EMiAiOi - p_est * pr * p_A * omega_TA -
    beta_est * (p_A * omega_TA - p_est * (1-p_B) * p_A * omega_TA)
)]
mean(df$c14_1)
cbar$C14 # ACTUAL COVARIANCE (EMPIRICAL)
cest$C14 # MEAN OF OUR ANALYTICAL COVARIANCE

df[, p := as.numeric(p)]
df[, test := truth * (1-p) / p * (p_A * omega_TAstar - p_A * (omega_TA**2 + omega_var))]

# Mean of estimated variance of W's
setnames(cest, c("C12", "C13", "C14", "C15", "C23", "C24", "C25", "C34", "C35", "C45"),
         c("C12E", "C13E", "C14E", "C15E", "C23E", "C24E", "C25E", "C34E", "C35E", "C45E"))

df_list <- list(wbar, west, vbar, vest, cbar, cest)

#merge all data frames in list
df_list <- df_list %>% reduce(full_join, by=id.vars)

df_list <- df_list[, c(id.vars,
            "EW1", "EW1E",
            "EW2", "EW2E",
            "EW3", "EW3E",
            "EW4", "EW4E",
            "EW5", "EW5E",
            "VW1", "VW1E",
            "VW2", "VW2E",
            "VW3", "VW3E",
            "VW4", "VW4E",
            "VW5", "VW5E",
            "C12", "C12E",
            "C13", "C13E",
            "C14", "C14E",
            "C15", "C15E",
            "C23", "C23E",
            "C24", "C24E",
            "C25", "C25E",
            "C34", "C34E",
            "C35", "C35E",
            "C45", "C45E"), with=F]

View(df_list)


# write.csv(df_s, "temp-variances.csv")
# df_s <- read.csv("temp-variances.csv") %>% data.table

df_s <- df[q == 1 & t_max == 2]
# df_s[, C14_test := n * (n_p/n) * (
#   eadj_est_est * (1-(n_p/n))/(n_p/n) * (omega_TAstar - omega_TA**2 - omega_TA_var) -
#   omega_TA * ((n_p/n) * n_r / n_p + beta_est * (1-(n_p/n)))
# )]

df_s[, C14] %>% mean
# df_s[, C14_test] %>% mean
# (df_s$C14 - df_s$C14_test) %>% mean


# NEW AS OF MONDAY
# df_s[, EMiAiOi]
# df_s[, EMiAiOi_test := truth * (1-0.29)/0.29 * (omega_TAstar - omega_TA**2 - omega_TA_var)]

# mean(df_s$EMiAiOi_test / df_s$EMiAiOi)
#
# with(df_s, plot(EMiAiOi ~ EMiAiOi_test))
# abline(a=0, b=1)

# Get the gamma parameters and baseline phi function
# params <- get.gamma.params(window=101/365.25, shadow=194/365.25)

# # Set up each type of phi function, will be overwritten
# phi.const <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
# phi.func <- phi.const
# ttime <- 2
# tval <- phi.none(ttime)
# phi.const <- function(t, ...) phi.none(t)*(t <= ttime) + tval*(t > ttime)
# phi.func <- phi.const

# This is \Omega_{T,A}
int1 <- integrate(phi.func, lower=0, upper=2, subdivisions=100)
int2 <- integrate(function(u) u * phi.func(u), lower=0, upper=2)
omega_TA_true <- int1$value - 0.5 * int2$value
omega_true <- int1$value

beta_true <- phi.func(4)

# This is \omega_{T,A}^*
int3 <- integrate(function(u) u**2 * phi.func(u), lower=0, upper=2)
omega_TAstar_true <- int1$value - 0.25 * int3$value

library(pracma)
myfun <- function(u, v){
  val <- (1-max(u, v)/2)*phi.func(u)*phi.func(v)
  return(val)
}
# This is \sigma^2_{\omega_{T,A}}
omega_TA2_true <- integral2(fun=myfun, xmin=0, xmax=2, ymin=0, ymax=2, reltol=1e-10)$Q
omega_var_true <- omega_TA2_true - omega_TA_true**2

# Compare them to the empirical ones
mean(df_s$omega_TA)
omega_TA_true

mean(df_s$omega_TAstar)
omega_TAstar_true

mean(df_s$omega_TA2)
omega_TA2_true

mean(df_s$omega_TA_var)
omega_var_true

df_s[, lamp := truth * (1-0.29)/0.29]

df_s[, EMiAiOi_analytical_true := (omega_TA_true * beta_true + lamp * (
  omega_TAstar_true - omega_TA2_true + (omega_true - beta_true*2)
))]

df_s[, EMiAiOi_analytical_est := (omega_TA*beta_est + eadj_est_est * (1-(n_p / n)) / (n_p/n) * (
  omega_TAstar - omega_TA2 + (omega_est - beta_est*2)
))]

plot(df_s$omega_TA ~ df_s$omega_TA2)

par(mfrow=c(1, 2))
mean(df_s$EMiAiOi_analytical_true / df_s$EMiAiOi)
with(df_s, plot(EMiAiOi_analytical_true ~ EMiAiOi))
abline(a=0, b=1)

mean(df_s$EMiAiOi_analytical_est / df_s$EMiAiOi)
with(df_s, plot(EMiAiOi_analytical_est ~ EMiAiOi))
abline(a=0, b=1)

# PLUG-IN VALUES FROM THE DERIVATION EXPRESSION
N <- 5000
p <- 0.29
nu <- 0.032 * (1-p) / p
bigT <- 2
# probability of recent -- the only thing we don't have "truth" for
pr <- mean(df_s$pr)

analytical_true <- N * p * (
  nu * (
    omega_TAstar_true - omega_TA2_true + (omega_true - beta_true * bigT)
  ) - omega_TA_true * p * ( pr - beta_true)
)

# this is very close to the estimated covariance. so our estimates are correct,
# our derivation of covariance is *not*.
analytical_true
df_s[, C14] %>% mean

# this is very close to the empirical covariance. so the issue
# is coming from the EMiAiOi analytical derivation being wrong
df_s[, analytical_modified := (
  n * p * (
    EMiAiOi - omega_TA_true * (p * pr + beta_true * (1-p))
  )
)]
mean(df_s$analytical_modified)
data.table(cbar)[q == 1 & t_max == 2, C14]

#
# df_s[, EMiAiOi_test2 := truth * (1-0.29)/0.29 * (omega_TAstar_true - omega_TA_true**2 - omega_TA_var)]
#
# mean(df_s$EMiAiOi_test2 / df_s$EMiAiOi)
#
# with(df_s, plot(EMiAiOi ~ EMiAiOi_test2))
# abline(a=0, b=1)

# df_s[, EMiAiOi_test3 := truth * (1-0.29)/0.29 *
#        (omega_TAstar_true - omega_TA_true**2 - omega_TA_var_true)]
# mean(df_s$EMiAiOi_test3 / df_s$EMiAiOi)
# with(df_s, plot(EMiAiOi ~ EMiAiOi_test3))
# abline(a=0, b=1)


# Do the same thing with Var(W4)
r_TA <- df_s[, r_TA] %>% mean
r_TAprime <- df_s[, r_TAprime] %>% mean
omega_TA_varE <- df_s[, omega_TA_var] %>% mean

analytical <- N * p * (
  r_TA + r_TAprime * (N * p - p) + omega_TA_true**2 * (1-p) + omega_TA_varE
)
analytical
vbar$VW4[5]
vest$VW4[5]






# Monday, Dec. 5th
df_list[q == 1 & t_max == 4, .(C45, C45E)]
df2 <- df[q == 1 & t_max == 4]
mean(df2$p_A)
mean(df2$p_B)
df2[, correct45 := - n * p^2 * beta_est * 0.5 * 0.5 * 3 * (1 - mean(df2$omega_TA))]
df2$correct45 %>% unique
mean(df2$C45)

df[, W4true := n * p * p_A * mu_TA - n * p * p_A * omega_TA]
df[, W5true := n * p * beta_est * p_B * mu_TB]


df2[, AiTi := n_p * p_A * mu_TA]
df2[, AiOi := n_p * p_A * omega_TA]
df2[, BiTi := n_p * p_B * mu_TB]

cov(df2$W4, df2$W5)
cov(df2$AiTi - df2$AiOi, df2$BiTi * df2$beta_est)
cov(df2$n_p, df2$n_p)


Ti <- runif(n=100000, min=0, max=2)
mean(Ti)
mean(Ti**2)

df2 <- df[q == 1 & t_max == 2]
df2$mu_TA %>% mean
(df2$mu_TA ** 2 + df2$var_TA) %>% mean
