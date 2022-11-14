rm(list=ls())
library(data.table)
library(reshape2)
library(tidyverse)
library(magrittr)
library(dplyr)
library(tidyr)

# Get the input and output directories
ISINFILE <- FALSE

if(ISINFILE){
  # in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/02-11-2022-11-40-45"
  in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/12-11-2022-14-52-04"

  # Read in files
  f <- list.files(in.dir, full.names=T)
  f <- f[!grepl("summary", f)]
  f <- f[!grepl("detail", f)]
  f <- f[!grepl("README.md", f)]
  df <- lapply(f, fread) %>% rbindlist(fill=T)
  df[, V1 := NULL]
}

id.vars <- c("truth", "n_sims", "seed", "n", "p", "inc", "tau", "bigT", "itype",
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

vbar <- df[, lapply(.SD, var), .SDcols=c("W1", "W2", "W3", "W4", "W5"), by=id.vars]
vest <- df[, lapply(.SD, mean), .SDcols=c("VW1", "VW2", "VW3", "VW4", "VW5"), by=id.vars]

# Empirical variance of W's
setnames(vbar, c("W1", "W2", "W3", "W4", "W5"),
         c("VW1", "VW2", "VW3", "VW4", "VW5"))
# Mean of estimated variance of W's
setnames(vest, c("VW1", "VW2", "VW3", "VW4", "VW5"),
         c("VW1E", "VW2E", "VW3E", "VW4E", "VW5E"))

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
    C14_A = cov(C14_1, C14_2),
    C14_B = cov(C14_2, C14_3),
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
cbar$C14
cest$C14

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

