# ----------------------------------------------------------------
# TABLES AND FIGURES FOR CROSS-SECTIONAL RECENCY ASSAY PERFORMANCE
# FEB 2021
# ----------------------------------------------------------------
# ----------------------------------------------------------------

rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(ggh4x)
library(RColorBrewer)
library(magrittr)
source("~/repos/XSRecency/R/phi-functions.R")
source("~/repos/XSRecency/R/data-generator.R")

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/enhanced/20-07-2022-10-39-02"
summ_all <- fread(paste0(version , "/summary.csv"))

for(q_val in unique(summ_all$q)){

  summ <- summ_all[q == q_val]

  # TABLE RESULTS ---------------------------------------------

  summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
  summ[, pname := "Constant"]
  summ[, sname := "1B"]

  summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

  summ <- summ[, .(q, q_eff, gamma, eta, xi, estimator_type, assay_vals, bias, se, mse)]
  summ[, bias := bias * 100]
  summ[, se := se * 100]
  summ[, mse := mse * 100]
  summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
  summ[, se := lapply(se, function(x) sprintf("%.3f", x))]
  summ[, mse := lapply(mse, function(x) sprintf("%.3f", x))]

  est <- summ[assay_vals == "est"]
  true <- summ[assay_vals == "true"]

  est.bias <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("bias"))
  true.bias <- dcast(true, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("bias"))

  est.se <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("se"))
  true.se <- dcast(true, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("se"))

  est.mse <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("mse"))
  true.mse <- dcast(true, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("mse"))

  EST <- cbind(est.bias, est.se[, -c(1:5)], est.mse[, -c(1:5)]) %>% data.table
  TRUTH <- cbind(true.bias, true.se[, -c(1:5)], true.mse[, -c(1:5)]) %>% data.table

  setorder(EST, q, eta, xi, gamma)
  setorder(TRUTH, q, eta, xi, gamma)

  addtorow <- list()
  addtorow$pos <- c(3, 9) %>% as.list
  addtorow$command <- rep("\\hline \n", length(addtorow$pos))

  tab <- xtable(TRUTH, align=rep("c", 12), digits=2)
  print(tab, include.rownames=FALSE,
        add.to.row = addtorow)

}
