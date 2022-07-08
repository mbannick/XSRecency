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

# MAIN VERSION, LAST POINT INTEGRATION + NEW PHI FUNCTION
version <- "~/Documents/FileZilla/xs-recent/enhanced/22-06-2022-13-51-46/"

summ <- fread(paste0(version , "summary.csv"))

# TABLE RESULTS ---------------------------------------------

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[, pname := "Constant"]
summ[, sname := "1B"]

summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

summ <- summ[, .(q_eff, gamma, eta, nu, xi, estimator_type, assay_vals, bias, se, mse)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
summ[, se := lapply(se, function(x) sprintf("%.3f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.3f", x))]

est <- summ[assay_vals == "est"]
true <- summ[assay_vals == "true"]

est.bias <- dcast(est, gamma + xi + eta + nu + q_eff ~ estimator_type, value.var=c("bias"))
true.bias <- dcast(true, gamma + xi + eta + nu + q_eff ~ estimator_type, value.var=c("bias"))

est.se <- dcast(est, gamma + xi + eta + nu + q_eff ~ estimator_type, value.var=c("se"))
true.se <- dcast(true, gamma + xi + eta + nu + q_eff ~ estimator_type, value.var=c("se"))

est.mse <- dcast(est, gamma + xi + eta + nu + q_eff ~ estimator_type, value.var=c("mse"))
true.mse <- dcast(true, gamma + xi + eta + nu + q_eff ~ estimator_type, value.var=c("mse"))

EST <- cbind(est.bias, est.se[, -c(1:5)], est.mse[, -c(1:5)]) %>% data.table
TRUTH <- cbind(true.bias, true.se[, -c(1:5)], true.mse[, -c(1:5)]) %>% data.table

setorder(EST, gamma, xi, eta, nu)
setorder(TRUTH, gamma, xi, eta, nu)

TRUTH1 <- TRUTH[(nu != 0.25) & (xi != 0.25) & (eta != 0.25)]
TRUTH25 <- TRUTH[(nu != 0.1) & (xi != 0.1) & (eta != 0.1)]

addtorow <- list()
addtorow$pos <- seq(5, nrow(TRUTH25), by=5) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(TRUTH25, align=rep("c", 12), digits=2)
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

