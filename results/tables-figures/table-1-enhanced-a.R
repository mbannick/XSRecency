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

version <- "~/Documents/FileZilla/xs-recent/enhanced/15-12-2022-17-23-00/"
summ <- fread(paste0(version , "summary.csv"))

detail <- fread(paste0(version, "/detail.csv"))

# DETAIL RESULTS FOR Q EFF ---------------------------------------------

detail <- detail[, .(q, mech2, q_eff)]
detail <- detail[, lapply(.SD, mean), .SDcols=c("q_eff"), by=c("q", "mech2")]

summ <- merge(summ, detail, by=c("q", "mech2"))

# TABLE RESULTS ---------------------------------------------

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[, pname := "Constant"]
summ[, sname := "1B"]
summ[mech2 == FALSE, mechtype := "Base"]
summ[mech2 == TRUE, mechtype := "Base + RI"]

summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

summ <- summ[, .(q, mechtype, q_eff, estimator_type, assay_vals, bias, se, mse)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
summ[, se := lapply(se, function(x) sprintf("%.3f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.3f", x))]

est <- summ[assay_vals == "est"]
true <- summ[assay_vals == "true"]

est.bias <- dcast(est, q + mechtype + q_eff ~ estimator_type, value.var=c("bias"))
true.bias <- dcast(true, q + mechtype + q_eff ~ estimator_type, value.var=c("bias"))

est.se <- dcast(est, q + mechtype + q_eff ~ estimator_type, value.var=c("se"))
true.se <- dcast(true, q + mechtype + q_eff ~ estimator_type, value.var=c("se"))

est.mse <- dcast(est, q + mechtype + q_eff ~ estimator_type, value.var=c("mse"))
true.mse <- dcast(true, q + mechtype + q_eff ~ estimator_type, value.var=c("mse"))

EST <- cbind(est.bias, est.se[, -c(1:3)], est.mse[, -c(1:3)]) %>% data.table
TRUTH <- cbind(true.bias, true.se[, -c(1:3)], true.mse[, -c(1:3)]) %>% data.table

addtorow <- list()
addtorow$pos <- seq(2, nrow(EST), by=2) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(EST, align=rep("c", 10), digits=2)
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

