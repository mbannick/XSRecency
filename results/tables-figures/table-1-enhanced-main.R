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
version <- "~/Documents/FileZilla/xs-recent/enhanced/15-06-22-12-2/"
version <- "~/Documents/FileZilla/xs-recent/enhanced/22-06-2022-12-14-23/"

summ <- fread(paste0(version , "summary.csv"))

# TABLE RESULTS ---------------------------------------------

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[, pname := "Constant"]
summ[, sname := "1B"]

summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

summ <- summ[, .(trange, q, estimator_type, assay_vals, bias, se, mse)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.5f", x))]
summ[, se := lapply(se, function(x) sprintf("%.5f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.5f", x))]

est <- summ[assay_vals == "est"]
true <- summ[assay_vals == "true"]

est.bias <- dcast(est, trange + q ~ estimator_type, value.var=c("bias"))
true.bias <- dcast(true, trange + q ~ estimator_type, value.var=c("bias"))

est.se <- dcast(est, trange + q ~ estimator_type, value.var=c("se"))
true.se <- dcast(true, trange + q ~ estimator_type, value.var=c("se"))

est.mse <- dcast(est, trange + q ~ estimator_type, value.var=c("mse"))
true.mse <- dcast(true, trange + q ~ estimator_type, value.var=c("mse"))

EST <- cbind(est.bias, est.se[, -c(1:2)], est.mse[, -c(1:2)]) %>% data.table
TRUTH <- cbind(true.bias, true.se[, -c(1:2)], true.mse[, -c(1:2)]) %>% data.table

addtorow <- list()
addtorow$pos <- seq(5, nrow(EST), by=5) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(EST, align=rep("c", 9), digits=2)
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

