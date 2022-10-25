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

version <- "~/Documents/FileZilla/xs-recent/enhanced/27-07-2022-21-53-56/"
summ <- fread(paste0(version , "/summary.csv"))

summ[, tname := ifelse(itype == "constant", "Constant", "Piecewise")]
summ[, pname := "Constant"]
summ[, sname := "1B"]
summ[exclude_pt_bigT == TRUE, estype := "only recent"]
summ[exclude_pt_bigT == FALSE, estype := "all tests"]
summ[, estype := factor(estype, levels=c("all tests", "only recent"))]

summ <- summ[, .(q, tname, t_max, estype, q_eff, estimator_type, assay_vals, bias, se, mse)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
summ[, se := lapply(se, function(x) sprintf("%.3f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.3f", x))]

est <- summ[assay_vals == "est"]
true <- summ[assay_vals == "true"]

est.bias <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("bias"))
true.bias <- dcast(true, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("bias"))

est.se <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("se"))
true.se <- dcast(true, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("se"))

est.mse <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("mse"))
true.mse <- dcast(true, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("mse"))

EST <- cbind(est.bias, est.se[, -c(1:5)], est.mse[, -c(1:5)]) %>% data.table
TRUTH <- cbind(true.bias, true.se[, -c(1:5)], true.mse[, -c(1:5)]) %>% data.table

addtorow <- list()
addtorow$pos <- seq(4, nrow(EST), by=4) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(EST, align=rep("c", 12), digits=2)
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

