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
detail <- fread(paste0(version, "/detail.csv"))

summ[, tname := ifelse(itype == "constant", "Constant", "Piecewise")]
summ[, pname := "Constant"]
summ[, sname := "1B"]
summ[exclude_pt_bigT == TRUE, estype := "only recent"]
summ[exclude_pt_bigT == FALSE, estype := "all tests"]
summ[, estype := factor(estype, levels=c("all tests", "only recent"))]

detail[, tname := ifelse(itype == "constant", "Constant", "Piecewise")]
detail[, pname := "Constant"]
detail[, sname := "1B"]
detail[exclude_pt_bigT == TRUE, estype := "only recent"]
detail[exclude_pt_bigT == FALSE, estype := "all tests"]
detail[, estype := factor(estype, levels=c("all tests", "only recent"))]

summ <- summ[, .(q, tname, t_max, estype, q_eff, estimator_type, assay_vals)]
detail <- detail[, .(q, tname, t_max, estype, estimator_type, assay_vals, estimate)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
summ[, se := lapply(se, function(x) sprintf("%.3f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.3f", x))]
