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

version <- "~/Documents/FileZilla/xs-recent/enhanced/15-12-2022-17-21-56/"
summ <- fread(paste0(version , "/summary.csv"))
detail <- fread(paste0(version, "/detail.csv"))

# DETAIL RESULTS FOR Q EFF ---------------------------------------------

detail <- detail[, .(q, itype, t_max, exclude_pt_bigT, q_eff)]
detail <- detail[, lapply(.SD, mean), .SDcols=c("q_eff"), by=c("q", "itype", "t_max", "exclude_pt_bigT")]

summ <- merge(summ, detail, by=c("q", "itype", "t_max", "exclude_pt_bigT"))
summ[, tname := ifelse(itype == "constant", "Constant", "Piecewise")]
summ[, pname := "Constant"]
summ[, sname := "1B"]
summ[exclude_pt_bigT == TRUE, estype := "only recent"]
summ[exclude_pt_bigT == FALSE, estype := "all tests"]
summ[, estype := factor(estype, levels=c("all tests", "only recent"))]

summ <- summ[, .(q, tname, t_max, estype, q_eff, estimator_type, assay_vals, bias, se, mse, cover_rob, cover_asm)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 1e4]
summ[, cover_rob := cover_rob * 100]
summ[, cover_asm := cover_asm * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.2f", x))]
summ[, se := lapply(se, function(x) sprintf("%.2f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.2f", x))]
summ[, cover_rob := lapply(cover_rob, function(x) sprintf("%.2f", x))]
summ[, cover_asm := lapply(cover_asm, function(x) sprintf("%.2f", x))]
summ[, q_eff := sapply(q_eff, function(x) sprintf("%.2f", x))]

est <- summ[assay_vals == "est"]

bias <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("bias"))
se <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("se"))
mse <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("mse"))
cover_rob <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("cover_rob"))
cover_asm <- dcast(est, q + tname + t_max + estype + q_eff ~ estimator_type, value.var=c("cover_asm"))

bias$estype <- as.character(bias$estype)

ADJ.C <- c(rep("--", 3),  bias[1,]$estype, "--", bias[1,]$adj, se[1,]$adj, mse[1,]$adj, cover_asm[1,]$adj)
ADJ.P <- c(rep("--", 3),  bias[5,]$estype, "--", bias[5,]$adj, se[5,]$adj, mse[5,]$adj, cover_asm[5,]$adj)

EST <- cbind(bias$q, bias$tname, bias$t_max, bias$estype, bias$q_eff, bias$eadj, se$eadj, mse$eadj, cover_rob$eadj)
DF <- rbind(ADJ.C, EST[1:4,], ADJ.P, EST[5:8,])
colnames(DF) <- c("q", "tname", "t_max", "estype", "q_eff", "Bias", "SE", "MSE", "Coverage")

addtorow <- list()
addtorow$pos <- c(1, 5, 6, 8) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(DF, align=rep("c", 10), digits=2)
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

