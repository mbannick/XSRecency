# ----------------------------------------------------------------
# TABLE FOR RECALL BIAS
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

version <- "~/Documents/FileZilla/xs-recent/enhanced/21-12-2022-10-49-35"
summ <- fread(paste0(version , "/summary.csv"))

detail <- fread(paste0(version, "/detail.csv"))

# TABLE RESULTS ---------------------------------------------

detail <- detail[, .(q, gamma, eta, xi, q_eff)]
detail <- detail[, lapply(.SD, mean), .SDcols=c("q_eff"), by=c("q", "gamma", "eta", "xi")]
summ <- merge(summ, detail, by=c("q", "gamma", "eta", "xi"))

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[, pname := "Constant"]
summ[, sname := "1B"]

summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

summ <- summ[, .(q, gamma, eta, xi, q_eff, estimator_type, assay_vals, bias, se, see_rob, see_asm, mse, cover_rob, cover_asm)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, mse := mse * 1e4]
summ[, see_rob := see_rob * 100]
summ[, see_asm := see_asm * 100]
summ[, cover_rob := cover_rob * 100]
summ[, cover_asm := cover_asm * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
summ[, se := lapply(se, function(x) sprintf("%.3f", x))]
summ[, mse := lapply(mse, function(x) sprintf("%.3f", x))]
summ[, see_rob := lapply(see_rob, function(x) sprintf("%.3f", x))]
summ[, see_asm := lapply(see_asm, function(x) sprintf("%.3f", x))]
summ[, cover_rob := lapply(cover_rob, function(x) sprintf("%.2f", x))]
summ[, cover_asm := lapply(cover_asm, function(x) sprintf("%.2f", x))]
summ[, q_eff := sapply(q_eff, function(x) sprintf("%.2f", x))]
summ[, eta := sapply(eta, function(x) sprintf("%.2f", x))]
summ[, xi := sapply(xi, function(x) sprintf("%.2f", x))]
summ[, gamma := sapply(gamma, function(x) sprintf("%.2f", x))]

est <- summ[assay_vals == "est"]

bias <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("bias"))
se <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("se"))
see_rob <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("see_rob"))
see_asm <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("see_asm"))
mse <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("mse"))
cover_rob <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("cover_rob"))
cover_asm <- dcast(est, q + eta + xi + gamma + q_eff ~ estimator_type, value.var=c("cover_asm"))

ADJ <- c(rep("--", 5),  bias[1,]$adj, se[1,]$adj, mse[1,]$adj, cover_rob[1,]$adj, "--")
EST <- cbind(bias$q, bias$eta, bias$xi, bias$gamma, bias$q_eff, bias$eadj, se$eadj, mse$eadj, cover_rob$eadj, cover_asm$eadj)
DF <- rbind(ADJ, EST)
colnames(DF) <- c("q", "eta", "xi", "gamma", "q_eff", "Bias", "SE", "MSE", "Coverage (Robust)", "Coverage")

tab <- xtable(DF, align=rep("c", 11), digits=2)
print(tab, include.rownames=FALSE)
