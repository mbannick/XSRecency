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

summ <- fread(paste0(version , "summary.csv"))

# TABLE RESULTS ---------------------------------------------

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[, pname := "Constant"]
summ[, sname := "2B"]

summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

summ <- summ[, .(trange, q, gamma, eta, nu, estimator_type, assay_vals, bias, se)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, bias := lapply(bias, function(x) sprintf("%.3f", x))]
summ[, se := lapply(se, function(x) sprintf("%.3f", x))]

est <- summ[assay_vals == "est"]
true <- summ[assay_vals == "true"]

est.bias <- dcast(est, trange + q + gamma + eta + nu ~ estimator_type, value.var=c("bias"))
true.bias <- dcast(true, trange + q + gamma + eta + nu ~ estimator_type, value.var=c("bias"))

est.se <- dcast(est, trange + q + gamma + eta + nu ~ estimator_type, value.var=c("se"))
true.se <- dcast(true, trange + q + gamma + eta + nu ~ estimator_type, value.var=c("se"))

EST <- cbind(est.bias, est.se[, -c(1:5)]) %>% data.table
TRUTH <- cbind(true.bias, true.se[, -c(1:5)]) %>% data.table

for(tr in unique(TRUTH$trange)){
  TRUTH.sub <- TRUTH[trange == tr]
  EST.sub <- EST[trange == tr]

  addtorow <- list()
  addtorow$pos <- seq(4, nrow(TRUTH.sub), by=4) %>% as.list
  addtorow$command <- rep("\\hline \n", length(addtorow$pos))

  # tab <- xtable(TRUTH.sub, align=rep("c", 10), digits=2,
  #               caption=paste0(tr, " TRUTH"))
  # print(tab, include.rownames=FALSE,
  #       add.to.row = addtorow)
  tab <- xtable(EST.sub, align=rep("c", 10), digits=2,
                caption=paste0(tr, " EST"))
  print(tab, include.rownames=FALSE,
        add.to.row = addtorow)

}

