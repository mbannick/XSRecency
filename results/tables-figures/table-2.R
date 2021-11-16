rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(ggh4x)
library(RColorBrewer)
library(magrittr)

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/17-04-21-15/"

summ <- fread(paste0(version , "summary.csv"))

setorder(summ, p, inc, estimator)
summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]

summ[window == 101, sname := "1 A-C"]
summ[window == 248, sname := "2 A-C"]

summ[pname == "Non-constant" & sname == "1 A-C", assay := "1C"]
summ[pname == "Non-constant" & sname == "2 A-C", assay := "2C"]

summ[ext_FRR == FALSE & is.na(max_FRR), dist := "Uniform 2-12"]
summ[ext_FRR == TRUE & is.na(max_FRR), dist := "Duong et al. 2015"]
summ[ext_FRR == TRUE & max_FRR == 5, dist := "Duong et al. 2015 (truncated)"]

summ <- summ[estimator == "adj_est"]

summ <- summ[, .(assay, dist, bias, se, see, cover)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, see := see * 100]
summ[, cover := cover * 100]

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{4}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Assay & Long Infected Distribution & \\multicolumn{4}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ")
tab <- xtable(summ, align=rep("c", 7), digits=c(2, 2, 3, rep(2, 4)), caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)
