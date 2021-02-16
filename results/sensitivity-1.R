# ----------------------------------------------------------------
# SENSITIVITY ANALYSIS -- COMPARING INC. PREV. SETTINGS
# FEB 2021
# ----------------------------------------------------------------
# ----------------------------------------------------------------

rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(ggh4x)
library(RColorBrewer)

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/15-02-21-12/"
detail <- fread(paste0(version, "detail.csv"))
summ <- fread(paste0(version , "summary.csv"))

summ <- summ[window == 71]
setorder(summ, p, inc, estimator)

exp <- summ[, .(p, inc, estimator_type, expected_bias)]
exp <- exp[!duplicated(exp)]

snap.exp <- exp[estimator_type == "snap"]
adj.exp <- exp[estimator_type == "adj"]

summ <- summ[, .(p, inc, estimator, bias, se, see, cover)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, see := see * 100]
summ[, cover := cover * 100]

bias <- reshape2::dcast(summ, p + inc ~ estimator, value.var="bias") %>% data.table
se <- reshape2::dcast(summ, p + inc ~ estimator, value.var="se") %>% data.table
see <- reshape2::dcast(summ, p + inc ~ estimator, value.var="see") %>% data.table
cover <- reshape2::dcast(summ, p + inc ~ estimator, value.var="cover") %>% data.table

df <- cbind(
  bias$p,
  bias$inc,

  snap.exp$expected_bias * 100,

  bias$snap_true,
  se$snap_true,
  see$snap_true,
  cover$snap_true,

  bias$snap_est,
  se$snap_est,
  see$snap_est,
  cover$snap_est,

  adj.exp$expected_bias * 100,

  bias$adj_true,
  se$adj_true,
  see$adj_true,
  cover$adj_true,

  bias$adj_est,
  se$adj_est,
  see$adj_est,
  cover$adj_est
)

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{9}{c}{Snapshot Estimator (1)} & \\multicolumn{9}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Prev. & Inc. & \\multicolumn{5}{c}{$\\mu$} & \\multicolumn{4}{c}{$\\hat{\\mu}$} & \\multicolumn{5}{c}{$\\Omega_{T^*}$, $\\beta_{T^*}$} & \\multicolumn{4}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & E[Bias] & Bias & SE & SEE & Cov & Bias & SE & SEE & Cov & E[Bias] & Bias & SE & SEE & Cov & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ")
tab <- xtable(df, align=rep("c", 21), digits=c(2, 2, 3, rep(2, 18)), caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)
