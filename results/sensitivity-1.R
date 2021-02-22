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
library(magrittr)

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/22-02-21-11/"
detail <- fread(paste0(version, "detail.csv"))
summ <- fread(paste0(version , "summary.csv"))

summ <- summ[window == 71]
setorder(summ, p, inc, estimator)

summ <- summ[, .(p, phi_tfrr, inc, estimator, bias, se, see, cover)]
summ[, bias := bias * 100]
summ[, se := se * 100]
summ[, see := see * 100]
summ[, cover := cover * 100]

snap <- summ[estimator == "snap_est" & is.na(phi_tfrr)]
adj <- summ[estimator == "adj_est" & !is.na(phi_tfrr)]

df <- cbind(
  snap[, .(p, inc, bias, se, see, cover)],
  adj[, .(bias, se, see, cover)]
)

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{4}{c}{Snapshot Estimator (1)} & \\multicolumn{4}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Prev. & Inc. & \\multicolumn{4}{c}{$\\hat{\\mu}$} & \\multicolumn{4}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & Bias & SE & SEE & Cov & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ")
tab <- xtable(df, align=rep("c", 11), digits=c(2, 2, 3, rep(2, 8)), caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)
