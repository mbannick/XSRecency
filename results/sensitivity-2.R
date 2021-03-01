# ----------------------------------------------------------------
# SENSITIVITY ANALYSIS -- COMPARING FRR MIXTURES SETTINGS
# FEB 2021
# ----------------------------------------------------------------
# ----------------------------------------------------------------

rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(ggh4x)
library(magrittr)
library(RColorBrewer)

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/28-02-21-16/"
summ <- fread(paste0(version , "summary.csv"))

summ[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]

summ[window == 71, sname := "Brookmeyer et al. 2013"]
summ[window == 248, sname := "Laeyendecker et al. 2018"]

summ <- summ[estimator == "adj_est"]
summ[is.na(frr_mix_start), mixtype := "Unif(2, 12)"]
summ[frr_mix_start == 2 & frr_mix_end == 5, mixtype := "Unif(2, 12) + Unif(2, 5)"]
summ[frr_mix_start == 2 & frr_mix_end == 9, mixtype := "Unif(2, 12) + Unif(2, 9)"]

setorder(summ, window, mixtype)

settings <- summ$sname %>% unique

data <- data.table()

for(setting in settings){
  for(mixtype in unique(summ$mixtype)){
    const <- 100

    name <- mixtype

    this.row <- summ[(sname == setting) & (mixtype == name)]

    snap.exp <- this.row[estimator == "snap_true", expected_bias]
    adj.exp <- this.row[estimator == "adj_true", expected_bias]

    row <- c(
      name,
      sprintf('%.2f', this.row[estimator == "adj_est", bias] * 100),
      sprintf('%.2f', this.row[estimator == "adj_est", se] * 100),
      sprintf('%.2f', this.row[estimator == "adj_est", see] * 100),
      sprintf('%05.2f', this.row[estimator == "adj_est", cover] * 100)
    )
    data <- rbind(data, t(row))
  }
}

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 3, 3, 3)
addtorow$command <- c("Setting & Bias & SE & SEE & Cover \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{5}{c}{Recency Assay 1C} \\\\\n",
                      "\\hline ",
                      "\\multicolumn{5}{c}{Recency Assay 2C} \\\\\n",
                      "\\hline ")
tab <- xtable(data, align=rep("c", 6), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)
