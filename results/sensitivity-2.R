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
library(RColorBrewer)

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/15-02-21-19/"
summ <- fread(paste0(version , "summary.csv"))

summ[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]

summ[window == 71, sname := "Brookmeyer et al. 2013"]
summ[window == 248, sname := "Laeyendecker et al. 2018"]

summ <- summ[estimator == "adj_est"]
summ[, mixture := !is.na(frr_mix_start)]
setorder(summ, window, mixture)

settings <- summ$sname %>% unique

data <- data.table()

for(setting in settings){
  for(mix in c(F, T)){
    const <- 100

    if(mix) name <- "Mixture" else name <- "Uniform"

    this.row <- summ[(sname == setting) & (mixture == mix)]

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
addtorow$pos <- list(0, 0, 0, 0, 2, 2, 2)
addtorow$command <- c("Setting & Bias & SE & SEE & Cover \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{5}{c}{Brookmeyer et al. 2013} \\\\\n",
                      "\\hline ",
                      "\\multicolumn{5}{c}{Laeyendecker et al. 2018} \\\\\n",
                      "\\hline ")
tab <- xtable(data, align=rep("c", 6), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)
