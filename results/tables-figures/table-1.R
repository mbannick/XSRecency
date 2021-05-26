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
version <- "~/Documents/FileZilla/xs-recent/16-05-21-10/"

detail <- fread(paste0(version, "detail.csv"))
summ <- fread(paste0(version , "summary.csv"))

# TABLE RESULTS ---------------------------------------------

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
detail[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]

summ[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Zero"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Constant"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Non-constant"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & !is.na(phi_pnorm_mu), pname := "Increasing"]
detail[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Zero"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Constant"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Non-constant"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & !is.na(phi_pnorm_mu), pname := "Increasing"]

summ[window == 101, sname := "1 A-D"]
summ[window == 248, sname := "2 A-D"]
detail[window == 101, sname := "1 A-D"]
detail[window == 248, sname := "2 A-D"]

summ[pname == "Zero" & sname == "1 A-D", assay := "1A"]
summ[pname == "Constant" & sname == "1 A-D", assay := "1B"]
summ[pname == "Non-constant" & sname == "1 A-D", assay := "1C"]
summ[pname == "Increasing" & sname == "1 A-D", assay := "1D"]
summ[pname == "Zero" & sname == "2 A-D", assay := "2A"]
summ[pname == "Constant" & sname == "2 A-D", assay := "2B"]
summ[pname == "Non-constant" & sname == "2 A-D", assay := "2C"]
summ[pname == "Increasing" & sname == "2 A-D", assay := "2D"]

detail[pname == "Zero" & sname == "1 A-D", assay := "1A"]
detail[pname == "Constant" & sname == "1 A-D", assay := "1B"]
detail[pname == "Non-constant" & sname == "1 A-D", assay := "1C"]
detail[pname == "Increasing" & sname == "1 A-D", assay := "1D"]
detail[pname == "Zero" & sname == "2 A-D", assay := "2A"]
detail[pname == "Constant" & sname == "2 A-D", assay := "2B"]
detail[pname == "Non-constant" & sname == "2 A-D", assay := "2C"]
detail[pname == "Increasing" & sname == "2 A-D", assay := "2D"]

settings <- summ$sname %>% unique
phis <- summ$pname %>% unique
# epis <- summ$tname %>% unique
epis <- c("Constant", "Linear", "Exponential")
phis <- c("Zero", "Constant", "Non-constant", "Increasing")

data <- data.table()

for(setting in settings){
  print(setting)
  for(phi in phis){
    print(phi)
    for(epi in epis){
      print(epi)

      const <- 100

      this.row <- summ[(sname == setting) & (tname == epi) & (pname == phi)]
      assay <- this.row[, assay][1]

      snap.exp <- this.row[estimator == "snap_true", expected_bias]
      adj.exp <- this.row[estimator == "adj_true", expected_bias]

      row <- c(
        epi,
        assay,

        "$\\times$",
        sprintf('%.2f', this.row[estimator == "snap_est", bias] * 100),
        sprintf('%.2f', this.row[estimator == "snap_est", se] * 100),
        sprintf('%.2f', this.row[estimator == "snap_est", see] * 100),
        sprintf('%05.2f', this.row[estimator == "snap_est", cover] * 100),

        "$\\times$",
        sprintf('%.2f', this.row[estimator == "adj_est", bias] * 100),
        sprintf('%.2f', this.row[estimator == "adj_est", se] * 100),
        sprintf('%.2f', this.row[estimator == "adj_est", see] * 100),
        sprintf('%05.2f', this.row[estimator == "adj_est", cover] * 100)
      )
      data <- rbind(data, t(row))
    }
  }
}

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 9, 9, 9, 12, 15)
addtorow$pos <- list(0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 9, 12, 12, 12, 15, 18, 21)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{5}{c}{Snapshot Estimator (1)} & \\multicolumn{5}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Incidence & Assay & \\multicolumn{5}{c}{$\\hat{\\mu}$} & \\multicolumn{5}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & Asm. & Bias & SE & SEE & Cov & Asm. & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{11}{c}{Recency Assay 1A-D} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{11}{c}{Recency Assay 2A-D} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ")
tab <- xtable(data, align=rep("c", 13), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE,
      sanitize.text.function=identity)

