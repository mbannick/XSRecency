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

# DUONG + GAMMA PHI FIX + N 5000 + SENSITIVITY FRR
version <- "~/Documents/FileZilla/xs-recent/02-04-21-10/"
summ <- fread(paste0(version , "summary.csv"))

summ <- summ[estimator_type == "adj"]

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
detail[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]

summ[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
summ[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]
detail[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
detail[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]

summ[window == 101, sname := "1 A-C"]
summ[window == 248, sname := "2 A-C"]
detail[window == 101, sname := "1 A-C"]
detail[window == 248, sname := "2 A-C"]

summ[pname == "Zero" & sname == "1 A-C", assay := "1A"]
summ[pname == "Constant" & sname == "1 A-C", assay := "1B"]
summ[pname == "Non-constant" & sname == "1 A-C", assay := "1C"]
summ[pname == "Zero" & sname == "2 A-C", assay := "2A"]
summ[pname == "Constant" & sname == "2 A-C", assay := "2B"]
summ[pname == "Non-constant" & sname == "2 A-C", assay := "2C"]

detail[pname == "Zero" & sname == "1 A-C", assay := "1A"]
detail[pname == "Constant" & sname == "1 A-C", assay := "1B"]
detail[pname == "Non-constant" & sname == "1 A-C", assay := "1C"]
detail[pname == "Zero" & sname == "2 A-C", assay := "2A"]
detail[pname == "Constant" & sname == "2 A-C", assay := "2B"]
detail[pname == "Non-constant" & sname == "2 A-C", assay := "2C"]

settings <- summ$sname %>% unique
phis <- summ$pname %>% unique
epis <- summ$tname %>% unique

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

      adj.exp <- this.row[estimator == "adj_true", expected_bias]

      row <- c(
        epi,
        assay,

        "$\\times$",
        # sprintf('%.2f', this.row[estimator == "adj_true", bias] * 100),
        # sprintf('%.2f', this.row[estimator == "adj_true", se] * 100),
        # sprintf('%.2f', this.row[estimator == "adj_true", see] * 100),
        # sprintf('%05.2f', this.row[estimator == "adj_true", cover] * 100),

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
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{5}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline ",
                      "Incidence & Assay & \\multicolumn{5}{c}{$\\hat{\\Omega}_{T^*}$, $\\hat{\\beta}_{T^*}$} \\\\\n",
                      "\\hline ",
                      "& & Asm. & Bias & SE & SEE & Cov \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{6}{c}{Recency Assay 1A-C} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ",
                      "\\multicolumn{6}{c}{Recency Assay 2A-C} \\\\\n",
                      "\\hline ",
                      "\\hline ",
                      "\\hline ")
tab <- xtable(data, align=rep("c", 8), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE,
      sanitize.text.function=identity)
