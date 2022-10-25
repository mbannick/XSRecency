# ----------------------------------------------------------------
# TABLES AND FIGURES FOR CROSS-SECTIONAL RECENCY ASSAY PERFORMANCE
# FEB 2021
# ----------------------------------------------------------------
# ----------------------------------------------------------------

rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(ggpattern)
library(ggh4x)
library(RColorBrewer)
library(magrittr)
source("~/repos/XSRecency/R/phi-functions.R")
source("~/repos/XSRecency/R/data-generator.R")

# READ IN VERSIONED RESULTS ---------------------------------

version <- "~/Documents/FileZilla/xs-recent/enhanced/20-07-2022-10-39-02"
summ <- fread(paste0(version , "/summary.csv"))
detail <- fread(paste0(version, "/detail.csv"))

# TABLE RESULTS ---------------------------------------------

TRUTH <- summ$truth %>% unique

summ[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
summ[, pname := "Constant"]
summ[, sname := "1B"]
summ[, trange := paste0("(", t_min, ", ", t_max, ")")]

detail[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
detail[, pname := "Constant"]
detail[, sname := "1B"]
detail[, trange := paste0("(", t_min, ", ", t_max, ")")]

vars <- c("q", "gamma", "eta", "xi", "estimator_type", "assay_vals")

detail <- detail[, c(vars, "estimate"), with=F]
summ <- summ[, c(vars, "q_eff"), with=F]
summ <- summ[gamma %in% c(0.0, 0.1, 0.5) & xi %in% c(0.0, 0.1) & eta %in% c(0.0, 0.1)]
detail <- detail[gamma %in% c(0.0, 0.1, 0.5) & xi %in% c(0.0, 0.1) & eta %in% c(0.0, 0.1)]
summ <- summ[assay_vals == "est"]
detail <- detail[assay_vals == "est"]

detail_ref <- detail[estimator_type == "adj" &
                       q == 1.0 & xi == 0.0 & gamma == 0.0 & eta == 0.0]
detail_ref[, q := 0.0]

summ_ref <- summ[estimator_type == "adj" &
                   q == 1.0 & xi == 0.0 & gamma == 0.0 & eta == 0.0]
summ_ref[, q := 0.0]

detail_enh <- detail[estimator_type == "eadj"]
summ_enh <- summ[estimator_type == "eadj"]

detail_plot <- rbind(detail_enh, detail_ref)
summ_plot <- rbind(summ_enh, summ_ref)

g1 <- "Correct and complete\nreporting of prior test results"
g2 <- "10% with positive prior\ntest say they've never been tested"
g3 <- "10% with positive prior\ntest say it was negative"

summ_plot[eta == 0.0 & xi == 0.0, group := g1]
summ_plot[eta == 0.1 & xi == 0.0, group := g3] # These are intentionally in this order
summ_plot[eta == 0.0 & xi == 0.1, group := g2]
detail_plot[eta == 0.0 & xi == 0.0, group := g1]
detail_plot[eta == 0.1 & xi == 0.0, group := g3]
detail_plot[eta == 0.0 & xi == 0.1, group := g2]
summ_plot[estimator_type == "adj", group := "Adjusted estimator"]
detail_plot[estimator_type == "adj", group := "Adjusted estimator"]

summ_plot[, group := factor(group,
                            levels=c("Adjusted estimator",
                                     g1,
                                     g2,
                                     g3))]
detail_plot[, group := factor(group,
                              levels=c("Adjusted estimator",
                                       g1,
                                       g2,
                                       g3))]

detail_plot[, q_group := factor(q,
                                levels=c(0, 0.5, 1.0),
                                labels=c("No testing", "50%", "100%"))]
detail_plot[estimator_type == "adj", gamma := -99]
detail_plot[, gamma_type := factor(gamma,
                                   levels=c(-99, 0, 0.1, 0.5),
                                   labels=c("No Tests",
                                            "No error",
                                            "1 month",
                                            "6 months"))]
detail_plot <- detail_plot[q != 1]
cols <- c("#000000", rev(brewer.pal(n=3,"Set1")))

pdf("~/repos/Recency-Algorithm-with-Prior-HIV-Testing/misspec-bias-b.pdf",
    height=6, width=10)
ggplot(detail_plot) +
  geom_hline(yintercept=TRUTH, color="black", linetype="dashed") +
  geom_boxplot(aes(x=group,
                   y=estimate,
                   color=gamma_type,
                   group=interaction(gamma_type, group, q_group)),
               position=position_dodge2(preserve="single",
                                        padding=0.3),
               outlier.size=1.5,
               outlier.alpha=0.5) +
  scale_color_manual(values=cols) +
  labs(y="Estimate",
       color="Measurement error in\ntiming of prior test (SD)") +
  theme(legend.position="top",
        axis.title.x=element_blank(),
        legend.box="vertical",
        legend.box.just = "left",
        legend.margin=margin(),
        legend.direction='horizontal',
        legend.justification='left') +
  guides(fill=guide_legend(nrow=2, byrow=TRUE, ncol=1))
dev.off()
