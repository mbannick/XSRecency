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

summ_plot[eta == 0.0 & xi == 0.0, group := "Correct and complete\nreporting of prior test results"]
summ_plot[eta == 0.1 & xi == 0.0, group := "10% with pos. prior\ntest report neg."]
summ_plot[eta == 0.0 & xi == 0.1, group := "10% with pos. prior\ntest do not report"]
detail_plot[eta == 0.0 & xi == 0.0, group := "Correct and complete\nreporting of prior test results"]
detail_plot[eta == 0.1 & xi == 0.0, group := "10% with pos. prior\ntest report neg."]
detail_plot[eta == 0.0 & xi == 0.1, group := "10% with pos. prior\ntest do not report"]
summ_plot[estimator_type == "adj", group := "Adjusted estimator"]
detail_plot[estimator_type == "adj", group := "Adjusted estimator"]

summ_plot[, group := factor(group,
                            levels=c("Adjusted estimator",
                                     "Correct and complete\nreporting of prior test results",
                                     "10% with pos. prior\ntest do not report",
                                     "10% with pos. prior\ntest report neg."))]
detail_plot[, group := factor(group,
                              levels=c("Adjusted estimator",
                                       "Correct and complete\nreporting of prior test results",
                                       "10% with pos. prior\ntest do not report",
                                       "10% with pos. prior\ntest report neg."))]

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
cols <- c("#000000", rev(brewer.pal(n=3,"Set1")))

pdf("~/repos/Recency-Algorithm-with-Prior-HIV-Testing/misspec-bias-b.pdf",
    height=6, width=10)
ggplot(detail_plot) +
  geom_hline(yintercept=TRUTH, color="black", linetype="dashed") +
  geom_boxplot_pattern(aes(x=group,
                   y=estimate,
                   color=gamma_type,
                   pattern=q_group,
                   group=interaction(gamma_type, group, q_group)),
               position=position_dodge2(preserve="single",
                                        padding=0.3),
               pattern_color="black",
               pattern_fill="black",
               pattern_angle=45,
               pattern_size=0.1,
               pattern_spacing=0.025,
               outlier.size=1.5,
               outlier.alpha=0.5) +
  scale_color_manual(values=cols) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch')) +
  labs(y="Estimate",
       color="Measurement error in\ntiming of prior test (SD)",
       pattern="Percent of people\nwith prior tests") +
  theme(legend.position="top",
        axis.title.x=element_blank(),
        legend.box="vertical",
        legend.box.just = "left",
        legend.margin=margin(),
        legend.direction='horizontal',
        legend.justification='left') +
  guides(color=guide_legend(override.aes=list(pattern=c('none', 'none', 'none', 'none'))),
         fill=guide_legend(nrow=2, byrow=TRUE, ncol=1))
dev.off()
