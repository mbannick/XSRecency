# TEST OUT NEW FIGURE RATHER THAN TABLE

rm(list=ls())

library(data.table)
library(xtable)
library(ggplot2)
library(viridis)
library(ggh4x)
library(RColorBrewer)
library(magrittr)
source("~/repos/XSRecency/R/phi-functions.R")
source("~/repos/XSRecency/R/data-generator.R")

# READ IN VERSIONED RESULTS ---------------------------------
version <- "~/Documents/FileZilla/xs-recent/enhanced/22-06-2022-12-14-23/"
df <- fread(paste0(version, "detail.csv"))

TRUTH <- df$truth %>% unique

df[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
df[, pname := "Constant"]
df[, sname := "1B"]

trange_levels <- c("0-4 Yrs", "1-3 Yrs", "0-2 Yrs")
df[, trange := paste0(t_min, "-", t_max, " Yrs")]
df[, trange := factor(trange, levels=trange_levels)]

ref_df <- df[estimator_type == "adj" & assay_vals == "est"
             & q == 0.2 & trange == "0-2 Yrs"]
ref_df[, q := 0.0]
ref_df[, trange := "No Testing"]

enh_df <- df[estimator_type == "eadj"]

plot_df <- rbind(enh_df, ref_df)
plot_df[, tlabs := factor(trange,
                          levels=c("No Testing", trange_levels),
                          ordered=TRUE)]
plot_df[, id := paste0(estimator_type, q, tlabs, sep="_")]
plot_df[, qlabs := factor(
  q, levels=c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
  labels=c("Adjusted\nEstimator", "20%", "40%", "60%", "80%", "100%"))
]

cols <- c("#000000", brewer.pal(n=3,"Set2"))

pdf("~/repos/Recency-Algorithm-with-Prior-HIV-Testing/correct-spec-bias.pdf",
    height=6, width=10)
ggplot(plot_df) +
  geom_hline(yintercept=TRUTH, color="black", linetype="dashed") +
  geom_boxplot(aes(x=qlabs,
                 y=estimate,
                 color=tlabs,
                 group=interaction(tlabs, q)),
               outlier.size=1.5,
               outlier.alpha=0.5,
               position=position_dodge2(preserve="single",
                                        padding=0.2)) +
  scale_color_manual(values=cols) +
  labs(x="Proportion with Tests Available",
       y="Estimate",
       color="Range of Prior Testing Times") +
  theme(legend.position="top")
dev.off()
