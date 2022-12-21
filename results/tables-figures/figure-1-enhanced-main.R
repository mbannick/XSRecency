# FIGURE 1 ENHANCED V. ADJUSTED ESTIMATOR

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
version <- "~/Documents/FileZilla/xs-recent/enhanced/15-12-2022-17-11-12/"
df <- fread(paste0(version, "detail.csv"))
summ <- fread(paste0(version, "summary.csv"))
adj.mse <- summ[estimator_type == "adj" & q == 1 & t_min == 0 & t_max == 2, mse]
summ[, rmse := (1 - (mse / adj.mse))*100]
summ <- summ[assay_vals == "est" & (estimator_type == "eadj" | (estimator_type == "adj" & t_min == 0 & t_max == 2 & q == 1))]

TRUTH <- df$truth %>% unique

df[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
df[, pname := "Constant"]
df[, sname := "1B"]

trange_levels <- c("2-4 Yrs", "0-4 Yrs", "0-2 Yrs")
df[, trange := paste0(t_min, "-", t_max, " Yrs")]
df[, trange := factor(trange, levels=trange_levels)]
summ[, trange := paste0(t_min, "-", t_max, " Yrs")]
summ[, trange := factor(trange, levels=trange_levels)]

# Get reference values for the adjusted estimator
ref_df <- df[estimator_type == "adj" & assay_vals == "est"
             & q == 0.2 & trange == "0-2 Yrs"]
ref_df[, q := 0.0]
ref_df[, trange := "No Testing"]

enh_df <- df[estimator == "eadj_est"]

plot_df <- rbind(enh_df, ref_df)
plot_df[, tlabs := factor(trange,
                          levels=c("No Testing", trange_levels),
                          ordered=TRUE)]
plot_df[, id := paste0(estimator_type, q, tlabs, sep="_")]
plot_df[, qlabs := factor(
  q, levels=c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
  labels=c("Adjusted\nEstimator", "20%", "40%", "60%", "80%", "100%"))
]
summ[, qlabs := factor(
  q, levels=c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
  labels=c("Adjusted\nEstimator", "20%", "40%", "60%", "80%", "100%"))
]
summ[, tlabs := factor(trange,
                          levels=c("No Testing", trange_levels),
                          ordered=TRUE)]
summ[estimator_type == "adj", qlabs := "Adjusted\nEstimator"]
summ[estimator_type == "adj", tlabs := "No Testing"]
summ[, rmse_labs := paste0(sprintf("%.0f", rmse), "%")]

cols <- c("#000000", brewer.pal(n=3,"Set2"))

pdf("~/repos/Recency-Algorithm-with-Prior-HIV-Testing/correct-spec-bias-2.pdf",
    height=7, width=10)
fig <- ggplot(plot_df) +
  geom_hline(yintercept=TRUTH, color="black", linetype="dashed") +
  geom_boxplot(aes(x=qlabs,
                 y=estimate,
                 color=tlabs,
                 group=interaction(tlabs, q)),
               outlier.size=1.5,
               outlier.alpha=0.5,
               position=position_dodge2(preserve="single",
                                        width=0.8, padding=0.2)) +
  scale_color_manual(values=cols) +
  labs(x="Proportion with Tests Available",
       y="Estimate",
       color="Range of Prior Testing Times") +
  theme(legend.position="top")

dt <- ggplot(summ) +
  geom_richtext(size=3.5, fontface="bold", aes(
            group=interaction(tlabs, q),
            x=qlabs, y=factor(""),
            color=tlabs,
            label=rmse_labs),
            show.legend=FALSE,
            inherit.aes=TRUE,
            position=position_dodge2(preserve="single", width=0.8, padding=0.2)
  ) +
  scale_color_manual(values=cols) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, vjust=0.3,
                                  size = 10, face = "italic")) +
  ggtitle("% Reduction in MSE Compared to Adjusted Estimator") +
  theme(panel.grid.major=element_blank(),
        panel.border=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  xlab(NULL) +
  ylab(NULL)

# Combine plot and table
plot_cmbd <- append_table(
  plot=fig,
  table=dt,
  extract.legend=FALSE
)

# Draw in RStudio viewer
grid.newpage()
grid.draw(plot_cmbd)

dev.off()
