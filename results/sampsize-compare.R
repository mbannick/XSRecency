rm(list=ls())

library(data.table)
library(ggplot2)
library(magrittr)

files <- c(
  "~/Documents/FileZilla/xs-recent/nsamps/18-05-21-16/summary.csv"
)

dfs <- lapply(files, fread) %>% rbindlist

dfs[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
dfs[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Zero"]
dfs[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Constant"]
dfs[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Non-constant"]
dfs[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & !is.na(phi_pnorm_mu), pname := "Increasing"]
dfs[window == 101, sname := "1 A-D"]
dfs[window == 248, sname := "2 A-D"]

dfs[pname == "Zero" & sname == "1 A-D", assay := "1A"]
dfs[pname == "Constant" & sname == "1 A-D", assay := "1B"]
dfs[pname == "Non-constant" & sname == "1 A-D", assay := "1C"]
dfs[pname == "Increasing" & sname == "1 A-D", assay := "1D"]
dfs[pname == "Zero" & sname == "2 A-D", assay := "2A"]
dfs[pname == "Constant" & sname == "2 A-D", assay := "2B"]
dfs[pname == "Non-constant" & sname == "2 A-D", assay := "2C"]
dfs[pname == "Increasing" & sname == "2 A-D", assay := "2D"]

id.cols <- c("assay", "estimator", "itype", "n")
bias.cols <- c(id.cols, "bias")
cover.cols <- c(id.cols, "cover")
see.cols <- c(id.cols, "see")
se.cols <- c(id.cols, "se")

bias <- dfs[, bias.cols, with=F]
cover <- dfs[, cover.cols, with=F]
se <- dfs[, se.cols, with=F]
see <- dfs[, see.cols, with=F]

cover <- cover[estimator %in% c("snap_est", "adj_est")]
cover[estimator == "snap_est", estimator := "Snapshot"]
cover[estimator == "adj_est", estimator := "Adjusted"]

pdf("coverage-samplesize.pdf", height=6, width=10)
ggplot(data=cover, aes(x=n, y=cover, color=estimator,
                       shape=itype)) +
  geom_hline(yintercept=0.95, alpha=0.5) +
  geom_line(linetype='dashed') + geom_point() +
  facet_wrap(~ assay, nrow=2) + theme_minimal() +
  labs(x="Trial Sample Size",
       y="Coverage",
       shape="Incidence",
       color="Estimator")
dev.off()

# pdf("bias-samplesize.pdf", height=10, width=7)
# ggplot(data=bias, aes(x=n, y=bias, color=estimator, group=estimator)) + geom_line() +
#   facet_grid(assay~itype)
# dev.off()
#
# pdf("se-samplesize.pdf", height=10, width=7)
# ggplot(data=se, aes(x=n, y=se, color=estimator, group=estimator)) + geom_line() +
#   facet_grid(assay~itype)
# dev.off()
#
# pdf("see-samplesize.pdf", height=10, width=7)
# ggplot(data=see, aes(x=n, y=see, color=estimator, group=estimator)) + geom_line() +
#   facet_grid(assay~itype)
# dev.off()
