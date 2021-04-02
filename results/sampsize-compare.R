rm(list=ls())

library(data.table)
library(ggplot2)

files <- c(
  "~/Documents/FileZilla/xs-recent/29-03-21-16-2000/summary.csv",
  "~/Documents/FileZilla/xs-recent/29-03-21-18-5000/summary.csv",
  "~/Documents/FileZilla/xs-recent/29-03-21-16-10000/summary.csv",
  "~/Documents/FileZilla/xs-recent/29-03-21-19-15000/summary.csv",
  "~/Documents/FileZilla/xs-recent/29-03-21-19-20000/summary.csv"
)

dfs <- lapply(files, fread) %>% rbindlist

dfs[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]
dfs[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu), pname := "Zero"]
dfs[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu), pname := "Constant"]
dfs[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu), pname := "Non-constant"]
dfs[window == 101, sname := "1 A-C"]
dfs[window == 248, sname := "2 A-C"]
dfs[pname == "Zero" & sname == "1 A-C", assay := "1A"]
dfs[pname == "Constant" & sname == "1 A-C", assay := "1B"]
dfs[pname == "Non-constant" & sname == "1 A-C", assay := "1C"]
dfs[pname == "Zero" & sname == "2 A-C", assay := "2A"]
dfs[pname == "Constant" & sname == "2 A-C", assay := "2B"]
dfs[pname == "Non-constant" & sname == "2 A-C", assay := "2C"]

id.cols <- c("assay", "estimator", "itype", "n")
bias.cols <- c(id.cols, "bias")
cover.cols <- c(id.cols, "cover")
see.cols <- c(id.cols, "see")
se.cols <- c(id.cols, "se")

bias <- dfs[, bias.cols, with=F]
cover <- dfs[, cover.cols, with=F]
se <- dfs[, se.cols, with=F]
see <- dfs[, see.cols, with=F]

pdf("coverage-samplesize.pdf", height=10, width=7)
ggplot(data=cover, aes(x=n, y=cover, color=estimator, group=estimator)) + geom_line() +
  facet_grid(assay~itype)
dev.off()

pdf("bias-samplesize.pdf", height=10, width=7)
ggplot(data=bias, aes(x=n, y=bias, color=estimator, group=estimator)) + geom_line() +
  facet_grid(assay~itype)
dev.off()

pdf("se-samplesize.pdf", height=10, width=7)
ggplot(data=se, aes(x=n, y=se, color=estimator, group=estimator)) + geom_line() +
  facet_grid(assay~itype)
dev.off()

pdf("see-samplesize.pdf", height=10, width=7)
ggplot(data=see, aes(x=n, y=see, color=estimator, group=estimator)) + geom_line() +
  facet_grid(assay~itype)
dev.off()
