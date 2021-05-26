# THIS HAS SENSITIVITY ANALYSES FOR INTEGRATING TO 8 VERSUS 12

rm(list=ls())
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbrace)

# Get the input and output directories
args <- commandArgs(trailingOnly=TRUE)
# in.dir <- "~/Documents/FileZilla/xs-recent/14-04-21-19/" # integrate to 12
# in.dir <- "~/Documents/FileZilla/xs-recent/14-04-21-20/" # integrate to 8
in.dir <- "~/Documents/FileZilla/xs-recent/12-05-21-10/"

# Read in files
f <- list.files(in.dir, full.names=T)
f <- f[!grepl("summary", f)]
f <- f[!grepl("detail", f)]
f <- f[!grepl("README.md", f)]
df <- lapply(f, fread) %>% rbindlist(fill=T)
df[, V1 := NULL]

id.vars <- c("truth", "n_sims", "sim", "seed", "n", "p", "inc", "tau", "bigT", "itype",
             "window", "shadow")

for(var in c("rho", "phi_frr", "phi_tfrr", "phi_norm_mu",
             "phi_norm_sd", "phi_norm_div", "frr_mix_start", "frr_mix_end",
             "ext_FRR", "duong_scale", "add_unif")){
  if(var %in% colnames(df)){
    id.vars <- c(id.vars, var)
  }
}

df[, tname := ifelse(itype == "constant", "Constant", ifelse(itype == "linear", "Linear", "Exponential"))]

df[is.na(phi_frr) & is.na(phi_tfrr) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Zero"]
df[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Constant"]
df[(!is.na(phi_frr) | !is.na(phi_tfrr)) & !is.na(phi_norm_mu) & is.na(phi_pnorm_mu), pname := "Non-constant"]
df[(!is.na(phi_frr) | !is.na(phi_tfrr)) & is.na(phi_norm_mu) & !is.na(phi_pnorm_mu), pname := "Increasing"]

df[window == 101, sname := "1 A-D"]
df[window == 248, sname := "2 A-D"]

df[pname == "Zero" & sname == "1 A-D", assay := "1A"]
df[pname == "Constant" & sname == "1 A-D", assay := "1B"]
df[pname == "Non-constant" & sname == "1 A-D", assay := "1C"]
df[pname == "Increasing" & sname == "1 A-D", assay := "1D"]
df[pname == "Zero" & sname == "2 A-D", assay := "2A"]
df[pname == "Constant" & sname == "2 A-D", assay := "2B"]
df[pname == "Non-constant" & sname == "2 A-D", assay := "2C"]
df[pname == "Increasing" & sname == "2 A-D", assay := "2D"]

df[pname == "Zero" & sname == "1 A-D", assay := "A"]
df[pname == "Constant" & sname == "1 A-D", assay := "B"]
df[pname == "Non-constant" & sname == "1 A-D", assay := "C"]
df[pname == "Increasing" & sname == "1 A-D", assay := "D"]
df[pname == "Zero" & sname == "2 A-D", assay := "A"]
df[pname == "Constant" & sname == "2 A-D", assay := "B"]
df[pname == "Non-constant" & sname == "2 A-D", assay := "C"]
df[pname == "Increasing" & sname == "2 A-D", assay := "D"]

df <- df[abs(adj_est_est) < 10]

pdf("adjusted-beta.pdf", height=6, width=10)
ggplot() + geom_point(data=df, aes(x=beta_est, y=adj_est_est), alpha=0.2) +
  facet_nested(sname ~ assay + tname) +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  geom_hline(yintercept=unique(df$truth),
             color='black', linetype='dashed') +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(color="Estimator") +
  ylab("Adjusted Estimate") + xlab("Beta")
dev.off()

pdf("adjusted-omega.pdf", height=6, width=10)
ggplot() + geom_point(data=df, aes(x=omega_est, y=adj_est_est), alpha=0.2) +
  facet_nested(sname ~ assay + tname) +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  geom_hline(yintercept=unique(df$truth),
             color='black', linetype='dashed') +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(color="Estimator") +
  ylab("Adjusted Estimate") + xlab("Omega")
dev.off()

pdf("adjusted-n-p.pdf", height=6, width=10)
ggplot() + geom_point(data=df, aes(x=beta_est, y=n_r), alpha=0.2) +
  facet_nested(sname ~ assay + tname) +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  geom_hline(yintercept=unique(df$truth),
             color='black', linetype='dashed') +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(color="Estimator") #
  # ylab("Adjusted Estimate") + xlab("Omega")
dev.off()

pdf("beta-hist.pdf", height=6, width=10)
ggplot() + geom_histogram(data=df, aes(y=beta_est), bins=40) +
  facet_nested(sname ~ assay + tname) +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(color="Estimator") +
  ylab("Beta Estimate") + xlab("")
dev.off()

pdf("mu-hist.pdf", height=6, width=10)
ggplot() + geom_histogram(data=df, aes(y=mu_est), bins=40) +
  facet_nested(sname ~ assay + tname) +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(color="Estimator") +
  ylab("Mu Estimate") + xlab("")
dev.off()

pdf("omega-hist.pdf", height=6, width=10)
ggplot() + geom_histogram(data=df, aes(y=omega_est), bins=40) +
  facet_nested(sname ~ assay + tname) +
  scale_colour_manual(values=c("#ffa75e", "#5ea4ff")) +
  scale_fill_manual(values=c("#ffa75e", "#5ea4ff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(color="Estimator") +
  ylab("Omega Estimate") + xlab("")
dev.off()
