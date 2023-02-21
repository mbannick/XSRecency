
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)

# Get the input and output directories
args <- commandArgs(trailingOnly=TRUE)
version <- args[1]
version <- "15-02-2023-14-27-56"
version <- "15-12-2022-17-11-12"
in.dir <- paste0("~/Documents/FileZilla/xs-recent/enhanced/", version)

# Read in files
f <- list.files(in.dir, full.names=T)
f <- f[!grepl("summary", f)]
f <- f[!grepl("detail", f)]
f <- f[!grepl("README.md", f)]
df <- lapply(f, fread) %>% rbindlist(fill=T)
df[, V1 := NULL]

id.vars <- c("truth", "n_sims", "sim", "n", "p", "inc", "tau", "bigT", "itype",
             "window", "shadow", "seed")

for(var in c("rho", "phi_frr", "phi_tfrr", "phi_norm_mu",
             "phi_norm_sd", "phi_norm_div",
             "phi_pnorm_mu", "phi_pnorm_sd", "phi_pnorm_div",
             "frr_mix_start", "frr_mix_end",
             "ext_FRR", "duong_scale", "max_FRR", "last_point",
             "pt", "t_min", "t_max", "q", "q_eff", "gamma", "eta", "nu",
             "xi", "mech2", "exclude_pt_bigT")){
  if(var %in% colnames(df)){
    id.vars <- c(id.vars, var)
  }
}

df[, adj_true_est := adj_est_est]
df[, eadj_true_est := eadj_est_est]


id.vars.est <- c("adj_true_est", "adj_est_est", "eadj_true_est", "eadj_est_est")

df2 <- df[, c(id.vars, id.vars.est), with=F]

df2 <- reshape2::melt(df2, id.vars=id.vars,
                      value.vars=c("adj_true_est", "adj_est_est", "eadj_true_est", "eadj_est_est"),
                      variable.name="estimator", value.name="estimate") %>% data.table

df2[, estimator := lapply(.SD, function(x) gsub("_est$", "", x)), .SDcols="estimator"]
df2[, bias := estimate - truth]

df2[, estimator_type := lapply(.SD, function(x) gsub("_est$", "", gsub("_true$", "", x))), .SDcols="estimator"]
df2[, assay_vals := lapply(.SD, function(x) ifelse(grepl("true", x), "true", "est")), .SDcols="estimator"]

write.csv(df2, paste0(in.dir, "/detail.csv"), row.names=F)

