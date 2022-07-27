
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)

# Get the input and output directories
args <- commandArgs(trailingOnly=TRUE)
in.dir <- args[1]

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
             "phi_norm_sd", "phi_norm_div",
             "phi_pnorm_mu", "phi_pnorm_sd", "phi_pnorm_div",
             "frr_mix_start", "frr_mix_end",
             "ext_FRR", "duong_scale", "max_FRR", "last_point",
             "pt", "t_min", "t_max", "q", "gamma", "eta", "nu",
             "xi", "mech2", "exclude_pt_bigT")){
  if(var %in% colnames(df)){
    id.vars <- c(id.vars, var)
  }
}

id.vars.nosim <- id.vars[!id.vars %in% c("sim")]
id.vars.nosim.est <- c(id.vars.nosim, "estimator")
id.vars.est <- c("adj_true_est", "adj_est_est", "eadj_true_est", "eadj_est_est")

df2 <- df[, c(id.vars, id.vars.est), with=F]
df3 <- df[, c(id.vars, "q_eff"), with=F]

df2 <- reshape2::melt(df2, id.vars=id.vars,
                           value.vars=c("adj_true_est", "adj_est_est", "eadj_true_est", "eadj_est_est"),
                           variable.name="estimator", value.name="estimate") %>% data.table

df2[, estimator := lapply(.SD, function(x) gsub("_est$", "", x)), .SDcols="estimator"]
df2[, bias := estimate - truth]

bias <- df2[, lapply(.SD, median), by=id.vars.nosim.est, .SDcols="bias"]
se <- df2[, lapply(.SD, function(x) var(x, na.rm=T) ** 0.5), by=id.vars.nosim.est, .SDcols="estimate"]

setnames(se, "estimate", "se")

results <- merge(bias, se, by=id.vars.nosim.est)
results[, mse := bias**2 + se**2]

results[, estimator_type := lapply(.SD, function(x) gsub("_est$", "", gsub("_true$", "", x))), .SDcols="estimator"]
results[, assay_vals := lapply(.SD, function(x) ifelse(grepl("true", x), "true", "est")), .SDcols="estimator"]

QEFF <- df3[, lapply(.SD, mean), by=id.vars.nosim, .SDcols="q_eff"]
results <- merge(results, QEFF, by=id.vars.nosim)

write.csv(results, paste0(in.dir, "/summary.csv"), row.names=F)

