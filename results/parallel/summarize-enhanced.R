rm(list=ls())
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(pbapply)

# Get the input and output directories
args <- commandArgs(trailingOnly=TRUE)
version <- args[1]
version <- "15-02-2023-14-27-56"
in.dir <- paste0("~/Documents/FileZilla/xs-recent/enhanced/", version)

# Read in files
f <- list.files(in.dir, full.names=T)
f <- f[!grepl("summary", f)]
f <- f[!grepl("detail", f)]
f <- f[!grepl("README.md", f)]
df <- pblapply(f, fread) %>% rbindlist(fill=T)
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

id.vars.nosim <- id.vars[!id.vars %in% c("sim", "seed")]
id.vars.nosim.est <- c(id.vars.nosim, "estimator")

df[, adj_true_var_rob := adj_est_var]
df[, adj_true_var_asm := adj_est_var]

df[, adj_est_var_rob := adj_est_var]
df[, adj_est_var_asm := adj_est_var]

df[, adj_true_est := adj_est_est]
df[, eadj_true_est := eadj_est_est]

df[, eadj_true_var_rob := eadj_est_var_rob]
df[, eadj_true_var_asm := eadj_est_var_asm]

est.variables <- c("adj_true_est", "adj_est_est", "eadj_true_est", "eadj_est_est")
var.rob.variables <- c("adj_true_var_rob", "adj_est_var_rob", "eadj_true_var_rob", "eadj_est_var_rob")
var.asm.variables <- c("adj_true_var_asm", "adj_est_var_asm", "eadj_true_var_asm", "eadj_est_var_asm")

df <- df[, c(id.vars, est.variables, var.rob.variables, var.asm.variables), with=FALSE]

estimate <- reshape2::melt(df[, c(id.vars, est.variables), with=F], id.vars=id.vars,
                           value.vars=est.variables,
                           variable.name="estimator", value.name="estimate") %>% data.table
variance.rob <- reshape2::melt(df[, c(id.vars, var.rob.variables), with=F], id.vars=id.vars,
                           value.vars=var.rob.variables,
                           variable.name="estimator", value.name="variance_rob") %>% data.table
variance.asm <- reshape2::melt(df[, c(id.vars, var.asm.variables), with=F], id.vars=id.vars,
                               value.vars=var.asm.variables,
                               variable.name="estimator", value.name="variance_asm") %>% data.table
estimate[, estimator := lapply(.SD, function(x) gsub("_est$", "", x)), .SDcols="estimator"]
variance.rob[, estimator := lapply(.SD, function(x) gsub("_var_rob$", "", x)), .SDcols="estimator"]
variance.asm[, estimator := lapply(.SD, function(x) gsub("_var_asm$", "", x)), .SDcols="estimator"]

df2 <- merge(estimate, variance.rob, by=c(id.vars, "estimator"))
df2 <- merge(df2, variance.asm, by=c(id.vars, "estimator"))
df2[, bias := estimate - truth]

df2[, width_rob := qnorm(0.975) * variance_rob ** 0.5]
df2[, lower_rob := estimate - width_rob]
df2[, upper_rob := estimate + width_rob]
df2[, cover_rob := (truth < upper_rob) & (truth > lower_rob)]

df2[, width_asm := qnorm(0.975) * variance_asm ** 0.5]
df2[, lower_asm := estimate - width_asm]
df2[, upper_asm := estimate + width_asm]
df2[, cover_asm := (truth < upper_asm) & (truth > lower_asm)]

bias <- df2[, lapply(.SD, median), by=id.vars.nosim.est, .SDcols="bias"]
se <- df2[, lapply(.SD, function(x) var(x, na.rm=T)**0.5), by=id.vars.nosim.est, .SDcols="estimate"]
see_rob <- df2[, lapply(.SD, function(x) median(x**0.5, na.rm=T)), by=id.vars.nosim.est, .SDcols="variance_rob"]
see_asm <- df2[, lapply(.SD, function(x) median(x**0.5, na.rm=T)), by=id.vars.nosim.est, .SDcols="variance_asm"]
cover_rob <- df2[, lapply(.SD, mean, na.rm=T), by=id.vars.nosim.est, .SDcols="cover_rob"]
cover_asm <- df2[, lapply(.SD, mean, na.rm=T), by=id.vars.nosim.est, .SDcols="cover_asm"]

setnames(se, "estimate", "se")
setnames(see_rob, "variance_rob", "see_rob")
setnames(see_asm, "variance_asm", "see_asm")

results <- merge(bias, se, by=id.vars.nosim.est)
results <- merge(results, see_rob, by=id.vars.nosim.est)
results <- merge(results, see_asm, by=id.vars.nosim.est)
results <- merge(results, cover_rob, by=id.vars.nosim.est)
results <- merge(results, cover_asm, by=id.vars.nosim.est)
results[, mse := bias**2 + se**2]

results[, estimator_type := lapply(.SD, function(x) gsub("_est$", "", gsub("_true$", "", x))), .SDcols="estimator"]
results[, assay_vals := lapply(.SD, function(x) ifelse(grepl("true", x), "true", "est")), .SDcols="estimator"]

write.csv(results, paste0(in.dir, "/summary.csv"), row.names=F)

