rm(list=ls())
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(pbapply)

# Get the input and output directories
args <- commandArgs(trailingOnly=TRUE)
in.dir <- args[1]
in.dir <- "~/Documents/FileZilla/xs-recent/enhanced/15-12-2022-17-11-12/"

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

df[, adj_true_var_rob := adj_true_var]
df[, adj_true_var_asm := adj_true_var]

df[, adj_est_var_rob := adj_est_var]
df[, adj_est_var_asm := adj_est_var]

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
df2[, width := qnorm(0.975) * variance ** 0.5]
df2[, lower := estimate - width]
df2[, upper := estimate + width]
df2[, cover := (truth < upper) & (truth > lower)]

df2[, width := qnorm(0.975) * variance**0.5]
df2[, lower := estimate - width]
df2[, upper := estimate + width]
df2[, cover := (truth < upper) & (truth > lower)]

bias <- df2[, lapply(.SD, median), by=id.vars.nosim.est, .SDcols="bias"]
se <- df2[, lapply(.SD, function(x) var(x, na.rm=T)**0.5), by=id.vars.nosim.est, .SDcols="estimate"]
see <- df2[, lapply(.SD, function(x) median(x**0.5, na.rm=T)), by=id.vars.nosim.est, .SDcols="variance"]
cover <- df2[, lapply(.SD, mean, na.rm=T), by=id.vars.nosim.est, .SDcols="cover"]

setnames(se, "estimate", "se")
setnames(see, "variance", "see")

results <- merge(bias, se, by=id.vars.nosim.est)
results <- merge(results, see, by=id.vars.nosim.est)
results <- merge(results, cover, by=id.vars.nosim.est)
results[, mse := bias**2 + se**2]

results[, estimator_type := lapply(.SD, function(x) gsub("_est$", "", gsub("_true$", "", x))), .SDcols="estimator"]
results[, assay_vals := lapply(.SD, function(x) ifelse(grepl("true", x), "true", "est")), .SDcols="estimator"]

write.csv(results, paste0(in.dir, "/summary.csv"), row.names=F)

