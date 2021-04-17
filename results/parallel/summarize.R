
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
             "phi_norm_sd", "phi_norm_div", "frr_mix_start", "frr_mix_end",
             "ext_FRR", "duong_scale", "max_FRR", "last_point")){
  if(var %in% colnames(df)){
    id.vars <- c(id.vars, var)
  }
}

id.vars.nosim <- id.vars[!id.vars %in% c("sim")]
id.vars.nosim.est <- c(id.vars.nosim, "estimator")

expected <- df[, c(id.vars, "snap_bias_exp", "adj_bias_exp"), with=F]
expected <- expected[, lapply(.SD, mean), .SDcols=c("snap_bias_exp", "adj_bias_exp"), by=id.vars.nosim]
expected <- reshape2::melt(expected, id.vars=id.vars.nosim, value.vars=c("snap_bias_exp", "adj_bias_exp"),
                           variable.name="estimator_type", value.name="expected_bias") %>% data.table
expected[, estimator_type := lapply(.SD, function(x) gsub("_bias_exp$", "", x)), .SDcols="estimator_type"]

df[, snap_bias_exp := NULL]
df[, adj_bias_exp := NULL]

estimate <- reshape2::melt(df, id.vars=id.vars,
                           value.vars=c("snap_true_est", "snap_est_est", "adj_true_est", "adj_est_est"),
                           variable.name="estimator", value.name="estimate") %>% data.table
variance <- reshape2::melt(df, id.vars=id.vars,
                           value.vars=c("snap_true_var", "snap_est_var", "adj_true_var", "adj_est_var"),
                           variable.name="estimator", value.name="variance") %>% data.table

estimate[, estimator := lapply(.SD, function(x) gsub("_est$", "", x)), .SDcols="estimator"]
variance[, estimator := lapply(.SD, function(x) gsub("_var$", "", x)), .SDcols="estimator"]

detail <- merge(estimate, variance, by=c(id.vars, "estimator"))
detail[, bias := estimate - truth]

detail[, width := qnorm(0.975) * variance ** 0.5]
detail[, lower := estimate - width]
detail[, upper := estimate + width]
detail[, cover := (truth < upper) & (truth > lower)]

bias <- detail[, lapply(.SD, median), by=id.vars.nosim.est, .SDcols="bias"]
se <- detail[, lapply(.SD, function(x) var(x) ** 0.5), by=id.vars.nosim.est, .SDcols="estimate"]
see <- detail[, lapply(.SD, function(x) mean(x**0.5)), by=id.vars.nosim.est, .SDcols="variance"]
cover <- detail[, lapply(.SD, mean), by=id.vars.nosim.est, .SDcols="cover"]

setnames(se, "estimate", "se")
setnames(see, "variance", "see")

results <- merge(bias, se, by=id.vars.nosim.est)
results <- merge(results, see, by=id.vars.nosim.est)
results <- merge(results, cover, by=id.vars.nosim.est)

results[, estimator_type := lapply(.SD, function(x) gsub("_est$", "", gsub("_true$", "", x))), .SDcols="estimator"]
results <- merge(results, expected, by=c(id.vars.nosim, "estimator_type"))

write.csv(results, paste0(in.dir, "/summary.csv"), row.names=F)
write.csv(detail, paste0(in.dir, "/detail.csv"), row.names=F)
