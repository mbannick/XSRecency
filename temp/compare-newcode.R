library(data.table)

old <- "15-12-2022-17-11-12/"
new <- "21-02-2023-11-41-02/"

folder <- "~/Documents/FileZilla/xs-recent/enhanced/"
file <- "results-1000_0________FALSE___0_0_0_0_FALSE_FALSE_5000_0.29_0.032_12_2_TRUE_constant_101_194_2_TRUE_0_2_0.2_100.csv"

df.old <- fread(paste0(folder, old, file))
df.new <- fread(paste0(folder, new, file))

par(mfrow=c(3, 4))
for(variable in c("omega_est", "beta_est", "n_n", "n_r", "n_p",
                  "n_r_pt", "num_beta", "den_omega", "den_beta", "q_eff")){
  cat(variable, "\n")
  oldvar <- df.old[[variable]]
  newvar <- df.new[[variable]]
  xmin <- min(oldvar, newvar)
  xmax <- max(oldvar, newvar)
  hist(oldvar, col=rgb(0,0,1,1/4), xlim=c(xmin, xmax), main=variable)  # first histogram
  hist(newvar, col=rgb(1,0,0,1/4), add=T)
}
