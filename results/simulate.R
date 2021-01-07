rm(list=ls())
library(XSRecency)
library(data.table)
library(xtable)
source("./results/sim-helpers.R")

# number of simulations
n_sims <- 2000

# number screened
n <- 5000

# time being considered
big_T <- 1

# prevalence
p <- 0.29

# baseline incidence
inc <- 0.032

for(type in c("constant", "linear", "exponential")){
  for(phi in c(1, 2, 3)){

    if(type == "constant"){
      inc.function <- c.incidence
      infection.function <- c.infections
      rho <- NA
    } else if(type == "linear"){
      inc.function <- l.incidence
      infection.function <- l.infections
      rho <- 0.0028
    } else {
      inc.function <- e.incidence
      infection.function <- e.infections
      rho <- 0.07
    }

    if(phi == 1){
      phi.func <- phi.character.1
      frr <- 0
      mdri <- 142
    } else if(phi == 2){
      phi.func <- phi.character.1
      frr <- 0.015
      mdri <- 142
    } else {
      phi.func <- phi.character.2
      frr <- 0.015
      mdri <- 142
    }

    sim <- simulate(n_sims=n_sims, n=n,
                    inc.function=inc.function,
                    infection.function=infection.function,
                    baseline_incidence=inc, prevalence=p, rho=rho,
                    phi.func=phi.func, frr=frr, mdri=mdri, big_T=big_T)
    assign(paste(type, phi, sep="_"), sim)
    summ <- summarize.simulation(sim)
    assign(paste(type, phi, "summary", sep="_"), summ)
    rm(sim)
    rm(summ)
  }
}

data <- data.table()

for(phi in c(1, 2, 3)){
  for(type in c("constant", "linear")){
    sim <- get(paste(type, phi, "summary", sep="_"))

    if(type == "constant") tname <- "Constant"
    if(type == "linear") tname <- "Non-constant"

    if(phi == 1) pname <- "Zero"
    if(phi == 2) pname <- "Constant"
    if(phi == 3) pname <- "Non-constant"

    row <- c(
      tname,
      pname,
      sprintf('%.4f', sim$snap_true$bias),
      sprintf('%.3f', sim$snap_true$cover),

      sprintf('%.4f', sim$snap_est$bias),
      sprintf('%.3f', sim$snap_est$cover),

      sprintf('%.4f', sim$adj_true$bias),
      sprintf('%.3f', sim$adj_true$cover),

      sprintf('%.4f', sim$adj_est$bias),
      sprintf('%.3f', sim$adj_est$cover)
    )
    data <- rbind(data, t(row))
  }
}

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 0, 0)
addtorow$command <- c("\\multicolumn{2}{c}{Setting} & \\multicolumn{4}{c}{Snapshot Estimator (1)} & \\multicolumn{4}{c}{Kassanjee Estimator (2)} \\\\\n",
                      "\\hline \\\\\n",
                      "Incidence & FRR & \\multicolumn{2}{c}{$\\mu$} & \\multicolumn{2}{c}{$\\hat{\\mu}$} & \\multicolumn{2}{c}{$\\Omega_T$, $\\beta_T$} & \\multicolumn{2}{c}{$\\hat{\\Omega_T}$, $\\hat{\\beta_T}$} \\\\\n",
                      "\\hline \\\\\n",
                      "& & Bias & Cov & Bias & Cov & Bias & Cov & Bias & Cov \\\\\n",
                      "\\hline \\\\\n")
tab <- xtable(data, align=rep("c", 11), digits=2, caption="Simulation results.")
print(tab, add.to.row = addtorow,
      include.colnames = FALSE, include.rownames = FALSE)

