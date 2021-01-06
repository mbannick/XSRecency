rm(list=ls())
library(XSRecency)
source("./results/sim-helpers.R")

# number of simulations
n_sims <- 100

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
  }
}
