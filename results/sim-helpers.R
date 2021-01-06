library(magrittr)

avg.incidence <- function(inc.function, big_T,
                          baseline_incidence, rho){
  ts <- -1 * seq(0, big_T, 0.001)
  incidence <- lapply(ts, inc.function,
                      lambda_0=baseline_incidence,
                      rho=rho) %>% unlist
  return(mean(incidence))
}

simulate <- function(n_sims, n, inc.function, infection.function, phi.func,
                     baseline_incidence, prevalence, rho, mdri, frr, big_T){

  assay <- assay.properties.sim(phi.func, frr=frr, mdri=mdri)
  avg_incidence <- avg.incidence(inc.function=inc.function, big_T=big_T,
                                 baseline_incidence=baseline_incidence, rho=rho)

  data <- generate.data(n=n, n_sims=n_sims,
                        infection.function=infection.function,
                        phi.func=phi.func,
                        baseline_incidence=baseline_incidence,
                        prevalence=prevalence, rho=rho, frr=frr, mdri=mdri)

  snap.true <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                            n=data$n, mu_sim=list(est=mdri/365.25, var=0))
  snap.est <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                           n=data$n, mu_sim=assay$mu)

  adj.true <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                           omega_sim=list(est=mdri/365.25, var=0),
                           beta_sim=list(est=frr, var=0), big_T=big_T)
  adj.est <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                          omega_sim=assay$omega,
                          beta_sim=assay$beta, big_T=big_T)

  return(list(truth=rep(avg_incidence, n_sims),
         snap_true_est=snap.true$est,
         snap_true_var=snap.true$var,
         snap_est_est=snap.est$est,
         snap_est_var=snap.est$var,
         adj_true_est=adj.true$est,
         adj_true_var=adj.true$var,
         adj_est_est=adj.est$est,
         adj_est_var=adj.est$var))
}

summarize <- function(truth, estimates, variance){

  # bias
  bias <- mean(estimates - truth)

  # standard error (true)
  se <- var(estimates) ** 0.5

  # standard error (estimated)
  see <- mean(variance ** 0.5)

  # coverage
  width <- qnorm(0.975) * variance ** 0.5

  cover <- ((estimates - width) < truth) &
    ((estimates + width) > truth)
  cover <- mean(cover)

  # truth should all be the same
  return(list(truth=truth[1], bias=bias, se=se, see=see, cover=cover))
}

summarize.simulation <- function(sim){

  snap_true <- summarize(sim$truth, sim$snap_true_est, sim$snap_true_var)
  snap_est <- summarize(sim$truth, sim$snap_est_est, sim$snap_est_var)
  adj_true <- summarize(sim$truth, sim$adj_true_est, sim$adj_true_var)
  adj_est <- summarize(sim$truth, sim$adj_est_est, sim$adj_est_var)

  return(list(
    snap_true=snap_true,
    snap_est=snap_est,
    adj_true=adj_true,
    adj_est=adj_est
  ))
}
