library(magrittr)

simulate <- function(n_sims, n, inc.function, infection.function, phi.func,
                     baseline_incidence, prevalence, rho, window, frr, shadow, big_T){

  assay <- assay.properties.nsim(n_sims, phi.func=phi.func, frr=frr, window=window, shadow=shadow)
  true_frr <- true.frr(phi.func=phi.func, window=window, frr=frr, shadow=shadow)
  true_mdri <- true.mdri(phi.func=phi.func, window=window, frr=frr, shadow=shadow)

  data <- generate.data(n=n, n_sims=n_sims,
                        infection.function=infection.function,
                        phi.func=phi.func,
                        baseline_incidence=baseline_incidence,
                        prevalence=prevalence, rho=rho, frr=frr, window=window,
                        shadow=shadow)

  snap.true <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                            n=data$n, mu=window/365.25, mu_var=0)
  snap.est <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                           n=data$n, mu=assay$mu_est, mu_var=assay$mu_var)

  adj.true <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                           omega=window/365.25, omega_var=0,
                           # omega=true_mdri, omega_var=0,
                           beta=true_frr, beta_var=0,
                           big_T=big_T)
  adj.est <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                          omega=assay$omega_est, omega_var=assay$omega_var,
                          beta=assay$beta_est, beta_var=assay$beta_var,
                          big_T=big_T)

  return(list(truth=rep(baseline_incidence, n_sims),
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
