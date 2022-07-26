setwd("~/repos/XSRecency/")
source("R/data-generator.R")
source("R/external-study-sim.R")
source("R/external-study-data.R")
source("R/estimators.R")
source("R/estimators-pt.R")

library(magrittr)

simulate <- function(n_sims, n, inc.function, infection.function, phi.func,
                     baseline_incidence, prevalence, rho, bigT, tau,
                     ext_FRR, ext_df=NULL, max_FRR=NULL, last_point=FALSE){

  # Generate trial data
  data <- generate.data(n=n, n_sims=n_sims,
                        infection.function=infection.function,
                        phi.func=phi.func,
                        baseline_incidence=baseline_incidence,
                        prevalence=prevalence, rho=rho)

  # Get assay parameters simulation based on external data simulation
  assay <- assay.properties.nsim(n_sims, phi.func=phi.func, bigT=bigT, tau=tau,
                                 ext_FRR=ext_FRR, ext_df=ext_df, max_FRR=max_FRR,
                                 last_point=last_point)

  # Calculate true assay parameters
  true_frr <- true.frr(phi.func=phi.func, bigT=bigT, tau=tau)
  true_mdri <- true.window.mdri(phi.func=phi.func, maxT=bigT)
  true_window <- true.window.mdri(phi.func=phi.func, maxT=tau)

  # Get expected bias of the estimators
  exp_bias_snap <- snap.bias(phi.func=phi.func, inc.func=inc.function,
                             inc.0=baseline_incidence, rho=rho, tau=tau)
  exp_bias_adj <- adj.bias(phi.func=phi.func, inc.func=inc.function,
                           inc.0=baseline_incidence, rho=rho, bigT=bigT, tau=tau)

  # Compute estimates for snapshot
  snap.true <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                            n=data$n, mu=true_window, mu_var=0)
  snap.est <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                           n=data$n, mu=assay$mu_est, mu_var=assay$mu_var)

  # Compute estimates for adjusted
  adj.true <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                           omega=true_mdri, omega_var=0,
                           beta=true_frr, beta_var=0,
                           big_T=bigT)
  adj.est <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                          omega=assay$omega_est, omega_var=assay$omega_var,
                          beta=assay$beta_est, beta_var=assay$beta_var,
                          big_T=bigT)

  return(list(truth=rep(baseline_incidence, n_sims),
         snap_bias_exp=rep(exp_bias_snap, n_sims),
         snap_true_est=snap.true$est,
         snap_true_var=snap.true$var,
         snap_est_est=snap.est$est,
         snap_est_var=snap.est$var,
         adj_bias_exp=rep(exp_bias_adj, n_sims),
         adj_true_est=adj.true$est,
         adj_true_var=adj.true$var,
         adj_est_est=adj.est$est,
         adj_est_var=adj.est$var,
         mu_est=assay$mu_est,
         mu_var=assay$mu_var,
         omega_est=assay$omega_est,
         omega_var=assay$omega_var,
         beta_est=assay$beta_est,
         beta_var=assay$beta_var,
         n_n=data$n_n,
         n_r=data$n_r,
         n_p=data$n_p,
         n=data$n
         ))
}

simulate.pt <- function(n_sims, n, infection.function, phi.func,
                     baseline_incidence, prevalence, rho, bigT, tau,
                     ext_FRR, ext_df=NULL, max_FRR=NULL, last_point=FALSE,
                     ptest.dist=NULL, ptest.prob=1.0,
                     t_range=NULL, t_noise=NULL,
                     d_misrep=0.0, q_misrep=0.0, p_misrep=0.0,
                     ptest.dist2=NULL,
                     exclude_pt_bigT=FALSE){

  # Generate trial data
  data <- generate.data(n=n, n_sims=n_sims,
                        infection.function=infection.function,
                        phi.func=phi.func,
                        baseline_incidence=baseline_incidence,
                        prevalence=prevalence, rho=rho,
                        bigT=bigT,
                        ptest.dist=ptest.dist, ptest.prob=ptest.prob,
                        t_range=t_range, t_noise=t_noise,
                        d_misrep=d_misrep, q_misrep=q_misrep,
                        p_misrep=p_misrep,
                        ptest.dist2=ptest.dist2,
                        exclude_pt_bigT=exclude_pt_bigT)

  # Get assay parameters simulation based on external data simulation
  assay <- assay.properties.nsim(n_sims, phi.func=phi.func, bigT=bigT, tau=tau,
                                 ext_FRR=ext_FRR, ext_df=ext_df, max_FRR=max_FRR,
                                 last_point=last_point)

  # # Calculate true assay parameters
  true_frr <- true.frr(phi.func=phi.func, bigT=bigT, tau=tau)
  true_mdri <- true.window.mdri(phi.func=phi.func, maxT=bigT)

  # Compute estimates for adjusted
  adj.true <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                           omega=true_mdri, omega_var=0,
                           beta=true_frr, beta_var=0,
                           big_T=bigT)
  adj.est <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p, n=data$n,
                          omega=assay$omega_est, omega_var=assay$omega_var,
                          beta=assay$beta_est, beta_var=assay$beta_var,
                          big_T=bigT)

  # Compute estimates for adjusted enhanced estimator
  eadj.true <- get.adjusted.pt(n_r=data$n_r_pt, n_n=data$n_n, n_p=data$n_p, n=data$n,
                           omega=true_mdri, omega_var=0,
                           beta=true_frr, beta_var=0,
                           big_T=bigT,
                           num_beta=data$num_beta,
                           den_omega=data$den_omega,
                           den_beta=data$den_beta)

  eadj.est <- get.adjusted.pt(n_r=data$n_r_pt, n_n=data$n_n, n_p=data$n_p, n=data$n,
                          omega=assay$omega_est, omega_var=assay$omega_var,
                          beta=assay$beta_est, beta_var=assay$beta_var,
                          big_T=bigT,
                          num_beta=data$num_beta,
                          den_omega=data$den_omega,
                          den_beta=data$den_beta)

  return(list(truth=rep(baseline_incidence, n_sims),
              adj_true_est=adj.true$est,
              adj_true_var=adj.true$var,
              adj_est_est=adj.est$est,
              adj_est_var=adj.est$var,
              eadj_true_est=eadj.true$est,
              eadj_true_var=eadj.true$var,
              eadj_est_est=eadj.est$est,
              eadj_est_var=eadj.est$var,
              mu_est=assay$mu_est,
              mu_var=assay$mu_var,
              omega_est=assay$omega_est,
              omega_var=assay$omega_var,
              beta_est=assay$beta_est,
              beta_var=assay$beta_var,
              n_n=data$n_n,
              n_r=data$n_r,
              n_p=data$n_p,
              n_r_pt=data$n_r_pt,
              num_beta=data$num_beta,
              den_omega=data$den_omega,
              den_beta=data$den_beta,
              q_eff=data$q_eff,
              n=data$n
  ))
}
