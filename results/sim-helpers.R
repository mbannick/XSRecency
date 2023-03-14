setwd("~/repos/XSRecency/")
source("R/data-generator.R")
source("R/data-summary.R")
source("R/external-study-sim.R")
source("R/external-study-data.R")
source("R/estimators.R")
source("R/estimators-pt.R")
source("R/utils.R")
source("R/phi-matrix.R")

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
                     t_min=NULL, t_max=NULL, t_noise=NULL,
                     d_misrep=0.0, q_misrep=0.0, p_misrep=0.0,
                     ptest.dist2=NULL,
                     exclude_pt_bigT=FALSE){

  # Generate assay simulations
  cat("Generating assay simulations\n")
  assay.nsim <- assay.nsim.pt(
    n_sims=n_sims, phi.func=phi.func, tau=tau, bigT=bigT,
    ext_FRR=ext_FRR, ext_df=ext_df,
    max_FRR=max_FRR
  )

  # Generate screening data at time 0 with n people
  cat("Generating screening data\n")
  sim <- sim.screening.generator(
    prevalence=prevalence,
    phi.func=phi.func,
    e.func=function(e) infection.function(e, t=0,
                                          p=prevalence,
                                          lambda_0=baseline_incidence,
                                          rho=rho)
  )
  dfs_screening <- replicate(n_sims, sim(n, return_n=FALSE), simplify=F)

  # Generate prior testing data
  cat("Generating prior testing data\n")
  sim.pt <- sim.pt.generator(
    ptest.dist=ptest.dist,
    ptest.prob=ptest.prob,
    ptest.dist2=ptest.dist2
  )

  pt.dfs <- lapply(dfs_screening, sim.pt)

  if(exclude_pt_bigT){
    if(is.null(t_range)){
      t_range.m <- c(0, bigT)
    } else {
      t_range.m <- c(min(t_range), min(max(t_range), bigT))
    }
  } else if(!is.null(t_min)){
    t_range.m <- c(t_min, t_max)
  }
  modify.pt <- modify.pt.generator(
    t_noise=t_noise,
    t_range=t_range.m,
    d_misrep=d_misrep,
    p_misrep=p_misrep
  )
  pt.dfs <- lapply(pt.dfs, modify.pt)

  ns <- rep(n, n_sims)
  n_p <- lapply(dfs_screening, function(x) nrow(x))

  rename <- function(x) setnames(x, c("id", "ui", "ri"))
  lapply(assay.nsim$studies, rename)

  cat("Applying enhanced estimator\n")
  eadj <- mapply(
    FUN=get.adjusted.pt,
    n=ns,
    n_p=n_p,
    ptdf=pt.dfs,
    beta=assay.nsim$beta_sim[1,],
    beta_var=assay.nsim$beta_sim[2,],
    phidat=assay.nsim$studies,
    use_geese=TRUE,
    formula="ri ~ poly(ui, 3, raw=TRUE)",
    family=replicate(n_sims, binomial(link="logit"), simplify=FALSE),
    big_T=bigT,
    plot_phi=FALSE
  )

  adj <- mapply(
    FUN=get.adjusted,
    n=ns,
    n_p=n_p,
    n_n=eadj["n_n",],
    n_r=eadj["n_r",],
    omega=eadj["omega",],
    omega_var=eadj["omega_var",],
    beta=assay.nsim$beta_sim[1,],
    beta_var=assay.nsim$beta_sim[2,],
    big_T=bigT
  )

  # # Generate trial data
  # data <- generate.data(
  #   n=n, n_sims=n_sims,
  #   infection.function=infection.function,
  #   phi.func=phi.func,
  #   baseline_incidence=baseline_incidence,
  #   prevalence=prevalence,
  #   rho=rho,
  #   bigT=bigT,
  #   ptest.dist=ptest.dist,
  #   ptest.prob=ptest.prob,
  #   t_range=t_range,
  #   t_noise=t_noise,
  #   d_misrep=d_misrep,
  #   q_misrep=q_misrep,
  #   p_misrep=p_misrep,
  #   ptest.dist2=ptest.dist2,
  #   exclude_pt_bigT=exclude_pt_bigT
  # )

  # Compute assay parameters based on external data simulation and
  # assay simulation from before
  # assay <- assay.properties.pt(
  #   studies=assay.nsim$studies,
  #   beta_sim=assay.nsim$beta_sim,
  #   bigT=bigT,
  #   tau=tau,
  #   last_point=last_point,
  #   ptest_times=data$ptest_times,
  #   ptest_delta=data$ptest_delta,
  #   ptest_avail=data$ptest_avail,
  #   ri=data$ri
  # )
  #
  # wvec <- get.w.vec(data, assay, bigT)

  # DEBUGGING FOR COV_14

  # Calculate true assay parameters
  # true_frr <- true.frr(phi.func=phi.func, bigT=bigT, tau=tau)
  # true_mdri <- true.window.mdri(phi.func=phi.func, maxT=bigT)

  # inputs <- data
  # inputs[["big_T"]] <- bigT

  # Carry over all of the assay properties
  # including those that have been estimated for prior test results
  # for(a in names(assay)){
  #   inputs[[a]] <- assay[[a]]
  # }

  # inputs_true <- inputs
  # inputs_est <- inputs

  # inputs_true[["omega"]] <- true_mdri
  # inputs_true[["omega_var"]] <- 0
  # inputs_true[["beta"]] <- true_frr
  # inputs_true[["beta_var"]] <- 0

  # inputs_est[["omega"]] <- assay$omega_est
  # inputs_est[["omega_var"]] <- assay$omega_var
  # inputs_est[["beta"]] <- assay$beta_est
  # inputs_est[["beta_var"]] <- assay$beta_var

  # Compute estimates for adjusted
  # adj.true <- R.utils::doCall(.fcn=get.adjusted, args=inputs_true, .ignoreUnusedArgs=TRUE)
  # adj.est <- R.utils::doCall(.fcn=get.adjusted, args=inputs_est, .ignoreUnusedArgs=TRUE)

  # Compute estimates for adjusted prior testing
  # eadj.true <- R.utils::doCall(.fcn=get.adjusted.pt, args=inputs_true, .ignoreUnusedArgs=TRUE)
  # eadj.est <- R.utils::doCall(.fcn=get.adjusted.pt, args=inputs_est, .ignoreUnusedArgs=TRUE)

  # The first component of get.adjusted.pt$var is the more robust variance (rob)
  # The second is assumption-based (asm)
  results <- list(
    truth=rep(baseline_incidence, n_sims),
    # adj_true_est=adj.true$est,
    # adj_true_var=adj.true$var,
    adj_est_est=unlist(adj["est",]),
    adj_est_var=unlist(adj["var",]),
    # eadj_true_est=eadj.true$est,
    # eadj_true_var_rob=eadj.true$var[[1]],
    # eadj_true_var_asm=eadj.true$var[[2]],
    eadj_est_est=unlist(eadj["est",]),
    eadj_est_var_rob=unlist(eadj["var",]),
    eadj_est_var_asm=unlist(eadj["var",]),
    omega_est=unlist(eadj["omega",]),
    omega_var=unlist(eadj["omega_var",]),
    beta_est=unlist(assay.nsim$beta_sim[1,]),
    beta_var=unlist(assay.nsim$beta_sim[2,]),
    n_n=unlist(eadj["n_n",]),
    n_r=unlist(eadj["n_r",]),
    n_p=unlist(eadj["n_p",]),
    n_r_pt=unlist(eadj["n_r_pt",]),
    q_eff=unlist(eadj["q_eff",]),
    n=unlist(eadj["n",])
  )
  for(elem in rownames(eadj)){
    results[[elem]] <- unlist(eadj[elem,])
  }
  # for(w in names(wvec)){
  #   results[[w]] <- wvec[[w]]
  # }
  # for(e in names(eadj.est$components_est)){
  #   e_rob <- paste0(e, "_rob")
  #   e_asm <- paste0(e, "_asm")
  #   results[[e_rob]] <- (eadj.est$components_est)[[e]][[1]]
  #   results[[e_asm]] <- (eadj.est$components_est)[[e]][[2]]
  # }
  # for(e in colnames(assay)){
  #   if(e %in% c("mu_est", "mu_var",
  #               "omega_est", "omega_var",
  #               "beta_est", "beta_var")) next
  #   results[[e]] <- assay[[e]]
  # }

  return(results)
}
