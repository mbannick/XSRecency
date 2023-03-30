library(simtools)

DEFAULTS <- list(

  seed=100,
  out_dir="./",

  n_sims=2,
  sim_blocksize=2,
  start_sim=1,

  n=5000,
  p=0.29,
  inc=0.032,
  window=101,
  shadow=194,
  itype="constant",
  rho=0,
  tau=12,
  bigT=2,

  phi_frr=NULL,
  phi_tfrr=2,
  phi_norm_mu=NULL,
  phi_norm_sd=NULL,
  phi_norm_div=NULL,
  phi_pnorm_mu=NULL,
  phi_pnorm_sd=NULL,
  phi_pnorm_div=NULL,

  ext_FRR=FALSE,
  duong_scale=NULL,
  max_FRR=NULL,
  last_point=FALSE,

  pt=TRUE,
  t_min=0,
  t_max=4,
  q=1,

  gamma=0.0, # variance for the Gaussian noise to add to prior test time
  eta=0.0, # the probability of incorrectly reporting negative test
  nu=0.0, # the probability of failing to report prior test result
  xi=0.0, # the probability of failing to report prior positive test results
  mech2=TRUE,
  t_min_exclude=NULL,
  exclude_pt_bigT=FALSE
)

get.startsims <- function(nsims, simblock) seq(1, nsims, simblock)

fill.overall <- function(params, overall){
  for(param in names(overall)){
    if(param == "out_dir") next()
    if(param %in% colnames(params)) next()
    params[[param]] <- OVERALL[[param]]
  }
  return(params)
}

get.paramlist.a <- function(startsims){

  params <- expand.grid(list(
    q=c(0.2, 0.4, 0.6, 0.8, 1.0),
    t_min=c(0, 2),
    t_max=c(2, 4),
    start_sim=startsims
  ))
  params <- data.table(params)
  params <- params[!(t_min == 2 & t_max == 2)]

  return(params)
}

get.paramlist.b <- function(startsims){

  params <- expand.grid(list(
    start_sim=startsims
  ))
  params <- data.table(params)
  params <- fill.overall(params, overall)

  return(params)
}

get.paramlist.c <- function(startsims){

  params <- expand.grid(list(
    q=0.5,
    gamma=c(0, 0.083, 0.5),
    eta=c(0, 0.1),
    xi=c(0, 0.1),
    start_sim=get.startsims(
      overall[["n_sims"]],
      overall[["sim_blocksize"]])
  ))
  params <- data.table(params)
  params <- params[!(eta == 0.1 & xi == 0.1)]

  return(params)
}

get.paramlist.d <- function(startsims){

  params <- expand.grid(list(
    q=0.5,
    t_min=0,
    t_max=c(4, 12),
    itype=c("constant", "piecewise"),
    rho=c(0, 0.0039),
    exclude_pt_bigT=c(FALSE, TRUE),
    start_sim=startsims
  ))
  params <- data.table(params)
  params <- params[!(itype == "constant" & rho != 0)]
  params <- params[!(itype == "piecewise" & rho == 0)]

  return(params)
}


