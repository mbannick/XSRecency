
get.phi <- function(window, shadow, tau,
                    phi_tfrr=NULL, phi_frr=NULL,
                    phi_norm_mu=NULL, phi_norm_sd=NULL, phi_norm_div=NULL){

  if(!is.null(phi_frr) & !is.null(phi_tfrr)){
    stop("Can't provide both frr and time for frr.")
  }
  if(!is.null(phi_norm_mu)){
    if(is.null(phi_norm_sd) | is.null(phi_norm_div)){
      stop("Need stdev and divided by params for normal.")
    }
  }

  pp <- get.gamma.params(window=window/365.25, shadow=shadow/365.25)

  # Set up each type of phi function, will be overwritten
  phi.none <- function(t) 1-pgamma(t, shape=pp[1], rate=params[2])
  phi.const <- function(t) 1-pgamma(t, shape=params[1], rate=params[2])
  phi.norm <- function(t) 1-pgamma(t, shape=params[1], rate=params[2])
  phit.pnorm <- function(t) 1-pgamma(t, shape=params[1], rate=params[2])

  phi.func <- phi.none

  # Get the phi function with constant FRR either past a certain time
  # or fixed after it hits some value.
  if(!is.null(phi_tfrr) | !is.null(phi_frr)){
    if(!is.null(phi_tfrr)){
      ttime <- phi_tfrr
      tval <- phi.none(ttime)
    }
    if(!is.null(phi_frr)){
      tval <- phi_frr
      ttime <- uniroot(function(t) phi.none(t) - tval, interval=c(0, tau))$root
    }
    phi.const <- function(t, ...) phi.none(t)*(t <= ttime) + tval*(t > ttime)
    phi.func <- phi.const
  }
  if(!is.null(phi_norm_mu)){
    phi.norm <- function(t) phi.const(t) + dnorm(t-phi_norm_mu, mean=0, sd=phi_norm_sd) / phi_norm_div
    phi.func <- phi.norm
  }
  if(!is.null(phi_pnorm_mu)){
    phi.pnorm <- function(t) phi.const(t) + pnorm(t-phi_pnorm_mu, mean=0, sd=phi_pnorm_sd) / phi_pnorm_div
    phi.func <- phi.pnorm
  }
}
