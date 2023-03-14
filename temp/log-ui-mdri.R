source("./R/phi-functions.R")
library(XSRecency)
library(data.table)
library(geepack)
library(pbapply)

params <- get.gamma.params(window=248/365.25, shadow=306/365.25)
phi.func <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

bigT <- 1
DT <- 0.0001

ts <- seq(DT, 10, by=DT)
plot(phi.func(ts) ~ ts, type='l')

true <- sum(phi.func(ts)[ts <= bigT]*DT)

studies <- simulate.studies(nsims=100, phi.func=phi.func)

fit.mod <- function(study){
  colnames(study) <- c("id", "ui", "ri")
  mod <- geese(ri ~ poly(log(ui), degree=4, raw=TRUE),
               family=binomial(link="logit"), data=study[study$ui > 0,])
  ts_plot <- seq(DT, 10, by=DT)
  idx <- max(which(ts_plot <= bigT))
  preds <- .predict.phi(ts=ts_plot,
                        model=mod, family=binomial(link="logit"),
                        varcov=FALSE)[["point"]]
  res <- sum(preds[1:idx]*DT)
  lines(preds ~ ts_plot, col='blue')
  return(res)
}

sapply(studies[1:10], fit.mod)

results <- pbsapply(studies, fit.mod)
mean(results)
