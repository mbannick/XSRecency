# PHI FIGURE ----------------------------------

pdf('RAplot-2.pdf',height=6,width=10)
par(mfrow=c(1,2))

for(i in 1:2){

  if(i == 1){
    mdri <- 101
    shadow <- 194
  } else {
    mdri <- 248
    shadow <- 306
  }

  t = seq(0, 12, 1e-3)
  params <- get.gamma.params(window=mdri/365.25, shadow=shadow/365.25)

  phit <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

  if(i == 1) ttime <- 2
  if(i == 1) tval <- phit(2)
  if(i == 1) name <- "Assay 1"

  if(i == 2) tval <- 0.02
  if(i == 2) ttime <- uniroot(function(t) phit(t) - tval, interval=c(0, 12))$root
  if(i == 2) name <- "Assay 2"

  phit.const <- function(t) phit(t)*(t <= ttime) + tval*(t > ttime)
  phit.const.dnorm <- function(t) phit.const(t) + dnorm(t-7, mean=0, sd=1) / 8
  phit.const.dnorm.ext <- function(t) phit.const(t) + pnorm(t, mean=10, sd=2) / 10
  # TRY THIS FUNCTION
  plot(t, phit(t), col='black', type='l', ylab=expression(phi(t)), xlab=expression(t))
  abline(h=tval, lty='dashed')

  ###### MDRI at 2 year 142/365, FRR 1.5%: 1.5% + Gamma distribution with mean 142/365.25-1.5%
  lines(t, phit.const(t), col='red')

  ##### Non-constant FRR -- peak at 7 years
  lines(t, phit.const.dnorm(t), col='blue')

  ##### Non-constant FRR -- increase after 7 years
  lines(t, phit.const.dnorm.ext(t), col='orange')

  if(i == 2){
    pgon <- c(seq(2, ttime, 0.01))
    pgon.phi <- 1-pgamma(pgon, shape = params[1], rate = params[2])
    polygon(x=c(pgon, rev(pgon)), y=c(rep(tval, length(pgon)), rev(pgon.phi)), col='lightgrey')
  }

  legend('topright', paste0(name, c("A", "B", "C", "D")),
         lty=rep(1, 4),col=c("black", "red", "blue", "orange"), cex=0.75)
  abline(v=2, lty='dashed')
}

dev.off()
