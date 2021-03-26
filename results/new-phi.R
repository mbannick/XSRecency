library(extraDistr)
library(XSRecency)

parms.1 <- get.gamma.params(71/365.25, 80/365.25)
parms.2 <- get.gamma.params(101/365.25, 194/365.25)

phi.gamma.1 <- function(t) 1 - pgamma(t, shape=parms.1[1], rate=parms.1[2])
phi.gamma.2 <- function(t) 1 - pgamma(t, shape=parms.2[1], rate=parms.2[2])

h.1 <- phi.gamma.1(2) # want this to be close to 0.026
md.1 <- true.window.mdri(phi.gamma.1, maxT=12)*365.25 # want this to be close to 0.273
s.1 <- true.shadow.snap(phi.gamma.1, tau=12)*365.25

h.2 <- phi.gamma.2(2) # want this to be close to 0.026
md.2 <- true.window.mdri(phi.gamma.2, maxT=12)*365.25 # want this to be close to 0.27
s.2 <- true.shadow.snap(phi.gamma.2, tau=12)*365.25

phi.hnorm <- function(t) 1 - phnorm(t, sigma=1)

h.3 <- phi.hnorm(2) # want this to be close to 0.026
md.3 <- true.window.mdri(phi.hnorm, maxT=12)*365.25 # want this to be close to 0.273
s.3 <- true.shadow.snap(phi.hnorm, tau=12)*365.25

phi.weib <- function(t) 1 - pweibull(t, shape=0.68, scale=0.25)
h.4 <- phi.weib(2)
md.4 <- true.window.mdri(phi.weib, maxT=12)*365.25
s.4 <- true.shadow.snap(phi.weib, tau=12)*365.25

pdf("test-distributions.pdf", height=8, width=8)
par(mfrow=c(2, 2))

t <- seq(0, 12, by=1e-3)

plot(phi.gamma.1(t) ~ t, type='l',
     main=paste0("Original Gamma CDF: ", round(h.1, 4), " at T=2 \n MDRI=", round(md.1, 4),
                 "\n Shadow=", round(s.1, 4)))
abline(v=2, lty='dashed')

plot(phi.gamma.2(t) ~ t, type='l',
     main=paste0("Gamma CDF: ", round(h.2, 4), " at T=2 \n MDRI=", round(md.2, 4),
                 "\n Shadow=", round(s.2, 4)))
abline(v=2, lty='dashed')

plot(phi.hnorm(t) ~ t, type='l',
     main=paste0("Half Normal CDF: ", round(h.3, 4), " at T=2 \n MDRI=", round(md.3, 4),
                 "\n Shadow=", round(s.3, 4)))
abline(v=2, lty='dashed')

plot(phi.weib(t) ~ t, type='l',
     main=paste0("Weibull CDF: ", round(h.4, 4), " at T=2 \n MDRI=", round(md.4, 4),
                 "\n Shadow=", round(s.4, 4)))
abline(v=2, lty='dashed')
dev.off()
