rm(list=ls())

window <- 71
shadow <- 80

# Divide by number of days
window <- window/365.25
shadow <- shadow/365.25

# Get gamma parameter based on this for the phi function
# From the overleaf, Appendix C
alpha <- window / (2 * shadow - window)
beta <- 1 / (2 * shadow - window)

phi.func <- function(t) 1-pgamma(t, shape=alpha, rate=beta)

# DELTA is the integration step
tau <- 12
DELTA <- 1e-3
taus <- seq(0, tau, DELTA)

# phi plot
par(mfrow=c(1,1))
plot(taus, phi.func(taus), type='l')
abline(v=2, lty='dashed')
abline(h=0)

# Get the mean window period
window.calc <- sum(phi.func(taus) * DELTA)
# This is slightly different, but not too bad
window.calc
window

# Calculate the shadow period
shadow.calc <- sum(phi.func(taus) * taus / window.calc * DELTA)
# **THIS is very different. Should be the same as shadow.
shadow.calc
shadow
