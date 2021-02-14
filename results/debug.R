rm(list=ls())

# Choose a setting
setting <- 2

# Which setting are we working in?
if(setting == 1){
  window <- 80
  shadow <- 25
} else {
  window <- 248
  shadow <- 306
}

# Divide by number of days
window <- window/365.25
shadow <- shadow/365.25

# Get gamma parameter based on this for the phi function
alpha <- window ** 2 / (2 * shadow - window ** 2)
beta <- window / (2 * shadow - window ** 2)

phi.func <- function(t) 1-pgamma(t, shape=alpha, rate=beta)

# T^* = 2, tau is the max value.
# DELTA is the integration step
t <- 2
tau <- 12
DELTA <- 1e-3
ts <- seq(0, t, DELTA)
ts.frr <- seq(t, tau, DELTA)
taus <- seq(0, tau, DELTA)

# phi plot
plot(taus, phi.func(taus), type='l')
abline(v=t, lty='dashed')
abline(h=0)

# Get the window, omega, and beta
window <- sum(phi.func(taus) * DELTA)
omega <- sum(phi.func(ts) * DELTA)
frr <- mean(phi.func(ts.frr))

# Calculate the shadow period for the snapshot
# And the adjusted estimators
shadow.snap <- sum(phi.func(taus) * taus / window * DELTA)
shadow.adj <- sum((phi.func(ts) - frr) * ts / (omega - 2 * frr) * DELTA)

# Incidence function
rho <- 0.0028
lambda_0 <- 0.032
inc <- function(t) lambda_0 - rho * t

# Compare the numerical bias to the analytical bias
# For the snapshot estimator
sum(phi.func(taus)/window * inc(-taus) * DELTA) - lambda_0
shadow.snap * rho

# Compare the numerical bias to the analytical bias
# For the adjusted estimator
sum((phi.func(ts) - frr)/(omega - 2 * frr) * inc(-ts) * DELTA) - lambda_0
shadow.adj * rho

