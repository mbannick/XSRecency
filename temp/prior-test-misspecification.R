# DISTRIBUTION OF PRIOR TESTS

# Model Misspecification

u <- seq(0, 10, by=0.01)

int <- 1
slope <- 0.5
expected <- int + slope * u

par(mfrow=c(1, 2))
plot(expected ~ u, type='l',
     ylab="Mean Prior Testing Time (Years)",
     xlab="Duration of Infection (Years)",
     main="Mean Prior Testing Time\nby Duration of Infection")

coefs <- c(-3, 0.03, 0.04, 0.001)
us <- cbind(1, u, u**2, u**3)
lprobs <- us %*% coefs
expit <- function(x) exp(x) / (1 + exp(x))
plot(expit(lprobs) ~ u, type='l',
     ylab="Probability of Prior Test",
     xlab="Duration of Infection (Years)",
     main="Probability of Prior Test\nby Duration of Infection",
     ylim=c(0, 1))
