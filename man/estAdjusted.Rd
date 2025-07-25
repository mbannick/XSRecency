% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{estAdjusted}
\alias{estAdjusted}
\title{Adjusted estimator and variance (Kassanjee et al. 2012)}
\usage{
estAdjusted(
  n_r,
  n_n,
  n_p,
  n,
  omega,
  omega_var,
  beta,
  beta_var,
  big_T,
  n_p_test = NULL,
  q_var = FALSE
)
}
\arguments{
\item{n_r}{Number of recent positives}

\item{n_n}{Number of negatives}

\item{n_p}{Number of positives}

\item{n}{Total number of observations (\eqn{n_n + n_p = n})}

\item{omega}{Mean duration of recent infection (MDRI) (in years, not days)}

\item{omega_var}{Variance of the estimator for MDRI (or 0 if MDRI known)}

\item{beta}{False recency rate (FRR)}

\item{beta_var}{Variance of the estimator for FRR (or 0 if FRR known)}

\item{big_T}{The \eqn{T^*} in the equation above}

\item{n_p_test}{Optional number of positives that were given a recency test. If left \code{NULL},
assumed to be equal to \code{n_p}.}
}
\value{
Returns a list of the estimate and the variance.
}
\description{
\deqn{
    \hat{\lambda} = \frac{N_{rec}/q - N_{pos} \hat{\beta}_{T^*}}{N_{neg} (\hat{\Omega}_{T^*} - \hat{\beta}_{T^*} T^*)}
  }
where \eqn{\hat{\Omega}_{T^*}} is an estimate of the mean duration of recent infection (MDRI),
\eqn{\hat{\beta}_{T^*}} is an estimate of the false recency rate, and \eqn{T^*} is the time cutoff for being
a recent infected or not.
}
\details{
\code{q} is an adjustment for the number of positives that are given recency tests.
}
\examples{
estAdjusted(n_r=2, n_n=50, n_p=10, n=60,
                  omega=0.36, omega_var=0, beta=0.02, beta_var=0, big_T=2, n_p_test=9)
estAdjusted(n_r=c(2, 3), n_n=c(50, 48), n_p=c(10, 12), n=c(60, 60),
                  omega=0.36, omega_var=0, beta=0.02, beta_var=0, big_T=2)
}
