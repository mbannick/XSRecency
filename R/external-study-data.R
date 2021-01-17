#' External study data parameters for estimating
#' omega, mu, and beta through simulation.
#'
#' @export
#'
PHI_PARAMS <- list(

  # what is the follow-up time
  FOLLOW_T=2,

  # tail probability for phi
  TAIL_PROB=0.016,

  # number of long infecteds for beta study
  N_LONG_INFECT=1388,

  # minimum and maximum infection times for long infecteds
  FRR_MIN=2,
  FRR_MAX=12,

  # cohorts
  coh=c(1, 2, 3),

  # number of people in each cohort
  coh_n=c(90, 162, 25),

  # mean number of samples per cohort
  mean_samp=c(6, 12, 4),

  # max number of samples per cohort
  max_samp=c(7, 20, 4),

  # multinomial probabilities buckets
  duration_years=list(c(0.0, 0.5),
                      c(0.5, 1.0),
                      c(1.0, 2.0),
                      c(2.0, 3.0),
                      c(3.0, 5.0),
                      c(5.0, 9.9)),

  # number of people in each duration of infection years category above
  # per cohort
  coh_buckets=list(
    c(159, 173, 088, 076, 022, 000),
    c(306, 262, 448, 105, 347, 371),
    c(042, 043, 000, 000, 000, 000)
  )
)
