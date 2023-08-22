get.w.vec <- function(data, assay, bigT){
  w1 <- data$n_r_pt - assay$beta_est * (data$n_p - assay$nB)
  w2 <- data$n_p
  w3 <- assay$omega_est - assay$beta_est * bigT
  w4 <- assay$nATO
  w5 <- assay$beta_est * assay$nBT
  return(list(
    W1=w1, W2=w2, W3=w3, W4=w4, W5=w5
  ))
}
