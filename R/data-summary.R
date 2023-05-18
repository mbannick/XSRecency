
# All elements that would be used in the estimators,
# defaults to zero for each.

ELEMENTS <- c(
  "n", "n_p", "n_n", "n_r",
  "n_r_pt", "q_eff",
  "p_A", "p_B",
  "num_beta", "den_beta", "den_omega",
  "omega_TA", "omega_TA_var", "omega_TAstar",
  "mu_TA", "var_TA", "mu_TB", "var_TB",
  "r_TA", "r_TAprime", "r_TAstar"
)

LISTKEY <- list()
for(elem in ELEMENTS){
  LISTKEY[[elem]] <- 0
}

# Function that summarizes a dataset
summarize.data <- function(df){

  s <- LISTKEY
  s[["n"]] <- nrow(df)
  s[["n_p"]] <- sum(df$di)
  s[["n_n"]] <- s[["n"]] - s[["n_p"]]
  s[["n_r"]] <- sum(df$ri, na.rm=TRUE)

  return(s)
}

.test.phi <- function(phidat, formula, family, use_geese=TRUE, ...){

  if(use_geese){
    modfunc <- geese
  } else {
    modfunc <- glm
  }
  mod <- modfunc(formula=as.formula(formula), family=family, ..., data=phidat)
  dat <- data.frame(
    ri=rnorm(10),
    ui=1:10
  )
  newdat <- model.matrix(as.formula(formula), data=dat)
  tryCatch(
    expr={
      if(use_geese){
        coefs <- mod$beta
      } else {
        coefs <- coef(mod)
      }
      preds <- newdat %*% coefs
    },
    error=function(e){
      msg <- paste0(
        "cannot create predictions from model specifications",
        "\noriginal error: ", e
      )
      stop(msg)
    }
  )
}

.get.ts <- function(minT, maxT, dt){
  index <- seq(minT/dt, maxT/dt, by=1)
  ts <- index*dt
  dat <- data.frame(index=index, ts=ts)
  return(dat)
}

# The algorithm to get recent classification
# based on recency assay + prior testing
.get.Mi <- function(ptdf, bigT){

  Np <- nrow(ptdf)

  Ai   <- rep(0, Np)
  Bi   <- rep(0, Np)
  BiTi <- rep(0, Np)

  Qi   <- as.logical(ptdf$qi)
  Di   <- ptdf$deltai
  Ti   <- ptdf$ti
  Ci   <- Ti <= bigT

  Ai[Ci & Qi]  <- 1
  Bi[!Ci & Qi] <- 1

  # This is the old recency indicator
  Ri <- ptdf$ri

  # This is the new recency indicator
  Mi <- Ri * (1 - Bi * Di) + (1 - Ri) * Ai * (1 - Di)
  Mi[!Qi] <- Ri[!Qi]

  return(data.frame(
    Mi=Mi,
    Ai=Ai,
    Bi=Bi,
    Ci=Ci
  ))
}

.map.to.tgrid <- function(ts_index, ptdf, bigT_index){

  closest <- ptdf[, c("ti", "qi", "Mi", "Ai", "Bi", "Ci")]
  closest <- data.table(closest)
  setnames(closest, c("ti", "qi"), c("ts_orig", "has_test"))

  # closest <- data.table(
  #   ts_orig=ptdf$ti,
  #   has_test=ptdf$qi,
  #   Mi=ptdf$Mi
  # )
  closest[, id := .I]
  closest[, index := NA]

  get.closest.index <- function(t){
    if(is.na(t)){
      idx <- NA
    } else {
      val <- which.min(abs(t - ts_index$ts))
      idx <- ts_index[val, "index"]
    }
  }
  closest_index <- unlist(sapply(closest$ts_orig, get.closest.index))
  closest[, index := closest_index]
  closest <- merge(closest, ts_index, by="index", all.x=T)

  # closest[, Ai := as.numeric((ts <= bigT_index[, "ts"]) & (has_test))]
  # closest[, Bi := as.numeric((ts > bigT_index[, "ts"]) & (has_test))]

  return(closest)
}

.pt.properties.est <- function(closest, bigTidx, cphi){
  elem <- list()

  # Calculate fraction of those with recent/non-recent prior tests
  p_A <- mean(closest$Ai)
  p_B <- mean(closest$Bi)

  elem[["p_A"]] <- p_A
  elem[["p_B"]] <- p_B

  # This is \sum 1 - I(T_i > T^*)
  elem[["num_beta"]] <- sum(!closest$Bi)

  if(p_A > 0){

    # Create point estimates and covariance matrix between
    # the phi estimates at different times

    cphi_point <- cphi$cphi
    cphi_covar <- cphi$csum

    # Get ids for those who have a recent prior test, and their
    # corresponding time index
    idmap <- closest[Ai == 1, .(id, index, ts, Mi)]
    # Get unique combinations of the ids
    # for the covariance calculation
    idmap_covar <- combn(idmap$id, 2) %>% t %>% data.table
    setnames(idmap_covar, c("idX", "idY"))
    idmap_covar <- merge(idmap_covar, idmap[, .(id, index)], by.x="idX", by.y="id")
    setnames(idmap_covar, "index", "indexX")
    idmap_covar <- merge(idmap_covar, idmap[, .(id, index)], by.x="idY", by.y="id")
    setnames(idmap_covar, "index", "indexY")

    # Calculate omega_T conditional expectation and variance
    omega_ta     <- merge(idmap, cphi_point, by="index", all.x=TRUE)
    elem[["omega_TA"]] <- mean(omega_ta$phi)
    elem[["omega_TA_var"]] <- var(omega_ta$phi)

    # Calculate omega_T * T conditional expectation
    omega_ta <- omega_ta[, tastar := ts * phi]
    omega_ta <- omega_ta[, taneg := ts - phi]
    elem[["omega_TAstar"]] <- mean(omega_ta$tastar)
    elem[["den_omega"]] <- sum(omega_ta$taneg)

    # Calculate r_TiTi conditional expectation and variance
    idmap2 <- idmap
    idmap2[, indexX := index]
    idmap2[, indexY := index]

    r_Tii  <- merge(idmap2, cphi_covar, by=c("indexX", "indexY"), all.x=T)
    elem[["r_TA"]]   <- mean(r_Tii$rho)
    elem[["var_TA"]] <- var(r_Tii$rho)

    # Calculate r_TiTj conditional expectation
    r_Tij <- merge(idmap_covar, cphi_covar, by=c("indexX", "indexY"), all.x=T)
    elem[["r_TAprime"]] <- mean(r_Tij$rho)

    # Calculate r_TiTstar conditional expectation
    idmap3 <- idmap
    idmap3[, indexX := index]
    idmap3[, indexY := bigTidx]

    r_Tistar <- merge(idmap3, cphi_covar, by=c("indexX", "indexY"), all.X=T)
    elem[["r_TAstar"]] <- mean(r_Tistar$rho)

    # Calculate the expected time of prior tests
    ta <- closest[Ai == 1, ts_orig]
    elem[["mu_TA"]] <- mean(ta)
    elem[["var_TA"]] <- var(ta)
  }

  if(p_B > 0){

    # Calculate the expected time of prior tests
    tb <- closest[Bi == 1, ts_orig]
    elem[["mu_TB"]] <- mean(tb)
    elem[["var_TB"]] <- var(tb)

    # This is I(T_i > T^*) * T_i
    elem[["den_beta"]] <- sum(tb)

  }
  return(elem)
}

#' Estimate a phi function, and output omega
#'
#' @param phidat Data frame used to estimate the phi function, needs to have
#'               columns ri (recency indicator), ui (duration), and id (if use_geese == TRUE)
#' @param maxT Maximum time to estimate phi until
#' @param min_dt Make the minimum time dt (necessary if doing a log transformation of ui)
#' @param bigT T^\* for estimating \Omega_{T^\*}
#' @param use_geese Indicator to fit a gee model using geese(), or a glm model with glm()
#' @param formula Formula for fitting the model to phidat.
#'                Formula argument must take in only ui as the newdata.
#'                For example, do not create a ui^2 variable. Use poly(ui, ...) function
#'                to fit polynomial terms.
#' @param family Family argument for glm or gee
#' @param dt An integration step size. Should be no more than 1 day.
#' @param plot_phi Whether or not to plot the estimated phi function
#' @param return_all Whether to return additional objects used in internal functions (will use more memory)
#' @param ... Additional arguments to glm or geese for fitting model
#' @export
#' @import formula.tools
#'
#' @returns List of results:
#' * `omega`: Estimate of \Omega_{T^\*}
#' * `omega_var`: Variance of \hat{\Omega_{T^\*}}
#' * `bigTidx`: Index of T^\* based on `ts_index`
#' * `ts_index`: Index mapping each time point to an integer
#' * `get.integral.est`: Function to give an estimate of \Omega_{T} for some T as arg to the function
#'
#' @examples
#' # Define phi function
#' params <- get.gamma.params(window=248/365.25, shadow=306/365.25)
#' phi.func <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
#'
#' # Simulate external study data
#' study <- simulate.studies(1, phi.func)[[1]]
#' setnames(study, c("id", "ui", "ri"))
#'
#' # Estimate omega based on 3rd degree polynomial
#' estimate.phi(phidat=study, maxT=5, bigT=1,
#'              formula="ri ~ poly(ui, raw=T, degree=4)", family=binomial(link="logit"),
#'              use_geese=TRUE, plot_phi=TRUE)
#'
#' # Estimate omega based on 3rd degree polynomial, log function
#' # NOTE: If you do this, need to not have any ui = 0
#' study.log <- study
#' study.log[study.log$ui == 0, "ui"] <- 0.01
#' estimate.phi(phidat=study.log, maxT=10, bigT=1, min_dt=TRUE,
#'              formula="ri ~ poly(log(ui), raw=T, degree=4)", family=binomial(link="logit"),
#'              use_geese=TRUE, plot_phi=TRUE)
#'
#' # Return additional function that can be used to
#' # quickly get other omega estimates
#' result <- estimate.phi(phidat=study, maxT=5, bigT=1,
#'                        formula="ri ~ poly(ui, raw=T, degree=3)", family=binomial(link="logit"),
#'                        use_geese=TRUE, plot_phi=FALSE, return_all=TRUE)
#'
#' # What is \int_0^T \phi(t) dt for ts=T?
#' result$get.integral.est(ts=2.1)
estimate.phi <- function(phidat, maxT, bigT, dt=1/365.25, min_dt=FALSE,
                         formula, family, use_geese,
                         return_all=FALSE, plot_phi=FALSE, ...){

  # Get map from time to integer
  minT <- ifelse(min_dt, dt, 0)
  ts_index <- .get.ts(minT=minT, maxT=maxT, dt=dt)

  # Make sure model works based on the user arguments
  # .test.phi(phidat=phidat, use_geese=use_geese, formula=formula, family=family, ...)

  if(use_geese){
    modfunc <- geese
  } else {
    modfunc <- glm
  }

  # Fit the model
  mod <- modfunc(formula=as.formula(formula), family=family, ..., data=phidat)

  start.time <- Sys.time()

  # Construct cumulative phi vector and cumulative covariance matrix
  cphi <- .matrix.phi(model=mod, family=family,
                      ts=ts_index$ts, ts_index=ts_index$index)

  c1 <- cphi$cphi
  c2 <- cphi$csum
  ts_index <- cphi$ts_index

  get.integral.est <- function(ts){
    # Convert bigT onto the grid
    val <- which.min(abs(ts - ts_index$ts))
    idx <- ts_index[val,]$index

    # Extract the cumulative phi at bigT (this is Omega est)
    # and the cumulative Cov(phi, phi) at (bigT, bigT) (this is Omega var)
    omega <- c1[c1$index == idx, phi]
    omega_var <- c2[c2$indexX == idx & c2$indexY == idx, rho]
    return(list(omega=omega, omega_var=omega_var, idx=idx))
  }

  est <- get.integral.est(bigT)

  results <- list(
    omega=est$omega,
    omega_var=est$omega_var
  )

  if(plot_phi){
    # PART 3B: Plot an estimated phi function and the prior testing data
    ts_plot <- seq(0, max(phidat$ui), by=dt)
    preds <- .predict.phi(ts=ts_plot,
                          model=mod, family=family, varcov=FALSE)[["point"]]
    plot(phidat$ri ~ phidat$ui,
         main=paste("Estimated Phi Function\n", as.character(formula),
                    "\n Omega:", round(est$omega*365.25), "days"),
         xlab="Infection Duration", ylab="Test Recent Probability")
    lines(preds ~ ts_plot, col="#4295f5", lwd=2)
    abline(v=bigT, lty="dashed")
  } else {
    end.time <- Sys.time()
    print(end.time - start.time)
  }

  if(return_all){
    results[["bigTidx"]] <- est$idx
    results[["ts_index"]] <- ts_index
    results[["get.integral.est"]] <- get.integral.est
    results[["cphi"]] <- cphi
  }

  return(results)
}

#' Function that summarizes a dataset + prior test results
#'
#' @param bigT The T^* parameter for defining recent infection.
#' @param use_geese Indicator to fit a gee model using geese(), or a glm model with glm()
#' @param formula Formula for fitting the model to phidat.
#'                Formula argument must take in only ui as the newdata.
#'                For example, do not create a ui^2 variable. Use poly(ui, ...) function
#'                to fit polynomial terms.
#' @param family Family argument for glm or gee
#' @param dt An integration step size. Should be no more than 1 day.
#' @param ptdf Data frame of prior testing data. These records should
#'             only be from people who are HIV positive. They need to have
#'             recency testing data as well as prior testing data.
#' @param n Required total sample size n if `ptdf` is used instead of `df`
#' @param n_p Required number of positives if `ptdf` is used instead of `df`
#' @param df Instead of passing ptdf, can pass a data frame with one row
#'           for each person, regardless of HIV status, and NA for all
#'           of the prior testing data for HIV negative people
#' @param phidat Data to fit a modfunc model to. Must have column names ri
#'               for recency indicator and ui for infection duration.
#'               Can have more columns, e.g., if an id is needed for geese.
#' @param ... Additional arguments to either glm or geese(..., data=phidat).
#' @export
summarizept.generator <- function(bigT, dt=1/365.25,
                                   use_geese=FALSE, formula, family,
                                   plot_phi=TRUE,
                                   ...){

  summarizept <- function(ptdf=NULL, n=NULL, n_p=NULL, df=NULL,
                           phidat=NULL){

    # PART 1: SUMMARIZE TRIAL DATA AND GET NEW RECENCY INDICATOR
    # ----------------------------------------------------------

    if(!is.null(ptdf)){

      if(is.null(n))   stop("must provide n with ptdf")
      if(is.null(n_p)) stop("must provide n_p with ptdr")
      if(!is.null(df)) stop("can only provide one of df or ptdf")

      s <- LISTKEY
      s[["n"]]   <- n
      s[["n_p"]] <- n_p
      s[["n_n"]] <- n - n_p

    } else {

      if(is.null(df))  stop("need to provide either df or ptdf")

      # First get the summary based on trial data
      s <- summarize.data(df)
      ptdf <- df[df$di == 1,]
    }

    s[["q"]] <- mean(!is.na(ptdf$ri))
    ptdf <- ptdf[!is.na(ptdf$ri),]

    # Generate recency indicators based on new algorithm
    newcols <- .get.Mi(ptdf, bigT)
    ptdf <- cbind(ptdf, newcols)

    # Calculate the new number of recents, and q effective
    # (proportion of positives w/ a prior test result)
    s[["n_r"]] <- sum(ptdf$ri)
    s[["n_r_pt"]] <- sum(ptdf$Mi)
    n_d <- sum(ptdf$qi)
    s[["q_eff"]] <- n_d / s[["n_p"]]

    # PART 2: GET AN ESTIMATED PHI FUNCTION AND COVARIANCE MATRIX
    # ALONG WITH ESTIMATED OMEGA AND OMEGA VAR
    # -----------------------------------------------------------

    maxT <- max(ptdf$ti, bigT, na.rm=TRUE)
    phi <- estimate.phi(
      phidat=phidat,
      maxT=maxT, bigT=bigT, dt=dt,
      formula=formula, family=family, use_geese=use_geese,
      plot_phi=FALSE, return_all=TRUE, ...
    )
    s[["omega"]] <- phi$omega
    s[["omega_var"]] <- phi$omega_var

    # PART 3: MAP TESTING TIMES TO A TIME GRID BASED ON DT
    # ----------------------------------------------------
    closest <- .map.to.tgrid(ts_index=phi$ts_index,
                             ptdf=ptdf,
                             bigT_index=phi$idx)

    # PART 4: ESTIMATE ALL OF THE OTHER THINGS THAT HAVE TO DO WITH PRIOR TESTING
    # ---------------------------------------------------------------------------
    ptlist <- .pt.properties.est(closest=closest, bigTidx=phi$bigTidx, cphi=phi$cphi)

    # Add all of the pt estimates to the main list
    for(el in names(ptlist)){
      s[[el]] <- ptlist[[el]]
    }

    return(s)
  }
  return(summarizept)
}

