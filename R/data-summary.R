
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
      idx <- ts_index[which.min(abs(t - ts_index$ts)), "index"]
    }
  }
  closest_index <- sapply(closest$ts_orig, get.closest.index)
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
summarize.pt.generator <- function(bigT, dt=1/365.25,
                                   use_geese=FALSE, formula, family,
                                   plot_phi=TRUE,
                                   ...){

  summarize.pt <- function(ptdf=NULL, n=NULL, n_p=NULL, df=NULL,
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

    # Generate recency indicators based on new algorithm
    newcols <- .get.Mi(ptdf, bigT)
    ptdf <- cbind(ptdf, newcols)

    # Calculate the new number of recents, and q effective
    # (proportion of positives w/ a prior test result)
    s[["n_r"]] <- sum(ptdf$ri)
    s[["n_r_pt"]] <- sum(ptdf$Mi)
    n_d <- sum(ptdf$qi)
    s[["q_eff"]] <- n_d / s[["n_p"]]

    # PART 2: MAP TESTING TIMES TO A TIME GRID BASED ON DT
    # ----------------------------------------------------

    ts_index <- .get.ts(minT=0, maxT=max(ptdf$ti, bigT, na.rm=TRUE), dt=dt)

    # Convert bigT onto the grid
    bigT_index <- ts_index[which.min(abs(bigT - ts_index$ts)),]

    closest <- .map.to.tgrid(ts_index=ts_index,
                             ptdf=ptdf,
                             bigT_index=bigT_index)

    # PART 3: GET AN ESTIMATED PHI FUNCTION AND COVARIANCE MATRIX
    # -----------------------------------------------------------

    # Make sure model works based on the user arguments
    .test.phi(phidat=phidat, use_geese=use_geese, formula=formula, family=family, ...)

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
    if(plot_phi){
      # PART 3B: Plot an estimated phi function and the prior testing data
      ts_plot <- seq(0, max(phidat$ui), by=dt)
      preds <- .predict.phi(ts=ts_plot,
                            model=mod, family=family, varcov=FALSE)[["point"]]
      plot(phidat$ri ~ phidat$ui, main="Estimated Phi Function",
           xlab="Infection Duration", ylab="Test Recent Probability")
      lines(preds ~ ts_plot, col="#4295f5", lwd=2)
    } else {
      end.time <- Sys.time()
      print(end.time - start.time)
    }

    # PART 4: ESTIMATE OMEGA_T^*
    # ----------------------
    bigTidx <- bigT_index$index

    # Extract the cumulative phi at bigT (this is Omega est)
    # and the cumulative Cov(phi, phi) at (bigT, bigT) (this is Omega var)
    c1 <- cphi$cphi
    c2 <- cphi$csum

    s[["omega"]] <- c1[c1$index == bigTidx, phi]
    s[["omega_var"]] <- c2[c2$indexX == bigTidx & c2$indexY == bigTidx, rho]

    # PART 5: ESTIMATE ALL OF THE OTHER THINGS THAT HAVE TO DO WITH PRIOR TESTING
    # ---------------------------------------------------------------------------
    ptlist <- .pt.properties.est(closest=closest, bigTidx=bigTidx, cphi=cphi)

    # Add all of the pt estimates to the main list
    for(el in names(ptlist)){
      s[[el]] <- ptlist[[el]]
    }

    return(s)
  }
  return(summarize.pt)
}

