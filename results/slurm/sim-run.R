rm(list=ls())

library(data.table)
library(magrittr)
library(stringr)
library(flexsurv)
library(R.utils)

# Load XSRecency
library(XSRecency)
library(simtools)

source("./paramlists.R")
source("./sim-helpers.R")
source("./phi-func.R")
source("./arg-construct.R")

if(is.parallel()){
  a <- commandArgs(trailingOnly=TRUE)
  OUTDIR <- a[1]
  TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  PARAMS <- read.csv(sprintf("%s/params.csv", OUTDIR))
} else {
  OUTDIR <- "."
  TASKID <- 1
  PARAMS <- data.table(do.call(data.table, DEFAULTS))
}

# Set up parameter grabber and set streaming seed for
# reproducibility
NTASKS <- nrow(PARAMS)
gp <- get.param.closure(TASKID, PARAMS)
set.stream.seed(gp("seed"), NTASKS, TASKID)

ARGS <- gp()

# Get the baseline phi function
ARGS[["phi.func"]] <- doCall(get.phi, args=gp())

# Get the epi functions
epi.funcs <- doCall(get.epi.funcs, args=gp())
ARGS[["infection.function"]] <- epi.funcs[["infection.function"]]
ARGS[["inc.function"]] <- epi.funcs[["inc.function"]]

# Get external data
if(!is.null(gp("duong_scale"))){
  df <- copy(XSRecency:::duong)
  df <- df[, days := days * as.numeric(gp("duong_scale"))]
  df <- df[, last.time := shift(days), by="id.key"]
  df <- df[, gap := days - last.time]
} else {
  df <- NULL
}
ARGS[["ext_df"]] <- df

ARGS[["prevalence"]] <- gp("p")
ARGS[["baseline_incidence"]] <- gp("inc")

if(!gp("pt")){

  sim <- doCall(simulate, args=ARGS)

} else {

  # Specific prior testing settings
  ARGS[["ptest.dist"]] <- function(u) runif(1, gp("t_min"), gp("t_max"))
  ARGS[["ptest.prob"]] <- function(u) gp("q")

  if(gp("mech2")){
    GAMMA_PARMS <- c(1.57243557, 1.45286770, -0.02105187)
    ptest.dist2 <- function(u) u - rgengamma(n=1,
                                             mu=GAMMA_PARMS[1],
                                             sigma=GAMMA_PARMS[2],
                                             Q=GAMMA_PARMS[3])
  } else {
    ptest.dist2 <- NULL
  }
  ARGS[["ptest.dist2"]] <- ptest.dist2

  if(!is.null(gp("gamma"))){
    t_noise <- function(t) max(0, t + rnorm(n=1, sd=gp("gamma")))
  } else {
    t_noise <- NULL
  }

  ARGS[["t_min"]] <- 0
  ARGS[["t_max"]] <- 100

  ARGS[["d_misrep"]] <- gp("eta")
  ARGS[["q_misrep"]] <- gp("nu")
  ARGS[["p_misrep"]] <- gp("xi")

  start <- Sys.time()
  sim <- doCall(simulate.pt, args=ARGS)
  end <- Sys.time()
  print(end - start)
}

df <- do.call(cbind, sim) %>% data.table
df[, sim := .I]
df[, simstart := gp("start_sim")]
df[, TASKID := TASKID]

filename <- paste0("results-", TASKID)
write.csv(df, file=paste0(OUTDIR, "/", filename, ".csv"))
