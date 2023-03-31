library(simtools)
library(R.utils)

source("./paramlists.R")

# Read in command line arguments -----------------
args <- commandArgs(
  trailingOnly=TRUE, asValues=TRUE,
  defaults=list(
    seed="1",
    n_sims="1000",
    sim_blocksize="10",
    scenario="a",
    desc="blank description"
  )
)

BASEDIR <- "~/xs-enhanced"

# Construct parameter grid based on the scenario ------

n_sims <- as.integer(args$n_sims)
blocksize <- as.integer(args$sim_blocksize)

OVERALL <- DEFAULTS
OVERALL[["n_sims"]] <- n_sims
OVERALL[["sim_blocksize"]] <- blocksize
OVERALL[["seed"]] <- as.integer(args$seed)

startsims <- get.startsims(n_sims, blocksize)

if(args$scenario == "a"){
  PARAMS <- get.paramlist.a(startsims)
} else if(args$scenario == "b"){
  PARAMS <- get.paramlist.b(startsims)
} else if(args$scenario == "c"){
  PARAMS <- get.paramlist.c(startsims)
} else if(args$scenario == "d"){
  PARAMS <- get.paramlist.d(startsims)
} else {
  stop("Unrecognized scenario")
}

PARAMS <- fill.overall(PARAMS, OVERALL)

# Setup directories ------------------------------
outdir <- setup.directories(BASEDIR, args$desc)

# Save parameter list -----------------------
njobs <- save.params(PARAMS, outdir)

# Submit job array
submit.jobs(njobs, "shell.sh", "simulate-parallel.R", outdir)
