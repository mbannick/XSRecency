# library(magrittr)
# library(data.table)
# library(R.utils)
#
# source("../parse-args.R")
#
# # GET COMMAND-LINE ARGUMENTS FOR THE DIFFERENT SCENARIOS OF
# # DATA GENERATION TO SUBMIT
# # ---------------------------------------------------------
#
# # Read in command line arguments, with the following defaults
# # Command-line arguments are for data generation only
# # because all estimators will be applied to each scenario
# args <- commandArgs(
#   trailingOnly=TRUE, asValues=TRUE,
#   defaults=DEFAULTS
# )
#
# # CONSTRUCT OUTPUT/ERROR DIRECTORIES,
# # WRITE DESCRIPTION AND GIT HASH TO A FILE
# # ----------------------------------------------
#
# # Output directory
# theday <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
# OUT_DIR <- paste0("~/g-calibration/", theday)
# if(dir.exists(OUT_DIR)) stop("Output directory already exists.")
# dir.create(OUT_DIR, recursive=TRUE)
#
# ERROR <- sprintf("%s/output", OUT_DIR)
# dir.create(ERROR, recursive=TRUE)
#
# RESULTS <- sprintf("%s/results", OUT_DIR)
# dir.create(RESULTS, recursive=TRUE)
#
# # Get the current git hash
# hash <- system("git rev-parse HEAD", intern=TRUE)
#
# # Install latest version of RobinCar and collect hash
# robin.git <- "--git-dir=~/repos/RobinCar/.git --work-tree=~/repos/RobinCar"
# if(args$update_robin == "yes"){
#   system(paste0("git ", robin.git, " pull"))
#   system("R-4.0.1 CMD INSTALL ~/repos/RobinCar")
# }
# hash.RobinCar <- system(paste0("git ", robin.git, " rev-parse HEAD"), intern=TRUE)
#
# # Write description to a text file
# fileConn <- file(sprintf("%s/DESCRIPTION.txt", OUT_DIR))
# writeLines(c(
#   "DESCRIPTION: ",
#   args$desc,
#   "GIT HASH: ",
#   hash,
#   "ROBINCAR HASH: ",
#   hash.RobinCar,
#   "DATE: ",
#   theday
# ),
# fileConn)
# close(fileConn)
#
# # GET THE PARAMETERS IN A GRID AND SAVE GRID TO A CSV
# # FOR THE JOB ARRAY TO READ IN LATER
# # ---------------------------------------------------
#
# # Number of simulations to run
# SEED <- as.integer(args$seed)
# N_SIMS <- as.integer(args$n_sims)
# SIMBLOCK <- as.integer(args$sim_blocksize)
#
# as.fraction <- function(x){
#   sapply(strsplit(x, split = "/"),
#          function(u) as.numeric(u[1]) / as.numeric(u[2]))
# }
#
# start_sims <- seq(1, N_SIMS, SIMBLOCK)
#
# # Parameter grid from arguments
# params <- list(
#   n             = parse.args(args$n, as.integer),
#   p             = parse.args(args$p, as.character),
#   case          = parse.args(args$case, as.integer),
#   car           = parse.args(args$car, as.character),
#   p_trt         = parse.args(args$p_trt, as.fraction),
#   ml            = parse.args(args$ml, as.character),
#   blocksize     = parse.args(args$blocksize, as.integer),
#   sim_blocksize = parse.args(args$sim_blocksize, as.integer),
#   start_sim     = start_sims
# )
#
# # Save parameter list and number of tasks
# param_grid <- expand.grid(params)
# param_grid <- data.table(param_grid)
#
# N_JOBS <- nrow(param_grid)
# write.csv(param_grid, file=sprintf("%s/params.csv", OUT_DIR))
#
# # SUBMIT JOB ARRAY
# # -------------------------------------------------
#
# command <- sprintf("sbatch --array=1-%s shell.sh simulate.R %s %s %s",
#                    N_JOBS, OUT_DIR, N_SIMS, SEED)
# print(command)
# system(command)
#
#
# OUT <- "/home/users/mnorwood/hutch/sgeout/"
# SHELL <- "/home/users/mnorwood/repos/XSRecency/results/parallel/shell.sh"
# SCRIPT <- "/home/users/mnorwood/repos/XSRecency/results/simulate-parallel.R"
# OUTPUT <- "/home/users/mnorwood/hutch/xs-recent/enhanced"
#
# CONSTANTS <- "-cwd -N ${1} -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"
#
# DT <- $(date '+%d-%m-%Y-%H-%M-%S')
#
# # Source shared parameters
# source ./shared-parameters.sh recency-a
# echo "${DT}: SCENARIO A with ${1} SIMS -- ${2}" >> "${OUTPUT}/RUN-LOG.txt"
#
#
#
# ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTDIR} -last_point"
# ASSAY="-window 101 -shadow 194"
# TYPE="-itype constant"
# TAIL="-phi_tfrr 2"
# PT="-pt"
#
# MECH2="-mech2"
# TMINEXC="-t_min_exclude 0.25"
#
# BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"
#
# for SEED in 100 2354 3482 8247 5893 48372 59547838 287484 13372 797683
# do
# for Q in 0.5
# do
# PTARGS="-t_min $TMIN -t_max $TMAX -q $Q"
# qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -seed ${SEED}
# qsub ${CONSTANTS} ${BASELINE} ${PTARGS} ${MECH2} -seed ${SEED}
# qsub ${CONSTANTS} ${BASELINE} ${PTARGS} ${MECH2} -seed ${SEED}
# done
# done
