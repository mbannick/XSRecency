#!/bin/sh

OUT="/home/users/mnorwood/hutch/sgeout/"
SHELL="/home/users/mnorwood/repos/XSRecency/results/parallel/shell.sh"
SCRIPT="/home/users/mnorwood/repos/XSRecency/results/simulate-parallel.R"
OUTPUT="/home/users/mnorwood/hutch/xs-recent/main"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

NSIMS=5000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT} -last_point"

SETTING1="-window 101 -shadow 194"
SETTING2="-window 248 -shadow 306"

NORMAL="-phi_norm_mu 7 -phi_norm_sd 1 -phi_norm_div 8"
PNORMAL="-phi_pnorm_mu 10 -phi_pnorm_sd 2 -phi_pnorm_div 10"

for t in constant linear exponential
do
  if [ ${t} = "constant" ]
  then
    RHO=""
  fi

  if [ ${t} = "linear" ]
  then
    RHO="-rho=0.0028"
  fi

  if [ ${t} = "exponential" ]
  then
    RHO="-rho=0.07"
  fi

  TYPE="-itype ${t}"

  BASELINE="${CONSTANTS} ${ARGS} ${RHO} ${TYPE}"

  # echo qsub ${BASELINE} ${SETTING1}
  qsub ${BASELINE} ${SETTING1}
  qsub ${BASELINE} ${SETTING1} -phi_tfrr 2
  qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL}
  qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${PNORMAL}

  qsub ${BASELINE} ${SETTING2}
  qsub ${BASELINE} ${SETTING2} -phi_frr 0.02
  qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL}
  qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${PNORMAL}
done
