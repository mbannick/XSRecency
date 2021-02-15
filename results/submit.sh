#!/bin/sh

OUT="/home/students/mnorwood/hutch/sgeout/"
SHELL="/home/students/mnorwood/repos/XSRecency/results/shell.sh"
SCRIPT="simulate-parallel.R"
OUTPUT="/home/students/mnorwood/hutch/xs-recent"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

NSIMS=1000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

ARGS="-n_sims ${NSIMS} -n ${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT}"

SETTING1="-window 071 -shadow 080"
SETTING2="-window 248 -shadow 306"

NORMAL="-phi_norm_mu 7 -phi_norm_sd 1 -phi_norm_div 8"

for t in "constant" "linear" "exponential"
do
  if [ t = "linear" ]
  then
    RHO=0.0028
    ARGS="${ARGS} -rho ${RHO}"
  fi
  if [ t = "exponential" ]
  then
    RHO=0.07
    ARGS="${ARGS} -rho ${RHO}"
  fi

  TYPE="-itype ${t}"

  BASELINE="${CONSTANTS} ${ARGS} ${TYPE}"

  # echo qsub ${BASELINE} ${SETTING1}
  qsub ${BASELINE} ${SETTING1}
  qsub ${BASELINE} ${SETTING1} -phi_tfrr 2
  qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL}

  qsub ${BASELINE} ${SETTING2}
  qsub ${BASELINE} ${SETTING2} -phi_frr 0.02
  qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL}
done
