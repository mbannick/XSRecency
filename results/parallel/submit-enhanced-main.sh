#!/bin/sh

OUT="/home/users/mnorwood/hutch/sgeout/"
SHELL="/home/users/mnorwood/repos/XSRecency/results/parallel/shell.sh"
SCRIPT="/home/users/mnorwood/repos/XSRecency/results/simulate-parallel.R"
OUTPUT="/home/users/mnorwood/hutch/xs-recent/enhanced"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

NSIMS=1000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT} -last_point"
ASSAY="-window 248 -shadow 306"
TYPE="-itype constant"
TAIL="-phi_tfrr 2"
PT="-pt"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"
OLDIFS=$IFS; IFS=','
for GAMMA in 0 0.1
do
  for Q in 0.2 0.4 0.6 0.8 1.0
  do
    for ETA in 0.0 0.25
    do
      for NU in 0.0 0.25
      do
        for TRANGE in 0,4 0,2 1,3
        do set -- $TRANGE;
          PTARGS="-t_min {$1} -t_max {$2} -q {$Q} -gamma {$GAMMA} -eta {$ETA} -nu {$NU}"
          IFS=$OLDIFS
          qsub ${CONSTANTS} ${BASELINE} ${PTARGS}
          # qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -qu_int 1 -qu_slope 1
          # qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -tu_int 1 -tu_slope 1
          # qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -qu_int 1 -qu_slope 1 -tu_int 1 -tu_slope 1
        done
      done
    done
  done
done
