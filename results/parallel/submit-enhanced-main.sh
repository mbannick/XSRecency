#!/bin/sh

# Source shared parameters
source ./shared-parameters.sh recency-main
echo "${DT}: MAIN with ${1} SIMS -- ${2}" >> "${OUTPUT}/RUN-LOG.txt"

OUTDIR="${OUTPUT}/${DT}"
mkdir ${OUTDIR}
NSIMS=$1
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTDIR} -last_point"
ASSAY="-window 101 -shadow 194"
TYPE="-itype constant"
TAIL="-phi_tfrr 2"
PT="-pt"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"

for SEED in 100 2354 3482 8247 5893
do
  for Q in 0.2 0.4 0.6 0.8 1.0
  do
    for TRANGE in "0 4" "0 2" "2 4"
    do
      set -- $TRANGE
      PTARGS="-t_min $1 -t_max $2 -q $Q"
      # echo "qsub ${CONSTANTS} ${BASELINE} ${PTARGS}"
      qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -seed ${SEED}
    done
  done
done
