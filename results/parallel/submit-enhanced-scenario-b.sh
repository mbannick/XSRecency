#!/bin/sh

# Source shared parameters
source ./shared-parameters.sh recency-b
echo "${DT}: SCENARIO B with ${1} SIMS -- ${2}" >> "${OUTPUT}/RUN-LOG.txt"

OUTDIR="${OUTPUT}/${DT}"
mkdir ${OUTDIR}
NSIMS=$1
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

TMIN=0
TMAX=4

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTDIR} -last_point"
ASSAY="-window 101 -shadow 194"
TYPE="-itype constant"
TAIL="-phi_tfrr 2"
PT="-pt"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"

for SEED in 100 2354 3482 8247 5893 48372 59547838 287484 13372 797683
do
  for Q in 0.5
  do
    for GAMMA in 0 0.083 0.5
    do
      for ETA in 0 0.1
      do
        PTARGS="-t_min $TMIN -t_max $TMAX -q $Q -gamma $GAMMA -eta $ETA -xi 0"
        qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -seed ${SEED}
      done
      for XI in 0.1
      do
        PTARGS="-t_min $TMIN -t_max $TMAX -q $Q -gamma $GAMMA -eta 0 -xi $XI"
        qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -seed ${SEED}
      done
    done
  done
done
