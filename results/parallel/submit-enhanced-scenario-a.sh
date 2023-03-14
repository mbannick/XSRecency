#!/bin/sh

# Source shared parameters
source ./shared-parameters.sh recency-a
echo "${DT}: SCENARIO A with ${1} SIMS -- ${2}" >> "${OUTPUT}/RUN-LOG.txt"

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

MECH2="-mech2"
TMINEXC="-t_min_exclude 0.25"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"

for SEED in 100 2354 3482 8247 5893 48372 59547838 287484 13372 797683
do
  for Q in 0.5
  do
    PTARGS="-t_min $TMIN -t_max $TMAX -q $Q"
    qsub ${CONSTANTS} ${BASELINE} ${PTARGS} -seed ${SEED}
    qsub ${CONSTANTS} ${BASELINE} ${PTARGS} ${MECH2} -seed ${SEED}
    qsub ${CONSTANTS} ${BASELINE} ${PTARGS} ${MECH2} -seed ${SEED}
  done
done
