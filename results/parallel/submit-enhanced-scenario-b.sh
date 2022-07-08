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

TMIN=1
TMAX=3
Q=1

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTDIR} -last_point"
ASSAY="-window 101 -shadow 194"
TYPE="-itype constant"
TAIL="-phi_tfrr 2"
PT="-pt"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"

for GAMMA in 0 0.1 0.5
do
  for NU in 0.0 0.1 0.25
  do
    PTARGS="-t_min $TMIN -t_max $TMAX -q $Q -gamma $GAMMA -nu $NU"
    qsub ${CONSTANTS} ${BASELINE} ${PTARGS}
  done
done

for GAMMA in 0 0.1 0.5
do
  for XI in 0.0 0.1 0.25
  do
    PTARGS="-t_min $TMIN -t_max $TMAX -q $Q -gamma $GAMMA -xi $XI"
    qsub ${CONSTANTS} ${BASELINE} ${PTARGS}
  done
done
