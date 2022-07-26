#!/bin/sh

# Source shared parameters
source ./shared-parameters.sh recency-c
echo "${DT}: SCENARIO C with ${1} SIMS -- ${2}" >> "${OUTPUT}/RUN-LOG.txt"

OUTDIR="${OUTPUT}/${DT}"
mkdir ${OUTDIR}
NSIMS=$1
N=5000
P=0.29
INC=0.034
TAU=12
BIGT=2
RHO=0.0039
TMIN=0
TMAX=4

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTDIR} -last_point"
ASSAY="-window 101 -shadow 194"
TYPE="-itype piecewise"
TAIL="-phi_tfrr 2"
PT="-pt"
TS="-t_min ${TMIN} -t_max ${TMAX}"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT} ${TS}"

for Q in 0.5 1.0
do
    qsub ${CONSTANTS} ${BASELINE} -q $Q
    qsub ${CONSTANTS} ${BASELINE} -q $Q -exclude_pt_bigT
done
