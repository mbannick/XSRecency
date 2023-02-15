#!/bin/sh

# Source shared parameters
source ./shared-parameters.sh recency-c
echo "${DT}: SCENARIO C with ${1} SIMS -- ${2}" >> "${OUTPUT}/RUN-LOG.txt"

OUTDIR="${OUTPUT}/${DT}"
mkdir ${OUTDIR}
NSIMS=$1
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2
RHO=0.0039
TMIN=0
TMAX=4

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTDIR} -last_point"
ASSAY="-window 101 -shadow 194"
TAIL="-phi_tfrr 2"
PT="-pt"

BASELINE="${ARGS} ${ASSAY} ${TAIL} ${PT}"

for SEED in 100 2354 3482 8247 5893
do
  for Q in 0.5
  do
    for TMAX in 4 $TAU
    do
      qsub ${CONSTANTS} ${BASELINE} -q $Q -itype constant -t_min $TMIN -t_max $TMAX -seed ${SEED}
      qsub ${CONSTANTS} ${BASELINE} -q $Q -itype constant -exclude_pt_bigT -t_min $TMIN -t_max $TMAX -seed ${SEED}
      qsub ${CONSTANTS} ${BASELINE} -q $Q -itype piecewise -rho ${RHO} -t_min $TMIN -t_max $TMAX -seed ${SEED}
      qsub ${CONSTANTS} ${BASELINE} -q $Q -itype piecewise -rho ${RHO} -exclude_pt_bigT -t_min $TMIN -t_max $TMAX -seed ${SEED}
    done
  done
done
