#!/bin/sh

sh shared-parameters.sh

NSIMS=1000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

TMIN=0
TMAX=$TAU

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT} -last_point"
ASSAY="-window 248 -shadow 306"
TYPE="-itype constant"
TAIL="-phi_tfrr 2"
PT="-pt"

MECH2="-mech2"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"

for Q in 0.2 0.4 0.6 0.8 1.0
do
  PTARGS="-t_min {$TMIN} -t_max {$TMAX} -q {$Q}"
  qsub ${CONSTANTS} ${BASELINE} ${PTARGS}
  qsub ${CONSTANTS} ${BASELINE} ${PTARGS} ${MECH2}
