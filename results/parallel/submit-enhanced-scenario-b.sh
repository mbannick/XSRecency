#!/bin/sh

sh shared-parameters.sh

NSIMS=1000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

TMIN=1
TMAX=3

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT} -last_point"
ASSAY="-window 248 -shadow 306"
TYPE="-itype constant"
TAIL="-phi_tfrr 2"
PT="-pt"

BASELINE="${ARGS} ${RHO} ${TYPE} ${ASSAY} ${TAIL} ${PT}"

for GAMMA in 0 0.1 0.5
do
  for ETA in 0.0 0.1 0.25
  do
    for NU in 0.0 0.1 0.25
    do
      for Q in 0.5 1.0
      do
        PTARGS="-t_min {$TMIN} -t_max {$TMAX} -q {$Q} -gamma {$GAMMA} -eta {$ETA} -nu {$NU}"
        qsub ${CONSTANTS} ${BASELINE} ${PTARGS}
