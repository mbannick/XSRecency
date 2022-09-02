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

for SD in 0715 0831 2022 1994 0410 0616 0320 9319 0810 0533 1850 1380 0926 1720 8919 3000 1000 4010 9884 0008
do
  for TRANGE in "0 2"
  do
    set -- $TRANGE
    PTARGS="-t_min $1 -t_max $2 -q $Q -seed $SD"
    # echo "qsub ${CONSTANTS} ${BASELINE} ${PTARGS}"
    qsub ${CONSTANTS} ${BASELINE} ${PTARGS}
  done
done
