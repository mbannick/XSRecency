#!/bin/sh

OUT="/home/students/mnorwood/hutch/sgeout/"
SHELL="/home/students/mnorwood/repos/XSRecency/results/shell.sh"
SCRIPT="simulate-parallel.R"
OUTPUT="/home/students/mnorwood/hutch/xs-recent"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

NSIMS=5000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT}"

SETTING1="-window 071 -shadow 080"
SETTING2="-window 248 -shadow 306"

TYPE="-itype constant"
NORMAL="-phi_norm_mu 7 -phi_norm_sd 1 -phi_norm_div 8"
FRR="-frr_mix_start 2 -frr_mix_end 5"

BASELINE="${CONSTANTS} ${ARGS} ${TYPE}"

qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL}
qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL}

qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL} ${FRR}
qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL} ${FRR}
