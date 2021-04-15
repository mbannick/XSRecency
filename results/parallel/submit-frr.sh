#!/bin/sh

OUT="/home/students/mnorwood/hutch/sgeout/"
SHELL="/home/students/mnorwood/repos/XSRecency/results/parallel/shell.sh"
SCRIPT="/home/students/mnorwood/repos/XSRecency/results/simulate-parallel.R"
OUTPUT="/home/students/mnorwood/hutch/xs-recent/sens-frr"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

NSIMS=5000
N=5000
P=0.29
INC=0.032
TAU=12
BIGT=2

ARGS="-n_sims=${NSIMS} -n=${N} -p ${P} -inc ${INC} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT}"

SETTING1="-window 101 -shadow 194"
SETTING2="-window 248 -shadow 306"

TYPE="-itype constant"
NORMAL="-phi_norm_mu 7 -phi_norm_sd 1 -phi_norm_div 8"

BASELINE="${CONSTANTS} ${ARGS} ${TYPE}"

qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL}
qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL}

qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL} -add_unif 5
qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL} -add_unif 5

qsub ${BASELINE} ${SETTING1} -phi_tfrr 2 ${NORMAL} "-ext_FRR"
qsub ${BASELINE} ${SETTING2} -phi_frr 0.02 ${NORMAL} "-ext_FRR"
