#!/bin/sh

OUT="/home/students/mnorwood/hutch/sgeout/"
SHELL="/home/students/mnorwood/repos/XSRecency/results/shell.sh"
SCRIPT="simulate-parallel.R"
OUTPUT="/home/students/mnorwood/hutch/xs-recent"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

NSIMS=10000
N=5000
TAU=12
BIGT=2

ARGS="-n_sims=${NSIMS} -n=${N} -tau ${TAU} -bigT ${BIGT} -out_dir ${OUTPUT}"

SETTING1="-window 071 -shadow 080"

BASELINE="${CONSTANTS} ${ARGS} ${RHO} ${TYPE}"

BASELINE="${CONSTANTS} ${ARGS} ${TYPE}"

qsub ${BASELINE} ${SETTING1} -inc 0.032 -p 0.29 -itype constant
qsub ${BASELINE} ${SETTING1} -inc 0.030 -p 0.30 -itype constant
qsub ${BASELINE} ${SETTING1} -inc 0.028 -p 0.32 -itype constant

qsub ${BASELINE} ${SETTING1} -inc 0.032 -p 0.29 -itype constant -phi_tfrr 2
qsub ${BASELINE} ${SETTING1} -inc 0.030 -p 0.30 -itype constant -phi_tfrr 2
qsub ${BASELINE} ${SETTING1} -inc 0.028 -p 0.32 -itype constant -phi_tfrr 2
