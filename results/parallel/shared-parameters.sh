#!/bin/sh

OUT="/home/users/mnorwood/hutch/sgeout/"
SHELL="/home/users/mnorwood/repos/XSRecency/results/parallel/shell.sh"
SCRIPT="/home/users/mnorwood/repos/XSRecency/results/simulate-parallel.R"
OUTPUT="/home/users/mnorwood/hutch/xs-recent/enhanced"

CONSTANTS="-cwd -N recency_sim -j y -o ${OUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

DT=$(date '+%d/%m/%Y %H:%M:%S')

NSIMS=$1
