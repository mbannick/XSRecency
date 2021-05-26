# Submit Parallel Simulations

## Bash Scripts

These scripts are used to submit simulations in parallel on an SGE-type cluster. They submit simulations that are parallelized by setting, rather than simulation number.

- `shell.sh`: A shell script used by all of the submit bash scripts below.
- `submit-main.sh`: This submits simulations that were used for the main tables and figures in the paper.
- `submit-frr.sh`: This submits simulations with settings for the sensitivity analysis of how to calculate FRR.
- `submit-nsamps.sh`: This is essentially identical to `submit-main.sh` but with an argument for sample size. To reproduce the sample size plot in the appendix, this was run with N = 2000, 5000, 10000, 15000, and 20000.

## Summarize Results

Each of the simulation settings will output one file in the designated folder.
To summarize the results in terms of bias and coverage, you can use the script `summarize.R`.

