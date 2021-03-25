# Submit Parallel Simulations

## Bash Scripts

These scripts are used to submit simulations in parallel on an SGE-type cluster. They submit simulations that are parallelized by setting, rather than simulation number.

- `shell.sh`: a shell script used by all of the submit bash scripts
- `submit-main.sh`: this submits simulations that were used for the main tables and figures in the paper
- `submit-frr.sh`: this submits simulations with settings for the sensitivity analysis of how to calculate frr
- `submit-epi.sh`: this submits simulations for varying incidence and prevalence, not included in the paper

## Summarize Results

Each of the simulation settings will output one file. To summarize the results in terms of bias and coverage, you can use the script `summarize.R`.
