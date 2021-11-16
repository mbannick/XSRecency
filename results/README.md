# Scripts to run the simulations and create tables and figures

The contents of this folder are:

- `parallel`: Folder of the bash scripts to submit simulations.
- `simulate-parallel.R`: Main simulation script that is called in parallel over many simulation settings. If you run this script
by itself (rather than calling it in parallel), it will run with whatever the default arguments are in `commandArgs`.
- `sim-helpers.R`: Functions called by `simulate-parallel.R` script, separated for organization purposes.
- `tables-figures`: Code to create the tables and figures that are output from running the `parallel` scripts, and then summarizing them with `parallel/summarize.R` within a folder. To create these tables and figures, you have to specify a version (i.e. a file path to a folder with the version of simulation results that you want to read in). Also includes the code to run the data analysis.
