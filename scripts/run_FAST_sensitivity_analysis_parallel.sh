#!/bin/bash
#PBS -N daesim2analysis
#PBS -l mem=24GB
#PBS -l ncpus=24
#PBS -l jobfs=24GB
#PBS -P xe2
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/xe2+gdata/v10+gdata/ka08
#PBS -q normal

wd=/home/272/ya6227/daesim2-analysis
source /g/data/xe2/ya6227/daesim2-analysis-env/bin/activate
cd $wd
n_processes=24
n_samples=600
dir_results=/g/data/xe2/ya6227/daesim2-analysis-data/FAST_results
path_df_forcing_1=/g/data/xe2/ya6227/daesim2-analysis-data/DAESim_forcing_data/Rutherglen_1971.csv
paths_df_forcing=("$path_df_forcing_1")

python3.9 notebooks/FAST_sensitivity_analysis_parallel.py \
  --n_processes $n_processes \
  --n_samples $n_samples \
  --dir_results $dir_results \
  --paths_df_forcing "$(IFS=,; echo "${paths_df_forcing[*]}")"
  
