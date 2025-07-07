#!/bin/bash
#PBS -N daesim2analysis
#PBS -l mem=24GB
#PBS -l ncpus=24
#PBS -l jobfs=24GB
#PBS -P xe2
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/xe2+gdata/v10+gdata/ka08
#PBS -q normal

wd=/g/data/xe2/ya6227/daesim2-analysis
source /g/data/xe2/ya6227/daesim2-analysis-env/bin/activate

cd $wd

n_processes=24
n_samples=100
dir_results=results
path_df_forcing_1=DAESIM_data/DAESim_forcing_data/DAESim_forcing_Milgadara_2018.csv
paths_df_forcing=("$path_df_forcing_1")
path_parameters_file=parameters/Fast1.json
crop=wheat

python3.9 notebooks/FAST.py \
  --crop $crop \
  --n_processes $n_processes \
  --n_samples $n_samples \
  --dir_results $dir_results \
  --paths_df_forcing "$(IFS=,; echo "${paths_df_forcing[*]}")" \
  --path_parameters_file $path_parameters_file
