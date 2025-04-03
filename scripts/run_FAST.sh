wd=/home/272/ya6227/daesim2-analysis
source /g/data/xe2/ya6227/daesim2-analysis-env/bin/activate
cd $wd

n_processes=4
n_samples=100
dir_results=results
path_df_forcing_1=DAESIM_Data/DAESim_forcing_data/Rutherglen_1971.csv
paths_df_forcing=("$path_df_forcing_1")
path_parameters_file='parameters/Fast1.json'

python3.9 notebooks/FAST_sensitivity_analysis_parallel.py \
  --n_processes $n_processes \
  --n_samples $n_samples \
  --dir_results $dir_results \
  --paths_df_forcing "$(IFS=,; echo "${paths_df_forcing[*]}")" \
  --path_parameters_file $path_parameters_file
