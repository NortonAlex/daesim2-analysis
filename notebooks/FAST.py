import numpy as np
## DAESIM2
from daesim.climate import *
from daesim.utils import ODEModelSolver
from daesim.utils import daesim_io_write_diag_to_nc
## DAESIM2 Analysis
from daesim2_analysis import fast_sensitivity as fastsa
from SALib.sample import fast_sampler
from SALib.analyze import fast
## Argument Parsing
import __main__
from argparse import ArgumentParser
## Parallelisation
from typing_extensions import Callable
from typing_extensions import Self
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Optional
from typing import List
from os import makedirs
from tqdm import tqdm
from time import time
from functools import partial

from daesim2_analysis.experiment import Experiment
from daesim2_analysis.experiment import is_interactive
from daesim2_analysis.parameters import Parameters
from daesim2_analysis.utils import *
from daesim2_analysis.forcing_data import ForcingData

experiment = Experiment.from_cli() if not is_interactive() else Experiment()

# %%
parameters: Parameters = experiment.parameters
param_values = parameters.sample(experiment.n_samples)

# %%
Mpx = []
iparamsets_per_run = get_iparamsets_per_run(param_values, experiment.n_processes)
n_runs = len(iparamsets_per_run)
Mpx_column_headers = "nparamset,W_P_peakW,W_L_peakW,W_R_peakW,W_S_peakW,W_S_spike0,W_S_anth0,GPP_int_seas,NPP_int_seas,Rml_int_seas,Rmr_int_seas,Rg_int_seas,trflux_int_seas,FCstem2grain_int_seas,NPP2grain_int_seas,E_int_seas,LAI_peakW,W_spike_anth1,GY_mature,Sdpot_mature,GN_mature" 

def evaluate_iparamset(iparamset: int):
    nparamset = iparamset + 1
    paramset = param_values[iparamset]
    # model_output = fastsa.update_and_run_model(paramset, PlantX, input_data, parameters_df, problem)
    model_output = fastsa.update_and_run_model(
        paramset,
        experiment.PlantX,
        experiment.input_data,
        parameters.df,
        parameters.problem
    )
    Mpxi, diagnostics = model_output[0], model_output[1]

    nsigfigures = len(str(np.shape(param_values)[0]))
    filename_write = f"FAST_results_{experiment.xsite}_paramset{nparamset:0{nsigfigures}}.nc"
    daesim_io_write_diag_to_nc(
      experiment.PlantX,
      diagnostics,
      experiment.dir_xsite_parameters + '/',
      filename_write,
      experiment.ForcingDataX.time_index,
      problem=parameters.problem,
      param_values=paramset,
      nc_attributes={'title': experiment.title, 'description': experiment.description}
    )
    return np.insert(Mpxi, 0, nparamset)

if __name__ == '__main__':
    for n_run in range(n_runs):
        with Pool(processes=experiment.n_processes) as pool:
            Mpx += pool.map(evaluate_iparamset, iparamsets_per_run[n_run])
