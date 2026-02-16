# %%
import numpy as np
import pandas as pd

## DAESIM2
from daesim.climate import *
from daesim.utils import ODEModelSolver
from daesim.utils import daesim_io_write_diag_to_nc

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
from time import time
from functools import partial
from pandas import Timestamp

## daesim2-analysis
from daesim2_analysis.experiment import Experiment
from daesim2_analysis.experiment import is_interactive
from daesim2_analysis.parameters import Parameters
from daesim2_analysis.utils import *
from daesim2_analysis.forcing_data import ForcingData
from daesim2_analysis.run import update_and_run_model

# %% [markdown]
# ### Define the `Experiment`

# %%
# experiment = Experiment.from_cli() if not is_interactive() else Experiment()
experiment = Experiment(
    df_forcing="../DAESIM_data/DAESim_forcing_data/Rutherglen_1971-1972_wheat_exp_v2.csv",
    parameters="../parameters/paramsRutherglen1971_wRootVcmaxOn_GYonly.json",
    daesim_config="../daesim_configs/DAESIM_site_Rutherglen_exp_wRootVcmaxOn_GYonly.json",
    dir_results="../results/",
    df_forcing_type='3',
    CLatDeg=-36.05,
    CLonDeg=146.50,
    sowing_dates=[Timestamp(year=1971,month=5,day=11), Timestamp(year=1972,month=5,day=19)],
    harvest_dates=[Timestamp(year=1971,month=12,day=23), Timestamp(year=1972,month=12,day=13)],
    # xsite="Rutherglen1971-72_lhs7",
    xsite="Rutherglen1971-72_tmp",
    n_samples=30,
    n_processes=2,
    # n_processes=1,
)

# %% [markdown]
# ### Define the `Parameters`

# %%
parameters: Parameters = experiment.parameters
#param_values = parameters.sample(experiment.n_samples)
param_values = parameters.sample(n=experiment.n_samples, method="lhs") #, seed=123)


# %% [markdown]
# ### Any other customized changes to the configuration

# %%
#experiment.PlantX.PlantCH2O.root_scale_vcmax = True

# %% [markdown]
# ### Customized version of `run_model_and_get_outputs` (distinct from the standard function defined in `run.py`)
#
# This allows you to customize the post-processing of standard model output so that you can output other metrics (e.g. objective function values) or diagnostics. 

# %%
def run_model_and_get_outputs(Plant, ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices):
    ## Define the callable calculator that defines the right-hand-side ODE function
    PlantCalc = Plant.calculate
    
    Model = ODEModelSolver(calculator=PlantCalc, states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_start=time_axis[0], log_diagnostics=True)

    Model.reset_diagnostics()
    
    ## Run the model solver
    res = Model.run(
        time_axis=time_axis,
        forcing_inputs=forcing_inputs,
        solver="euler",
        zero_crossing_indices=zero_crossing_indices,
        reset_days=reset_days,
    )

    # Convert the defaultdict to a regular dictionary
    _diagnostics = dict(Model.diagnostics)
    # Convert each list in the dictionary to a NumPy array
    diagnostics = {key: np.array(value) for key, value in _diagnostics.items()}

    # Convert the array to a numeric type, handling mixed int and float types
    diagnostics['idevphase_numeric'] = np.array(diagnostics['idevphase'],dtype=np.float64)
    
    # In the model idevphase can equal None but that is not useable in post-processing, so we set None values to np.nan
    diagnostics["idevphase_numeric"][diagnostics["idevphase"] == None] = np.nan

    ## Conversion notes: When _E units are mol m-2 s-1, multiply by molar mass H2O to get g m-2 s-1, divide by 1000 to get kg m-2 s-1, multiply by 60*60*24 to get kg m-2 d-1, and 1 kg m-2 d-1 = 1 mm d-1. 
    ## Noting that 1 kg of water is equivalent to 1 liter (L) of water (because the density of water is 1000 kg/m³), and 1 liter of water spread over 1 square meter results in a depth of 1 mm
    diagnostics["E_mmd"] = diagnostics["E"]*18.015/1000*(60*60*24)
    
    # Turnover rates per pool
    _tr_Leaf = np.zeros(diagnostics['t'].size)
    _tr_Root = np.zeros(diagnostics['t'].size)
    _tr_Stem = np.zeros(diagnostics['t'].size)
    _tr_Seed = np.zeros(diagnostics['t'].size)
    for it, t in enumerate(diagnostics['t']):
        if np.isnan(diagnostics['idevphase_numeric'][it]):
            tr_ = Plant.PlantDev.turnover_rates[-1]
            _tr_Leaf[it] = tr_[Plant.PlantDev.ileaf]
            _tr_Root[it] = tr_[Plant.PlantDev.iroot]
            _tr_Stem[it] = tr_[Plant.PlantDev.istem]
            _tr_Seed[it] = tr_[Plant.PlantDev.iseed]
        else:
            tr_ = Plant.PlantDev.turnover_rates[diagnostics['idevphase'][it]]
            _tr_Leaf[it] = tr_[Plant.PlantDev.ileaf]
            _tr_Root[it] = tr_[Plant.PlantDev.iroot]
            _tr_Stem[it] = tr_[Plant.PlantDev.istem]
            _tr_Seed[it] = tr_[Plant.PlantDev.iseed]
    
    diagnostics['tr_Leaf'] = _tr_Leaf
    diagnostics['tr_Root'] = _tr_Root
    diagnostics['tr_Stem'] = _tr_Stem
    diagnostics['tr_Seed'] = _tr_Seed

    # Add np.nan to the end of each array in the dictionary to represent the last time point in the time_axis (corresponds to the last time point of the state vector)
    for key in diagnostics:
        if key == "t":
            diagnostics[key] = np.append(diagnostics[key], res["t"][-1])
        else:
            diagnostics[key] = np.append(diagnostics[key], np.nan)

    
    # Add state variables to the diagnostics dictionary
    diagnostics["Cleaf"] = res["y"][0,:]
    diagnostics["Cstem"] = res["y"][1,:]
    diagnostics["Croot"] = res["y"][2,:]
    diagnostics["Cseed"] = res["y"][3,:]
    diagnostics["Bio_time"] = res["y"][4,:]
    diagnostics["VRN_time"] = res["y"][5,:]
    diagnostics["Cstate"] = res["y"][7,:]
    diagnostics["Cseedbed"] = res["y"][8,:]

    # Add forcing inputs to diagnostics dictionary
    for i,f in enumerate(forcing_inputs):
        ni = i+1
        if f(time_axis[0]).size == 1:
            fstr = f"forcing {ni:02}"
            diagnostics[fstr] = f(time_axis)
        elif f(time_axis[0]).size > 1:
            # this forcing input has levels/layers (e.g. multilayer soil moisture)
            nz = f(time_axis[0]).size
            for iz in range(nz):
                fstr = f"forcing {ni:02} z{iz}"
                diagnostics[fstr] = f(time_axis)[:,iz]


    xdoy = np.floor(forcing_inputs[-2](time_axis))
    xyear = np.array(forcing_inputs[-1](time_axis), dtype=int)
    time_index = pd.to_datetime(xyear.astype(str), format='%Y') + pd.to_timedelta(xdoy - 1, unit='D')

    ## First season
    iseason = 0

    # Time indexing for model output data, to determine outputs at specific times in the growing season
    itax_sowing, itax_mature, itax_harvest, itax_phase_transitions = Plant.Site.time_index_growing_season(time_index, diagnostics['idevphase_numeric'], Plant.Management, Plant.PlantDev, iseason=iseason)

    # Total plant carbon at maturity
    # - includes the accumulation of turned over leafs and stems (excludes roots and grain)
    total_plant_carbon_at_maturity_s0 = (res["y"][0,itax_mature] + res["y"][1,itax_mature] + res["y"][3,itax_mature]) + np.nansum(diagnostics["tr_Leaf"][itax_sowing:itax_mature+1]*res["y"][0,itax_sowing:itax_mature+1] + diagnostics["tr_Stem"][itax_sowing:itax_mature+1]*res["y"][1,itax_sowing:itax_mature+1])

    # Root:Shoot ratio and LAI at selected points in the growing season
    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('vegetative'))[0][0]
    root_shoot_ratio_at_emergence_s0 = diagnostics["Croot"][itax_phase_transitions[ip]+1] / (diagnostics["Cleaf"][itax_phase_transitions[ip]+1] + diagnostics["Cstem"][itax_phase_transitions[ip]+1])

    ip0 = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('vegetative'))[0][0]
    ip1 = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('spike'))[0][0]
    it = int(itax_phase_transitions[ip0] + 0.5*(itax_phase_transitions[ip1] - itax_phase_transitions[ip0]))
    root_shoot_ratio_at_mid_veg_s0 = diagnostics["Croot"][it] / (diagnostics["Cleaf"][it] + diagnostics["Cstem"][it])
    lai_at_mid_veg_s0 = diagnostics['LAI'][it]

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('spike'))[0][0]
    root_shoot_ratio_at_spike_s0 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])
    lai_at_spike_s0 = diagnostics['LAI'][itax_phase_transitions[ip]]

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('anthesis'))[0][0]
    root_shoot_ratio_at_anthesis_s0 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('grainfill'))[0][0]
    root_shoot_ratio_at_grainfill_s0 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])
    lai_at_grainfill_s0 = diagnostics['LAI'][itax_phase_transitions[ip]]

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('maturity'))[0][0]
    root_shoot_ratio_at_maturity_s0 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])

    # Grain yield variables
    grain_yield_at_harvest_s0 = diagnostics["Cseed"][itax_harvest]/100 * (1/Plant.PlantCH2O.f_C)   ## convert gC m-2 to t dry biomass ha-1
    grain_number_harvest_s0 = diagnostics["Cseed"][itax_harvest]/Plant.PlantCH2O.f_C/Plant.W_seedTKW0    # Actual grain number at harvest


    ## Second season
    iseason = 1

    # Time indexing for model output data, to determine outputs at specific times in the growing season
    itax_sowing, itax_mature, itax_harvest, itax_phase_transitions = Plant.Site.time_index_growing_season(time_index, diagnostics['idevphase_numeric'], Plant.Management, Plant.PlantDev, iseason=iseason)

    # Total plant carbon at maturity
    # - includes the accumulation of turned over leafs and stems (excludes roots and grain)
    total_plant_carbon_at_maturity_s1 = (res["y"][0,itax_mature] + res["y"][1,itax_mature] + res["y"][3,itax_mature]) + np.nansum(diagnostics["tr_Leaf"][itax_sowing:itax_mature+1]*res["y"][0,itax_sowing:itax_mature+1] + diagnostics["tr_Stem"][itax_sowing:itax_mature+1]*res["y"][1,itax_sowing:itax_mature+1])

    # Root:Shoot ratio and LAI at selected points in the growing season
    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('vegetative'))[0][0]
    root_shoot_ratio_at_emergence_s1 = diagnostics["Croot"][itax_phase_transitions[ip]+1] / (diagnostics["Cleaf"][itax_phase_transitions[ip]+1] + diagnostics["Cstem"][itax_phase_transitions[ip]+1])

    ip0 = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('vegetative'))[0][0]
    ip1 = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('spike'))[0][0]
    it = int(itax_phase_transitions[ip0] + 0.5*(itax_phase_transitions[ip1] - itax_phase_transitions[ip0]))
    root_shoot_ratio_at_mid_veg_s1 = diagnostics["Croot"][it] / (diagnostics["Cleaf"][it] + diagnostics["Cstem"][it])
    lai_at_mid_veg_s1 = diagnostics['LAI'][it]

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('spike'))[0][0]
    root_shoot_ratio_at_spike_s1 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])
    lai_at_spike_s1 = diagnostics['LAI'][itax_phase_transitions[ip]]

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('anthesis'))[0][0]
    root_shoot_ratio_at_anthesis_s1 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('grainfill'))[0][0]
    root_shoot_ratio_at_grainfill_s1 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])
    lai_at_grainfill_s1 = diagnostics['LAI'][itax_phase_transitions[ip]]

    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('maturity'))[0][0]
    root_shoot_ratio_at_maturity_s1 = diagnostics["Croot"][itax_phase_transitions[ip]] / (diagnostics["Cleaf"][itax_phase_transitions[ip]] + diagnostics["Cstem"][itax_phase_transitions[ip]])

    # Grain yield variables
    grain_yield_at_harvest_s1 = diagnostics["Cseed"][itax_harvest]/100 * (1/Plant.PlantCH2O.f_C)   ## convert gC m-2 to t dry biomass ha-1
    grain_number_harvest_s1 = diagnostics["Cseed"][itax_harvest]/Plant.PlantCH2O.f_C/Plant.W_seedTKW0    # Actual grain number at harvest

    
    # Model output given the parameter vector p
    # - this is the model output that we compare to observations and use to calibrate the parameters
    M_p = np.array([
        total_plant_carbon_at_maturity_s0, 
        lai_at_mid_veg_s0,
        lai_at_spike_s0,
        lai_at_grainfill_s0,
        root_shoot_ratio_at_emergence_s0, 
        root_shoot_ratio_at_mid_veg_s0,
        root_shoot_ratio_at_spike_s0,
        root_shoot_ratio_at_anthesis_s0,
        root_shoot_ratio_at_grainfill_s0,
        root_shoot_ratio_at_maturity_s0,
        grain_yield_at_harvest_s0,
        grain_number_harvest_s0,
        total_plant_carbon_at_maturity_s1, 
        lai_at_mid_veg_s1,
        lai_at_spike_s1,
        lai_at_grainfill_s1,
        root_shoot_ratio_at_emergence_s1, 
        root_shoot_ratio_at_mid_veg_s1,
        root_shoot_ratio_at_spike_s1,
        root_shoot_ratio_at_anthesis_s1,
        root_shoot_ratio_at_grainfill_s1,
        root_shoot_ratio_at_maturity_s1,
        grain_yield_at_harvest_s1,
        grain_number_harvest_s1,
    ])

    return {
        "metrics": M_p,              # 1D np.array
        "diagnostics": diagnostics   # dict[str, time-series]
    }

# %% [markdown]
# ### Run the parameter sampling

# %%
Mpx = []
iparamsets_per_run = get_iparamsets_per_run(param_values, experiment.n_processes)
n_runs = len(iparamsets_per_run)
Mpx_column_headers = "nparamset,W_P_peakW,W_L_peakW,W_R_peakW,W_S_peakW,W_S_spike0,W_S_anth0,GPP_int_seas,NPP_int_seas,Rml_int_seas,Rmr_int_seas,Rg_int_seas,trflux_int_seas,FCstem2grain_int_seas,NPP2grain_int_seas,E_int_seas,LAI_peakW,W_spike_anth1,GY_mature,Sdpot_mature,GN_mature" 

def evaluate_iparamset(iparamset: int):
    nparamset = iparamset + 1
    paramset = param_values[iparamset]
    # Update and run model, with customised 'run_model_and_get_outputs' function passed in
    model_output = update_and_run_model(
        paramset,
        experiment.PlantX,
        experiment.input_data,
        parameters.df,
        parameters.problem,
        run_fn=run_model_and_get_outputs,
    )
    Mpxi = model_output['metrics']
    diagnostics = model_output['diagnostics']

    nsigfigures = len(str(np.shape(param_values)[0]))
    filename_write = f"DAESIM2_results_{experiment.xsite}_paramset{nparamset:0{nsigfigures}}.nc"
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
