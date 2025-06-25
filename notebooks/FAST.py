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

from daesim2_analysis.args import Args
from daesim2_analysis.args import is_interactive
from daesim2_analysis.parameters import Parameters
from daesim2_analysis.utils import *
from daesim2_analysis.forcing_data import ForcingData

args = Args.from_cli() if not is_interactive() else Args()

# %%
parameters: Parameters = Parameters.__from_file__(args.path_parameters_file)
param_values = parameters.sample(args.n_samples)

# %%
df_forcing = load_df_forcing(args.paths_df_forcing)
SiteX = ClimateModule(CLatDeg=args.CLatDeg,CLonDeg=args.CLonDeg,timezone=args.tz)
forcing_data = ForcingData(
    SiteX=SiteX,
    sowing_dates=args.sowing_dates,
    harvest_dates=args.harvest_dates,
    df=df_forcing
)
ManagementX = args.management.create(
    cropType=args.crop_type,
    sowingDays=forcing_data.sowing_days,
    harvestDays=forcing_data.harvest_dates,
    sowingYears=forcing_data.sowing_years,
    harvestYears=forcing_data.harvest_years,
)

PlantDevX = args.plant_growth_phases.create()
BoundaryLayerX = args.boundary_layer.create(Site=SiteX)
LeafX = args.leaf_exchange.create(Site=SiteX)
CanopyX = args.canopy.create(nlevmlcan=3)
CanopyRadX = args.canopy_rad.create(Canopy=CanopyX)
CanopyGasExchangeX = args.canopy_gas_exchange.create(Leaf=LeafX, Canopy=CanopyX, CanopyRad=CanopyRadX)
SoilLayersX = args.soil_layers.create(nlevmlsoil=6)
PlantCH2OX = args.plant_ch2o.create(
    Site=SiteX,
    SoilLayers=SoilLayersX,
    CanopyGasExchange=CanopyGasExchangeX,
    BoundaryLayer=BoundaryLayerX
)
PlantAllocX = args.plant_optimal_allocation.create(Plant=PlantCH2OX)
PlantX = args.plant_module_calculator.create(
    Site=SiteX,
    Management=ManagementX,
    PlantDev=PlantDevX,
    PlantCH2O = PlantCH2OX,
    PlantAlloc=PlantAllocX
)
PlantXCalc = PlantX.calculate
Model = ODEModelSolver(
    calculator=PlantXCalc,
    states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    time_start=forcing_data.time_axis[0],
    log_diagnostics=True
)

input_data = [
    ODEModelSolver,
    forcing_data.time_axis,
    forcing_data.time_index,
    forcing_data.inputs,
    forcing_data.reset_days,
    forcing_data.zero_crossing_indices
]
Mpx = []
iparamsets_per_run = get_iparamsets_per_run(param_values, args.n_processes)
n_runs = len(iparamsets_per_run)
Mpx_column_headers = "nparamset,W_P_peakW,W_L_peakW,W_R_peakW,W_S_peakW,W_S_spike0,W_S_anth0,GPP_int_seas,NPP_int_seas,Rml_int_seas,Rmr_int_seas,Rg_int_seas,trflux_int_seas,FCstem2grain_int_seas,NPP2grain_int_seas,E_int_seas,LAI_peakW,W_spike_anth1,GY_mature,Sdpot_mature,GN_mature" 

def evaluate_iparamset(iparamset: int):
    nparamset = iparamset + 1
    paramset = param_values[iparamset]
    # model_output = fastsa.update_and_run_model(paramset, PlantX, input_data, parameters_df, problem)
    model_output = fastsa.update_and_run_model(
        paramset,
        PlantX,
        input_data,
        parameters.df,
        parameters.problem
    )
    Mpxi, diagnostics = model_output[0], model_output[1]

    nsigfigures = len(str(np.shape(param_values)[0]))
    filename_write = f"FAST_results_{args.xsite}_paramset{nparamset:0{nsigfigures}}.nc"
    daesim_io_write_diag_to_nc(
      PlantX,
      diagnostics,
      args.dir_xsite_parameters + '/',
      filename_write,
      forcing_data.time_index,
      problem=parameters.problem,
      param_values=paramset,
      nc_attributes={'title': args.title, 'description': args.description}
    )
    return np.insert(Mpxi, 0, nparamset)

if __name__ == '__main__':
    print(args.dir_xsite_FAST_results)
    print(args.dir_xsite_parameters)
    for n_run in range(n_runs):
        with Pool(processes=args.n_processes) as pool:
            Mpx+= = pool.map(evaluate_iparamset, iparamsets_per_run[n_run])
            print(result)
            break
    



