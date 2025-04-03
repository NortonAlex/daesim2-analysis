from SALib.sample import fast_sampler
from SALib.analyze import fast

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## DAESIM2
from daesim.plant_1000 import PlantModuleCalculator
from daesim.climate import *
from daesim.plantgrowthphases import PlantGrowthPhases
from daesim.management import ManagementModule
from daesim.soillayers import SoilLayers
from daesim.canopylayers import CanopyLayers
from daesim.canopyradiation import CanopyRadiation
from daesim.boundarylayer import BoundaryLayerModule
from daesim.leafgasexchange import LeafGasExchangeModule
from daesim.leafgasexchange2 import LeafGasExchangeModule2
from daesim.canopygasexchange import CanopyGasExchange
from daesim.plantcarbonwater import PlantModel as PlantCH2O
from daesim.plantallocoptimal import PlantOptimalAllocation
from daesim.utils import ODEModelSolver
from daesim.utils import daesim_io_write_diag_to_nc

## DAESIM2 Analysis
from daesim2_analysis import fast_sensitivity as fastsa

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

from daesim2_analysis.args import Args
from daesim2_analysis.args import is_interactive
from daesim2_analysis.parameters import Parameters

args = Args.from_cli() if not is_interactive() else Args()

# %%
parameters = Parameters.__from_file__(args.path_parameters_file)

# # Generate samples using the FAST method
# num_samples = 300  # Number of samples to be generated
# param_values = fast_sampler.sample(problem, num_samples, seed=0)

# # param_values will contain the sampled input sets for the parameters
# param_values.shape
# print(param_values.shape)

# # %%
# def load_df_forcing(paths_df_forcing: List[str])->pd.DataFrame:
#   dfs: List[pd.DataFrame] = []
#   for path_df_forcing in paths_df_forcing:
#     df = pd.read_csv(path_df_forcing)
#     df['Date'] = pd.to_datetime(df['Date'])
#     dfs += [df]
#   df_forcing = pd.concat(dfs)
#   df_forcing = df_forcing.drop_duplicates(subset='Date')
#   df_forcing = df_forcing.sort_values(by='Date')
#   df_forcing = df_forcing.reset_index(drop=True)
#   df_forcing['DOY'] = df_forcing['Date'].dt.dayofyear
#   df_forcing['Year'] = df_forcing['Date'].dt.year
#   return df_forcing

# df_forcing = load_df_forcing(paths_df_forcing=paths_df_forcing)

# SiteX = ClimateModule(CLatDeg=-36.05,CLonDeg=146.5,timezone=10)
# start_doy_f = df_forcing["DOY"].values[0]
# start_year_f = df_forcing["Year"].values[0]
# nrundays_f = df_forcing.index.size

# ## Time discretisation
# time_nday_f, time_doy_f, time_year_f = SiteX.time_discretisation(start_doy_f, start_year_f, nrundays=nrundays_f)
# ## Adjust daily time-step to represent midday on each day
# time_doy_f = [time_doy_f[i]+0.5 for i in range(len(time_doy_f))]

# ## Time discretization for forcing data
# time_index_f = pd.to_datetime(df_forcing["Date"].values)

# # Define lists of sowing and harvest dates
# sowing_dates = [
#     pd.Timestamp(year=1971, month=5, day=11)
# ]

# harvest_dates = [
#     pd.Timestamp(year=1971, month=12, day=23)
# ]

# _soilTheta_z = np.column_stack((
#     df_forcing["Soil moisture 5 cm"].values,
#     df_forcing["Soil moisture 8 cm"].values,
#     df_forcing["Soil moisture 14 cm"].values,
#     df_forcing["Soil moisture 22 cm"].values,
#     df_forcing["Soil moisture 34 cm"].values,
#     df_forcing["Soil moisture 52 cm"].values,))

# Climate_doy_f = interp_forcing(time_nday_f, time_doy_f, kind="pconst") #, fill_value=(time_doy[0],time_doy[-1]))
# Climate_year_f = interp_forcing(time_nday_f, time_year_f, kind="pconst") #, fill_value=(time_year[0],time_year[-1]))
# Climate_airTempCMin_f = interp1d(time_nday_f, df_forcing["Minimum temperature"].values)
# Climate_airTempCMax_f = interp1d(time_nday_f, df_forcing["Maximum temperature"].values)
# Climate_airTempC_f = interp1d(time_nday_f, (df_forcing["Minimum temperature"].values+df_forcing["Maximum temperature"].values)/2)
# Climate_solRadswskyb_f = interp1d(time_nday_f, 10*(df_forcing["Global Radiation"].values-df_forcing["Diffuse Radiation"].values))
# Climate_solRadswskyd_f = interp1d(time_nday_f, 10*df_forcing["Diffuse Radiation"].values)
# Climate_airPressure_f = interp1d(time_nday_f, 100*df_forcing["Pressure"].values)
# Climate_airRH_f = interp1d(time_nday_f, df_forcing["Relative Humidity"].values)
# Climate_airU_f = interp1d(time_nday_f, df_forcing["Uavg"].values)
# Climate_airCO2_f = interp1d(time_nday_f, df_forcing["Atmospheric CO2 Concentration (bar)"].values)
# Climate_airO2_f = interp1d(time_nday_f, df_forcing["Atmospheric O2 Concentration (bar)"].values)
# Climate_soilTheta_z_f = interp1d(time_nday_f, _soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
# Climate_nday_f = interp1d(time_nday_f, time_nday_f)   ## nday represents the ordinal day-of-year plus each simulation day (e.g. a model run starting on Jan 30 and going for 2 years will have nday=30+np.arange(2*365))

# # %% [markdown]
# # ### Step 4. Initialise the Model

# # %%
# ## Check that the first sowing date and last harvest date are within the available forcing data period
# # Apply validation to sowing and harvest dates
# SiteX.validate_event_dates(sowing_dates, time_index_f, event_name="Sowing")
# SiteX.validate_event_dates(harvest_dates, time_index_f, event_name="Harvest")

# # Find steps for all sowing and harvest dates in the forcing data
# sowing_steps_f = SiteX.find_event_steps(sowing_dates, time_index_f)
# harvest_steps_f = SiteX.find_event_steps(harvest_dates, time_index_f)

# # 3. Generate time axis for ODE Solver, which is set to start at the first sowing event and end just after the last harvest event
# time_axis = time_nday_f[sowing_steps_f[0]:harvest_steps_f[-1]+2]
# # Corresponding date time index
# time_index = time_index_f[sowing_steps_f[0]:harvest_steps_f[-1]+2]

# ## Determine the ODE solver reset days
# sowing_steps_itax = SiteX.find_event_steps(sowing_dates, time_index)
# harvest_steps_itax = SiteX.find_event_steps(harvest_dates, time_index)
# reset_days_itax = SiteX.find_event_steps(sowing_dates+harvest_dates, time_index)
# reset_days_itax.sort()
# reset_days_tax = list(time_axis[reset_days_itax])

# ## Determine the sowing and harvest days and years for the DAESIM2 Management module
# sowingDays = np.floor(Climate_doy_f(time_axis[sowing_steps_itax]))
# sowingYears = np.floor(Climate_year_f(time_axis[sowing_steps_itax]))
# harvestDays = np.floor(Climate_doy_f(time_axis[harvest_steps_itax]))
# harvestYears = np.floor(Climate_year_f(time_axis[harvest_steps_itax]))

# # %%
# # time_axis = np.arange(135, 391, 1)   ## Note: time_axis represents the simulation day (_nday) and must be the same x-axis upon which the forcing data was interpolated on
# # sowing_date = 135
# # harvest_date = None

# ManagementX = ManagementModule(cropType="Wheat",sowingDays=sowingDays,harvestDays=harvestDays,sowingYears=sowingYears,harvestYears=harvestYears)

# PlantDevX = PlantGrowthPhases(
#     phases=["germination", "vegetative", "spike", "anthesis", "grainfill", "maturity"],
#     # gdd_requirements=[120,800,250,200,300,200],
#     gdd_requirements=[50,800,280,150,300,300],
#     # vd_requirements=[0, 40, 0, 0, 0, 0],
#     vd_requirements=[0, 30, 0, 0, 0, 0],
#     allocation_coeffs = [
#         [0.2, 0.1, 0.7, 0.0, 0.0],
#         [0.5, 0.1, 0.4, 0.0, 0.0],
#         [0.30, 0.4, 0.30, 0.0, 0.0],
#         [0.30, 0.4, 0.30, 0.0, 0.0],
#         [0.1, 0.1, 0.1, 0.7, 0.0],
#         [0.1, 0.1, 0.1, 0.7, 0.0]
#     ],
#     turnover_rates = [[0.001,  0.001, 0.001, 0.0, 0.0],
#                       [0.01, 0.002, 0.008, 0.0, 0.0],
#                       [0.01, 0.002, 0.008, 0.0, 0.0],
#                       [0.01, 0.002, 0.008, 0.0, 0.0],
#                       [0.033, 0.016, 0.033, 0.0002, 0.0],
#                       [0.10, 0.033, 0.10, 0.0002, 0.0]])

# BoundLayerX = BoundaryLayerModule(Site=SiteX)
# LeafX = LeafGasExchangeModule2(Site=SiteX,Jmax_opt_rVcmax=0.89,Jmax_opt_rVcmax_method="log")
# CanopyX = CanopyLayers(nlevmlcan=3)
# CanopyRadX = CanopyRadiation(Canopy=CanopyX)
# CanopyGasExchangeX = CanopyGasExchange(Leaf=LeafX,Canopy=CanopyX,CanopyRad=CanopyRadX)
# # SoilLayersX = SoilLayers(nlevmlsoil=2,z_max=2.0)
# SoilLayersX = SoilLayers(nlevmlsoil=6,z_max=0.66,z_top=0.10,discretise_method="horizon",
#                          z_horizon=[0.06, 0.06, 0.06, 0.10, 0.10, 0.28],
#                         Psi_e=[-1.38E-03, -1.38E-03, -1.38E-03, -1.32E-03, -2.58E-03, -0.960E-03],
#                         b_soil = [4.74, 4.74, 4.74, 6.77, 8.17, 10.73],
#                          K_sat = [29.7, 29.7, 29.7, 25.2, 13.9, 40.9],
#                          soilThetaMax = [0.12, 0.12, 0.12, 0.20, 0.3, 0.4])
# # PlantCH2OX = PlantCH2O(Site=SiteX,SoilLayers=SoilLayersX,CanopyGasExchange=CanopyGasExchangeX,BoundaryLayer=BoundLayerX,maxLAI=6.0,ksr_coeff=1000,SLA=0.030)
# PlantCH2OX = PlantCH2O(Site=SiteX,SoilLayers=SoilLayersX,CanopyGasExchange=CanopyGasExchangeX,BoundaryLayer=BoundLayerX,maxLAI=6.5,ksr_coeff=1500,SLA=0.02,sf=1.0,Psi_f=-5.0)
# PlantAllocX = PlantOptimalAllocation(Plant=PlantCH2OX)
# PlantX = PlantModuleCalculator(
#     Site=SiteX,
#     Management=ManagementX,
#     PlantDev=PlantDevX,
#     PlantCH2O=PlantCH2OX,
#     PlantAlloc=PlantAllocX,
#     GDD_method="nonlinear",
#     GDD_Tbase=0.0,
#     GDD_Topt=22.5,
#     GDD_Tupp=35.0,
#     hc_max_GDDindex=sum(PlantDevX.gdd_requirements[0:2])/PlantDevX.totalgdd,
#     d_r_max=2.0,
#     Vmaxremob=3.0,
#     Kmremob=0.5,
#     remob_phase=["grainfill","maturity"],
#     specified_phase="spike",
#     grainfill_phase=["grainfill","maturity"],
# )

# # %%
# ## Define the callable calculator that defines the right-hand-side ODE function
# PlantXCalc = PlantX.calculate

# Model = ODEModelSolver(calculator=PlantXCalc, states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_start=time_axis[0], log_diagnostics=True)

# forcing_inputs = [Climate_solRadswskyb_f,
#                   Climate_solRadswskyd_f,
#                   Climate_airTempCMin_f,
#                   Climate_airTempCMax_f,
#                   Climate_airPressure_f,
#                   Climate_airRH_f,
#                   Climate_airCO2_f,
#                   Climate_airO2_f,
#                   Climate_airU_f,
#                   Climate_soilTheta_z_f,
#                   Climate_doy_f,
#                   Climate_year_f]

# reset_days = reset_days_tax

# zero_crossing_indices = [4,5,6]

# # %%

# write_to_nc = True

# iparamsets_per_run = [list(range(i, i+n_processes)) for i in range(0, len(param_values), n_processes)]
# input_data = [ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices]
# Mpx_column_headers = "nparamset,W_P_peakW,W_L_peakW,W_R_peakW,W_S_peakW,W_S_spike0,W_S_anth0,GPP_int_seas,NPP_int_seas,Rml_int_seas,Rmr_int_seas,Rg_int_seas,trflux_int_seas,FCstem2grain_int_seas,NPP2grain_int_seas,E_int_seas,LAI_peakW,W_spike_anth1,GY_mature,Sdpot_mature,GN_mature"   # to save as header in csv file. N.B. must match the custom output in the update_and_run_model function
# Mpx = []
# n_runs = len(iparamsets_per_run)

# def evaluate_paramset(iparamset: int):
#   nparamset = iparamset + 1
#   paramset = param_values[iparamset]
#   # model_output = fastsa.update_and_run_model(paramset, PlantX, input_data, parameters_df, problem)
#   model_output = fastsa.update_and_run_model(paramset, PlantX, input_data, parameters_df, problem)
  
#   Mpxi, diagnostics = model_output[0], model_output[1]
#   if write_to_nc:
#     nsigfigures = len(str(np.shape(param_values)[0]))
#     filename_write = f"FAST_results_{xsite}_paramset{nparamset:0{nsigfigures}}.nc"
#     daesim_io_write_diag_to_nc(
#       PlantX,
#       diagnostics,
#       dir_xsite_parameters,
#       filename_write,
#       time_index,
#       problem=problem,
#       param_values=paramset,
#       nc_attributes={'title': title, 'description': description}
#     )
#   return np.insert(Mpxi, 0, nparamset)

# # %%

# for n_run in range(n_runs):
#    with Pool(processes=n_processes) as pool:
#       print(iparamsets_per_run[n_run])
#       results = pool.map(evaluate_paramset, iparamsets_per_run[n_run])
#       Mpx += results
  
# np.save(path_Mpx, np.array(Mpx))