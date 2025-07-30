# %%
"""
Forward run of the DAESIM2 model using the daesim2_analysis helper functions
"""

import numpy as np
## DAESIM2
from daesim.climate import *
from daesim.utils import ODEModelSolver
from daesim.utils import daesim_io_write_diag_to_nc
## DAESIM2 Analysis
# from daesim2_analysis import fast_sensitivity as fastsa
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
# from tqdm import tqdm
from time import time
from functools import partial

from daesim2_analysis.experiment import Experiment
from daesim2_analysis.experiment import is_interactive
from daesim2_analysis.parameters import Parameters
from daesim2_analysis.utils import *
from daesim2_analysis.forcing_data import ForcingData
from daesim2_analysis.run import *

# %%
# experiment = Experiment.from_cli() if not is_interactive() else Experiment()
experiment = Experiment(
    df_forcing="../DAESIM_data/DAESim_forcing_data/DAESim_forcing_Milgadara_2018.csv",
    parameters="../parameters/Fast1.json",
    daesim_config="../daesim_configs/DAESIM1.json",
    dir_results="../results/",
)

# %%
parameters: Parameters = experiment.parameters
param_values = parameters.sample(experiment.n_samples)

# %%
# Call the function that updates parameters, runs the model and returns selected outputs
model_output = update_and_run_model(
    parameters.init, 
    experiment.PlantX,
    experiment.input_data,
    parameters.df,
    parameters.problem)

# %% [markdown]
# ## Save Model Output to File

# %%
filename_write = f"DAESIM2-Plant_{experiment.xsite}_{experiment.crop_type}.nc"
daesim_io_write_diag_to_nc(
    experiment.PlantX,
    model_output,
    experiment.dir_results,
    filename_write,
    experiment.ForcingDataX.time_index,
    nc_attributes={'title': experiment.title, 'description': experiment.description}
    )

# %% [markdown]
# ## Plot Results

# %%
import matplotlib.pyplot as plt

# %%
## Forcing data mapping of input variable names to output variable names
d_fd_mapping = {
    'Climate_solRadswskyb_f': 'forcing 01',
    'Climate_solRadswskyd_f': 'forcing 02',
    'Climate_airTempCMin_f': 'forcing 03',
    'Climate_airTempCMax_f': 'forcing 04',
    'Climate_airPressure_f': 'forcing 05',
    'Climate_airRH_f': 'forcing 06',
    'Climate_airCO2_f': 'forcing 07',
    'Climate_airO2_f': 'forcing 08',
    'Climate_airU_f': 'forcing 09',
    'Climate_soilTheta_z_f': 'forcing 10',
    'Climate_doy_f': 'forcing 11',
    'Climate_year_f': 'forcing 12'}

# %%
fig, axes = plt.subplots(4,1,figsize=(8,8),sharex=True)

axes[0].plot(model_output[d_fd_mapping['Climate_doy_f']], model_output[d_fd_mapping['Climate_solRadswskyb_f']] + model_output[d_fd_mapping['Climate_solRadswskyd_f']], c='0.4', label="Global")
axes[0].plot(model_output[d_fd_mapping['Climate_doy_f']], model_output[d_fd_mapping['Climate_solRadswskyb_f']], c='goldenrod', alpha=0.5, label="Direct")
axes[0].plot(model_output[d_fd_mapping['Climate_doy_f']], model_output[d_fd_mapping['Climate_solRadswskyd_f']], c='C0', alpha=0.5, label="Diffuse")
axes[0].set_ylabel("Solar radiation\n"+r"($\rm W \; m^{-2}$)")
axes[0].legend(loc=1,handlelength=0.75)
# axes[0].tick_params(axis='x', labelrotation=45)
axes[0].annotate("Solar radiation", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
axes[0].set_ylim([0,400])

axes[1].plot(model_output[d_fd_mapping['Climate_doy_f']], model_output[d_fd_mapping['Climate_airTempCMin_f']], c="lightsteelblue", label="Min")
axes[1].plot(model_output[d_fd_mapping['Climate_doy_f']], model_output[d_fd_mapping['Climate_airTempCMax_f']], c="indianred", label="Max")
leafTempC = experiment.PlantX.Site.compute_skin_temp((model_output[d_fd_mapping['Climate_airTempCMin_f']] + model_output[d_fd_mapping['Climate_airTempCMax_f']])/2, model_output[d_fd_mapping['Climate_solRadswskyb_f']] + model_output[d_fd_mapping['Climate_solRadswskyd_f']])
axes[1].plot(model_output[d_fd_mapping['Climate_doy_f']], leafTempC, c="darkgreen", label="Leaf")
axes[1].plot(model_output[d_fd_mapping['Climate_doy_f']], (model_output[d_fd_mapping['Climate_airTempCMin_f']] + model_output[d_fd_mapping['Climate_airTempCMax_f']])/2, c="0.5", label="Air")
axes[1].set_ylabel("Air Temperature\n"+r"($\rm ^{\circ}C$)")
# axes[1].tick_params(axis='x', labelrotation=45)
axes[1].legend(loc=1,handlelength=0.75)
axes[1].annotate("Air temperature", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
axes[1].set_ylim([-5,45])

axes[2].plot(model_output[d_fd_mapping['Climate_doy_f']], model_output[d_fd_mapping['Climate_airRH_f']], c="0.4")
axes[2].set_ylabel("Relative humidity\n"+r"(%)")
# axes[2].tick_params(axis='x', labelrotation=45)
axes[2].annotate("Relative humidity", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)

xcolors = np.linspace(0.9,0.1,experiment.PlantX.PlantCH2O.SoilLayers.nlevmlsoil).astype(str)
for iz in range(experiment.PlantX.PlantCH2O.SoilLayers.nlevmlsoil):
    zlayer_forcing_label = d_fd_mapping['Climate_soilTheta_z_f'] + f' z{iz}'
    axes[3].plot(model_output[d_fd_mapping['Climate_doy_f']], 100*model_output[zlayer_forcing_label], c=xcolors[iz])
axes[3].set_ylabel("Soil moisture\n"+r"(%)")
# axes[3].tick_params(axis='x', labelrotation=45)
axes[3].annotate("Soil moisture", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
# axes[3].set_ylim([20,50])
axes[3].set_xlabel("Time (day of year)")

ax2 = axes[3].twinx()
# i0, i1 = time_axis[0]-1, time_axis[-1]
ax2.bar(model_output[d_fd_mapping['Climate_doy_f']], experiment.df_forcing["Precipitation"].values, color="0.4")
ax2.set_ylabel("Daily Precipitation\n(mm)")
axes[3].annotate("Precipitation", (0.98,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='right', fontsize=12)

axes[0].set_xlim([experiment.PlantX.Management.sowingDays[0],model_output[d_fd_mapping['Climate_doy_f']][-1]])

plt.tight_layout()
# plt.savefig("/Users/alexandernorton/ANU/Projects/DAESim/DAESIM/results/DAESIM2_%s_climate.png" % site_filename,dpi=300,bbox_inches='tight')

# %%

# %%
