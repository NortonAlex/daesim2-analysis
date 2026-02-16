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
from pandas import Timestamp

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
    sowing_dates=[Timestamp(year=2018,month=5,day=15)],
    harvest_dates=[Timestamp(year=2018,month=12,day=31)],
    xsite="Milgadara_2018",
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
## Create datetime array from doy and year outputs
xyear = np.array(model_output[d_fd_mapping['Climate_year_f']], dtype=int)
xdoy = model_output[d_fd_mapping['Climate_doy_f']]
model_output['Date'] = pd.to_datetime(xyear.astype(str), format='%Y') + pd.to_timedelta(xdoy - 1, unit='D')

# %%
fig, axes = plt.subplots(4,1,figsize=(8,8),sharex=True)

axes[0].plot(model_output['Date'], model_output[d_fd_mapping['Climate_solRadswskyb_f']] + model_output[d_fd_mapping['Climate_solRadswskyd_f']], c='0.4', label="Global")
axes[0].plot(model_output['Date'], model_output[d_fd_mapping['Climate_solRadswskyb_f']], c='goldenrod', alpha=0.5, label="Direct")
axes[0].plot(model_output['Date'], model_output[d_fd_mapping['Climate_solRadswskyd_f']], c='C0', alpha=0.5, label="Diffuse")
axes[0].set_ylabel("Solar radiation\n"+r"($\rm W \; m^{-2}$)")
axes[0].legend(loc=1,handlelength=0.75)
axes[0].annotate("Solar radiation", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
axes[0].set_ylim([0,400])

axes[1].plot(model_output['Date'], model_output[d_fd_mapping['Climate_airTempCMin_f']], c="lightsteelblue", label="Min")
axes[1].plot(model_output['Date'], model_output[d_fd_mapping['Climate_airTempCMax_f']], c="indianred", label="Max")
leafTempC = experiment.PlantX.Site.compute_skin_temp((model_output[d_fd_mapping['Climate_airTempCMin_f']] + model_output[d_fd_mapping['Climate_airTempCMax_f']])/2, model_output[d_fd_mapping['Climate_solRadswskyb_f']] + model_output[d_fd_mapping['Climate_solRadswskyd_f']])
axes[1].plot(model_output['Date'], leafTempC, c="darkgreen", label="Leaf")
axes[1].plot(model_output['Date'], (model_output[d_fd_mapping['Climate_airTempCMin_f']] + model_output[d_fd_mapping['Climate_airTempCMax_f']])/2, c="0.5", label="Air")
axes[1].set_ylabel("Air Temperature\n"+r"($\rm ^{\circ}C$)")
axes[1].legend(loc=1,handlelength=0.75)
axes[1].annotate("Air temperature", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
axes[1].set_ylim([-5,45])

axes[2].plot(model_output['Date'], model_output[d_fd_mapping['Climate_airRH_f']], c="0.4")
axes[2].set_ylabel("Relative humidity\n"+r"(%)")
axes[2].annotate("Relative humidity", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)

xcolors = np.linspace(0.9,0.1,experiment.PlantX.PlantCH2O.SoilLayers.nlevmlsoil).astype(str)
for iz in range(experiment.PlantX.PlantCH2O.SoilLayers.nlevmlsoil):
    zlayer_forcing_label = d_fd_mapping['Climate_soilTheta_z_f'] + f' z{iz}'
    axes[3].plot(model_output['Date'], 100*model_output[zlayer_forcing_label], c=xcolors[iz])
axes[3].set_ylabel("Soil moisture\n"+r"(%)")
# axes[3].tick_params(axis='x', labelrotation=45)
axes[3].annotate("Soil moisture", (0.02,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)

if "Precipitation" in experiment.df_forcing.columns:
    ax2 = axes[3].twinx()
    i0, i1 = experiment.ForcingDataX.time_axis[0]-1, experiment.ForcingDataX.time_axis[-1]
    ax2.bar(model_output['Date'], experiment.df_forcing["Precipitation"].values[i0:i1], color="0.4")
    ax2.set_ylabel("Daily Precipitation\n(mm)")
    axes[3].annotate("Precipitation", (0.98,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='right', fontsize=12)

# axes[0].set_xlim([experiment.PlantX.Management.sowingDays[0],model_output[d_fd_mapping['Climate_doy_f']][-1]])

plt.tight_layout()

# %%
fig, axes = plt.subplots(5,1,figsize=(8,10),sharex=True)

axes[0].plot(model_output['Date'], model_output["LAI"])
axes[0].set_ylabel("LAI\n"+r"($\rm m^2 \; m^{-2}$)")
axes[0].tick_params(axis='x', labelrotation=45)
axes[0].annotate("Leaf area index", (0.01,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
# axes[0].set_ylim([0,6.5])

axes[1].plot(model_output['Date'], model_output["GPP"])
axes[1].set_ylabel("GPP\n"+r"($\rm g C \; m^{-2} \; d^{-1}$)")
axes[1].tick_params(axis='x', labelrotation=45)
axes[1].annotate("Photosynthesis", (0.01,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
# axes[1].set_ylim([0,30])

axes[2].plot(model_output['Date'], model_output["E_mmd"])
axes[2].set_ylabel(r"$\rm E$"+"\n"+r"($\rm mm \; d^{-1}$)")
axes[2].tick_params(axis='x', labelrotation=45)
axes[2].annotate("Transpiration Rate", (0.01,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)
# axes[2].set_ylim([0,6])

axes[3].plot(model_output['Date'], model_output["Bio_time"])
axes[3].set_ylabel("Thermal Time\n"+r"($\rm ^{\circ}$C d)")
axes[3].annotate("Growing Degree Days", (0.01,0.93), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=12)


alp = 0.6
axes[4].plot(model_output['Date'], model_output["Cleaf"]+model_output["Croot"]+model_output["Cstem"]+model_output["Cseed"],c='k',label="Plant", alpha=alp)
axes[4].plot(model_output['Date'], model_output["Cleaf"],label="Leaf", alpha=alp)
axes[4].plot(model_output['Date'], model_output["Cstem"],label="Stem", alpha=alp)
axes[4].plot(model_output['Date'], model_output["Croot"],label="Root", alpha=alp)
axes[4].plot(model_output['Date'], model_output["Cseed"],label="Seed", alpha=alp)
axes[4].set_ylabel("Carbon Pool Size\n"+r"(g C $\rm m^{-2}$)")
axes[4].legend(fontsize=9,handlelength=0.8)


xlimmin, xlimmax = experiment.ForcingDataX.sowing_dates[0], experiment.ForcingDataX.harvest_dates[-1]

for iseason in range(len(experiment.ForcingDataX.sowing_dates)):
    # Time indexing for model output data, to determine outputs at specific times in the growing season
    itax_sowing, itax_mature, itax_harvest, itax_phase_transitions = experiment.PlantX.Site.time_index_growing_season(experiment.ForcingDataX.time_index, model_output['idevphase_numeric'], experiment.PlantX.Management, experiment.PlantX.PlantDev, iseason=iseason)
    harvest_index_maturity = model_output["Cseed"][itax_harvest] / (model_output["Cleaf"][itax_mature]+model_output["Croot"][itax_mature]+model_output["Cstem"][itax_mature])
    yield_from_seed_Cpool = model_output["Cseed"][itax_harvest]/100 * (1/experiment.PlantX.PlantCH2O.f_C)   ## convert gC m-2 to t dry biomass ha-1
    xyear = int(model_output[d_fd_mapping['Climate_year_f']][itax_harvest])
    axes[4].annotate(f"{xyear} Yield = %1.2f t/ha" % (yield_from_seed_Cpool), (0.01,0.96-0.1*iseason), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=9)
    # axes[4].annotate(f"{xyear} Harvest index = %1.2f" % (harvest_index_maturity), (0.01,0.5-0.1*iseason), xycoords='axes fraction', verticalalignment='top', horizontalalignment='left', fontsize=9)

    # Add annotations for developmental phases
    itax_sowing, itax_mature, itax_harvest, itax_phase_transitions = experiment.PlantX.Site.time_index_growing_season(experiment.ForcingDataX.time_index, model_output['idevphase_numeric'], experiment.PlantX.Management, experiment.PlantX.PlantDev, iseason=iseason)
    # for ax in axes:
    ax = axes[3]
    ylimmin, ylimmax = 0, ax.get_ylim()[1]
    for itime in itax_phase_transitions:
        ax.vlines(x=model_output['Date'][itime], ymin=ylimmin, ymax=ylimmax, color='0.5',linestyle="--")
        text_x = model_output['Date'][itime] + pd.Timedelta(days=1)
        text_y = 0.5 * ylimmax
        if (text_x < xlimmin) or (text_x > xlimmax):
            continue
        elif ~np.isnan(model_output['idevphase_numeric'][itime]):
            idevphase = int(model_output['idevphase'][itime])
            phase = experiment.PlantX.PlantDev.phases[idevphase]
            ax.text(text_x, text_y, phase, horizontalalignment='left', verticalalignment='center',
                    fontsize=8, alpha=0.7, rotation=90)
        elif np.isnan(model_output['idevphase_numeric'][itime]):
            phase = "mature"
            ax.text(text_x, text_y, phase, horizontalalignment='left', verticalalignment='center',
                    fontsize=8, alpha=0.7, rotation=90)
    ax.set_ylim([ylimmin, ylimmax])

plt.tight_layout()

# %%

# %%

# %%

# %%
