# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Sobol' Sensitivity Analysis
#
# **Description**: The Sobol’ method decomposes the variance of the model output into contributions from individual input parameters and their interactions. It provides first-order indices (measuring the contribution of a single parameter) and total-order indices (which account for both individual effects and interactions with other parameters). The underlying assumption is that the total variance in model output can be decomposed (partitioned) into components of the variability in inputs. 
#
# **Sampling**: Typically uses Monte Carlo sampling, though quasi-random sequences like Sobol sequences (low-discrepancy sequences) are often employed to improve efficiency. Other efficient sampling approaches can be employed such as a so-called "radial design" (Campolongo et al., 2011, doi: 10.1016/j.cpc.2010.12.039). 
#
# **Advantages**: Provides a comprehensive picture of both individual parameter effects and higher-order interactions.
#
# **Disadvantages**: Requires many model runs, particularly for models with many parameters. In addition, as noted by Peteers et al. (2014, doi:10.5194/hess-18-3777-2014) "the main drawback of variance-based methods is that it assumes that the entire effect of a parameter can be summarized by the variance (Borgonovo, 2007; Borgonovo et al., 2011). Variance-based sensitivity indices will therefore be less reliable if the response to a parameter has a skewed or multi-modal distribution."
#
# **Use case**: One of the most widely used methods in GSA for models with complex interactions between parameters.
#
# ### Step 1.  Install Required Libraries
#
# The main Python package for performing Sobol' sensitivity analysis is SALib (Sensitivity Analysis Library). If you don’t have it installed, you can install it via pip or conda. In addition, we use Sobol' low discrepancy sequences for a more efficiency sampling of the parameter space, which requires the Scipy package qmc. 

# %%
import scipy.stats.qmc as qmc
from SALib.analyze import sobol

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
from daesim.utils import ODEModelSolver, daesim_io_write_diag_to_nc

## DAESIM2 Analysis
from daesim2_analysis import fast_sensitivity as fastsa

# %% [markdown]
# ### Step 2. Define the model problem
#
# Need to define the number of parameters, their names, and the ranges (bounds) over which the sensitivity analysis will be performed. For example, if your model has three parameters (e.g., param1, param2, param3), you would define the problem in the following format:
#
# This includes both model parameters (e.g. here this includes Vcmax_opt, Jmax_opt, g1) and the forcing data ranges (e.g. here this includes T, Q, fgsw).

# %%
## Parameters
parameter_modulepath = ["PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""]
parameter_module = ["Leaf", "Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""]
parameter_names  = ["Vcmax_opt", "g1", "SLA", "maxLAI", "ksr_coeff", "Psi_f", "sf", "gdd_requirements", "gdd_requirements", "GY_FE", "GY_SDW_50", "CI", "d_r_max"]
parameter_units  = ["mol CO2 m-2 s-1", "kPa^0.5", "m2 g d.wt-1", "m2 m-2", "g d.wt-1 m-1", "MPa", "MPa-1", "deg C d", "deg C d", "thsnd grains g d.wt spike-1", "g d.wt m-2", "-", "m"]
parameter_init   = [60e-6, 3, 0.03, 6, 1000, -3.5, 3.5, 900, 650, 0.1, 100, 0.75, 0.5]
parameter_min    = [30e-6, 1, 0.015, 5, 300, -8.0, 1.5, 600, 350, 0.08, 80, 0.5, 0.15]
parameter_max    = [120e-6, 6, 0.035, 7, 5000, -1.0, 7.0, 1800, 700, 0.21, 150, 1.0, 0.66]
parameter_phase_specific = [False, False, False, False, False, False, False, True, True, False, False, False, False]
parameter_phase = [None, None, None, None, None, None, None, "vegetative", "grainfill", None, None, None, None]

# %%
# Check if all parameter vectors have the same length
lengths = [len(parameter_modulepath), len(parameter_module), len(parameter_names), len(parameter_units), len(parameter_init), len(parameter_min), len(parameter_max)]
# Print result of the length check
if all(length == lengths[0] for length in lengths):
    print("All parameter vectors are of the same length.")
else:
    print("The parameter vectors are not of the same length. Lengths found:", lengths)
    
# Create a dataframe to combine the parameter information into one data structure
parameters_df = pd.DataFrame({
    "Module Path": parameter_modulepath,
    "Module": parameter_module,
    "Phase Specific": parameter_phase_specific,
    "Phase": parameter_phase,
    "Name": parameter_names,
    "Units": parameter_units,
    "Initial Value": parameter_init,
    "Min": parameter_min,
    "Max": parameter_max
})

# Dictionary required as input to SALib FAST sampler
problem = {
    "num_vars": len(parameters_df),    # Number of input parameters
    "names": parameters_df["Name"].values,   # Parameter names
    "bounds": [[row['Min'], row['Max']] for _, row in parameters_df.iterrows()],    # Parameter ranges
}

# %% [markdown]
# ### Step 3. Generate the samples using Sobol' low discrepancy sequences
#
# As with any sensitivity analysis method, you need to generate the sampling points for the input parameters. Here, we use the Scipy qmc function which implements Sobol' low discrepancy sequences. Here, you specify the total number of samples you want which must be a power of 2 for the Sobol' method (it is important to choose a sufficiently large number to ensure the analysis is robust). 
#
# N.B. For a typical setup:
#
# - num_samples refers to the number of frequencies used to sample each parameter.
# - The actual number of parameter sets generated depends on num_vars and the sampling process.
#
# The total number of samples is typically num_samples * num_vars, and thus, Y should have this size. You can adjust the setup to reflect this. 

# %%
# Generate samples using the Sobol' low discrepancy sequences method
num_samples = 1024  # Number of samples to be generated

# Create a reproducible NumPy random number generator
# rng = np.random.default_rng(42)

# Generate parameter samples
# sobol_sampler = qmc.Sobol(d=problem["num_vars"], scramble=True, rng=rng)  # Scrambled Sobol' with fixed rng
sobol_sampler = qmc.Sobol(d=problem["num_vars"], scramble=True, seed=0)  # Scrambled Sobol' with fixed seed
sample = sobol_sampler.random(num_samples)

# Scale sample to parameter bounds
lower_bounds = np.array([b[0] for b in problem['bounds']])
upper_bounds = np.array([b[1] for b in problem['bounds']])
param_values = lower_bounds + sample * (upper_bounds - lower_bounds)  # Scale Sobol' samples

# param_values will contain the sampled input sets for the parameters
param_values.shape

# %% [markdown]
# ### Step 4. Initialise the study site and import forcing data

# %% [markdown]
# #### a. Milgadara site 2018-2019

# %%
# Import two years of forcing files and combine then
file = "/Users/alexandernorton/ANU/Projects/DAESIM/daesim/data/DAESim_forcing_Milgadara_2018.csv"
df_forcing1 = pd.read_csv(file)

file = "/Users/alexandernorton/ANU/Projects/DAESIM/daesim/data/DAESim_forcing_Milgadara_2019.csv"
df_forcing2 = pd.read_csv(file)

df_forcing1['Date'] = pd.to_datetime(df_forcing1['Date'])
df_forcing2['Date'] = pd.to_datetime(df_forcing2['Date'])

# Append the new dataframe and drop duplicates based on the "Date" column
df_forcing = pd.concat([df_forcing1, df_forcing2]).drop_duplicates(subset="Date").sort_values(by="Date")

# Reset the index for the combined dataframe (optional)
df_forcing.reset_index(drop=True, inplace=True)

# Add ordinal day of year (DOY) and Year variables
df_forcing["DOY"] = df_forcing["Date"].dt.dayofyear
df_forcing["Year"] = df_forcing["Date"].dt.year


## Interpolate discrete soil moisture data
df_forcing["Soil moisture interp"] = df_forcing["Soil moisture"].interpolate('quadratic')

## Assume that the forcing data (units: mm) can be equated to relative changes in volumetric soil moisture between two arbitrary minimum and maximum values
f_soilTheta_min = 0.25
f_soilTheta_max = 0.40

f_soilTheta_min_mm = df_forcing["Soil moisture interp"].min()
f_soilTheta_max_mm = df_forcing["Soil moisture interp"].max()

f_soilTheta_norm_mm = (df_forcing["Soil moisture interp"].values - f_soilTheta_min_mm)/(f_soilTheta_max_mm - f_soilTheta_min_mm)
f_soilTheta_norm = f_soilTheta_min + f_soilTheta_norm_mm * (f_soilTheta_max - f_soilTheta_min)


## Milgadara site location-34.38904277303204, 148.46949938279096
SiteX = ClimateModule(CLatDeg=-34.389,CLonDeg=148.469,timezone=10)
start_doy_f = df_forcing["DOY"].values[0]
start_year_f = df_forcing["Year"].values[0]
nrundays_f = df_forcing.index.size

## Time discretisation
time_nday_f, time_doy_f, time_year_f = SiteX.time_discretisation(start_doy_f, start_year_f, nrundays=nrundays_f)
## Adjust daily time-step to represent midday on each day
time_doy_f = [time_doy_f[i]+0.5 for i in range(len(time_doy_f))]

## Time discretization for forcing data
time_index_f = pd.to_datetime(df_forcing["Date"].values)

# Define lists of sowing and harvest dates
sowing_dates = [
    pd.Timestamp(year=2018, month=7, day=20),  # First season
    pd.Timestamp(year=2019, month=7, day=15)   # Second season
]

harvest_dates = [
    pd.Timestamp(year=2018, month=12, day=10),  # First season
    pd.Timestamp(year=2019, month=12, day=5)    # Second season
]


## Make some assumption about the fraction of diffuse radiation
diffuse_fraction = 0.2

## Shortwave radiation at surface (convert MJ m-2 d-1 to W m-2)
_Rsb_Wm2 = (1-diffuse_fraction) * df_forcing["SRAD"].values * 1e6 / (60*60*24)
_Rsd_Wm2 = diffuse_fraction * df_forcing["SRAD"].values * 1e6 / (60*60*24)

## Create synthetic data for other forcing variables
_p = 101325*np.ones(nrundays_f)
_es = SiteX.compute_sat_vapor_pressure_daily(df_forcing["Minimum temperature"].values,df_forcing["Maximum temperature"].values)
_RH = SiteX.compute_relative_humidity(df_forcing["VPeff"].values/10,_es/1000)
_RH[_RH > 100] = 100
_CO2 = 400*(_p/1e5)*1e-6     ## carbon dioxide partial pressure (bar)
_O2 = 209000*(_p/1e5)*1e-6   ## oxygen partial pressure (bar)
_soilTheta =  0.35*np.ones(nrundays_f)   ## volumetric soil moisture content (m3 m-3)
_soilTheta = f_soilTheta_norm

## Create a multi-layer soil moisture forcing dataset
## Option 1: Same soil moisture across all layers
nlevmlsoil = 2
_soilTheta_z = np.repeat(_soilTheta[:, np.newaxis], nlevmlsoil, axis=1)
## Option 2: Adjust soil moisture in each layer
_soilTheta_z0 = _soilTheta-0.06
_soilTheta_z1 = _soilTheta+0.02
_soilTheta_z = np.column_stack((_soilTheta_z0, _soilTheta_z1))


Climate_doy_f = interp_forcing(time_nday_f, time_doy_f, kind="pconst", fill_value=(time_doy_f[0],time_doy_f[-1]))
Climate_year_f = interp_forcing(time_nday_f, time_year_f, kind="pconst", fill_value=(time_year_f[0],time_year_f[-1]))
Climate_airTempCMin_f = interp1d(time_nday_f, df_forcing["Minimum temperature"].values)
Climate_airTempCMax_f = interp1d(time_nday_f, df_forcing["Maximum temperature"].values)
Climate_airTempC_f = interp1d(time_nday_f, (df_forcing["Minimum temperature"].values+df_forcing["Maximum temperature"].values)/2)
Climate_solRadswskyb_f = interp1d(time_nday_f, _Rsb_Wm2)
Climate_solRadswskyd_f = interp1d(time_nday_f, _Rsd_Wm2)
Climate_airPressure_f = interp1d(time_nday_f, _p)
Climate_airRH_f = interp1d(time_nday_f, _RH)
Climate_airU_f = interp1d(time_nday_f, df_forcing["Uavg"].values)
Climate_airCO2_f = interp1d(time_nday_f, _CO2)
Climate_airO2_f = interp1d(time_nday_f, _O2)
Climate_soilTheta_f = interp1d(time_nday_f, _soilTheta)
Climate_soilTheta_z_f = interp1d(time_nday_f, _soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
Climate_nday_f = interp1d(time_nday_f, time_nday_f)   ## nday represents the ordinal day-of-year plus each simulation day (e.g. a model run starting on Jan 30 and going for 2 years will have nday=30+np.arange(2*365))

# %% [markdown]
# #### b. Rutherglen Site 1971

# %%
df_forcing = pd.read_csv("/Users/alexandernorton/ANU/Projects/DAESim/DAESIM/data/DAESim_forcing_Rutherglen_1971.csv")

SiteX = ClimateModule(CLatDeg=-36.05,CLonDeg=146.5,timezone=10)
start_doy_f = df_forcing["DOY"].values[0]
start_year_f = df_forcing["Year"].values[0]
nrundays_f = df_forcing.index.size

## Time discretisation
time_nday_f, time_doy_f, time_year_f = SiteX.time_discretisation(start_doy_f, start_year_f, nrundays=nrundays_f)
## Adjust daily time-step to represent midday on each day
time_doy_f = [time_doy_f[i]+0.5 for i in range(len(time_doy_f))]

## Time discretization for forcing data
time_index_f = pd.to_datetime(df_forcing["Date"].values)

# Define lists of sowing and harvest dates
sowing_dates = [
    pd.Timestamp(year=1971, month=5, day=11)
]

harvest_dates = [
    pd.Timestamp(year=1971, month=12, day=23)
]

_soilTheta_z = np.column_stack((
    df_forcing["Soil moisture 5 cm"].values,
    df_forcing["Soil moisture 8 cm"].values,
    df_forcing["Soil moisture 14 cm"].values,
    df_forcing["Soil moisture 22 cm"].values,
    df_forcing["Soil moisture 34 cm"].values,
    df_forcing["Soil moisture 52 cm"].values,))

Climate_doy_f = interp_forcing(time_nday_f, time_doy_f, kind="pconst") #, fill_value=(time_doy[0],time_doy[-1]))
Climate_year_f = interp_forcing(time_nday_f, time_year_f, kind="pconst") #, fill_value=(time_year[0],time_year[-1]))
Climate_airTempCMin_f = interp1d(time_nday_f, df_forcing["Minimum temperature"].values)
Climate_airTempCMax_f = interp1d(time_nday_f, df_forcing["Maximum temperature"].values)
Climate_airTempC_f = interp1d(time_nday_f, (df_forcing["Minimum temperature"].values+df_forcing["Maximum temperature"].values)/2)
Climate_solRadswskyb_f = interp1d(time_nday_f, 10*(df_forcing["Global Radiation"].values-df_forcing["Diffuse Radiation"].values))
Climate_solRadswskyd_f = interp1d(time_nday_f, 10*df_forcing["Diffuse Radiation"].values)
Climate_airPressure_f = interp1d(time_nday_f, 100*df_forcing["Pressure"].values)
Climate_airRH_f = interp1d(time_nday_f, df_forcing["Relative Humidity"].values)
Climate_airU_f = interp1d(time_nday_f, df_forcing["Uavg"].values)
Climate_airCO2_f = interp1d(time_nday_f, df_forcing["Atmospheric CO2 Concentration (bar)"].values)
Climate_airO2_f = interp1d(time_nday_f, df_forcing["Atmospheric O2 Concentration (bar)"].values)
Climate_soilTheta_z_f = interp1d(time_nday_f, _soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
Climate_nday_f = interp1d(time_nday_f, time_nday_f)   ## nday represents the ordinal day-of-year plus each simulation day (e.g. a model run starting on Jan 30 and going for 2 years will have nday=30+np.arange(2*365))

# %% [markdown]
# #### c. CSIRO Harden site

# %%
# Step 1: Download data from the two files

## Local site station data
df_site = pd.read_csv("/Users/alexandernorton/ANU/Projects/Kirkegaard_John_18_Mar_2024/data/HardenLongTermExperiment/HardenClimate.csv", skiprows=list(range(21))+[22], encoding='ISO-8859-1')

# Convert "year" and "day" into datetime
df_site["Date"] = pd.to_datetime(df_site['year'].astype(str) + df_site['day'].astype(str).str.zfill(3), format='%Y%j')

df_site["_doy"] = df_site["year"].values + df_site["day"].values/365


## Gridded data
file = "/Users/alexandernorton/ANU/Projects/DAESIM/daesim/data/DAESim_forcing_Harden_2000-2019.csv"

df_gridded = pd.read_csv(file)

df_gridded = df_gridded.rename(columns={"date":"Date"})

df_gridded['Date'] = pd.to_datetime(df_gridded['Date'])

# Add ordinal day of year (DOY) and Year variables
df_gridded["DOY"] = df_gridded["Date"].dt.dayofyear
df_gridded["Year"] = df_gridded["Date"].dt.year


# Step 2: Define the overlapping time window
start_date = max(df_site['Date'].min(), df_gridded['Date'].min())
end_date = min(df_site['Date'].max(), df_gridded['Date'].max())

# Step 3: Filter both DataFrames by the time window
df_site_tsubset = df_site[(df_site['Date'] >= start_date) & (df_site['Date'] <= end_date)]
df_gridded_tsubset = df_gridded[(df_gridded['Date'] >= start_date) & (df_gridded['Date'] <= end_date)]


# Step 4: Merge the DataFrames on the 'Date' column
merged_df = pd.merge(df_site_tsubset[['Date', 'radn', 'maxt', 'mint', 'rain', 'vp']], df_gridded_tsubset[['Date', 'Year', 'DOY', 'Soil moisture', 'Uavg']], on='Date')

# Rename columns
df_forcing_all = merged_df.rename(columns={"radn":"SRAD", "maxt":"Maximum temperature", "mint":"Minimum temperature", "rain":"Precipitation", "vp": "VPeff"})

# Interpolate sparse soil moisture values to daily
df_forcing_all["Soil moisture interp"] = df_forcing_all["Soil moisture"].interpolate('quadratic')

f = "/Users/alexandernorton/ANU/Resources/Gladish et al. (2021) Supplementary Material - tabulated - Cluster 4 Mean.csv"
df_soil = pd.read_csv(f, skiprows=2)

depths_z = df_soil["Depth (mm)"].values
LL15_z = df_soil["LL15 (mm/mm)"].values
CLL_z = df_soil["CLL (mm/mm)"].values
DUL_z = df_soil["DUL (mm/mm)"].values

soiltheta = df_forcing_all["Soil moisture interp"].values

# Estimate min/max soil moisture (could use empirical percentiles)
soiltheta_min = np.nanmin(soiltheta)
soiltheta_max = np.nanmax(soiltheta)
#print("Soil moisture range: %1.1f to %1.1f mm" % (soiltheta_min, soiltheta_max))

# Compute relative wetness index (W)
W = (soiltheta - soiltheta_min) / (soiltheta_max - soiltheta_min)
W = np.clip(W, 0, 1)   # Ensure values stay within [0,1]

# Initialize an array to store the time-varying soil moisture profile
theta_profile = np.zeros((len(soiltheta), len(depths_z)))

for i, w_t in enumerate(W):
    theta_profile[i, :] = CLL_z + w_t * (DUL_z - CLL_z)
    # theta_profile[i, :] = LL15_z + w_t * (DUL_z - LL15_z)

# Add the soil moisture profile data into the forcing dataframe
df_forcing_all["Soil moisture interp 0-100 mm"] = theta_profile[:,0]
df_forcing_all["Soil moisture interp 100-200 mm"] = theta_profile[:,1]
df_forcing_all["Soil moisture interp 200-300 mm"] = theta_profile[:,2]
df_forcing_all["Soil moisture interp 300-400 mm"] = theta_profile[:,3]
df_forcing_all["Soil moisture interp 400-500 mm"] = theta_profile[:,4]
df_forcing_all["Soil moisture interp 500-600 mm"] = theta_profile[:,5]
df_forcing_all["Soil moisture interp 600-700 mm"] = theta_profile[:,6]
df_forcing_all["Soil moisture interp 700-800 mm"] = theta_profile[:,7]
df_forcing_all["Soil moisture interp 800-900 mm"] = theta_profile[:,8]
df_forcing_all["Soil moisture interp 900-1000 mm"] = theta_profile[:,9]
df_forcing_all["Soil moisture interp 1000-1100 mm"] = theta_profile[:,10]
df_forcing_all["Soil moisture interp 1100-1200 mm"] = theta_profile[:,11]
df_forcing_all["Soil moisture interp 1200-1300 mm"] = theta_profile[:,12]
df_forcing_all["Soil moisture interp 1300-1400 mm"] = theta_profile[:,13]
df_forcing_all["Soil moisture interp 1400-1500 mm"] = theta_profile[:,14]
df_forcing_all["Soil moisture interp 1500-1600 mm"] = theta_profile[:,15]
df_forcing_all["Soil moisture interp 1600-1700 mm"] = theta_profile[:,16]
df_forcing_all["Soil moisture interp 1700-1800 mm"] = theta_profile[:,17]
df_forcing_all["Soil moisture interp 1800-1900 mm"] = theta_profile[:,18]
df_forcing_all["Soil moisture interp 1900-2000 mm"] = theta_profile[:,19]

# Select one year of data only
# df_forcing = df_forcing_all.loc[(df_forcing_all["Date"] >= "2008-01-01") & (df_forcing_all["Date"] <= "2008-12-31")]
df_forcing = df_forcing_all.loc[(df_forcing_all["Date"] >= "2012-01-01") & (df_forcing_all["Date"] <= "2012-12-31")]

## Harden CSIRO site location-34.52194, 148.30472
## NOTES:
SiteX = ClimateModule(CLatDeg=-34.52194,CLonDeg=148.30472,timezone=10)
start_doy_f = df_forcing["DOY"].values[0]
start_year_f = df_forcing["Year"].values[0]
nrundays_f = df_forcing.index.size

## Time discretisation
time_nday_f, time_doy_f, time_year_f = SiteX.time_discretisation(start_doy_f, start_year_f, nrundays=nrundays_f)
## Adjust daily time-step to represent midday on each day
time_doy_f = [time_doy_f[i]+0.5 for i in range(len(time_doy_f))]

## Time discretization for forcing data
time_index_f = pd.to_datetime(df_forcing["Date"].values)

# Define lists of sowing and harvest dates
sowing_dates = [
    pd.Timestamp(year=2012, month=5, day=20),  # First season
    # pd.Timestamp(year=2008, month=5, day=14),
]

harvest_dates = [
    pd.Timestamp(year=2012, month=12, day=10),  # First season
    # pd.Timestamp(year=2008, month=11, day=17),
]

## Make some assumption about the fraction of diffuse radiation
diffuse_fraction = 0.2

## Shortwave radiation at surface (convert MJ m-2 d-1 to W m-2)
_Rsb_Wm2 = (1-diffuse_fraction) * df_forcing["SRAD"].values * 1e6 / (60*60*24)
_Rsd_Wm2 = diffuse_fraction * df_forcing["SRAD"].values * 1e6 / (60*60*24)

## Create synthetic data for other forcing variables
_p = 101325*np.ones(nrundays_f)
_es = SiteX.compute_sat_vapor_pressure_daily(df_forcing["Minimum temperature"].values,df_forcing["Maximum temperature"].values)
_RH = SiteX.compute_relative_humidity(df_forcing["VPeff"].values/10,_es/1000)
_RH[_RH > 100] = 100
_CO2 = 400*(_p/1e5)*1e-6     ## carbon dioxide partial pressure (bar)
_O2 = 209000*(_p/1e5)*1e-6   ## oxygen partial pressure (bar)
_soilTheta =  0.35*np.ones(nrundays_f)   ## volumetric soil moisture content (m3 m-3)
# _soilTheta = f_soilTheta_norm

## Create a multi-layer soil moisture forcing dataset by using the scaled soil moisture following the a Gladish et al. (2021) profile
_soilTheta_z = np.column_stack((
    df_forcing["Soil moisture interp 0-100 mm"].values,
    df_forcing["Soil moisture interp 100-200 mm"].values,
    df_forcing["Soil moisture interp 200-300 mm"].values,
    df_forcing["Soil moisture interp 300-400 mm"].values,
    df_forcing["Soil moisture interp 400-500 mm"].values,
    df_forcing["Soil moisture interp 500-600 mm"].values,
    df_forcing["Soil moisture interp 600-700 mm"].values,
    df_forcing["Soil moisture interp 700-800 mm"].values,
    df_forcing["Soil moisture interp 800-900 mm"].values,
    df_forcing["Soil moisture interp 900-1000 mm"].values,
    df_forcing["Soil moisture interp 1000-1100 mm"].values,
    df_forcing["Soil moisture interp 1100-1200 mm"].values,
    df_forcing["Soil moisture interp 1200-1300 mm"].values,
    df_forcing["Soil moisture interp 1300-1400 mm"].values,
    df_forcing["Soil moisture interp 1400-1500 mm"].values,
    df_forcing["Soil moisture interp 1500-1600 mm"].values,
    df_forcing["Soil moisture interp 1600-1700 mm"].values,
    df_forcing["Soil moisture interp 1700-1800 mm"].values,
    df_forcing["Soil moisture interp 1800-1900 mm"].values,
    df_forcing["Soil moisture interp 1900-2000 mm"].values,
    ))

## Convert discrete to continuous forcing data
Climate_doy_f = interp_forcing(time_nday_f, time_doy_f, kind="pconst") #, fill_value=(time_doy[0],time_doy[-1]))
Climate_year_f = interp_forcing(time_nday_f, time_year_f, kind="pconst") #, fill_value=(time_year[0],time_year[-1]))
Climate_airTempCMin_f = interp1d(time_nday_f, df_forcing["Minimum temperature"].values)
Climate_airTempCMax_f = interp1d(time_nday_f, df_forcing["Maximum temperature"].values)
Climate_airTempC_f = interp1d(time_nday_f, (df_forcing["Minimum temperature"].values+df_forcing["Maximum temperature"].values)/2)
Climate_solRadswskyb_f = interp1d(time_nday_f, _Rsb_Wm2)
Climate_solRadswskyd_f = interp1d(time_nday_f, _Rsd_Wm2)
Climate_airPressure_f = interp1d(time_nday_f, _p)
Climate_airRH_f = interp1d(time_nday_f, _RH)
Climate_airU_f = interp1d(time_nday_f, df_forcing["Uavg"].values)
Climate_airCO2_f = interp1d(time_nday_f, _CO2)
Climate_airO2_f = interp1d(time_nday_f, _O2)
Climate_soilTheta_f = interp1d(time_nday_f, _soilTheta)
Climate_soilTheta_z_f = interp1d(time_nday_f, _soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
Climate_nday_f = interp1d(time_nday_f, time_nday_f)   ## nday represents the ordinal day-of-year plus each simulation day (e.g. a model run starting on Jan 30 and going for 2 years will have nday=30+np.arange(2*365))

# %% [markdown]
# ### Step 4. Initialise the Model

# %%
## Check that the first sowing date and last harvest date are within the available forcing data period
# Apply validation to sowing and harvest dates
SiteX.validate_event_dates(sowing_dates, time_index_f, event_name="Sowing")
SiteX.validate_event_dates(harvest_dates, time_index_f, event_name="Harvest")

# Find steps for all sowing and harvest dates in the forcing data
sowing_steps_f = SiteX.find_event_steps(sowing_dates, time_index_f)
harvest_steps_f = SiteX.find_event_steps(harvest_dates, time_index_f)

# 3. Generate time axis for ODE Solver, which is set to start at the first sowing event and end just after the last harvest event
time_axis = time_nday_f[sowing_steps_f[0]:harvest_steps_f[-1]+2]
# Corresponding date time index
time_index = time_index_f[sowing_steps_f[0]:harvest_steps_f[-1]+2]

## Determine the ODE solver reset days
sowing_steps_itax = SiteX.find_event_steps(sowing_dates, time_index)
harvest_steps_itax = SiteX.find_event_steps(harvest_dates, time_index)
reset_days_itax = SiteX.find_event_steps(sowing_dates+harvest_dates, time_index)
reset_days_itax.sort()
reset_days_tax = list(time_axis[reset_days_itax])

## Determine the sowing and harvest days and years for the DAESIM2 Management module
sowingDays = np.floor(Climate_doy_f(time_axis[sowing_steps_itax]))
sowingYears = np.floor(Climate_year_f(time_axis[sowing_steps_itax]))
harvestDays = np.floor(Climate_doy_f(time_axis[harvest_steps_itax]))
harvestYears = np.floor(Climate_year_f(time_axis[harvest_steps_itax]))

# %%
ManagementX = ManagementModule(cropType="Wheat",sowingDays=sowingDays,harvestDays=harvestDays,sowingYears=sowingYears,harvestYears=harvestYears)

PlantDevX = PlantGrowthPhases(
    phases=["germination", "vegetative", "spike", "anthesis", "grainfill", "maturity"],
    # gdd_requirements=[120,800,250,200,300,200],
    gdd_requirements=[50,800,280,150,300,300],
    # vd_requirements=[0, 40, 0, 0, 0, 0],
    vd_requirements=[0, 30, 0, 0, 0, 0],
    allocation_coeffs = [
        [0.2, 0.1, 0.7, 0.0, 0.0],
        [0.5, 0.1, 0.4, 0.0, 0.0],
        [0.30, 0.4, 0.30, 0.0, 0.0],
        [0.30, 0.4, 0.30, 0.0, 0.0],
        [0.1, 0.1, 0.1, 0.7, 0.0],
        [0.1, 0.1, 0.1, 0.7, 0.0]
    ],
    turnover_rates = [[0.001,  0.001, 0.001, 0.0, 0.0],
                      [0.01, 0.002, 0.008, 0.0, 0.0],
                      [0.01, 0.002, 0.008, 0.0, 0.0],
                      [0.01, 0.002, 0.008, 0.0, 0.0],
                      [0.033, 0.016, 0.033, 0.0002, 0.0],
                      [0.10, 0.033, 0.10, 0.0002, 0.0]])

BoundLayerX = BoundaryLayerModule(Site=SiteX)
LeafX = LeafGasExchangeModule2(Site=SiteX,Jmax_opt_rVcmax=0.89,Jmax_opt_rVcmax_method="log")
CanopyX = CanopyLayers(nlevmlcan=3)
CanopyRadX = CanopyRadiation(Canopy=CanopyX)
CanopyGasExchangeX = CanopyGasExchange(Leaf=LeafX,Canopy=CanopyX,CanopyRad=CanopyRadX)
# SoilLayersX = SoilLayers(nlevmlsoil=2,z_max=2.0)
# SoilLayersX = SoilLayers(nlevmlsoil=6,z_max=0.66,z_top=0.10,discretise_method="horizon",
#                          z_horizon=[0.06, 0.06, 0.06, 0.10, 0.10, 0.28],
#                         Psi_e=[-1.38E-03, -1.38E-03, -1.38E-03, -1.32E-03, -2.58E-03, -0.960E-03],
#                         b_soil = [4.74, 4.74, 4.74, 6.77, 8.17, 10.73],
#                          K_sat = [29.7, 29.7, 29.7, 25.2, 13.9, 40.9],
#                          soilThetaMax = [0.12, 0.12, 0.12, 0.20, 0.3, 0.4])
SoilLayersX = SoilLayers(nlevmlsoil=20,z_max=2.0,discretise_method="horizon",
                         z_horizon=0.1*np.ones(20),
                         b_soil = 4.74, Psi_e = -1.38E-03, K_sat = 29.7,  ## Sandy loam soil type (Duursma et al., 2008)
                         soilThetaMax=df_soil["SAT (mm/mm) "].values)

# PlantCH2OX = PlantCH2O(Site=SiteX,SoilLayers=SoilLayersX,CanopyGasExchange=CanopyGasExchangeX,BoundaryLayer=BoundLayerX,maxLAI=6.0,ksr_coeff=1000,SLA=0.030)
PlantCH2OX = PlantCH2O(Site=SiteX,SoilLayers=SoilLayersX,CanopyGasExchange=CanopyGasExchangeX,BoundaryLayer=BoundLayerX,maxLAI=6.5,ksr_coeff=1500,SLA=0.02,sf=1.0,Psi_f=-5.0)
PlantAllocX = PlantOptimalAllocation(Plant=PlantCH2OX)
PlantX = PlantModuleCalculator(
    Site=SiteX,
    Management=ManagementX,
    PlantDev=PlantDevX,
    PlantCH2O=PlantCH2OX,
    PlantAlloc=PlantAllocX,
    GDD_method="nonlinear",
    GDD_Tbase=0.0,
    GDD_Topt=22.5,
    GDD_Tupp=35.0,
    hc_max_GDDindex=sum(PlantDevX.gdd_requirements[0:2])/PlantDevX.totalgdd,
    d_r_max=1.5,
    Vmaxremob=3.0,
    Kmremob=0.5,
    remob_phase=["grainfill","maturity"],
    specified_phase="spike",
    grainfill_phase=["grainfill","maturity"],
)

# %%
## Define the callable calculator that defines the right-hand-side ODE function
PlantXCalc = PlantX.calculate

Model = ODEModelSolver(calculator=PlantXCalc, states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_start=time_axis[0], log_diagnostics=True)

forcing_inputs = [Climate_solRadswskyb_f,
                  Climate_solRadswskyd_f,
                  Climate_airTempCMin_f,
                  Climate_airTempCMax_f,
                  Climate_airPressure_f,
                  Climate_airRH_f,
                  Climate_airCO2_f,
                  Climate_airO2_f,
                  Climate_airU_f,
                  Climate_soilTheta_z_f,
                  Climate_doy_f,
                  Climate_year_f]

reset_days = reset_days_tax

zero_crossing_indices = [4,5,6]

# %% [markdown]
# ### Step 5. Re-Initialise and Run the Model Using the Sampled Input Sets
#
# Once you have the samples, you can evaluate your model at each of the sampled input points. This typically involves running a loop where each set of sampled inputs is passed through the model, and the output (or outputs) is recorded.

# %%
# write_to_nc = True

# # # Location/site of the simulations
# # xsite = "Milgadara_2021_test_single"
# # xsite = "Rutherglen_1971_test_single"
# xsite = "Harden2012_test_single"
# title = "DAESIM2-Plant Sobol Sensitivity Analysis: "+xsite
# description = "DAESIM2-Plant Sobol sensitivity analysis for the site " + xsite

# # Path for writing outputs to file
# filepath_write = "/Users/alexandernorton/ANU/Projects/DAESim/daesim2-analysis/results/Sobol/"

# # Create input_data for model run
# input_data = [ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices]

# # Create output array for target variables
# Mpx = []

# # Sub-sample the FAST samples
# nsamples = 100
# isamples = np.array([1]) #np.arange(0,param_values.shape[0],param_values.shape[0]/nsamples,dtype=int)

# iparamset = 0
# # for iparamset in isamples:
# param_set = param_values[iparamset]
# # Select parameter set from FAST samples
# nparamset = iparamset+1

# # Call the function that updates parameters, runs the model and returns selected outputs
# model_output = fastsa.update_and_run_model(param_set, PlantX, input_data, parameters_df, problem)

# # Separate model output into FAST target variables and model diagnostics
# Mpxi, diagnostics = model_output[0], model_output[1]

# # Append target variables to output
# Mpx.append(np.insert(Mpxi, 0, nparamset))

# if write_to_nc:
#     # Write diagnostics to file
#     nsigfigures = len(str(np.shape(param_values)[0]))   # number of significant figures of parameter sets to insert zero-padded nparamset into filename 
#     # File name for writing outputs to file
#     filename_write = f"FAST_results_{xsite}_paramset{nparamset:0{nsigfigures}}.nc"
#     paramset = param_values[iparamset]
#     daesim_io_write_diag_to_nc(PlantX, diagnostics, 
#                                    filepath_write, filename_write, 
#                                    time_index,
#                                    problem=problem,
#                                    param_values=paramset,
#                                    nc_attributes={"title": title, "description": description})


# %%

# %%
write_to_nc = True

# Location/site of the simulations
# xsite = "Milgadara_2018_test"
# xsite = "Rutherglen_1971_test2"
xsite = "Harden_2012_test"
title = "DAESIM2-Plant Sobol Sensitivity Analysis: "+xsite
description = "DAESIM2-Plant Sobol sensitivity analysis for the site " + xsite

# Path for writing outputs to file
filepath_write = "/Users/alexandernorton/ANU/Projects/DAESim/daesim2-analysis/results/Sobol/"

# Create input_data for model run
input_data = [ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices]

# Create output array for target variables
Mpx_column_headers = "nparamset,W_P_peakW,W_L_peakW,W_R_peakW,W_S_peakW,W_S_spike0,W_S_anth0,GPP_int_seas,NPP_int_seas,Rml_int_seas,Rmr_int_seas,Rg_int_seas,trflux_int_seas,FCstem2grain_int_seas,NPP2grain_int_seas,E_int_seas,LAI_peakW,W_spike_anth1,GY_mature,Sdpot_mature,GN_mature"   # to save as header in csv file. N.B. must match the custom output in the update_and_run_model function
Mpx = []

# Option: Run a sub-sample of the FAST samples
# nsamples = 300
# isamples = np.arange(0,param_values.shape[0],param_values.shape[0]/nsamples,dtype=int)

# Option: Run the full set of FAST samples
nsamples = param_values.shape[0]

for iparamset in range(nsamples):
    param_set = param_values[iparamset]
    # Select parameter set from FAST samples
    nparamset = iparamset+1

    # Call the function that updates parameters, runs the model and returns selected outputs
    model_output = fastsa.update_and_run_model(param_set, PlantX, input_data, parameters_df, problem)
    
    # Separate model output into FAST target variables and model diagnostics
    Mpxi, diagnostics = model_output[0], model_output[1]

    # Append target variables to output, inserting nparamset as the first element
    Mpx.append(np.insert(Mpxi, 0, nparamset))
    
    if write_to_nc:
        # Write diagnostics to file
        nsigfigures = len(str(np.shape(param_values)[0]))   # number of significant figures of parameter sets to insert zero-padded nparamset into filename 
        # File name for writing outputs to file
        filename_write = f"FAST_results_{xsite}_paramset{nparamset:0{nsigfigures}}.nc"
        paramset = param_values[iparamset]
        daesim_io_write_diag_to_nc(PlantX, diagnostics, 
                                   filepath_write, filename_write, 
                                   time_index,
                                   problem=problem,
                                   param_values=paramset,
                                   nc_attributes={"title": title, "description": description})

# Write the target variables data to a csv file
fname_target_variables = f"{filepath_write}FAST_results_{xsite}_target_variables.csv"
np.savetxt(fname_target_variables, np.array(Mpx), delimiter=",")


# %%

# %%
param_set = [Mpx[ipar][0] for ipar in range(len(Mpx))]
W_P_peakW = [Mpx[ipar][1] for ipar in range(len(Mpx))]
W_L_peakW = [Mpx[ipar][2] for ipar in range(len(Mpx))]
W_R_peakW = [Mpx[ipar][3] for ipar in range(len(Mpx))]
W_S_peakW = [Mpx[ipar][4] for ipar in range(len(Mpx))]
W_S_spike0 = [Mpx[ipar][5] for ipar in range(len(Mpx))]
W_S_anth0 = [Mpx[ipar][6] for ipar in range(len(Mpx))]
GPP_int_seas = [Mpx[ipar][7] for ipar in range(len(Mpx))]
NPP_int_seas = [Mpx[ipar][8] for ipar in range(len(Mpx))]
Rml_int_seas = [Mpx[ipar][9] for ipar in range(len(Mpx))]
Rmr_int_seas = [Mpx[ipar][10] for ipar in range(len(Mpx))]
Rg_int_seas = [Mpx[ipar][11] for ipar in range(len(Mpx))]
trflux_int_seas = [Mpx[ipar][12] for ipar in range(len(Mpx))]
FCstem2grain_int_seas = [Mpx[ipar][13] for ipar in range(len(Mpx))]
NPP2grain_int_seas = [Mpx[ipar][14] for ipar in range(len(Mpx))]
E_int_seas = [Mpx[ipar][15] for ipar in range(len(Mpx))]
LAI_peakW = [Mpx[ipar][16] for ipar in range(len(Mpx))]
W_spike_anth1 = [Mpx[ipar][17] for ipar in range(len(Mpx))]
GY_mature = [Mpx[ipar][18] for ipar in range(len(Mpx))]
Sdpot_mature = [Mpx[ipar][19] for ipar in range(len(Mpx))]
GN_mature = [Mpx[ipar][20] for ipar in range(len(Mpx))]

# %%

# %% [markdown]
# ## Evaluate Sobol' Results

# %%
file = "/Users/alexandernorton/ANU/Projects/DAESim/daesim2-analysis/results/Rutherglen_1971_test_Mpx.npy"

Mpx = np.load(file)

# %%
# Define column names
column_names = [
    "Simulation_Number", "W_P_peakW", "W_L_peakW", "W_R_peakW", "W_S_peakW",
    "W_S_spike0", "W_S_anth0", "GPP_int_seas", "NPP_int_seas", "Rml_int_seas",
    "Rmr_int_seas", "Rg_int_seas", "trflux_int_seas", "FCstem2grain_int_seas",
    "NPP2grain_int_seas", "E_int_seas", "LAI_peakW", "W_spike_anth1",
    "GY_mature", "Sdpot_mature", "GN_mature"
]

# Create Pandas DataFrame
df = pd.DataFrame(Mpx, columns=column_names)

# Sort the DataFrame by 'Simulation_Number'
df_sorted = df.sort_values(by="Simulation_Number").reset_index(drop=True)

# If needed, convert back to a NumPy array
sorted_arr = df_sorted.to_numpy()

# Display first few rows
print(df_sorted.head())



# %% [markdown]
# ### ***Temporary fix of code bug***

# %%
df_sorted = df_sorted.rename(columns={
    "W_R_peakW": "W_S_peakW", 
    "W_S_peakW": "W_R_peakW"})

# %% [markdown]
# ### Analyse Sensitivity Indices
#
# When you run fast.analyze(problem, Y), the output is a dictionary Si containing various sensitivity indices. Here’s a breakdown of the key components:
#
# - S1 (First-order sensitivity index): This value represents the direct effect of each parameter on the model output, ignoring interactions with other parameters. It measures the fraction of output variance that can be attributed to each parameter alone.
#
# Interpretation: A high S1 value (close to 1) for a parameter indicates that this parameter has a strong influence on the model output, independent of other parameters.
#
#
# - ST (Total-order sensitivity index): This represents the total effect of each parameter on the model output, including interactions with other parameters. It shows how much the output variance would be reduced if the parameter were fixed at any value within its range.
#
# Interpretation: A high ST value (close to 1) indicates that the parameter not only has a significant individual effect but also contributes to interactions with other parameters.
#
# Names: Corresponds to the parameter names in the problem['names'] list, so you can map each index to a parameter.

# %%
color_S1 = 'lightsteelblue'
color_ST = '0.25'

# %%
# Perform the Sobol' sensitivity analysis

## Selected variable
variable_key = "GY_mature"
fast_output_array = df_sorted[variable_key].values

# Step 5: Perform Sobol' Sensitivity Analysis
Si1 = sobol.analyze(problem, fast_output_array, print_to_console=True)

# Print the first-order and total-order sensitivity indices
print("--%s--" % variable_key)
print("First-order sensitivity indices (S1):", Si1['S1'])
print("Total-order sensitivity indices (ST):", Si1['ST'])
print()


# Create a bar plot for S1 (first-order) and ST (total-order) sensitivity indices
parameters = problem['names']
S1 = Si1['S1']
ST = Si1['ST']

# Bar width
bar_width = 0.35

# Set up the figure
fig, ax = plt.subplots(1,1,figsize=(4,3))

# Plot first-order indices

ax.bar(np.arange(len(parameters)), S1, bar_width, label='First-order', color=color_S1)
ax.bar([p + bar_width for p in np.arange(len(parameters))], ST, bar_width, label='Total-order', color=color_ST)
# Add labels and title
ax.set_ylabel('Sensitivity Index')
ax.set_title('First-Order Sensitivity: %s' % variable_key)
ax.legend()
ax.set_xticks(np.arange(len(parameters))+0.5*bar_width, parameters, rotation=90)

# Show plot
plt.tight_layout()
plt.show()

# %%
# Perform the FAST sensitivity analysis

## Selected variable
variable_key = "Root:Shoot Ratio"
fast_output_array = df_sorted["W_R_peakW"].values/(df_sorted["W_L_peakW"].values+df_sorted["W_S_peakW"].values)
Si1 = fast.analyze(problem, fast_output_array)

# Print the first-order and total-order sensitivity indices
print("--%s--" % variable_key)
print("First-order sensitivity indices (S1):", Si1['S1'])
print("Total-order sensitivity indices (ST):", Si1['ST'])
print()


# Create a bar plot for S1 (first-order) and ST (total-order) sensitivity indices
parameters = problem['names']
S1 = Si1['S1']
ST = Si1['ST']

# Bar width
bar_width = 0.35

# Set up the figure
fig, ax = plt.subplots(1,1,figsize=(4,3))

# Plot first-order indices

ax.bar(np.arange(len(parameters)), S1, bar_width, label='First-order', color=color_S1)
ax.bar([p + bar_width for p in np.arange(len(parameters))], ST, bar_width, label='Total-order', color=color_ST)
# Add labels and title
ax.set_ylabel('Sensitivity Index')
ax.set_title('First-Order Sensitivity: %s' % variable_key)
ax.legend()
ax.set_xticks(np.arange(len(parameters))+0.5*bar_width, parameters, rotation=90)

# Show plot
plt.tight_layout()
plt.show()

# %%
# Example list of variable keys
variable_keys = ["NPP_int_seas", "LAI_peakW", "GY_mature"]  # Replace with your actual variables

# Create a figure with subplots
num_vars = len(variable_keys)
fig, axes = plt.subplots(nrows=num_vars, ncols=1, figsize=(6, 3 * num_vars), sharex=True)

# Ensure axes is iterable even for a single subplot
if num_vars == 1:
    axes = [axes]

# Loop over each variable key
i = 0
for ax, variable_key in zip(axes, variable_keys):
    # Extract data
    fast_output_array = df_sorted[variable_key].values
    Si = fast.analyze(problem, fast_output_array)

    # Print sensitivity indices
    print(f"--{variable_key}--")
    print("First-order sensitivity indices (S1):", Si['S1'])
    print("Total-order sensitivity indices (ST):", Si['ST'])
    print()

    # Bar plot data
    parameters = problem['names']
    S1 = Si['S1']
    ST = Si['ST']
    bar_width = 0.35
    x = np.arange(len(parameters))

    # Plot first-order and total-order sensitivity indices
    ax.bar(x, S1, bar_width, label='First-order', color=color_S1)
    ax.bar(x + bar_width, ST, bar_width, label='Total-order', color=color_ST)

    # Formatting
    ax.set_ylabel('Sensitivity Index')
    ax.set_title(f'Sensitivity Analysis: {variable_key}')
    if i == 0:
        ax.legend()
    ax.set_xticks(x + 0.5 * bar_width)
    ax.set_xticklabels(parameters, rotation=90)

    i += 1

# Adjust layout and show the figure
plt.tight_layout()
plt.show()

# %%

# %%
from scipy.stats import pearsonr

# Function to compute partial correlation
def partial_corr(x, y, covar):
    if covar.shape[1] > 0:  # Ensure covariates exist
        try:
            # Fit linear models: regress x and y on covariates
            beta_x = np.linalg.lstsq(covar, x, rcond=None)[0]
            beta_y = np.linalg.lstsq(covar, y, rcond=None)[0]

            # Compute residuals
            x_resid = x - covar @ beta_x
            y_resid = y - covar @ beta_y
        except np.linalg.LinAlgError:
            print(f"Warning: SVD did not converge for parameter. Skipping it.")
            return np.nan  # Return NaN for failed cases
    else:
        x_resid, y_resid = x, y  # No covariates, use original values
    
    # Compute Pearson correlation on residuals
    return pearsonr(x_resid, y_resid)[0]


# %%

# %%
# Convert data to Pandas DataFrame
df = pd.DataFrame(param_values, columns=[f'Param_{i+1}' for i in range(param_values.shape[1])])
df['Output'] = df_sorted["GY_mature"].values

# Compute PCC for each parameter
pcc_results = {}
for param in df.columns[:-1]:
    covariates = df.drop(columns=[param, 'Output']).values
    pcc_value = partial_corr(df[param].values, df['Output'].values, covariates)
    pcc_results[param] = pcc_value
    print(f"{param}: PCC = {pcc_value:.3f}")

# Plot PCC values
plt.figure(figsize=(10, 5))
# pcc_df.sort_values(by='PCC', ascending=True).plot(kind='barh', legend=False)
pcc_df.plot(kind='barh', legend=False)
plt.xlabel('Partial Correlation Coefficient')
plt.title('Partial Correlation Coefficients (PCC)')
plt.yticks(ticks=np.arange(pcc_df.index.size), labels=parameter_names) 
plt.show()


# %%

# %%

# %%

# %%

# %%

# %%

# %%
plt.plot(param_values[:,0])

# %%

# %%

# %%
model_variable = "GY_mature"

iparam = 9

index0, index1 = (0+iparam)*num_samples, (1+iparam)*num_samples

xparam_vector = param_values[index0:index1, iparam]
xmodeloutput_vector = df_sorted[model_variable].values[index0:index1]

plt.scatter(xparam_vector, xmodeloutput_vector)
plt.title(parameter_names[iparam])
plt.show()

nsubplots = 12
p = np.arange(13)
mask = (p == iparam)
p_plot = p[~mask]
fig, axes = plt.subplots(3,4,figsize=(3*3,4*2))

for i, ax in enumerate(axes.flatten()):
    ip = p_plot[i]
    ax.scatter(param_values[index0:index1, ip], df_sorted[model_variable].values[index0:index1], marker='o', alpha=0.5, s=10)
    ax.set_title(parameter_names[ip])

plt.tight_layout()

# %%
model_variable = "GY_mature"

iparam = 12

index0, index1 = (0+iparam)*num_samples, (1+iparam)*num_samples

xparam_vector = param_values[index0:index1, iparam]
xmodeloutput_vector = df_sorted[model_variable].values[index0:index1]

plt.scatter(xparam_vector, xmodeloutput_vector)
plt.title(parameter_names[iparam])
plt.show()


index0, index1 = (0+iparam)*num_samples, (1+iparam)*num_samples

xparam_vector = param_values[index0:index1, 0]
xmodeloutput_vector = df_sorted[model_variable].values[index0:index1]

plt.scatter(xparam_vector, xmodeloutput_vector)
plt.title(parameter_names[iparam])
plt.show()


# %%

# %%

# %%
# Convert data to Pandas DataFrame
df_pcc = df_sorted.copy(deep=True)
df_pcc['Output'] = np.zeros(df_pcc.index.size)

# Compute PCC for each parameter (ignoring the first column which should represent "Simulation number")
for param in df_pcc.columns[7:12]:
    print(param)
    covariates = df_pcc[[p for p in df.columns if p not in [param, 'Output']]].values
    pcc_value = partial_corr(df_pcc[param].values, df_pcc['Output'].values, covariates)
    print(f"{param}: PCC = {pcc_value:.3f}")

# %%
