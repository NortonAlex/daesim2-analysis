# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pandas import Timestamp
from daesim.utils import ODEModelSolver

from daesim.climate import *
from daesim.plantgrowthphases import PlantGrowthPhases
from daesim.management import ManagementModule
from daesim.plant_1000_thermaltime import PlantModuleCalculator

from daesim2_analysis.parameters import Parameters
from daesim2_analysis.forcing_data import ForcingData

# %% [markdown]
# ## Calibrate Growing-Degree-Days Module
#
# Model Simulation Steps:
#
# 1. Set up the model inputs for model forward run, including forcing data time-series and parameter vector.
# 2. Initialise the model with model inputs.
# 3. Run the model give the inputs.
# 4. Collate the key model outputs that we compare with observations in a format that can be easily compared (ready for a calibration/optimzation algorithm).                                                                                                            
#
# Calibration Routine:
#
# 1. Define parameter vector (free parameters) including their respective min and max ranges (or their probability density functions).
# 2. Define the observation vector including uncertainties.
# 3. Define the optimization algorithm to do the parameter optimization, including the objective function to be minimised. 
# 4. Set up the interface between the optimization algorithm and the model simulation. 

# %% [markdown]
# ### 1. Define the model analysis problem (sensitivity, optimization or uncertainty)
#
# Need to define the number of parameters, their names, and the ranges (bounds) over which the analysis will be performed. For example, if your model has three parameters (e.g., param1, param2, param3), you would define the problem in the following format:
#
# problem = {
#     "num_vars": 6,    # Number of input parameters
#     "names": ["Vcmax_opt", "Jmax_opt", "g1"],   # Parameter names
#     "bounds": [[10e-6, 200e-6], [20e-6, 400e-6], [2, 8]]    # Parameter ranges
# }
#

# %%
parameters = Parameters.__from_file__("../parameters/PlantDevCalibration.json")

# %% [markdown]
# ### 2. Initialise the study site and import forcing data

# %%
f = "../DAESIM_data/DAESim_forcing_data/DAESim_forcing_Milgadara_2019_v1.csv"

df_forcing = pd.read_csv(f)

# %%
SiteX = ClimateModule(CLatDeg=-34.389,CLonDeg=148.469,timezone=10)

ForcingDataX = ForcingData(
    SiteX=SiteX,
    sowing_dates=[Timestamp(year=2019,month=4,day=10)],
    harvest_dates=[Timestamp(year=2019,month=12,day=26)],
    df=df_forcing,
    df_type="3",
    zero_crossing_indices=[0,1],
        )

# %% [markdown]
# ### 3. Initialise the Model

# %%
# Insert the parameters into the model class
ManagementX = ManagementModule(cropType="Canola", sowingDays=ForcingDataX.sowing_days, harvestDays=ForcingDataX.harvest_days, sowingYears=ForcingDataX.sowing_years, harvestYears=ForcingDataX.harvest_years)
PlantDevX = PlantGrowthPhases(
    phases=["germination", "vegetative", "anthesis", "grainfill", "maturity"],
    gdd_requirements=[120, 500, 200, 350, 200],
    vd_requirements=[0, 25, 0, 0, 0],
    allocation_coeffs=[
        [0.2, 0.1, 0.7, 0.0, 0.0],   # Phase 1
        [0.5, 0.1, 0.4, 0.0, 0.0],   # Phase 2
        [0.25, 0.5, 0.25, 0.0, 0.0], # Phase 3
        [0.1, 0.1, 0.1, 0.7, 0.0],   # Phase 4
        [0.1, 0.1, 0.1, 0.7, 0.0]    # Phase 5
    ],
    turnover_rates = [
        [0.001, 0.001, 0.001, 0.0, 0.0],  # Phase 1
        [0.01,  0.002, 0.01,  0.0, 0.0],  # Phase 2
        [0.02,  0.002, 0.04,  0.0, 0.0],  # Phase 3
        [0.10,  0.008, 0.10,  0.0, 0.0],  # Phase 4
        [0.50,  0.017, 0.50,  0.0, 0.0]   # Phase 5
    ]    ## Turnover rates per pool and developmental phase (days-1))
)

PlantX = PlantModuleCalculator(
    Site=SiteX,
    Management=ManagementX,
    PlantDev=PlantDevX,
    GDD_method="linear1",
    GDD_Tbase=0.0,
    GDD_Tupp=25.0,
)

# %%
## Define the callable calculator that defines the right-hand-side ODE function
PlantXCalc = PlantX.calculate

Model = ODEModelSolver(calculator=PlantXCalc, states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_start=ForcingDataX.time_axis[0], log_diagnostics=True)


# %% [markdown]
# ### 4. Define the observations

# %%
## Observations

## Example 1: Milgadara, Rocky East, Canola TT Y44, 2021
observables_names = ["start of vegetative", "start of flowering", "start of grainfill", "harvest"]
observables_units = ["ordinal day of year", "ordinal day of year", "ordinal day of year", "ordinal day of year"]
observables_values = [143, 234, 262, 339]
observables_uncertainty = [5, 5, 5, 5]

## Example 2: Milgadara, Horse Paddock Big, Canola, 2019
observables_names = ["start of vegetative", "start of flowering", "start of grainfill", "harvest"]
observables_units = ["ordinal day of year", "ordinal day of year", "ordinal day of year", "ordinal day of year"]
observables_values = [167, 258, 265, 307]    # 16-Jun, 3-Nov, 15-Sep, 22-Sep (2021)
observables_uncertainty = [5, 5, 5, 5]

## Put these into a combined pandas dataframe
target_df = pd.DataFrame({
    "Name": observables_names,
    "Units": observables_units,
    "Values": observables_values,
    "Uncertainty": observables_uncertainty,
})

# Put the observations vectors into arrays
y = target_df["Values"].values
U_y = target_df["Uncertainty"].values

# %% [markdown]
# ### 5. Set up the interface between the optimisation algorithm and the model

# %%
from scipy.optimize import differential_evolution

from daesim2_analysis.run import update_attribute, update_attribute_in_phase


# %%
def run_model_and_get_outputs(Plant, ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices):
    ## Define the callable calculator that defines the right-hand-side ODE function
    PlantCalc = Plant.calculate
    
    Model = ODEModelSolver(calculator=PlantCalc, states_init=[0.0, 0.0], time_start=time_axis[0], log_diagnostics=True)
    
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

    # Add np.nan to the end of each array in the dictionary to represent the last time point in the time_axis (corresponds to the last time point of the state vector)
    for key in diagnostics:
        if key == "t":
            diagnostics[key] = np.append(diagnostics[key], res["t"][-1])
        else:
            diagnostics[key] = np.append(diagnostics[key], np.nan)

    # Add state variables to the diagnostics dictionary
    diagnostics["GDD"] = res["y"][0,:]
    diagnostics["VD"] = res["y"][1,:]

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
    
    # Observation Operator
    # Calculate model-equivalent observations from model run output

    # # Diagnose time indexes when developmental phase transitions occur
    ngrowing_seasons = (len(Plant.Management.sowingDays) if (isinstance(Plant.Management.sowingDays, int) == False) else 1)
    if ngrowing_seasons > 1:
        # print("Multiple sowing and harvest events occur. Only returning results for first growing season.")
        ## ignore any time steps before first sowing event and after last harvest event
        it_sowing = np.where(time_axis == reset_days[0])[0][0]  #sowing_steps_itax[0]
        
        if Plant.Management.harvestDays is not None:
            it_harvest = np.where(time_axis == reset_days[1])[0][0]  #harvest_steps_itax[0]   # np.where(np.floor(Climate_doy_f(time_axis)) == Plant.Management.harvestDay)[0][0]
        else:
            it_harvest = -1   # if there is no harvest day specified, we just take the last day of the simulation. 
    else:
        # print("Just one sowing event and one harvest event occurs. Returning results for first (and only) growing season.")
        ## ignore any time steps before first sowing event and after last harvest event
        it_sowing = np.where(time_axis == reset_days[0])[0][0]  #sowing_steps_itax[0]
        
        if Plant.Management.harvestDays is not None:
            it_harvest = np.where(time_axis == reset_days[1])[0][0]  #harvest_steps_itax[0]   # np.where(np.floor(Climate_doy_f(time_axis)) == Plant.Management.harvestDay)[0][0]
        else:
            it_harvest = -1   # if there is no harvest day specified, we just take the last day of the simulation. 

    # Diagnose time indexes when developmental phase transitions occur

    # Convert the array to a numeric type, handling mixed int and float types
    idevphase = diagnostics["idevphase_numeric"]   #[it_sowing:it_harvest+1]
    valid_mask = ~np.isnan(idevphase)
    
    # Identify all transitions (number-to-NaN, NaN-to-number, or number-to-different-number)
    it_phase_transitions = np.where(
        ~valid_mask[:-1] & valid_mask[1:] |  # NaN-to-number
        valid_mask[:-1] & ~valid_mask[1:] |  # Number-to-NaN
        (valid_mask[:-1] & valid_mask[1:] & (np.diff(idevphase) != 0))  # Number-to-different-number
    )[0] + 1
    
    # Time index for the end of the maturity phase
    if PlantX.PlantDev.phases.index('maturity') in idevphase:
        it_mature = np.where(idevphase == PlantX.PlantDev.phases.index('maturity'))[0][-1]    # Index for end of maturity phase
    elif PlantX.Management.harvestDays is not None: 
        it_mature = it_harvest    # Maturity developmental phase not completed, so take harvest as the end of growing season
    else:
        it_mature = -1    # if there is no harvest day specified, we just take the last day of the simulation. 

    # it_sowing = np.where(time_axis == Plant.Management.sowingDay)[0][0]
    # if Plant.Management.harvestDay is not None:
    #     it_harvest = np.where(time_axis == Plant.Management.harvestDay)[0][0]
    # else:
    #     it_harvest = -1   # if there is no harvest day specified, we just take the last day of the simulation. 

    # # Convert the array to a numeric type, handling mixed int and float types
    # idevphase = diagnostics["idevphase_numeric"]
    # valid_mask = ~np.isnan(idevphase)
    
    # # Identify all transitions (number-to-NaN, NaN-to-number, or number-to-different-number)
    # it_phase_transitions = np.where(
    #     ~valid_mask[:-1] & valid_mask[1:] |  # NaN-to-number
    #     valid_mask[:-1] & ~valid_mask[1:] |  # Number-to-NaN
    #     (valid_mask[:-1] & valid_mask[1:] & (np.diff(idevphase) != 0))  # Number-to-different-number
    # )[0] + 1
    
    # # Time index for the end of the maturity phase
    # if Plant.PlantDev.phases.index('maturity') in idevphase:
    #     it_mature = np.where(idevphase == Plant.PlantDev.phases.index('maturity'))[0][-1]    # Index for end of maturity phase
    # elif Plant.Management.harvestDay is not None: 
    #     it_mature = it_harvest    # Maturity developmental phase not completed, so take harvest as the end of growing season
    # else:
    #     it_mature = -1    # if there is no harvest day specified, we just take the last day of the simulation. 

    # import pdb; pdb.set_trace()
    # Filter out transitions that occur on or before the sowing day
    # it_phase_transitions = [t for t in it_phase_transitions if time_axis[t] > time_axis[it_sowing+1]]
    it_phase_transitions = [t for t in it_phase_transitions if t > int(it_sowing+1)]
    # Filter out transitions that occur after the maturity or harvest day
    # it_phase_transitions = [t for t in it_phase_transitions if time_axis[t] <= time_axis[it_mature]]
    it_phase_transitions = [t for t in it_phase_transitions if t <= it_mature]

    # Developmental phase indexes
    igermination = Plant.PlantDev.phases.index("germination")
    ivegetative = Plant.PlantDev.phases.index("vegetative")
    if Plant.Management.cropType == "Wheat":
        ispike = Plant.PlantDev.phases.index("spike")
    ianthesis = Plant.PlantDev.phases.index("anthesis")
    igrainfill = Plant.PlantDev.phases.index("grainfill")
    imaturity = Plant.PlantDev.phases.index("maturity")

    ip = np.where(diagnostics['idevphase'][it_phase_transitions] == Plant.PlantDev.phases.index('vegetative'))[0][0]
    tdoy_vegetative = time_axis[it_phase_transitions[ip]]   # ordinal day-of-year at transition point into vegetative phase
    if Plant.PlantDev.phases.index('anthesis') in idevphase[it_sowing+1:it_harvest+1]:
        ip = np.where(diagnostics['idevphase'][it_phase_transitions] == Plant.PlantDev.phases.index('anthesis'))[0][0]
        tdoy_anth0 = time_axis[it_phase_transitions[ip]]   # ordinal day-of-year at transition point into anthesis phase
    else:
        tdoy_anth0 = time_axis[it_harvest]
    if Plant.PlantDev.phases.index('grainfill') in idevphase[it_sowing+1:it_harvest+1]:
        ip = np.where(diagnostics['idevphase'][it_phase_transitions] == Plant.PlantDev.phases.index('grainfill'))[0][0]
        tdoy_anth1 = time_axis[it_phase_transitions[ip]]   # ordinal day-of-year at transition point into grainfill stage (out of anthesis phase)
    else:
        tdoy_anth1 = time_axis[it_harvest]
    tdoy_harvest = time_axis[it_harvest]   # ordinal day-of-year at harvest
    
    # import pdb; pdb.set_trace()
    # Model output (of observables) given the parameter vector p
    # - this is the model output that we compare to observations and use to calibrate the parameters
    M_p = np.array([
        tdoy_vegetative, 
        tdoy_anth0, 
        tdoy_anth1,
        tdoy_harvest,
    ])

    return M_p

def model_function(params, model_instance, input_data, param_info):
    # Update the model class with the new parameters
    # Plant.PlantDev.gdd_requirements = [params[ip_GDD_germ],params[ip_GDD_veg],params[ip_GDD_ant],params[ip_GDD_gf],params[ip_GDD_mat]]
    # Plant.Management.sowingDay = params[ip_sowingDay]
    # Plant.Management.harvestDay = params[ip_harvestDay]

    # Collate input data to pass to model run function
    ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices, time_nday_f, time_doy_f, time_year_f = input_data
    
    for idx, value in enumerate(params):
        param_name = param_info["Name"].values[idx]
        param_path = param_info["Module Path"].values[idx]
        full_path = f"{param_path}.{param_name}"
        phase_specific = param_info["Phase Specific"].values[idx]
        
        if phase_specific:
            # Handle phase-specific parameters
            phase = param_info["Phase"].values[idx]
            update_attribute_in_phase(model_instance, full_path, value, phase)
        else:
            if (param_name == "sowingDays") or (param_name == "harvestDays"):
                # Update parameters that must be defined as a list type
                update_attribute(model_instance, full_path, [value])
            else:
                # Update regular parameters
                update_attribute(model_instance, full_path, value)

        # Make sure the solver knows about the sowing and harvest dates as well (to reset the state variables like GDD and VD)
        if (param_name == "sowingDays") or (param_name == "harvestDays"):
            # Find value of time_nday_f where time_doy_f == sowingDay and time_year_f == sowingYear.
            sowingDay, sowingYear = model_instance.Management.sowingDays, model_instance.Management.sowingYears
            sowing_nday = time_nday_f[(np.floor(time_doy_f) == sowingDay) & (np.array(time_year_f) == sowingYear)]
            
            # Find value of time_nday_f where time_doy_f == sowingDay and time_year_f == sowingYear.
            harvestDay, harvestYear = model_instance.Management.harvestDays, model_instance.Management.harvestYears
            harvest_nday = time_nday_f[(np.floor(time_doy_f) == harvestDay) & (np.array(time_year_f) == harvestYear)]
            
            # Set reset_days to be the updated sowing and harvest nday
            reset_days = [sowing_nday[0], harvest_nday[0]]
     
    model_output = run_model_and_get_outputs(model_instance, ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices)

    return model_output

# Define an objective function to minimize
def objective_function_mse(params, observations, Plant, input_data, param_info):
    """
    Objective function using mean squared error (MSE). 

    Notes
    -----
    """
    # Round the parameters off to integers as these parameters must be integers representing day-of-year
    int_params = np.round(params).astype(int)
    # Calculate model outputs
    model_outputs = model_function(int_params, Plant, input_data, param_info)
    # Calculate the error (e.g., mean squared error) TODO: Include model-obs uncertainties here too
    error = np.mean((model_outputs - observations) ** 2)
    return error

# Define an objective function to minimize
def objective_function_wls(params, observations, observation_unc_sigma, Plant, input_data, param_info):
    """
    Objective function using weighted least squares (WLS). 

    Notes
    -----
    The cost function, $J$, is defined using a weighted least-squares as follows: 
    
    $J = \sum_{i=1}^n \frac{(M_i(p) - y_i)^2}{\sigma^2}$
    
    Where $y$ is the vector of observations, $M$ is the vector model predicted observables 
    given parameter set $p$, and $\sigma$ is the observation errors (assumed to include 
    structural model errors). Note that this formulation ignores the priors. 
    """
    # Round the parameters off to integers as these parameters must be integers representing day-of-year
    int_params = np.round(params).astype(int)
    # Calculate model outputs
    model_outputs = model_function(int_params, Plant, input_data, param_info)
    # Calculate the error as the weighted Least Squares: 
    # 
    # Error is the model - observed difference squared, normalised by the uncertainty (as a variance), and summed over all obs
    error = np.mean( ((model_outputs - observations) ** 2) / (observation_unc_sigma**2))
    print(f"Current error: {error}")
    return error



# %% [markdown]
# #### - Test the wrapper functions used to run the model

# %%
input_data = [ODEModelSolver, ForcingDataX.time_axis, ForcingDataX.inputs, ForcingDataX.reset_days, ForcingDataX.zero_crossing_indices, ForcingDataX.time_nday_f, ForcingDataX.time_doy_f, ForcingDataX.time_year_f]

param_info = parameters.df
params = parameters.df["Initial Value"].values

model_function(params, PlantX, input_data, param_info)

# %% [markdown]
# ### 6. Run the optimisation

# %%
input_data = [ODEModelSolver, ForcingDataX.time_axis, ForcingDataX.inputs, ForcingDataX.reset_days, ForcingDataX.zero_crossing_indices, ForcingDataX.time_nday_f, ForcingDataX.time_doy_f, ForcingDataX.time_year_f]

# Creating a list of (min, max) tuples
param_bounds = list(zip(parameters.df["Min"].values, parameters.df["Max"].values))

# Run the differential evolution algorithm
result = differential_evolution(objective_function_wls, bounds=param_bounds, args=(y, U_y, PlantX, input_data, parameters.df), popsize=5, tol=0.05, maxiter=500)

# Display results
print("Optimal Parameters:", result.x)
print("Minimum Objective Function Value:", result.fun)

# %% [markdown]
# ### 7. Display optimised parameters and model-observed comparison

# %%

# %%
np.around(np.array(result.x)).astype(int)

# %%
optimised_parameters = np.around(np.array(result.x)).astype(int)

for i, p in enumerate(parameters.df["Name"].values):
    print(i,p,' initial =',parameters.df["Initial Value"].values[i],', optimised =',optimised_parameters[i])


# %%
# Add new column to parameters class for optimised parameters
df = parameters.to_dataframe().copy()
df['Optimised Value'] = optimised_parameters
object.__setattr__(parameters, 'df', df)

# Write to new JSON file
parameters.to_file("../parameters/PlantDevCalibrated.json")

# %% [markdown]
# #### - Run the model with the optimised parameters

# %%
PlantX.Management.sowingDays = parameters.df.loc[parameters.df['Name'] == 'sowingDays']['Optimised Value'].values
PlantX.Management.harvestDays = parameters.df.loc[parameters.df['Name'] == 'harvestDays']['Optimised Value'].values

# %%
input_data = [ODEModelSolver, ForcingDataX.time_axis, ForcingDataX.inputs, ForcingDataX.reset_days, ForcingDataX.zero_crossing_indices, ForcingDataX.time_nday_f, ForcingDataX.time_doy_f, ForcingDataX.time_year_f]

## Run with optimised parameters
M_x = model_function(optimised_parameters, PlantX, input_data, parameters.df)

## Run with initial values for comparison
M_x0 = model_function(parameters.df["Initial Value"].values, PlantX, input_data, parameters.df)

# %%
fig, axes = plt.subplots(1,1,figsize=(4,3))

axes.scatter(np.arange(len(target_df)), M_x0, label="Initial Values", color="C0", marker="o", alpha=0.5)
axes.scatter(np.arange(len(target_df)), M_x, label="Optimized", c="C1", marker="o", alpha=0.5)
axes.scatter(np.arange(len(target_df)), observables_values, c='k', marker='x', alpha=0.75, label="Observed")
axes.legend()
axes.set_xticks(np.arange(M_x.size)) 
axes.set_xticklabels(observables_names,rotation=45)
axes.set_ylabel("Ordinal Day of Year")

# %%
num_params = len(parameters.df)
fig, axes = plt.subplots(1, num_params, figsize=(num_params * 2, 3))

# If only one parameter, axes won't be a list, so handle that case
if num_params == 1:
    axes = [axes]

# Loop over each parameter and create a subplot
for i, ax in enumerate(axes):
    param_name = parameters.df.loc[i, "Name"]
    param_unit = parameters.df.loc[i, "Unit"]
    param_min = parameters.df.loc[i, "Min"]
    param_max = parameters.df.loc[i, "Max"]
    param_init = parameters.df.loc[i, "Initial Value"]
    solution_value = result.x[i]
    param_population = result.population[:,i]

    # Plot the initial value as a scatter point
    ax.scatter([1], [param_init], color='0.5', s=50, marker='x')
    # Plot the solution value as a scatter point
    ax.scatter([1], [solution_value], color='red', s=50, marker='x')
    # Plot the box and whisker plot for the parameter "population"
    bplot = ax.boxplot(param_population, widths=0.6, patch_artist=True, boxprops=dict(facecolor='r', color='k', alpha=0.25))
    
    # Customize the subplot
    ax.set_title(param_name)
    ax.set_ylabel(f'{param_unit}')
    ax.set_ylim(param_min, param_max)
    ax.set_xticks([])  # Remove x-ticks for a cleaner look

# Adjust layout to prevent overlapping
plt.tight_layout()

# Display the plot
plt.show()

# %%

# %% [markdown]
# ### Run Full Biophysical Model
#
# **Here you can now use the optimised parameters (assigned to the Management and PlantDev modules) as input to the full DAESIM2 biophysical model and run the full model to get simulated carbon, water, yield, etc.**

# %%
