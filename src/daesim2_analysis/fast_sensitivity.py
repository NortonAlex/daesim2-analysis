"""
Analysis helper functions to support DAESIM2 analysis, sensitivity, calibration, etc.
"""

from typing import Any
from netCDF4 import date2num, Dataset
from datetime import datetime, timedelta
import time
import subprocess
import numpy as np
import pandas as pd

# from daesim2_analysis.run import update_attribute, update_attribute_in_phase, update_and_run_model

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

    return M_p, diagnostics


def update_attribute(obj: Any, path: str, new_value: Any) -> None:
    """
    Update the attribute specified by the path with a new value.

    Parameters:
    obj (Any): The main class instance to update.
    path (str): Dot-separated path to the attribute (e.g., "Site.temperature").
    new_value (Any): The new value to set for the attribute.
    """
    attributes = path.split('.')
    # Traverse the path except for the last attribute
    if attributes[0] == "":
        # attribute is in the parent class (i.e. not in a sub-class of the obj)
        setattr(obj, attributes[-1], new_value)
    else:
        for attr in attributes[:-1]:
            obj = getattr(obj, attr)
            # Set the new value to the last attribute in the path
        setattr(obj, attributes[-1], new_value)

def update_attribute_in_phase(obj: Any, path: str, new_value: Any, phase: str) -> None:
    """
    Update the attribute specified by the path and a specific developmental phase with a new value.

    Parameters:
    obj (Any): The main class instance to update.
    path (str): Dot-separated path to the attribute (e.g., "Site.temperature").
    new_value (Any): The new value to set for the attribute.
    phase (str): The plant developmental phase of the attribute
    """
    attributes = path.split('.')
    # Traverse the path except for the last attribute
    for attr in attributes[:-1]:
        obj = getattr(obj, attr)
    # Make a copy of the phase-specific values
    new_phase_values = getattr(obj, attributes[-1])
    # Get the phase index for the value, assuming the obj is the correct module that includes the "phases" attribute
    phase_i = obj.phases.index(phase)
    # Update only the phase-specific value
    new_phase_values[phase_i] = new_value
    # Set the new phase value list to the last attribute in the path
    setattr(obj, attributes[-1], new_phase_values)

def update_and_run_model(param_values, model_instance, input_data, param_info, salib_problem):
    """
    Updates model parameters, runs the biophysical model, and returns the outputs.

    Parameters
    ----------
    param_values : array_like
        Vector of parameter values to be updated in the model.
    model_instance : object
        Instance of the model class that includes methods for running simulations.
    input_data : tuple
        A tuple containing input data required for the model run, structured as:
        (ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices, time_nday_f, time_doy_f, time_year_f).
    param_info : pd.DataFrame
        DataFrame with parameter metadata, including:
        - "Name": Parameter name as used in the model.
        - "Module Path": Path to the model attribute to update (e.g., "Plant.Parameters").
        - "Phase Specific": Boolean indicating whether the parameter is specific to a growth phase.
        - "Phase": (Optional) Index of the growth phase if "Phase Specific" is True.
    salib_problem : dict
        Dictionary required for SALib sensitivity analysis. Includes:
        - "num_vars": Number of parameters.
        - "names": List of parameter names.
        - "bounds": List of parameter bounds.

    Returns
    -------
    array_like
        Model outputs after running the simulation.
    """

    # Unpack input data to pass to model run function
    ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices, time_nday_f, time_doy_f, time_year_f = input_data

    # Update the model instance with the provided parameters
    for idx, value in enumerate(param_values):
        param_name = salib_problem["names"][idx]
        param_path = param_info.loc[idx]["Module Path"]
        full_path = f"{param_path}.{param_name}"

        if param_info.loc[idx]["Phase Specific"]:
            # Handle phase-specific parameters
            phase = param_info.loc[idx]['Phase']
            update_attribute_in_phase(model_instance, full_path, value, phase)
        else:
            if (param_name == "sowingDays") or (param_name == "harvestDays"):
                # Update parameters that must be defined as a list type
                if isinstance(value, list):
                    update_attribute(model_instance, full_path, value)
                else:
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

    # Run the model and get the outputs
    model_outputs = run_model_and_get_outputs(
        model_instance, ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices
    )

    return model_outputs


# def run_model_and_get_outputs(Plant, ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices):
#     ## Define the callable calculator that defines the right-hand-side ODE function
#     PlantCalc = Plant.calculate
    
#     Model = ODEModelSolver(calculator=PlantCalc, states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_start=time_axis[0], log_diagnostics=True)

#     Model.reset_diagnostics()
    
#     ## Run the model solver
#     res = Model.run(
#         time_axis=time_axis,
#         forcing_inputs=forcing_inputs,
#         solver="euler",
#         zero_crossing_indices=zero_crossing_indices,
#         reset_days=reset_days,
#     )

#     # Convert the defaultdict to a regular dictionary
#     _diagnostics = dict(Model.diagnostics)
#     # Convert each list in the dictionary to a NumPy array
#     diagnostics = {key: np.array(value) for key, value in _diagnostics.items()}

#     # Convert the array to a numeric type, handling mixed int and float types
#     diagnostics['idevphase_numeric'] = np.array(diagnostics['idevphase'],dtype=np.float64)
    
#     # In the model idevphase can equal None but that is not useable in post-processing, so we set None values to np.nan
#     diagnostics["idevphase_numeric"][diagnostics["idevphase"] == None] = np.nan

#     ## Conversion notes: When _E units are mol m-2 s-1, multiply by molar mass H2O to get g m-2 s-1, divide by 1000 to get kg m-2 s-1, multiply by 60*60*24 to get kg m-2 d-1, and 1 kg m-2 d-1 = 1 mm d-1. 
#     ## Noting that 1 kg of water is equivalent to 1 liter (L) of water (because the density of water is 1000 kg/m³), and 1 liter of water spread over 1 square meter results in a depth of 1 mm
#     diagnostics["E_mmd"] = diagnostics["E"]*18.015/1000*(60*60*24)
    
#     # Turnover rates per pool
#     _tr_Leaf = np.zeros(diagnostics['t'].size)
#     _tr_Root = np.zeros(diagnostics['t'].size)
#     _tr_Stem = np.zeros(diagnostics['t'].size)
#     _tr_Seed = np.zeros(diagnostics['t'].size)
#     for it, t in enumerate(diagnostics['t']):
#         if np.isnan(diagnostics['idevphase_numeric'][it]):
#             tr_ = Plant.PlantDev.turnover_rates[-1]
#             _tr_Leaf[it] = tr_[Plant.PlantDev.ileaf]
#             _tr_Root[it] = tr_[Plant.PlantDev.iroot]
#             _tr_Stem[it] = tr_[Plant.PlantDev.istem]
#             _tr_Seed[it] = tr_[Plant.PlantDev.iseed]
#         else:
#             tr_ = Plant.PlantDev.turnover_rates[diagnostics['idevphase'][it]]
#             _tr_Leaf[it] = tr_[Plant.PlantDev.ileaf]
#             _tr_Root[it] = tr_[Plant.PlantDev.iroot]
#             _tr_Stem[it] = tr_[Plant.PlantDev.istem]
#             _tr_Seed[it] = tr_[Plant.PlantDev.iseed]
    
#     diagnostics['tr_Leaf'] = _tr_Leaf
#     diagnostics['tr_Root'] = _tr_Root
#     diagnostics['tr_Stem'] = _tr_Stem
#     diagnostics['tr_Seed'] = _tr_Seed

#     # Add np.nan to the end of each array in the dictionary to represent the last time point in the time_axis (corresponds to the last time point of the state vector)
#     for key in diagnostics:
#         if key == "t":
#             diagnostics[key] = np.append(diagnostics[key], res["t"][-1])
#         else:
#             diagnostics[key] = np.append(diagnostics[key], np.nan)

    
#     # Add state variables to the diagnostics dictionary
#     diagnostics["Cleaf"] = res["y"][0,:]
#     diagnostics["Cstem"] = res["y"][1,:]
#     diagnostics["Croot"] = res["y"][2,:]
#     diagnostics["Cseed"] = res["y"][3,:]
#     diagnostics["Bio_time"] = res["y"][4,:]
#     diagnostics["VRN_time"] = res["y"][5,:]
#     diagnostics["Cstate"] = res["y"][7,:]
#     diagnostics["Cseedbed"] = res["y"][8,:]

#     # Add forcing inputs to diagnostics dictionary
#     for i,f in enumerate(forcing_inputs):
#         ni = i+1
#         if f(time_axis[0]).size == 1:
#             fstr = f"forcing {ni:02}"
#             diagnostics[fstr] = f(time_axis)
#         elif f(time_axis[0]).size > 1:
#             # this forcing input has levels/layers (e.g. multilayer soil moisture)
#             nz = f(time_axis[0]).size
#             for iz in range(nz):
#                 fstr = f"forcing {ni:02} z{iz}"
#                 diagnostics[fstr] = f(time_axis)[:,iz]

    # Time indexing for model output data, to determine outputs at specific times in the growing season
    itax_sowing, itax_mature, itax_harvest, itax_phase_transitions = Plant.Site.time_index_growing_season(time_index, diagnostics['idevphase_numeric'], Plant.Management, Plant.PlantDev)

    # Time of peak biomass (with and without the seed pool)
    total_carbon_t = res["y"][Plant.PlantDev.ileaf,:] + res["y"][Plant.PlantDev.istem,:] + res["y"][Plant.PlantDev.iroot,:] + res["y"][Plant.PlantDev.iseed,:]
    total_carbon_exclseed_t = res["y"][Plant.PlantDev.ileaf,:] + res["y"][Plant.PlantDev.istem,:] + res["y"][Plant.PlantDev.iroot,:]

    itax_peakbiomass = np.argmax(total_carbon_t[itax_sowing:itax_harvest+1]) + itax_sowing
    itax_peakbiomass_exclseed = np.argmax(total_carbon_exclseed_t[itax_sowing:itax_harvest+1]) + itax_sowing

    # Carbon and Water
    W_P_peakW = total_carbon_t[itax_peakbiomass_exclseed]/Plant.PlantCH2O.f_C    # Total dry biomass at peak biomass
    W_L_peakW = res["y"][Plant.PlantDev.ileaf,itax_peakbiomass_exclseed]/Plant.PlantCH2O.f_C    # Leaf dry biomass at peak biomass
    W_S_peakW = res["y"][Plant.PlantDev.istem,itax_peakbiomass_exclseed]/Plant.PlantCH2O.f_C    # Stem dry biomass at peak biomass
    W_R_peakW = res["y"][Plant.PlantDev.iroot,itax_peakbiomass_exclseed]/Plant.PlantCH2O.f_C    # Root dry biomass at peak biomass
    if Plant.Management.cropType == "Wheat":
        ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('spike'))[0][0]
        W_S_spike0 = res["y"][Plant.PlantDev.istem,itax_phase_transitions[ip]]/Plant.PlantCH2O.f_C    # Stem dry biomass at start of spike
    ip = np.where(diagnostics['idevphase'][itax_phase_transitions] == Plant.PlantDev.phases.index('anthesis'))[0][0]
    W_S_anth0 = res["y"][Plant.PlantDev.istem,itax_phase_transitions[ip]]/Plant.PlantCH2O.f_C    # Stem dry biomass at start of anthesis
    GPP_int_seas = np.sum(diagnostics['GPP'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal GPP
    NPP_int_seas = np.sum(diagnostics['NPP'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal NPP
    Rml_int_seas = np.sum(diagnostics['Rml'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal Rml
    Rmr_int_seas = np.sum(diagnostics['Rmr'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal Rmr
    Rg_int_seas = np.sum(diagnostics['Rg'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal Rg
    trflux_int_seas = np.sum(diagnostics['trflux_total'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal turnover losses
    FCstem2grain_int_seas = np.sum(diagnostics['F_C_stem2grain'][itax_sowing:itax_mature+1])    # Total (integrated) remobilisation to grain
    _Cflux_NPP2grain = diagnostics['u_Seed'] * diagnostics['NPP']    # NPP carbon allocation flux to grain
    NPP2grain_int_seas = np.sum(_Cflux_NPP2grain[itax_sowing:itax_mature+1])    # Total (integrated) NPP carbon allocation to grain
    E_int_seas = np.nansum(diagnostics['E_mmd'][itax_sowing:itax_mature+1])    # Total (integrated) seasonal transpiration
    LAI_peakW = diagnostics['LAI'][itax_peakbiomass_exclseed]    # Leaf area index at peak biomass
    
    # Grain Production
    W_spike_anth1 = res["y"][7,itax_mature]/Plant.PlantCH2O.f_C    # Spike dry biomass at anthesis
    GY_mature = res["y"][Plant.PlantDev.iseed,itax_mature]/Plant.PlantCH2O.f_C    # Grain yield at maturity
    GY_harvest = res["y"][Plant.PlantDev.iseed,itax_harvest]/Plant.PlantCH2O.f_C    # Grain yield at harvest
    Sdpot_mature = diagnostics['S_d_pot'][itax_mature]    # Potential seed density (grain number density) at maturity
    Sdpot_harvest = diagnostics['S_d_pot'][itax_harvest]    # Potential seed density (grain number density) at harvest
    if Plant.Management.cropType == "Wheat":
        GN_mature = res["y"][Plant.PlantDev.iseed,itax_mature]/Plant.PlantCH2O.f_C/Plant.W_seedTKW0    # Actual grain number at maturity
        GN_harvest = res["y"][Plant.PlantDev.iseed,itax_harvest]/Plant.PlantCH2O.f_C/Plant.W_seedTKW0    # Actual grain number at harvest
    
    # Model output (of observables) given the parameter vector p
    # - this is the model output that we compare to observations and use to calibrate the parameters
    M_p = np.array([
        W_P_peakW, 
        W_L_peakW,
        W_R_peakW,
        W_S_peakW,
        W_S_spike0,
        W_S_anth0,
        GPP_int_seas,
        NPP_int_seas,
        Rml_int_seas,
        Rmr_int_seas,
        Rg_int_seas,
        trflux_int_seas,
        FCstem2grain_int_seas,
        NPP2grain_int_seas,
        E_int_seas,
        LAI_peakW,
        W_spike_anth1,
        GY_mature,
        Sdpot_mature,
        GN_mature,
    ])

    return M_p, diagnostics

# def update_attribute(obj: Any, path: str, new_value: Any) -> None:
#     """
#     Update the attribute specified by the path with a new value.

#     Parameters:
#     obj (Any): The main class instance to update.
#     path (str): Dot-separated path to the attribute (e.g., "Site.temperature").
#     new_value (Any): The new value to set for the attribute.
#     """
#     attributes = path.split('.')
#     # Traverse the path except for the last attribute
#     if attributes[0] == "":
#         # attribute is in the parent class (i.e. not in a sub-class of the obj)
#         setattr(obj, attributes[-1], new_value)
#     else:
#         for attr in attributes[:-1]:
#             obj = getattr(obj, attr)
#             # Set the new value to the last attribute in the path
#         setattr(obj, attributes[-1], new_value)

# def update_attribute_in_phase(obj: Any, path: str, new_value: Any, phase: str) -> None:
#     """
#     Update the attribute specified by the path and a specific developmental phase with a new value.

#     Parameters:
#     obj (Any): The main class instance to update.
#     path (str): Dot-separated path to the attribute (e.g., "Site.temperature").
#     new_value (Any): The new value to set for the attribute.
#     phase (str): The plant developmental phase of the attribute
#     """
#     attributes = path.split('.')
#     # Traverse the path except for the last attribute
#     for attr in attributes[:-1]:
#         obj = getattr(obj, attr)
#     # Make a copy of the phase-specific values
#     new_phase_values = getattr(obj, attributes[-1])
#     # Get the phase index for the value, assuming the obj is the correct module that includes the "phases" attribute
#     phase_i = obj.phases.index(phase)
#     # Update only the phase-specific value
#     new_phase_values[phase_i] = new_value
#     # Set the new phase value list to the last attribute in the path
#     setattr(obj, attributes[-1], new_phase_values)

# def model_function(params, Plant, input_data, params_info, problem):
#     """
#     Function to update the model parameters, run the model and return the outputs

#     Parameters
#     ----------
#     params : array_like
#         Vector of parameter values
#     Plant : Callable
#         Model class that is updated and run
#     input_data : list
#         List of input data for running the model (e.g. forcing_inputs, ODEModelSolver inputs)
#     params_info : pd.DataFrame()
#         Pandas dataframe with parameter information including parameter name ("Name") and model class path ("Module Path")
#     problem : dict
#         Dictionary required as input to SALib FAST sampler (includes keys: "num_vars", "names", "bounds")

#     Returns
#     -------
#     array_like
#         Output from model run
#     """
#     # Update the model class with the new parameters
#     p_x = params
#     for j, p_xi in enumerate(p_x):
#         # print("Parameter:",problem["names"][j])
#         # print("Module path:",parameters_df.loc[parameters_df["Name"] == problem["names"][j]]["Module Path"].values)
#         # print("Full module path:",parameters_df.loc[parameters_df["Name"] == problem["names"][j]]["Module Path"].values + "." + str(problem["names"][j])
#         # print("Value:",p_xi)
#         path = params_info.loc[params_info["Name"] == problem["names"][j]]["Module Path"].values[0] + "." + str(problem["names"][j])
#         if params_info.loc[params_info["Name"] == problem["names"][j]]["Phase Specific"].values[0] == True:
#             # if the parameter is specific to a phase in PlantDev then we handle it differently, as the attribute is a list of values corresponding to each phase but we only want to update one value in that list
#             phase = params_info.loc[params_info["Name"] == problem["names"][j]]["Phase"].values[0]
#             update_attribute_in_phase(Plant, path, p_xi, phase)
#         else:
#             update_attribute(Plant, path, p_xi)

#     # Collate input data to pass to model run function
#     ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices = input_data
#     # N.B. If any of the parameters are sowingDay or harvestDay, we must make sure the solver knows about the sowing and harvest dates as well (to reset the state variables like GDD and VD)
#     # reset_days = [Plant.Management.sowingDay, Plant.Management.harvestDay]
    
#     model_output = run_model_and_get_outputs(Plant, ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices)

#     return model_output

# def update_and_run_model(param_values, model_instance, input_data, param_info, salib_problem):
#     """
#     Updates model parameters, runs the biophysical model, and returns the outputs.

#     Parameters
#     ----------
#     param_values : array_like
#         Vector of parameter values to be updated in the model.
#     model_instance : object
#         Instance of the model class that includes methods for running simulations.
#     input_data : tuple
#         A tuple containing input data required for the model run, structured as:
#         (ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices).
#     param_info : pd.DataFrame
#         DataFrame with parameter metadata, including:
#         - "Name": Parameter name as used in the model.
#         - "Module Path": Path to the model attribute to update (e.g., "Plant.Parameters").
#         - "Phase Specific": Boolean indicating whether the parameter is specific to a growth phase.
#         - "Phase": (Optional) Index of the growth phase if "Phase Specific" is True.
#     salib_problem : dict
#         Dictionary required for SALib sensitivity analysis. Includes:
#         - "num_vars": Number of parameters.
#         - "names": List of parameter names.
#         - "bounds": List of parameter bounds.

#     Returns
#     -------
#     array_like
#         Model outputs after running the simulation.
#     """
#     # Update the model instance with the provided parameters
#     for idx, value in enumerate(param_values):
#         param_name = salib_problem["names"][idx]
#         param_path = param_info.loc[param_info["Name"] == param_name, "Module Path"].values[0]
#         full_path = f"{param_path}.{param_name}"

#         if param_info.loc[param_info["Name"] == param_name, "Phase Specific"].values[0]:
#             # Handle phase-specific parameters
#             phase = param_info.loc[param_info["Name"] == param_name, "Phase"].values[0]
#             update_attribute_in_phase(model_instance, full_path, value, phase)
#         else:
#             # Update regular parameters
#             update_attribute(model_instance, full_path, value)

#     # Unpack input data
#     ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices = input_data

#     # Run the model and get the outputs
#     model_outputs = run_model_and_get_outputs(
#         model_instance, ODEModelSolver, time_axis, time_index, forcing_inputs, reset_days, zero_crossing_indices
#     )

#     return model_outputs
